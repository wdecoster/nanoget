import logging
from functools import reduce
import nanoget.utils as ut
import pandas as pd
import sys
import pysam
import re
from Bio import SeqIO
import concurrent.futures as cfutures
from itertools import repeat


def process_summary(summaryfile, **kwargs):
    """Extracting information from an albacore summary file.

    Only reads which have a >0 length are returned.

    The fields below may or may not exist, depending on the type of sequencing performed.
    Fields 1-14 are for 1D sequencing.
    Fields 1-23 for 2D sequencing.
    Fields 24-27, 2-5, 22-23 for 1D^2 (1D2) sequencing
    Fields 28-38 for barcoded workflows
     1  filename
     2  read_id
     3  run_id
     4  channel
     5  start_time
     6  duration
     7  num_events
     8  template_start
     9  num_events_template
    10  template_duration
    11  num_called_template
    12  sequence_length_template
    13  mean_qscore_template
    14  strand_score_template
    15  complement_start
    16    num_events_complement
    17    complement_duration
    18    num_called_complement
    19    sequence_length_complement
    20    mean_qscore_complement
    21    strand_score_complement
    22    sequence_length_2d
    23    mean_qscore_2d
    24    filename1
    25    filename2
    26    read_id1
    27    read_id2
    28    barcode_arrangement
    29    barcode_score
    30    barcode_full_arrangement
    31    front_score
    32    rear_score
    33    front_begin_index
    34    front_foundseq_length
    35    rear_end_index
    36    rear_foundseq_length
    37    kit
    38    variant
    """
    logging.info(
        f"Nanoget: Collecting metrics from summary file {summaryfile} for {kwargs['readtype']} sequencing"
    )
    ut.check_existance(summaryfile)
    if kwargs["readtype"] == "1D":
        cols = [
            "channel",
            "start_time",
            "duration",
            "sequence_length_template",
            "mean_qscore_template",
        ]
    elif kwargs["readtype"] in ["2D", "1D2"]:
        cols = ["channel", "start_time", "duration", "sequence_length_2d", "mean_qscore_2d"]
    if kwargs["barcoded"]:
        cols.append("barcode_arrangement")
        logging.info("Nanoget: Extracting metrics per barcode.")
    try:
        datadf = pd.read_csv(
            filepath_or_buffer=summaryfile,
            sep="\t",
            usecols=cols,
        )
    except ValueError:
        logging.error(
            "Nanoget: did not find expected columns in summary file {}:\n {}".format(
                summaryfile, ", ".join(cols)
            )
        )
        sys.exit(
            "ERROR: expected columns in summary file {} not found:\n {}".format(
                summaryfile, ", ".join(cols)
            )
        )
    if kwargs["barcoded"]:
        datadf.columns = ["channelIDs", "time", "duration", "lengths", "quals", "barcode"]
    else:
        datadf.columns = ["channelIDs", "time", "duration", "lengths", "quals"]
    logging.info("Nanoget: Finished collecting statistics from summary file {}".format(summaryfile))
    return ut.reduce_memory_usage(datadf.loc[datadf["lengths"] != 0].copy())


def check_bam(bam, samtype="bam"):
    """Check if bam file is valid.

    Bam file should:
    - exists
    - has an index (create if necessary)
    - is sorted by coordinate
    - has at least one mapped read
    """
    ut.check_existance(bam)
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile.has_index():
        pysam.index(bam)
        samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
        logging.info("Nanoget: No index for bam file could be found, created index.")
    if not samfile.header["HD"]["SO"] == "coordinate":
        logging.error("Nanoget: Bam file {} not sorted by coordinate!.".format(bam))
        sys.exit("Please use a bam file sorted by coordinate.")
    if samtype == "bam":
        logging.info(
            "Nanoget: Bam file {} contains {} mapped and {} unmapped reads.".format(
                bam, samfile.mapped, samfile.unmapped
            )
        )
        if samfile.mapped == 0:
            logging.error("Nanoget: Bam file {} does not contain aligned reads.".format(bam))
            sys.exit("FATAL: not a single read was mapped in bam file {}".format(bam))
    return samfile


def process_ubam(bam, **kwargs):
    """Extracting metrics from unaligned bam format
    Extracting lengths
    """
    logging.info("Nanoget: Starting to collect statistics from ubam file {}.".format(bam))
    samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
    if not samfile.has_index():
        pysam.index(bam)
        # Need to reload the samfile after creating index
        samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
        logging.info("Nanoget: No index for bam file could be found, created index.")
    datadf = (
        pd.DataFrame(
            data=[
                (read.query_name, ut.ave_qual(read.query_qualities), read.query_length)
                for read in samfile.fetch(until_eof=True)
            ],
            columns=["readIDs", "quals", "lengths"],
        )
        .dropna(axis="columns", how="all")
        .dropna(axis="index", how="any")
    )
    logging.info("Nanoget: ubam {} contains {} reads.".format(bam, datadf["lengths"].size))
    return ut.reduce_memory_usage(datadf)


def process_bam(bam, **kwargs):
    """Combines metrics from bam after extraction.

    Processing function: calls pool of worker functions
    to extract from a bam file the following metrics:
    -lengths
    -aligned lengths
    -qualities
    -aligned qualities
    -mapping qualities
    -edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    logging.info("Nanoget: Starting to collect statistics from bam file {}.".format(bam))
    samfile = check_bam(bam)
    chromosomes = samfile.references
    if len(chromosomes) > 200 or kwargs["huge"]:
        logging.info("Nanoget: lots of contigs (>200) or --huge, not running in separate processes")
        datadf = (
            pd.DataFrame(
                data=extract_from_bam(bam, None, kwargs["keep_supp"]),
                columns=[
                    "readIDs",
                    "quals",
                    "aligned_quals",
                    "lengths",
                    "aligned_lengths",
                    "mapQ",
                    "percentIdentity",
                ],
            )
            .dropna(axis="columns", how="all")
            .dropna(axis="index", how="any")
        )

    else:
        unit = chromosomes
        with cfutures.ProcessPoolExecutor(max_workers=kwargs["threads"]) as executor:
            datadf = (
                pd.DataFrame(
                    data=[
                        res
                        for sublist in executor.map(
                            extract_from_bam, repeat(bam), unit, repeat(kwargs["keep_supp"])
                        )
                        for res in sublist
                    ],
                    columns=[
                        "readIDs",
                        "quals",
                        "aligned_quals",
                        "lengths",
                        "aligned_lengths",
                        "mapQ",
                        "percentIdentity",
                    ],
                )
                .dropna(axis="columns", how="all")
                .dropna(axis="index", how="any")
            )
    logging.info(f"Nanoget: bam {bam} contains {datadf['lengths'].size} primary alignments.")
    return ut.reduce_memory_usage(datadf)


def process_cram(cram, **kwargs):
    """Combines metrics from cram after extraction.

    Processing function: calls pool of worker functions
    to extract from a cram file the following metrics:
    -lengths
    -aligned lengths
    -qualities
    -aligned qualities
    -mapping qualities
    -edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    logging.info("Nanoget: Starting to collect statistics from cram file {}.".format(cram))
    samfile = check_bam(cram, samtype="cram")
    chromosomes = samfile.references
    if len(chromosomes) > 100:
        unit = [None]
        logging.info("Nanoget: lots of contigs (>100), not running in separate processes")
    else:
        unit = chromosomes
    with cfutures.ProcessPoolExecutor(max_workers=kwargs["threads"]) as executor:
        datadf = (
            pd.DataFrame(
                data=[
                    res
                    for sublist in executor.map(
                        extract_from_bam, repeat(cram), unit, repeat(kwargs["keep_supp"])
                    )
                    for res in sublist
                ],
                columns=[
                    "readIDs",
                    "quals",
                    "aligned_quals",
                    "lengths",
                    "aligned_lengths",
                    "mapQ",
                    "percentIdentity",
                ],
            )
            .dropna(axis="columns", how="all")
            .dropna(axis="index", how="any")
        )
    logging.info(f"Nanoget: cram {cram} contains {datadf['lengths'].size} primary alignments.")
    return ut.reduce_memory_usage(datadf)


def extract_from_bam(bam, chromosome, keep_supplementary=True):
    """Extracts metrics from bam.

    Worker function per chromosome
    loop over a bam file and create list with tuples containing metrics:
    -qualities
    -aligned qualities
    -lengths
    -aligned lengths
    -mapping qualities
    -edit distances to the reference genome scaled by read length
    """
    samfile = pysam.AlignmentFile(bam, "rb")
    if keep_supplementary:
        return [
            (
                read.query_name,
                ut.ave_qual(read.query_qualities),
                ut.ave_qual(read.query_alignment_qualities),
                read.query_length,
                read.query_alignment_length,
                read.mapping_quality,
                get_pID(read),
            )
            for read in samfile.fetch(reference=chromosome, multiple_iterators=True)
            if not read.is_secondary and not read.is_unmapped
        ]
    else:
        return [
            (
                read.query_name,
                ut.ave_qual(read.query_qualities),
                ut.ave_qual(read.query_alignment_qualities),
                read.query_length,
                read.query_alignment_length,
                read.mapping_quality,
                get_pID(read),
            )
            for read in samfile.fetch(reference=chromosome, multiple_iterators=True)
            if not read.is_secondary and not read.is_unmapped and not read.is_supplementary
        ]


def get_pID(read):
    """Return the percent identity of a read.

    based on the NM tag if present,
    if not calculate from MD tag and CIGAR string

    read.query_alignment_length can be zero in the case of ultra long reads aligned with minimap2 -L
    """
    match = reduce(lambda x, y: x + y[1] if y[0] in (0, 7, 8) else x, read.cigartuples, 0)
    ins = reduce(lambda x, y: x + y[1] if y[0] == 1 else x, read.cigartuples, 0)
    delt = reduce(lambda x, y: x + y[1] if y[0] == 2 else x, read.cigartuples, 0)
    alignment_length = match + ins + delt
    try:
        return (1 - read.get_tag("NM") / alignment_length) * 100
    except KeyError:
        try:
            return 100 * (
                1
                - (parse_MD(read.get_tag("MD")) + parse_CIGAR(read.cigartuples)) / alignment_length
            )
        except KeyError:
            return None
    except ZeroDivisionError:
        return None


def parse_MD(MDlist):
    """Parse MD string to get number of mismatches and deletions."""
    return sum([len(item) for item in re.split("[0-9^]", MDlist)])


def parse_CIGAR(cigartuples):
    """Count the insertions in the read using the CIGAR string."""
    return sum([item[1] for item in cigartuples if item[0] == 1])


def handle_compressed_input(inputfq, file_type="fastq"):
    """Return handles from compressed files according to extension.

    Check for which fastq input is presented and open a handle accordingly
    Can read from compressed files (gz, bz2, bgz) or uncompressed
    Relies on file extensions to recognize compression
    """
    ut.check_existance(inputfq)
    if inputfq.endswith((".gz", "bgz")):
        import gzip

        logging.info("Nanoget: Decompressing gzipped {} {}".format(file_type, inputfq))
        return gzip.open(inputfq, "rt")
    elif inputfq.endswith(".bz2"):
        import bz2

        logging.info("Nanoget: Decompressing bz2 compressed {} {}".format(file_type, inputfq))
        return bz2.open(inputfq, "rt")
    elif inputfq.endswith((".fastq", ".fq", "fasta", ".fa", ".fas")):
        return open(inputfq, "r")
    else:
        logging.error("INPUT ERROR: Unrecognized file extension {}".format(inputfq))
        sys.exit(
            "INPUT ERROR:\nUnrecognized file extension in {}\n"
            "Supported are gz, bz2, bgz, fastq, fq, fasta, fa and fas".format(inputfq)
        )


def process_fasta(fasta, **kwargs):
    """Combine metrics extracted from a fasta file."""
    logging.info("Nanoget: Starting to collect statistics from a fasta file.")
    inputfasta = handle_compressed_input(fasta, file_type="fasta")
    return ut.reduce_memory_usage(
        pd.DataFrame(
            data=[len(rec) for rec in SeqIO.parse(inputfasta, "fasta")], columns=["lengths"]
        ).dropna()
    )


def process_fastq_plain(fastq, **kwargs):
    """Combine metrics extracted from a fastq file."""
    logging.info("Nanoget: Starting to collect statistics from plain fastq file.")
    inputfastq = handle_compressed_input(fastq)
    return ut.reduce_memory_usage(
        pd.DataFrame(
            data=[res for res in extract_from_fastq(inputfastq) if res],
            columns=["quals", "lengths"],
        ).dropna()
    )


def extract_from_fastq(fq):
    """Extract metrics from a fastq file.

    Return average quality and read length
    """
    for rec in SeqIO.parse(fq, "fastq"):
        yield ut.ave_qual(rec.letter_annotations["phred_quality"]), len(rec)


def stream_fastq_full(fastq, threads):
    """Generator for returning metrics extracted from fastq.

    Extract from a fastq file:
    -readname
    -average and median quality
    -read_lenght
    """
    logging.info("Nanoget: Starting to collect full metrics from plain fastq file.")
    inputfastq = handle_compressed_input(fastq)
    with cfutures.ProcessPoolExecutor(max_workers=threads) as executor:
        for results in executor.map(extract_all_from_fastq, SeqIO.parse(inputfastq, "fastq")):
            yield results
    logging.info("Nanoget: Finished collecting statistics from plain fastq file.")


def extract_all_from_fastq(rec):
    """Extract metrics from a fastq file.

    Return identifier, read length, average quality and median quality
    """
    return (rec.id, len(rec), ut.ave_qual(rec.letter_annotations["phred_quality"]), None)


def info_to_dict(info):
    """Get the key-value pairs from the albacore/minknow fastq description and return dict"""
    return {field.split("=")[0]: field.split("=")[1] for field in info.split(" ")[1:]}


def process_fastq_rich(fastq, **kwargs):
    """Extract metrics from a richer fastq file.

    Extract information from fastq files generated by albacore or MinKNOW,
    containing richer information in the header (key-value pairs)
    read=<int> [72]
    ch=<int> [159]
    start_time=<timestamp> [2016-07-15T14:23:22Z]  # UTC ISO 8601 ISO 3339 timestamp
    Z indicates UTC time, T is the delimiter between date expression and time expression
    dateutil.parser.parse("2016-07-15T14:23:22Z") imported as dparse
    -> datetime.datetime(2016, 7, 15, 14, 23, 22, tzinfo=tzutc())
    """
    logging.info("Nanoget: Starting to collect statistics from rich fastq file.")
    inputfastq = handle_compressed_input(fastq)
    res = []
    for record in SeqIO.parse(inputfastq, "fastq"):
        try:
            read_info = info_to_dict(record.description)
            res.append(
                (
                    ut.ave_qual(record.letter_annotations["phred_quality"]),
                    len(record),
                    read_info["ch"],
                    read_info["start_time"],
                    read_info["runid"],
                )
            )
        except KeyError:
            logging.error("Nanoget: keyerror when processing record {}".format(record.description))
            sys.exit(
                "Unexpected fastq identifier:\n{}\n\n \
            missing one or more of expected fields 'ch', 'start_time' or 'runid'".format(
                    record.description
                )
            )
    df = pd.DataFrame(
        data=res, columns=["quals", "lengths", "channelIDs", "timestamp", "runIDs"]
    ).dropna()
    df["timestamp"] = df["timestamp"].astype("datetime64[ns]")
    df["channelIDs"] = df["channelIDs"].astype("int64")
    return ut.reduce_memory_usage(df)


def readfq(fp):
    """Generator function adapted from https://github.com/lh3/readfq."""
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in ">@":  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in "@+>":
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != "+":  # this is a fasta record
            yield name, "".join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = "".join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, "".join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def fq_minimal(fq):
    """Minimal fastq metrics extractor.

    Quickly parse a fasta/fastq file - but makes expectations on the file format
    There will be dragons if unexpected format is used
    Expects a fastq_rich format, but extracts only timestamp and length
    """
    try:
        while True:
            time = next(fq)[1:].split(" ")[4][11:-1]
            length = len(next(fq))
            next(fq)
            next(fq)
            yield time, length
    except StopIteration:
        yield None


def process_fastq_minimal(fastq, **kwargs):
    """Swiftly extract minimal features (length and timestamp) from a rich fastq file"""
    infastq = handle_compressed_input(fastq)
    try:
        df = pd.DataFrame(
            data=[rec for rec in fq_minimal(infastq) if rec], columns=["timestamp", "lengths"]
        )
    except IndexError:
        logging.error("Fatal: Incorrect file structure for fastq_minimal")
        sys.exit("Error: file does not match expected structure for fastq_minimal")
    return ut.reduce_memory_usage(df)
