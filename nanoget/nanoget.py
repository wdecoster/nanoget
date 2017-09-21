# wdecoster

from __future__ import division
import sys
import os
import logging
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
import dateutil.parser
import pysam
import nanomath


def checkExistance(f):
    '''
    Check if the file supplied as input exists
    '''
    if not os.path.isfile(f):
        logging.error("Nanoget: File provided doesn't exist or the path is incorrect: {}".format(f))
        sys.exit("File provided doesn't exist or the path is incorrect: {}".format(f))


def processSummary(summaryfile, readtype):
    '''
    Extracting information from an albacore summary file.
    Only reads which have a >0 length are returned.

    The fields below may or may not exist, depending on the type of sequencing performed.
    Fields 1-14 are for 1D sequencing.
    Fields 1-23 for 2D sequencing.
    Fields 24-27, 2-5, 22-23 for 1D^2 (1D2) sequencing
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
    11  num_called_templatexticks
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
    '''
    logging.info("Nanoget: Staring to collect statistics from summary file.")
    checkExistance(summaryfile)
    logging.info("Collecting statistics for {} sequencing".format(readtype))
    if readtype == "1D":
        cols = ["read_id", "run_id", "channel", "start_time",
                "sequence_length_template", "mean_qscore_template"]
    elif readtype in ["2D", "1D2"]:
        cols = ["read_id", "run_id", "channel", "start_time",
                "sequence_length_2d", "mean_qscore_2d"]
    try:
        datadf = pd.read_csv(
            filepath_or_buffer=summaryfile,
            sep="\t",
            usecols=cols,
        )
    except ValueError:
        logging.error(
            "Nanoget: did not find expected columns in summary file:\n {}.".format(', '.join(cols)))
        sys.exit("ERROR: did not find expected columns in summary file:\n {}".format(', '.join(cols)))
    datadf.columns = ["readIDs", "runIDs", "channelIDs", "time", "lengths", "quals"]
    a_time_stamps = np.array(datadf["time"], dtype='datetime64[s]')
    datadf["start_time"] = a_time_stamps - np.amin(a_time_stamps)
    logging.info("Nanoget: Finished collecting statistics from summary file.")
    return datadf[datadf["lengths"] != 0]


def checkBam(bam):
    '''
    Check if bam file
    - exists
    - has an index (create if necessary)
    - is sorted by coordinate
    - has at least one mapped read
    '''
    checkExistance(bam)
    samfile = pysam.AlignmentFile(bam, "rb")
    if not samfile.has_index():
        pysam.index(bam)
        samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
        logging.info("Nanoget: No index for bam file could be found, created index.")
    if not samfile.header['HD']['SO'] == 'coordinate':
        logging.error("Nanoget: Bam file not sorted by coordinate!.")
        sys.exit("Please use a bam file sorted by coordinate.")
    logging.info("Nanoget: Bam file contains {} mapped and {} unmapped reads.".format(
        samfile.mapped, samfile.unmapped))
    if samfile.mapped == 0:
        logging.error("Nanoget: Bam file does not contain aligned reads.")
        sys.exit("FATAL: not a single read was mapped in the bam file.")
    return samfile


def processBam(bam, threads):
    '''
    Processing function: calls pool of worker functions
    to extract from a bam file the following metrics:
    -lengths
    -aligned lengths
    -qualities
    -aligned qualities
    -mapping qualities
    -edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    '''
    logging.info("Nanoget: Starting to collect statistics from bam file.")
    samfile = checkBam(bam)
    chromosomes = samfile.references
    pool = Pool(processes=threads)
    params = zip([bam] * len(chromosomes), chromosomes)
    try:
        output = [results for results in pool.imap(extractFromBam, params)]
    except KeyboardInterrupt:
        sys.stderr.write("Terminating worker threads")
        pool.terminate()
        pool.join()
        sys.exit()
    # Output contains a tuple per worker
    # Each tuple contains lists per metric
    # Unpacked by following nested list comprehensions
    datadf = pd.DataFrame(data={
        "lengths": np.array([x for y in [elem[0] for elem in output] for x in y]),
        "aligned_lengths": np.array([x for y in [elem[1] for elem in output] for x in y]),
        "quals": np.array([x for y in [elem[2] for elem in output] for x in y]),
        "aligned_quals": np.array([x for y in [elem[3] for elem in output] for x in y]),
        "mapQ": np.array([x for y in [elem[4] for elem in output] for x in y]),
        "percentIdentity": np.array([x for y in [elem[5] for elem in output] for x in y]),
    })
    logging.info("Nanoget: bam contains {} primary alignments.".format(datadf["lengths"].size))
    logging.info("Nanoget: Finished collecting statistics from bam file.")
    return datadf


def extractFromBam(params):
    '''
    Worker function per chromosome
    loop over a bam file and create tuple with lists containing metrics:
    -lengths
    -aligned lengths
    -qualities
    -aligned qualities
    -mapping qualities
    -edit distances to the reference genome scaled by read length
    '''
    bam, chromosome = params
    samfile = pysam.AlignmentFile(bam, "rb")
    lengths = []
    alignedLengths = []
    quals = []
    alignedQuals = []
    mapQ = []
    pID = []
    for read in samfile.fetch(reference=chromosome, multiple_iterators=True):
        if not read.is_secondary:
            quals.append(nanomath.aveQual(read.query_qualities))
            alignedQuals.append(nanomath.aveQual(read.query_alignment_qualities))
            lengths.append(read.query_length)
            alignedLengths.append(read.query_alignment_length)
            mapQ.append(read.mapping_quality)
            try:
                pID.append((1 - read.get_tag("NM") / read.query_alignment_length) * 100)
            except KeyError:
                pID.append((1 -
                            (parseMD(read.get_tag("MD")) + parseCIGAR(read.cigartuples))
                            / read.query_alignment_length) * 100)
    return (lengths, alignedLengths, quals, alignedQuals, mapQ, pID)


def parseMD(MDlist):
    return sum([len(item) for item in re.split('[0-9^]', MDlist)])


def parseCIGAR(cigartuples):
    return sum([item[1] for item in cigartuples if item[0] == 1])


def handlecompressedFastq(inputfq):
    '''
    Check for which fastq input is presented and open a handle accordingly
    Can read from stdin, compressed files (gz, bz2, bgz) or uncompressed
    Relies on file extensions to recognize compression
    '''
    if inputfq == 'stdin':
        logging.info("Nanoget: Reading from stdin.")
        return sys.stdin
    else:
        checkExistance(inputfq)
        if inputfq.endswith('.gz'):
            import gzip
            logging.info("Nanoget: Decompressing gzipped fastq.")
            return gzip.open(inputfq, 'rt')
        elif inputfq.endswith('.bz2'):
            import bz2
            logging.info("Nanoget: Decompressing bz2 compressed fastq.")
            return bz2.BZ2File(inputfq, 'rt')
        elif inputfq.endswith(('.fastq', '.fq', '.bgz')):
            return open(inputfq, 'r')
        else:
            logging.error("INPUT ERROR: Unrecognized file extension")
            sys.exit('''INPUT ERROR: Unrecognized file extension\n,
                        supported formats for --fastq are .gz, .bz2, .bgz, .fastq and .fq''')


def processFastqPlain(fastq, threads):
    '''
    Processing function
    Iterate over a fastq file and extract metrics
    Parallelized, although there is not much to gain with using many threads
    Saturation already starts at threads=3
    '''
    logging.info("Nanoget: Starting to collect statistics from plain fastq file.")
    inputfastq = handlecompressedFastq(fastq)
    pool = Pool(processes=threads)
    output = [results for results in pool.imap(extractFromFastq, SeqIO.parse(inputfastq, "fastq"))]
    logging.info("Nanoget: Finished collecting statistics from plain fastq file.")
    return pd.DataFrame(data={"lengths": np.array([item[0] for item in output]), "quals": np.array([item[0] for item in output])})


def extractFromFastq(rec):
    try:
        return (nanomath.aveQual(rec.letter_annotations["phred_quality"]), len(rec))
    except ZeroDivisionError:  # If length 0, nanomath.aveQual will throw a ZeroDivisionError
        pass


def info2dict(info):
    '''
    Get the key-value pairs from the albacore/minknow fastq description and return as dictionary
    '''
    return {field.split('=')[0]: field.split('=')[1] for field in info.split(' ')[1:]}


def processFastq_rich(fastq):
    '''
    Extract information from fastq files generated by albacore or MinKNOW,
    containing richer information in the header (key-value pairs)
    read=<int> [72]
    ch=<int> [159]
    start_time=<timestamp> [2016-07-15T14:23:22Z]  # UTC ISO 8601 ISO 3339 timestamp
    Z indicates UTC time, T is the delimiter between date expression and time expression
    dateutil.parser.parse("2016-07-15T14:23:22Z")
    -> datetime.datetime(2016, 7, 15, 14, 23, 22, tzinfo=tzutc())
    '''
    logging.info("Nanoget: Starting to collect statistics from rich fastq file.")
    inputfastq = handlecompressedFastq(fastq)
    lengths = []
    quals = []
    channels = []
    time_stamps = []
    runids = []
    for record in SeqIO.parse(inputfastq, "fastq"):
        try:
            quals.append(nanomath.aveQual(record.letter_annotations["phred_quality"]))
            lengths.append(len(record))
            read_info = info2dict(record.description)
            channels.append(read_info["ch"])
            time_stamps.append(dateutil.parser.parse(read_info["start_time"]))
            runids.append(read_info["runid"])
        except ZeroDivisionError:  # If length 0, nanomath.aveQual will throw a ZeroDivisionError
            pass
        except KeyError:
            logging.error("Nanoget: keyerror when processing record {}".format(record.description))
            sys.exit("Unexpected fastq identifier:\n{}\n\n \
            missing one or more of expected fields 'ch', 'start_time' or 'runid'".format(
                record.description))
    a_time_stamps = np.array(time_stamps, dtype='datetime64[s]')
    datadf = pd.DataFrame(data={
        "lengths": np.array(lengths),
        "quals": np.array(quals),
        "channelIDs": np.int16(channels),
        "runIDs": np.array(runids),
        "start_time": a_time_stamps - np.amin(a_time_stamps)
    })
    logging.info("Nanoget: Finished collecting statistics from rich fastq file.")
    return datadf
