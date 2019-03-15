def process_ubam(bam, **kwargs):
    """Extracting metrics from unaligned bam format
    Extracting lengths
    """
    logging.info("Nanoget: Starting to collect statistics from ubam file {}.".format(bam))
    samfile = pysam.AlignmentFile(bam, "rb", check_sq=False)
    if not samfile.has_index():
        pysam.index(bam)
        # Need to reload the samfile after creating index
        samfile = pysam.AlignmentFile(bam, "rb")
        logging.info("Nanoget: No index for bam file could be found, created index.")
    datadf = pd.DataFrame(
        data=[(read.query_name, nanomath.ave_qual(read.query_qualities), read.query_length)
              for read in samfile.fetch(until_eof=True)],
        columns=["readIDs", "quals", "lengths"]) \
        .dropna(axis='columns', how='all') \
        .dropna(axis='index', how='any')
    logging.info("Nanoget: ubam {} contains {} reads.".format(
        bam, datadf["lengths"].size))
    return ut.reduce_memory_usage(datadf)

