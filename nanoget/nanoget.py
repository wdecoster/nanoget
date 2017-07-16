# wdecoster

from __future__ import division
import sys
import os
import time
import logging
import re
import pandas as pd
import numpy as np
from Bio import SeqIO
from multiprocessing import Pool
import dateutil.parser
import pysam
import nanomath
from .version import __version__


def checkExistance(f):
	'''
	Check if the file supplied as input exists
	'''
	if not os.path.isfile(f):
		logging.error("Nanoget: File provided doesn't exist or the path is incorrect: {}".format(f))
		sys.exit("File provided doesn't exist or the path is incorrect: {}".format(f))


def processSummary(summaryfile):
	'''
	Extracting information from an albacore summary file with the following header:
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
	'''
	logging.info("Nanoget: Staring to collect statistics from summary file.")
	checkExistance(summaryfile)
	cols = ["read_id", "run_id", "channel", "start_time", "sequence_length_template", "mean_qscore_template"]
	try:
		datadf = pd.read_csv(
			filepath_or_buffer=summaryfile,
			sep="\t",
			usecols=cols,
			)
	except ValueError:
		logging.error("Nanoget: did not find expected columns in summary file.")
		sys.exit("ERROR: did not find expected columns in summary file:\n {}".format(', '.join(cols)))
	datadf.columns = ["readIDs", "runIDs", "channelIDs", "time", "lengths", "quals"]
	a_time_stamps = np.array(datadf["time"], dtype='datetime64[s]')
	datadf["start_time"] = a_time_stamps - np.amin(a_time_stamps)
	logging.info("Nanoget: Finished collecting statistics from summary file.")
	return datadf[datadf["lengths"] != 0]


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
	logging.info("Nanoget: Staring to collect statistics from bam file.")
	checkExistance(bam)
	samfile = pysam.AlignmentFile(bam, "rb")
	if not samfile.has_index():
		pysam.index(bam)
		samfile = pysam.AlignmentFile(bam, "rb")  # Need to reload the samfile after creating index
		logging.info("Nanoget: No index for bam file could be found, created index.")
	if not samfile.header['HD']['SO'] == 'coordinate':
		logging.info("Nanoget: Bam file not sorted by coordinate!.")
		sys.exit("Please use a bam file sorted by coordinate.")
	NumberOfmappedReads = samfile.mapped
	NumberOfunmappedReads = samfile.unmapped
	logging.info("Nanoget: Bam file contains {} mapped and {} unmapped reads.".format(NumberOfmappedReads, NumberOfunmappedReads))
	if NumberOfmappedReads == 0:
		sys.exit("FATAL: not a single read was mapped in the bam file.")
	chromosomes = samfile.references
	datadf = pd.DataFrame()
	pool = Pool(processes=threads)
	params = zip([bam]*len(chromosomes), chromosomes)
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
	datadf["lengths"] = np.array([x for y in [elem[0] for elem in output] for x in y])
	datadf["aligned_lengths"] = np.array([x for y in [elem[1] for elem in output] for x in y])
	datadf["quals"] = np.array([x for y in [elem[2] for elem in output] for x in y])
	datadf["aligned_quals"] = np.array([x for y in [elem[3] for elem in output] for x in y])
	datadf["mapQ"] = np.array([x for y in [elem[4] for elem in output] for x in y])
	datadf["percentIdentity"] = np.array([x for y in [elem[5] for elem in output] for x in y])
	assert datadf["lengths"].size == NumberOfmappedReads, "Unexpected difference in length of entries in datadict"
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
		lengths.append(read.query_length)
		alignedLengths.append(read.query_alignment_length)
		quals.append(nanomath.aveQual(read.query_qualities))
		alignedQuals.append(nanomath.aveQual(read.query_alignment_qualities))
		mapQ.append(read.mapping_quality)
		try:
			pID.append((1- read.get_tag("NM")/read.query_alignment_length)*100)
		except KeyError:
			pID.append((1 -
				( parseMD(read.get_tag("MD")) + parseCIGAR(read.cigartuples))
				/read.query_alignment_length)*100)
	return (lengths, alignedLengths, quals, alignedQuals, mapQ, pID)


def parseMD(MDlist):
	return sum([len(item) for item in re.split('[0-9^]', MDlist )])


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
		get_compression_type(inputfq)
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


def get_compression_type(filename):
	"""
	Attempts to guess the compression (if any) on a file using the first few bytes.
	Based on http://stackoverflow.com/questions/13044562 and https://github.com/rrwick/Porechop/blob/master/porechop/misc.py#L68
	"""
	magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
					'bz2': (b'\x42', b'\x5a', b'\x68'),
					'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
	max_len = max(len(x) for x in magic_dict.values())

	unknown_file = open(filename, 'rb')
	file_start = unknown_file.read(max_len)
	unknown_file.close()
	compression_type = 'plain'
	for filetype, magic_bytes in magic_dict.items():
		if file_start.startswith(magic_bytes):
			compression_type = filetype
	if compression_type == 'bz2':
		sys.exit('Error: cannot use bzip2 format - use gzip instead')
	if compression_type == 'zip':
		sys.exit('Error: cannot use zip format - use gzip instead')
	return compression_type


def processFastqPlain(fastq):
	'''
	Processing function
	Iterate over a fastq file and extract metrics
	'''
	logging.info("Nanoget: Starting to collect statistics from plain fastq file.")
	inputfastq = handlecompressedFastq(fastq)
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		try:
			quals.append(nanomath.aveQual(record.letter_annotations["phred_quality"]))
			lengths.append(len(record))
		except ZeroDivisionError:
			pass
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	logging.info("Nanoget: Finished collecting statistics from plain fastq file.")
	return datadf


def processFastq_rich(fastq):
	'''
	Extract information from fastq files generated by albacore or MinKNOW, containing richer information in the header
	containing key-value pairs
	read=<int> [72]
	ch=<int> [159]
	start_time=<timestamp> [2016-07-15T14:23:22Z]  # UTC ISO 8601 ISO 3339 timestamp
	Z indicates UTC time, T is the delimiter between date expression and time expression
	dateutil.parser.parse("2016-07-15T14:23:22Z") # -> datetime.datetime(2016, 7, 15, 14, 23, 22, tzinfo=tzutc())
	'''
	logging.info("Nanoget: Starting to collect statistics from rich fastq file.")
	inputfastq = handlecompressedFastq(fastq)
	datadf = pd.DataFrame()
	lengths = []
	quals = []
	channels = []
	time_stamps = []
	runids = []
	for record in SeqIO.parse(inputfastq, "fastq"):
		try:
			quals.append(nanomath.aveQual(record.letter_annotations["phred_quality"]))
			lengths.append(len(record))
			r_info = {field.split('=')[0] : field.split('=')[1] for field in record.description.split(' ')[1:]}
			channels.append(r_info["ch"])
			time_stamps.append(dateutil.parser.parse(r_info["start_time"]))
			runids.append(r_info["runid"])
		except ZeroDivisionError:  # For reads with length 0, nanomath.aveQual will throw a ZeroDivisionError
			pass
		except KeyError:
			logging.error("Nanoget: keyerror when processing record {}".format(record.description))
			sys.exit("Unexpected fastq identifier:\n{}\n\nmissing one or more of expected fields 'ch', 'start_time' or 'runid'".format(record.description))
	datadf["lengths"] = np.array(lengths)
	datadf["quals"] = np.array(quals)
	datadf["channelIDs"] = np.int16(channels)
	datadf["runIDs"] = np.array(runids)
	a_time_stamps = np.array(time_stamps, dtype='datetime64[s]')
	datadf["start_time"] = a_time_stamps - np.amin(a_time_stamps)
	logging.info("Nanoget: Finished collecting statistics from rich fastq file.")
	return datadf
