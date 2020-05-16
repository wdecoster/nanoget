import nanoget


def run_tests():
    """Test functions using testdata from the nanotest repo."""
    nanoget.get_input("bam", ["nanotest/alignment.bam"])
    nanoget.get_input("bam", ["nanotest/alignment.bam"], keep_supp=False)
    nanoget.get_input("fastq_rich", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("summary", ["nanotest/sequencing_summary.txt"], combine="track")
    nanoget.get_input("fastq_minimal", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("fastq", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("fasta", ["nanotest/reads.fa.gz"])


if __name__ == '__main__':
    run_tests()
