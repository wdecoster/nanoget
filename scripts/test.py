import os
import tempfile

import nanoget


def run_tests():
    """Test functions using testdata from the nanotest repo."""
    nanoget.get_input("bam", ["nanotest/alignment.bam"])
    nanoget.get_input("bam", ["nanotest/alignment.bam"], keep_supp=False)
    nanoget.get_input("fastq_rich", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("fastq_rich", ["nanotest/reads-mixed-timestamp.fastq"])
    nanoget.get_input("summary", ["nanotest/sequencing_summary.txt"], combine="track")
    nanoget.get_input("fastq_minimal", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("fastq", ["nanotest/reads.fastq.gz"])
    nanoget.get_input("fasta", ["nanotest/reads.fa.gz"])


HEADER = ["channel", "start_time", "duration", "sequence_length_template", "mean_qscore_template"]


def write_summary(path, extra_columns, extra_values):
    """Write a minimal summary file with 10 reads, cycling through extra_values."""
    with open(path, "w") as fh:
        fh.write("\t".join(HEADER + extra_columns) + "\n")
        for i in range(10):
            row = [str(i + 1), str(i * 2.5), "1.5", str(500 + i), "12.5"]
            fh.write("\t".join(row + list(extra_values[i % len(extra_values)])) + "\n")


def barcodes(path):
    df = nanoget.get_input("summary", [path], barcoded=True)
    assert "alias" not in df.columns, "alias column leaked into the output"
    return sorted(df["barcode"].unique())


def run_barcode_tests():
    """Test which column barcodes are taken from, see NanoPlot issue 440."""
    with tempfile.TemporaryDirectory() as tmpdir:
        def path(name):
            return os.path.join(tmpdir, name)

        # barcode_arrangement is authoritative when it carries barcodes, even though
        # a sample sheet puts different values in alias
        write_summary(
            path("sample_sheet.txt"),
            ["alias", "barcode_arrangement"],
            [("sample_A", "barcode01"), ("sample_B", "barcode02")],
        )
        assert barcodes(path("sample_sheet.txt")) == ["barcode01", "barcode02"]

        # dorado 2.1.0 and 2.1.1 left barcode_arrangement "unclassified" for every read
        # while writing the real barcode to alias, fall back to alias in that case
        write_summary(
            path("dorado211.txt"),
            ["alias", "barcode_arrangement"],
            [("barcode01", "unclassified"), ("barcode02", "unclassified")],
        )
        assert barcodes(path("dorado211.txt")) == ["barcode01", "barcode02"]

        # a genuinely unclassified run has "unclassified" in alias too, so it stays that way
        write_summary(
            path("unclassified.txt"),
            ["alias", "barcode_arrangement"],
            [("unclassified", "unclassified")],
        )
        assert barcodes(path("unclassified.txt")) == ["unclassified"]

        # older summary files have no alias column at all
        write_summary(
            path("legacy.txt"),
            ["barcode_arrangement"],
            [("barcode01",), ("barcode02",)],
        )
        assert barcodes(path("legacy.txt")) == ["barcode01", "barcode02"]

        # columns are renamed by name, not by position, so any column order works
        write_summary(
            path("reordered.txt"),
            ["barcode_arrangement"],
            [("barcode01",), ("barcode02",)],
        )
        with open(path("reordered.txt")) as fh:
            rows = [line.rstrip("\n").split("\t") for line in fh]
        order = sorted(range(len(rows[0])), key=lambda i: rows[0][i])
        with open(path("reordered.txt"), "w") as fh:
            for row in rows:
                fh.write("\t".join(row[i] for i in order) + "\n")
        df = nanoget.get_input("summary", [path("reordered.txt")], barcoded=True)
        assert sorted(df["barcode"].unique()) == ["barcode01", "barcode02"]
        assert sorted(df["channelIDs"].unique()) == list(range(1, 11))
        assert sorted(df["lengths"].unique()) == list(range(500, 510))


if __name__ == "__main__":
    run_tests()
    run_barcode_tests()
