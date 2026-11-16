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


def run_sam_style_header_tests():
    """SAM-style headers (MinKNOW/Dorado >= 26.01), including tab-separated, see issue 38."""
    from nanoget.extraction_functions import info_to_dict

    legacy = "id runid=abc read=1 ch=5 start_time=2016-07-15T14:23:22Z"
    assert info_to_dict(legacy) == {"runid": "abc", "read": "1", "ch": "5",
                                    "start_time": "2016-07-15T14:23:22Z"}
    spaced = "id qs:f:28.1 ch:i:353 st:Z:2026-11-09T17:10:18+01:00 RG:Z:run1_model@v5_bc01"
    expected = {"ch": "353", "start_time": "2026-11-09T17:10:18+01:00", "runid": "run1"}
    assert info_to_dict(spaced) == expected
    assert info_to_dict(spaced.replace(" ", "\t")) == expected
    # no RG tag (raw dorado output): runid falls back to "."
    no_rg = "id\tqs:f:28.1\tmx:i:2\tch:i:353\trn:i:1\tst:Z:2026-11-09T17:10:18.5+01:00"
    assert info_to_dict(no_rg) == {"ch": "353", "start_time": "2026-11-09T17:10:18.5+01:00",
                                   "runid": "."}
    with tempfile.TemporaryDirectory() as tmpdir:
        fq = os.path.join(tmpdir, "reads.fastq")
        with open(fq, "w") as fh:
            fh.write("@r1\tqs:f:20\tch:i:7\tst:Z:2026-11-09T17:10:18+01:00\nACGT\n+\nIIII\n")
        df = nanoget.get_input("fastq_rich", [fq])
        assert list(df["channelIDs"]) == [7] and list(df["runIDs"]) == ["."]


if __name__ == "__main__":
    run_tests()
    run_barcode_tests()
    run_sam_style_header_tests()
