"""
Tests for sequencing-read loading, quality trimming and alignment (no GUI).

The .ab1 files are simulated with tools/abif_writer.py, so the tests don't
need real sequencing data.
"""
import random

import pytest

from core import seqtools as st
from tools.abif_writer import simulate_read, write_ab1
from conftest import ROOT

EXAMPLES = ROOT / "examples"
MUTATION = 1000  # the demo reads come from a clone with a point mutation here


@pytest.fixture
def ref(demo_record):
    return st.record_sequence(demo_record)


# ------------------------------------------------------------------ trimming
def test_quality_trim():
    assert st.quality_trim([5] * 10 + [40] * 50 + [5] * 10) == (10, 60)
    assert st.quality_trim([5] * 50) == (0, 0)          # nothing worth keeping
    assert st.quality_trim([40] * 30) == (0, 30)        # keep all
    assert st.quality_trim(None) == (0, 0)


# ------------------------------------------------------------------ loading
def test_load_ab1(tmp_path, ref):
    seq, quals, traces, peaks = simulate_read(ref[100:400], "test", seed=7)
    path = tmp_path / "test.ab1"
    write_ab1(path, "test", seq, quals, traces, peaks)
    read = st.load_read(path)
    assert read["name"] == "test" and read["seq"] == seq and read["quals"] == quals
    assert set(read["trace"]["channels"]) == set("ACGT")
    assert read["trace"]["peaks"] == peaks


def test_load_fastq_fasta_and_plain_text(tmp_path):
    (tmp_path / "r.fastq").write_text("@r\nACGTACGT\n+\nIIIIIIII\n")
    (tmp_path / "r.fasta").write_text(">r\nACGT\nACGT\n")
    (tmp_path / "r.seq").write_text("acgt acgt\nAC\n")
    fastq = st.load_read(tmp_path / "r.fastq")
    assert fastq["seq"] == "ACGTACGT" and fastq["quals"] == [40] * 8 and fastq["trace"] is None
    assert st.load_read(tmp_path / "r.fasta")["quals"] is None
    assert st.load_read(tmp_path / "r.seq")["seq"] == "ACGTACGTAC"


# ------------------------------------------------------------------ demo reads
# Where the mutation sits in each demo read, counted in the read as sequenced
# (see tools/make_demo_reads.py: GFP_fwd reads clone[640:1460] forwards,
# GFP_rev reads clone[620:1480] backwards).
@pytest.mark.parametrize("name, strand, index", [("GFP_fwd", 1, MUTATION - 1 - 640),
                                                 ("GFP_rev", -1, 1480 - MUTATION)])
def test_demo_reads_find_the_mutation(ref, name, strand, index):
    read = st.load_read(EXAMPLES / f"pDEMO_{name}.ab1")
    result = st.align_read(ref, read["seq"], read["quals"], circular=True)
    assert result["strand"] == strand
    real = [d for d in result["differences"] if not d["low_quality"]]
    assert [(d["kind"], d["ref_pos"]) for d in real] == [("mismatch", MUTATION)]
    difference = real[0]
    assert difference["quality"] >= 40
    # read_index points at the base in the read as sequenced (used for the chromatogram)
    assert difference["read_index"] == index
    assert difference["quality"] == read["quals"][index]
    called = read["seq"][index]
    assert called == (difference["read"] if strand == 1 else st.reverse_complement(difference["read"]))


def test_clean_demo_read(ref):
    read = st.load_read(EXAMPLES / "pDEMO_ori_fwd.ab1")
    result = st.align_read(ref, read["seq"], read["quals"], circular=True)
    assert result["identity"] == 100 and result["differences"] == []


# ------------------------------------------------------------------ alignment
def test_insertion_and_deletion(ref):
    clone = ref[:1100] + ref[1103:1200] + "GG" + ref[1200:]  # delete 1101-1103, insert GG after 1200
    result = st.align_read(ref, clone[900:1500], None, circular=True)
    found = [(d["kind"], d["ref_pos"], d["ref"], d["read"]) for d in result["differences"]]
    assert found == [("deletion", 1101, ref[1100:1103], "-"), ("insertion", 1200, "-", "GG")]


def test_reverse_read_across_origin(ref):
    read = st.reverse_complement(ref[-300:] + ref[:300])
    result = st.align_read(ref, read, None, circular=True)
    assert (result["strand"], result["start"], result["end"]) == (-1, 3701, 300)
    assert result["identity"] == 100


def test_unrelated_read_gives_none(ref):
    rng = random.Random(5)
    junk = "".join(rng.choice("ACGT") for _ in range(500))
    assert st.align_read(ref, junk, None, circular=True) is None


def test_long_reference_uses_seeding():
    rng = random.Random(9)
    big = "".join(rng.choice("ACGT") for _ in range(400_000))
    read = big[250_000:250_800]
    changed = {"A": "C", "C": "G", "G": "T", "T": "A"}[read[400]]
    read = read[:400] + changed + read[401:]
    result = st.align_read(big, read, None)
    assert (result["start"], result["end"]) == (250_001, 250_800)
    assert [(d["ref_pos"], d["read"]) for d in result["differences"]] == [(250_401, changed)]


def test_align_reads_reports_unreadable_files(tmp_path, ref):
    bad = tmp_path / "broken.ab1"
    bad.write_bytes(b"not a chromatogram")
    good = EXAMPLES / "pDEMO_ori_fwd.ab1"
    results = st.align_reads(ref, [good, bad], circular=True)
    assert results[0]["alignment"] and results[0]["error"] is None
    assert results[1]["alignment"] is None and "Could not read" in results[1]["error"]


def test_format_alignment(ref):
    result = st.align_read(ref, ref[100:250], None, circular=True)
    text = st.format_alignment(result["rows"], width=60)
    assert text.startswith("Ref      101  ")
    assert "|" * 60 in text
