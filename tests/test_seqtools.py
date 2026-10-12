"""
Tests for core/seqtools.py: the sequence logic, no GUI involved.
These are fast and catch most bugs, so most new tests should go here.
"""
import warnings

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from core import seqtools as st


# ------------------------------------------------------------------ type detection
@pytest.mark.parametrize("sequence, expected", [
    ("ACGT", "DNA"),
    ("acgu nnn", "DNA"),          # RNA, lowercase and spaces are fine
    ("MKVLAAGIVG*", "Protein"),
    ("ACGT123", "Unknown"),
    ("", "Unknown"),
])
def test_detect_seq_type(sequence, expected):
    assert st.detect_seq_type(sequence) == expected


def test_gc_percent():
    assert st.gc_percent("GGCCAATT") == 50
    assert st.gc_percent("") == 0


# ------------------------------------------------------------------ translation
def test_translate_ignores_leftover_bases():
    assert st.translate_dna("ATGGCTTAAGG") == "MA*"


def test_translate_codon_start():
    assert st.translate_dna("AATGTAA", codon_start=2) == "M*"


def test_translation_settings_from_qualifiers():
    feature = st.make_feature("CDS", "x", 1, 30, 1, 100, False)
    assert st.feature_translation_settings(feature) == (1, 1)
    feature.qualifiers["transl_table"] = ["11"]
    feature.qualifiers["codon_start"] = ["2"]
    assert st.feature_translation_settings(feature) == (11, 2)


# ------------------------------------------------------------------ locations
def test_simple_location_is_zero_based_end_exclusive():
    location = st.make_location(101, 200, 1, 4000, circular=True)
    assert (int(location.start), int(location.end)) == (100, 200)


@pytest.mark.parametrize("strand", [1, -1])
def test_feature_across_origin(strand):
    feature = st.make_feature("CDS", "x", 3901, 60, strand, 4000, circular=True)
    assert st.wraps_origin(feature, 4000)
    assert st.feature_span(feature) == (3901, 60)
    assert st.feature_strand(feature) == strand


def test_feature_across_origin_refused_on_linear():
    with pytest.raises(ValueError, match="circular"):
        st.make_location(300, 100, 1, 4000, circular=False)


def test_location_out_of_range():
    with pytest.raises(ValueError):
        st.make_location(0, 100, 1, 4000, circular=True)


def test_reverse_feature_sequence_across_origin():
    sequence = "AAAACCCCGGGGTTTT" * 250  # 4,000 bp
    feature = st.make_feature("CDS", "x", 3901, 60, -1, 4000, circular=True)
    expected = st.reverse_complement(sequence[3900:] + sequence[:60])
    assert st.feature_sequence(sequence, feature) == expected


# ------------------------------------------------------------------ demo plasmid
def test_demo_genes_translate_without_internal_stops(demo_record):
    sequence = st.record_sequence(demo_record)
    for feature in st.display_features(demo_record):
        if feature.type == "CDS":
            protein = st.translate_dna(st.feature_sequence(sequence, feature))
            assert protein.startswith("M") and protein.endswith("*")
            assert "*" not in protein[:-1], st.feature_label(feature)


def test_demo_mcs_enzymes_cut_once(demo_record):
    sites = st.find_restriction_sites(st.record_sequence(demo_record), circular=True)
    for enzyme in ["EcoRI", "SacI", "KpnI", "SmaI", "BamHI",
                   "XbaI", "SalI", "PstI", "SphI", "HindIII"]:
        assert len(sites[enzyme]) == 1, enzyme
    assert sites["EcoRI"] == [602]


# ------------------------------------------------------------------ saving
def test_genbank_save_shortens_long_names(tmp_path):
    record = SeqRecord(Seq("ATGC" * 50), id="a_really_long_name_from_fasta",
                       name="a_really_long_name_from_fasta")
    record.features.append(st.make_feature("CDS", "gene X", 11, 40, 1, 200, False))
    path = tmp_path / "out.gb"

    with warnings.catch_warnings():
        warnings.simplefilter("error")  # fail on any Biopython warning, too
        renamed = st.save_records([record], path, "genbank")

    assert renamed == [("a_really_long_name_from_fasta", "a_really_long_na")]
    back = st.load_records(path)[0][0]
    assert back.name == "a_really_long_na"
    assert len(back.features) == 1
    assert record.name == "a_really_long_name_from_fasta"  # original untouched


def test_genbank_round_trip_keeps_topology_and_features(demo_record, tmp_path):
    path = tmp_path / "copy.gb"
    st.save_records([demo_record], path, "genbank")
    back = st.load_records(path)[0][0]
    assert st.is_circular(back)
    assert ([st.feature_label(f) for f in st.display_features(back)]
            == [st.feature_label(f) for f in st.display_features(demo_record)])
