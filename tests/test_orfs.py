"""
Tests for the ORF finder (seqtools.find_orfs).

Besides hand-made cases, this compares the finder against a deliberately
simple (slow, obvious) implementation on random DNA: if both agree, the fast
one is very likely right.
"""
import random

import pytest

from core import seqtools as st

FILLER = "C" * 30
STOPS = {"TAA", "TAG", "TGA"}


def orf_dna(n_aa):
    """ATG + (n_aa - 1) alanine codons + stop: an ORF of exactly n_aa amino acids."""
    return "ATG" + "GCT" * (n_aa - 1) + "TAA"


def protein_of(sequence, orf):
    return st.translate_dna(st.feature_sequence(sequence, orf))


def test_forward_orf():
    sequence = FILLER + orf_dna(120) + FILLER
    [orf] = st.find_orfs(sequence, min_aa=100)
    assert st.feature_span(orf) == (31, 393)
    assert protein_of(sequence, orf) == "M" + "A" * 119 + "*"


def test_reverse_orf():
    sequence = FILLER + st.reverse_complement(orf_dna(150)) + FILLER
    [orf] = st.find_orfs(sequence, min_aa=100)
    assert st.feature_strand(orf) == -1
    assert protein_of(sequence, orf) == "M" + "A" * 149 + "*"


def test_minimum_length():
    sequence = FILLER + orf_dna(120) + FILLER
    assert st.find_orfs(sequence, min_aa=120)
    assert not st.find_orfs(sequence, min_aa=121)


def test_nested_atgs_give_one_orf():
    sequence = FILLER + "ATG" + "GCT" * 50 + "ATG" + "GCT" * 120 + "TAA" + FILLER
    [orf] = st.find_orfs(sequence, min_aa=100)
    assert orf.qualifiers["aa"][0] == 172  # from the first ATG


@pytest.mark.parametrize("reverse", [False, True])
def test_orf_across_origin(reverse):
    dna = orf_dna(130)
    sequence = dna[200:] + "C" * 60 + dna[:200]  # ORF wraps around position 0
    if reverse:
        sequence = st.reverse_complement(sequence)
    assert st.find_orfs(sequence, min_aa=100, circular=False) == []
    [orf] = st.find_orfs(sequence, min_aa=100, circular=True)
    assert st.wraps_origin(orf, len(sequence))
    assert protein_of(sequence, orf) == "M" + "A" * 129 + "*"


def simple_orf_finder(seq, min_aa, circular):
    """Slow but obvious: for every stop codon, walk back to the furthest ATG
    before the previous stop. Returns {(strand, stop end, aa length)}."""
    length = len(seq)
    found = set()
    for strand, s in ((1, seq), (-1, st.reverse_complement(seq))):
        def codon(j):
            if circular:
                return "".join(s[(j + k) % length] for k in range(3))
            return s[j:j + 3]
        for stop in range(length if circular else length - 2):
            if codon(stop) not in STOPS:
                continue
            start, j = None, stop - 3
            for _ in range(length // 3 + 1):
                if not circular and j < 0:
                    break
                c = codon(j % length if circular else j)
                if c in STOPS:
                    break
                if c == "ATG":
                    start = j
                j -= 3
            if start is not None and stop - start + 3 <= length:
                aa = (stop - start) // 3
                if aa >= min_aa:
                    end = (stop + 3) % length if circular else stop + 3
                    found.add((strand, end, aa))
    return found


@pytest.mark.parametrize("circular", [False, True])
def test_matches_simple_implementation(circular):
    rng = random.Random(11)
    for _ in range(25):
        length = rng.randint(200, 600)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        mine = set()
        for orf in st.find_orfs(seq, min_aa=15, circular=circular):
            parts = st.feature_parts(orf)
            strand = st.feature_strand(orf)
            stop_end = parts[-1][1] if strand == 1 else length - parts[-1][0]
            mine.add((strand, stop_end % length if circular else stop_end, orf.qualifiers["aa"][0]))
        assert mine == simple_orf_finder(seq, 15, circular)


def test_every_orf_is_valid_on_random_dna():
    rng = random.Random(3)
    for _ in range(20):
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(300, 900)))
        for circular in (False, True):
            for orf in st.find_orfs(seq, min_aa=20, circular=circular):
                protein = protein_of(seq, orf)
                assert protein[0] == "M" and protein[-1] == "*"
                assert "*" not in protein[:-1]
