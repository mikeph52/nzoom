"""
Tests for primer search, Tm and virtual digests (seqtools, no GUI).
"""
import pytest

from core import seqtools as st


@pytest.fixture
def demo_seq(demo_record):
    return st.record_sequence(demo_record)


# ------------------------------------------------------------------ primers
def test_forward_primer(demo_seq):
    primer = demo_seq[700:720]  # first 20 bases of GFP
    [site] = st.find_primer_sites(demo_seq, primer, circular=True)
    assert (site["strand"], site["start"], site["end"]) == (1, 701, 720)
    assert site["anneal"] == 20 and site["tail"] == 0


def test_reverse_primer(demo_seq):
    primer = st.reverse_complement(demo_seq[1397:1417])  # last 20 bases of GFP
    [site] = st.find_primer_sites(demo_seq, primer, circular=True)
    assert (site["strand"], site["start"], site["end"]) == (-1, 1398, 1417)


def test_primer_with_5_prime_tail(demo_seq):
    tail = "TTTTTTTTTT"  # can't anneal next to GFP's start (base 700 isn't T)
    assert demo_seq[699] != "T"
    [site] = st.find_primer_sites(demo_seq, tail + demo_seq[700:720], circular=True)
    assert site["anneal"] == 20 and site["tail"] == 10
    assert (site["start"], site["end"]) == (701, 720)


@pytest.mark.parametrize("reverse", [False, True])
def test_primer_across_origin(demo_seq, reverse):
    primer = demo_seq[-10:] + demo_seq[:12]
    if reverse:
        primer = st.reverse_complement(primer)
    [site] = st.find_primer_sites(demo_seq, primer, circular=True)
    assert (site["start"], site["end"], site["anneal"]) == (3991, 12, 22)
    assert st.find_primer_sites(demo_seq, primer, circular=False) == []


def test_primer_binding_twice():
    template = "GGGG" + "ACGTACGTTGCAAGCT" * 1 + "CCCCCCCC" + "ACGTACGTTGCAAGCT" + "GGGG"
    sites = st.find_primer_sites(template, "ACGTACGTTGCAAGCT")
    assert [s["start"] for s in sites] == [5, 29]


def test_3_prime_mismatch_means_no_site(demo_seq):
    primer = demo_seq[700:719] + ("A" if demo_seq[719] != "A" else "C")
    assert st.find_primer_sites(demo_seq, primer, circular=True) == []


def test_invalid_primer():
    with pytest.raises(ValueError, match="only contain"):
        st.clean_primer("ACGTNNACGT")
    assert st.clean_primer("acg u t") == "ACGTT"


def test_primer_tm_is_plausible():
    tm = st.primer_tm("GTAAAACGACGGCCAGTGAG")  # a common 20-mer, ~55% GC
    assert 52 < tm < 64
    assert st.primer_tm("GCGCGCGCGCGCGCGCGCGC") > st.primer_tm("ATATATATATATATATATAT")


# ------------------------------------------------------------------ digests
def test_digest_linear(demo_seq):
    cuts, fragments = st.digest(demo_seq, ["EcoRI", "HindIII"], circular=False)
    assert cuts == [(602, "EcoRI"), (653, "HindIII")]
    assert [(f["start"], f["end"], f["length"]) for f in fragments] == [
        (1, 601, 601), (602, 652, 51), (653, 4000, 3348)]
    assert fragments[0]["left"] is None and fragments[-1]["right"] is None


def test_digest_circular(demo_seq):
    _, fragments = st.digest(demo_seq, ["EcoRI", "HindIII"], circular=True)
    assert [(f["start"], f["end"], f["length"]) for f in fragments] == [
        (602, 652, 51), (653, 601, 3949)]  # the big one runs across the origin


def test_single_cut_linearises_a_circle(demo_seq):
    _, [fragment] = st.digest(demo_seq, ["EcoRI"], circular=True)
    assert fragment["length"] == 4000 and not fragment["uncut"]


def test_no_cuts(demo_seq):
    _, [fragment] = st.digest(demo_seq, [], circular=True)
    assert fragment["uncut"] and fragment["length"] == 4000


@pytest.mark.parametrize("circular", [False, True])
def test_fragments_add_up(demo_seq, circular):
    sites = st.find_restriction_sites(demo_seq, circular=circular)
    three_cutters = [name for name, positions in sites.items() if len(positions) == 3][:5]
    for enzyme in three_cutters:
        _, fragments = st.digest(demo_seq, [enzyme], circular)
        assert sum(f["length"] for f in fragments) == len(demo_seq)
        assert len(fragments) == (3 if circular else 4)


def test_same_position_enzymes_are_merged():
    # SpeI and BcuI recognise A^CTAGT and cut at the same place: one cut, both names.
    # Acc65I (G^GTACC) and KpnI (GGTAC^C) share a site but cut at different places.
    sequence = "AAAAAAAAAA" + "ACTAGT" + "A" * 20 + "GGTACC" + "A" * 20
    cuts, fragments = st.digest(sequence, ["SpeI", "BcuI", "Acc65I", "KpnI"], circular=False)
    assert cuts == [(12, "BcuI/SpeI"), (38, "Acc65I"), (42, "KpnI")]
    assert [f["length"] for f in fragments] == [11, 26, 4, len(sequence) - 41]
