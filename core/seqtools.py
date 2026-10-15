"""
seqtools.py — sequence logic for N-Zoom (no GUI code in here).

Type detection, file loading/saving, GenBank feature helpers,
restriction-site search and NCBI fetching.
"""
import copy
import re
from pathlib import Path

from Bio import SeqIO, Entrez
from Bio.Seq import Seq
from Bio.SeqFeature import CompoundLocation, SeqFeature, SimpleLocation

# ---------------------------------------------------------------- constants
DNA_CHARS = set("ACGTUNRYKMSWBDHV-")
PROTEIN_CHARS = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-")
MAX_COLORED_LENGTH = 50_000  # above this the sequence view shows plain text

GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank", ".gbff"}
GENBANK_NAME_LIMIT = 16      # max length of the name on a GenBank LOCUS line
DEFAULT_MIN_ORF_AA = 100     # default minimum ORF length, in amino acids

DNA_COLORS = {
    "A": "#2e7d32",
    "T": "#c62828",
    "U": "#c62828",
    "G": "#f9a825",
    "C": "#1565c0",
    "N": "#616161",
}

PROTEIN_COLORS = {
    "A": "#388e3c", "V": "#388e3c", "L": "#388e3c", "I": "#388e3c", "M": "#388e3c",
    "F": "#6a1b9a", "W": "#6a1b9a", "Y": "#6a1b9a",
    "K": "#1565c0", "R": "#1565c0", "H": "#1565c0",
    "D": "#c62828", "E": "#c62828",
    "S": "#00838f", "T": "#00838f", "N": "#00838f", "Q": "#00838f",
    "G": "#5d4037", "P": "#5d4037", "C": "#ef6c00",
    "*": "#000000",
}

FEATURE_COLORS = {
    "CDS": "#1e88e5",
    "gene": "#43a047",
    "mRNA": "#3949ab",
    "promoter": "#fb8c00",
    "terminator": "#e53935",
    "rep_origin": "#fdd835",
    "primer_bind": "#00897b",
    "protein_bind": "#6d4c41",
    "RBS": "#8e24aa",
    "regulatory": "#f4511e",
    "misc_feature": "#757575",
    "ORF": "#9ccc65",        # found by the ORF finder, not saved in the file
    "read": "#78909c",       # aligned sequencing reads, not saved in the file
}

# feature types offered in the "Add feature" dialog (any text is allowed too)
COMMON_FEATURE_TYPES = [
    "CDS", "gene", "promoter", "terminator", "rep_origin", "primer_bind",
    "protein_bind", "RBS", "regulatory", "misc_feature", "mRNA", "polyA_signal",
    "enhancer", "LTR", "sig_peptide", "mat_peptide",
]
DEFAULT_FEATURE_COLOR = "#90a4ae"
SKIP_FEATURE_TYPES = {"source"}  # covers the whole sequence, not useful to draw


# ---------------------------------------------------------------- sequences
def normalize_seq(seq):
    """Remove whitespace and force uppercase."""
    return "".join(seq.split()).upper()


def detect_seq_type(seq):
    """Return 'DNA', 'Protein' or 'Unknown'. RNA counts as DNA."""
    chars = set(normalize_seq(seq))
    if not chars:
        return "Unknown"
    if chars <= DNA_CHARS:
        return "DNA"
    if chars <= PROTEIN_CHARS:
        return "Protein"
    return "Unknown"


def record_sequence(record):
    """The record's sequence as a clean string ('' if the file has none)."""
    try:
        return normalize_seq(str(record.seq))
    except Exception:  # e.g. GenBank CONTIG records without sequence data
        return ""


# ---------------------------------------------------------------- files
def guess_format(path):
    return "genbank" if Path(path).suffix.lower() in GENBANK_SUFFIXES else "fasta"


def load_records(path):
    """Read every record from a FASTA or GenBank file. Returns (records, format)."""
    fmt = guess_format(path)
    return list(SeqIO.parse(path, fmt)), fmt


def genbank_locus_name(record):
    """A name that fits the GenBank LOCUS line: no spaces, at most 16 characters."""
    name = record.name if record.name and record.name != "<unknown name>" else record.id
    name = re.sub(r"\s+", "_", name or "unnamed")
    return name[:GENBANK_NAME_LIMIT]


def save_records(records, path, fmt):
    """Write records as 'genbank' or 'fasta'.

    The loaded records are not changed. Returns a list of (old name, new name)
    for records whose name had to be shortened to fit the GenBank format.
    """
    renamed = []
    if fmt == "genbank":
        prepared = []
        for record in records:
            out = copy.copy(record)  # shallow copy: same sequence and features
            out.annotations = dict(record.annotations)
            # GenBank output needs a molecule type; FASTA-loaded records lack one
            if "molecule_type" not in out.annotations:
                is_protein = detect_seq_type(record_sequence(record)) == "Protein"
                out.annotations["molecule_type"] = "protein" if is_protein else "DNA"
            safe = genbank_locus_name(record)
            if safe != record.name:
                if len(record.name or "") > GENBANK_NAME_LIMIT or " " in (record.name or ""):
                    renamed.append((record.name, safe))
                out.name = safe
            prepared.append(out)
        records = prepared
    SeqIO.write(records, path, fmt)
    return renamed


# ---------------------------------------------------------------- features
def is_circular(record):
    return str(record.annotations.get("topology", "")).lower() == "circular"


def topology(record):
    return str(record.annotations.get("topology", "")).lower() or "unknown"


def display_features(record):
    """Features worth showing (everything except 'source')."""
    return [
        f for f in record.features
        if f.type not in SKIP_FEATURE_TYPES and f.location is not None
    ]


def feature_label(feature):
    for key in ("label", "gene", "product", "locus_tag", "note"):
        values = feature.qualifiers.get(key)
        if values:
            return str(values[0])
    return feature.type


def feature_color(feature):
    return FEATURE_COLORS.get(feature.type, DEFAULT_FEATURE_COLOR)


def feature_parts(feature):
    """List of (start, end) — 0-based, end exclusive — one per location part.

    Parts are in biological order, so a feature that crosses the origin of a
    circular plasmid comes back as two parts.
    """
    if feature.location is None:
        return []
    return [(int(part.start), int(part.end)) for part in feature.location.parts]


def feature_strand(feature):
    """1, -1, or 0 for unknown / mixed."""
    strand = feature.location.strand if feature.location is not None else None
    return strand if strand in (1, -1) else 0


def wraps_origin(feature, length):
    """True for a two-part feature that runs over position 0 of a circular sequence."""
    parts = sorted(feature_parts(feature))
    return len(parts) == 2 and parts[0][0] == 0 and parts[1][1] == length


def has_simple_location(feature, length):
    """One continuous range, or one range that crosses the origin.
    (Features made of several pieces, such as exons, are not 'simple'.)"""
    return len(feature_parts(feature)) == 1 or wraps_origin(feature, length)


def make_location(start, end, strand, length, circular):
    """Build a location from 1-based, inclusive start/end (as people type them).

    On a circular sequence start > end means the feature crosses the origin.
    strand: 1, -1, or 0/None for no strand.
    """
    if not (1 <= start <= length and 1 <= end <= length):
        raise ValueError(f"Start and end must be between 1 and {length:,}.")
    strand = strand if strand in (1, -1) else None
    if start <= end:
        return SimpleLocation(start - 1, end, strand)
    if not circular:
        raise ValueError(
            "Start is after end. Only features on circular sequences can cross the origin."
        )
    tail = SimpleLocation(start - 1, length, strand)  # start .. end of sequence
    head = SimpleLocation(0, end, strand)             # 1 .. end
    # parts are stored in biological (reading) order
    return CompoundLocation([head, tail] if strand == -1 else [tail, head])


def new_feature(feature_type, label, location, note=""):
    qualifiers = {"label": [label]} if label else {}
    if note:
        qualifiers["note"] = [note]
    return SeqFeature(location, type=feature_type, qualifiers=qualifiers)


def make_feature(feature_type, label, start, end, strand, length, circular, note=""):
    """Feature from 1-based, inclusive coordinates (see make_location)."""
    return new_feature(feature_type, label,
                       make_location(start, end, strand, length, circular), note)


def update_feature(feature, feature_type, label, location=None, note=None):
    """Change a feature in place."""
    feature.type = feature_type
    if label:
        feature.qualifiers["label"] = [label]
    else:
        feature.qualifiers.pop("label", None)
    if note is not None:
        if note:
            feature.qualifiers["note"] = [note]
        else:
            feature.qualifiers.pop("note", None)
    if location is not None:
        feature.location = location


def feature_note(feature):
    values = feature.qualifiers.get("note")
    return str(values[0]) if values else ""


# ---------------------------------------------------------------- regions
def feature_sequence(sequence, feature):
    """The feature's sequence read 5'->3' on its own strand (joins handled)."""
    return str(feature.location.extract(Seq(sequence)))


def reverse_complement(dna):
    return str(Seq(dna).reverse_complement())


def gc_percent(dna):
    return (dna.count("G") + dna.count("C")) / len(dna) * 100 if dna else 0.0


def translate_dna(dna, table=1, codon_start=1):
    """Translate from base `codon_start` (1-3). Leftover bases are ignored."""
    clean = normalize_seq(dna).replace("-", "").replace("U", "T")[codon_start - 1:]
    usable = len(clean) - len(clean) % 3
    return str(Seq(clean[:usable]).translate(table=table)) if usable else ""


def feature_translation_settings(feature):
    """(table, codon_start) from the feature's qualifiers, defaulting to (1, 1)."""
    def first_int(key, default):
        try:
            return int(feature.qualifiers.get(key, [default])[0])
        except (TypeError, ValueError):
            return default
    codon_start = first_int("codon_start", 1)
    return first_int("transl_table", 1), codon_start if codon_start in (1, 2, 3) else 1


# ---------------------------------------------------------------- ORFs
def find_orfs(sequence, min_aa=DEFAULT_MIN_ORF_AA, circular=False, table=1):
    """Open reading frames (ATG ... stop) in all six frames.

    Returns SeqFeatures of type "ORF", sorted by position. Each runs from the
    first ATG after a stop to the next stop codon (included). min_aa counts
    amino acids without the stop. On circular sequences ORFs may cross the
    origin.
    """
    from Bio.Data import CodonTable

    stops = set(CodonTable.unambiguous_dna_by_id[table].stop_codons)
    seq = normalize_seq(sequence).replace("U", "T")
    length = len(seq)
    if length < 6:
        return []

    best = {}  # (strand, stop position) -> (aa, start, end) in that strand's coordinates
    for strand, strand_seq in ((1, seq), (-1, reverse_complement(seq))):
        # scanning the sequence twice in a row finds ORFs that run over the origin
        scan = strand_seq + strand_seq if circular else strand_seq
        for frame in range(3):
            orf_start = None
            for i in range(frame, len(scan) - 2, 3):
                if orf_start is None:
                    if i >= length:
                        break  # no new ORFs start in the second copy
                    if scan[i:i + 3] == "ATG":
                        orf_start = i
                    continue
                if i - orf_start >= length:
                    orf_start = None  # went all the way round without a stop
                    continue
                if scan[i:i + 3] in stops:
                    end = i + 3
                    aa = (end - orf_start) // 3 - 1
                    key = (strand, end % length)
                    if aa >= min_aa and aa > best.get(key, (0,))[0]:
                        best[key] = (aa, orf_start, end)
                    orf_start = None

    orfs = []
    for (strand, _), (aa, a, b) in best.items():
        if strand == 1:
            if b <= length:
                location = SimpleLocation(a, b, 1)
            else:
                location = CompoundLocation([SimpleLocation(a, length, 1),
                                             SimpleLocation(0, b - length, 1)])
            frame = f"+{a % 3 + 1}"
        else:
            # reverse-complement coordinate x is forward coordinate length - x
            if b <= length:
                location = SimpleLocation(length - b, length - a, -1)
            else:
                location = CompoundLocation([SimpleLocation(0, length - a, -1),
                                             SimpleLocation(2 * length - b, length, -1)])
            frame = f"-{a % 3 + 1}"
        orfs.append(SeqFeature(location, type="ORF",
                               qualifiers={"note": [f"frame {frame}, {aa} aa"], "aa": [aa]}))

    orfs.sort(key=lambda f: min(s for s, _ in feature_parts(f)))
    for number, orf in enumerate(orfs, 1):
        orf.qualifiers["label"] = [f"ORF{number} ({orf.qualifiers['aa'][0]} aa)"]
    return orfs


# ---------------------------------------------------------------- motif search
IUPAC_PATTERNS = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[CG]", "W": "[AT]", "K": "[GT]", "M": "[AC]",
    "B": "[CGT]", "D": "[AGT]", "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}


def motif_pattern(motif, dna=True):
    """Regular expression for a search motif.
    DNA: IUPAC codes allowed (N = any base). Protein: letters, X = any residue."""
    if dna:
        bad = sorted(set(motif) - set(IUPAC_PATTERNS))
        if bad:
            raise ValueError(f"Not a DNA letter or IUPAC code: {', '.join(bad)}")
        return "".join(IUPAC_PATTERNS[c] for c in motif)
    bad = sorted(set(motif) - set(PROTEIN_CHARS) - {"X"})
    if bad:
        raise ValueError(f"Not an amino-acid letter: {', '.join(bad)}")
    return "".join("." if c == "X" else re.escape(c) for c in motif)


def find_motif(sequence, motif, circular=False, both_strands=True, dna=True, max_hits=None):
    """Every occurrence of `motif`, including overlapping ones.

    DNA motifs may use IUPAC codes and are searched on both strands (unless
    both_strands is False); a palindromic motif such as GAATTC is reported
    once per site, with strand 0. On circular sequences matches may run over
    the origin. Protein motifs are searched as text, with X matching anything.

    Returns a list of dicts sorted by position: strand (1, -1 or 0),
    start, end (1-based, inclusive; start > end across the origin), location.
    """
    motif = normalize_seq(motif)
    if dna:
        motif = motif.replace("U", "T")
    template = normalize_seq(sequence)
    if dna:
        template = template.replace("U", "T")
    length, size = len(template), len(motif)
    if not motif or size > length:
        return []

    searches = [(1, motif)]
    if dna and both_strands:
        reverse = reverse_complement(motif)  # Biopython complements IUPAC codes too
        if reverse == motif:
            searches = [(0, motif)]  # palindrome: same sites on both strands
        else:
            searches.append((-1, reverse))
    searchable = template + template[:size - 1] if circular else template

    hits = []
    for strand, pattern_motif in searches:
        regex = re.compile(f"(?=({motif_pattern(pattern_motif, dna)}))")
        for match in regex.finditer(searchable):
            i = match.start()
            if i >= length:
                break
            start1, end1 = i + 1, (i + size - 1) % length + 1
            hits.append({"strand": strand, "start": start1, "end": end1,
                         "location": make_location(start1, end1, strand, length, circular)})
            if max_hits and len(hits) >= max_hits:
                break
    hits.sort(key=lambda hit: (hit["start"], -hit["strand"]))
    return hits[:max_hits] if max_hits else hits


# ---------------------------------------------------------------- sequencing reads
READ_SUFFIXES = {".ab1": "abi", ".abi": "abi", ".fastq": "fastq", ".fq": "fastq",
                 ".fasta": "fasta", ".fa": "fasta", ".fas": "fasta", ".fna": "fasta",
                 ".seq": "text", ".txt": "text"}
TRIM_ERROR_LIMIT = 0.05      # quality trimming: keep bases better than 5% error
LOW_QUALITY = 20             # differences at bases below this are flagged
MIN_ALIGNED = 30             # shorter alignments count as "no match"
SEED_SIZE = 14               # k-mer size used to locate reads in long sequences


def load_read(path):
    """Read one sequencing read. Returns a dict:
    name, seq, quals (list or None), trace (dict or None, .ab1 only).

    .ab1/.abi: chromatogram files from ABI sequencers.  .fastq: bases with
    qualities.  .fasta and .seq/.txt: bases only."""
    path = Path(path)
    kind = READ_SUFFIXES.get(path.suffix.lower(), "fasta")
    name = path.stem
    if kind == "text":
        text = path.read_text(errors="replace")
        lines = [line for line in text.splitlines() if not line.startswith(">")]
        seq = normalize_seq("".join(lines))
        return {"name": name, "seq": seq, "quals": None, "trace": None, "path": str(path)}
    record = SeqIO.read(str(path), kind)
    quals = record.letter_annotations.get("phred_quality")
    trace = None
    if kind == "abi":
        raw = record.annotations.get("abif_raw", {})
        order = raw.get("FWO_1", b"GATC")
        order = order.decode() if isinstance(order, bytes) else str(order)
        channels = {base: list(raw.get(f"DATA{9 + i}") or ()) for i, base in enumerate(order[:4])}
        peaks = raw.get("PLOC2") or raw.get("PLOC1")
        if all(channels.values()) and peaks:
            trace = {"channels": channels, "peaks": list(peaks)}
    return {"name": name, "seq": normalize_seq(str(record.seq)),
            "quals": list(quals) if quals else None, "trace": trace, "path": str(path)}


def quality_trim(quals, limit=TRIM_ERROR_LIMIT):
    """(start, end) of the region worth keeping, by the "modified Mott"
    method used by phred: each base scores limit - P(error), and the
    best-scoring stretch is kept. Without qualities, keep everything."""
    if not quals:
        return 0, 0
    best_sum, best = 0.0, (0, 0)
    running, start = 0.0, 0
    for i, q in enumerate(quals):
        running += limit - 10 ** (-q / 10)
        if running <= 0:
            running, start = 0.0, i + 1
        elif running > best_sum:
            best_sum, best = running, (start, i + 1)
    return best


def _seed_offset(reference, read):
    """Most common (reference - read) offset of shared k-mers, with its votes."""
    votes = {}
    step = max(1, SEED_SIZE // 2)
    for i in range(0, len(read) - SEED_SIZE + 1, step):
        kmer = read[i:i + SEED_SIZE]
        if "N" in kmer:
            continue
        j = reference.find(kmer)
        hits = 0
        while j != -1 and hits < 20:  # ignore k-mers in highly repetitive sequence
            votes[j - i] = votes.get(j - i, 0) + 1
            j = reference.find(kmer, j + 1)
            hits += 1
    if not votes:
        return None, 0
    offset = max(votes, key=votes.get)
    return offset, votes[offset]


def align_read(reference, read, quals=None, circular=False):
    """Align a sequencing read to the reference, trying both orientations.

    The read is quality-trimmed first (if it has qualities). Returns a dict
    describing the best local alignment, or None if the read doesn't match:
      strand          1 = read runs along the top strand, -1 = reverse
      start, end      1-based, inclusive, reference coordinates covered
                      (start > end across the origin of a circular sequence)
      location        the same as a Biopython location
      identity        % of alignment columns that match
      columns, matches
      differences     list of dicts: kind ("mismatch", "insertion",
                      "deletion"), ref_pos (1-based), ref, read (bases),
                      quality (or None), read_index (0-based, in the read as
                      sequenced), low_quality (True if quality < 20)
      rows            (reference row, match row, read row, first ref
                      position) for printing the alignment
      trim            (start, end) of the read that was used
    """
    from Bio import Align

    reference = normalize_seq(reference).replace("U", "T")
    read = normalize_seq(read).replace("U", "T")
    length = len(reference)
    trim_start, trim_end = quality_trim(quals) if quals else (0, len(read))
    if trim_end - trim_start < MIN_ALIGNED:
        return None
    used = read[trim_start:trim_end]
    used_quals = quals[trim_start:trim_end] if quals else None

    target = reference + reference[:len(used)] if circular else reference
    aligner = Align.PairwiseAligner(mode="local", match_score=2, mismatch_score=-3,
                                    open_gap_score=-5, extend_gap_score=-2)
    best = None
    for strand, query in ((1, used), (-1, reverse_complement(used))):
        window_start = 0
        window = target
        if len(target) > 50_000:  # long references: align only near the k-mer hits
            offset, votes = _seed_offset(target, query)
            if offset is None or votes < 2:
                continue
            window_start = max(0, offset - 500)
            window = target[window_start:offset + len(query) + 500]
        score = aligner.score(window, query)
        if best is None or score > best[0]:
            best = (score, strand, query, window, window_start)
    if best is None or best[0] <= 0:
        return None
    score, strand, query, window, window_start = best
    alignment = aligner.align(window, query)[0]

    # walk the alignment block by block
    oriented_quals = (used_quals if strand == 1 else used_quals[::-1]) if used_quals else None

    def read_index(q):  # position in the read as the sequencer produced it
        return trim_start + q if strand == 1 else trim_end - 1 - q

    def quality(q):
        return oriented_quals[q] if oriented_quals else None

    coords = alignment.coordinates
    ref_row, read_row, mid_row = [], [], []
    differences = []
    matches = columns = 0
    for k in range(coords.shape[1] - 1):
        t0, t1 = int(coords[0][k]), int(coords[0][k + 1])
        q0, q1 = int(coords[1][k]), int(coords[1][k + 1])
        if t1 > t0 and q1 > q0:  # aligned block
            for t, q in zip(range(t0, t1), range(q0, q1)):
                r_base, q_base = window[t], query[q]
                columns += 1
                ref_row.append(r_base)
                read_row.append(q_base)
                if r_base == q_base:
                    matches += 1
                    mid_row.append("|")
                else:
                    mid_row.append(" ")
                    differences.append({
                        "kind": "mismatch", "ref_pos": (window_start + t) % length + 1,
                        "ref": r_base, "read": q_base, "quality": quality(q),
                        "read_index": read_index(q)})
        elif t1 > t0:  # bases missing from the read
            columns += t1 - t0
            ref_row.append(window[t0:t1])
            read_row.append("-" * (t1 - t0))
            mid_row.append(" " * (t1 - t0))
            neighbours = [x for x in (quality(q0 - 1) if q0 > 0 else None,
                                      quality(q0) if q0 < len(query) else None) if x is not None]
            differences.append({
                "kind": "deletion", "ref_pos": (window_start + t0) % length + 1,
                "ref": window[t0:t1], "read": "-", "quality": min(neighbours) if neighbours else None,
                "read_index": read_index(min(q0, len(query) - 1))})
        elif q1 > q0:  # extra bases in the read
            columns += q1 - q0
            ref_row.append("-" * (q1 - q0))
            read_row.append(query[q0:q1])
            mid_row.append(" " * (q1 - q0))
            qs = [quality(q) for q in range(q0, q1) if quality(q) is not None]
            differences.append({
                "kind": "insertion", "ref_pos": (window_start + t0 - 1) % length + 1,
                "ref": "-", "read": query[q0:q1], "quality": min(qs) if qs else None,
                "read_index": read_index(q0)})
    for difference in differences:
        q = difference["quality"]
        difference["low_quality"] = q is not None and q < LOW_QUALITY

    if columns < MIN_ALIGNED:
        return None
    first = int(coords[0][0]) + window_start
    last = int(coords[0][-1]) + window_start  # exclusive
    start1 = first % length + 1
    end1 = (last - 1) % length + 1
    return {
        "strand": strand, "start": start1, "end": end1,
        "location": make_location(start1, end1, strand, length, circular),
        "identity": 100.0 * matches / columns, "columns": columns, "matches": matches,
        "differences": differences, "score": score,
        "rows": ("".join(ref_row), "".join(mid_row), "".join(read_row), start1),
        "trim": (trim_start, trim_end),
    }


def align_reads(reference, paths, circular=False):
    """Load and align several read files. Each result is load_read()'s dict
    plus "alignment" (align_read's dict, or None) and "error" (or None)."""
    results = []
    for path in paths:
        try:
            read = load_read(path)
        except Exception as exc:  # unreadable file: report it, carry on with the rest
            results.append({"name": Path(path).stem, "seq": "", "quals": None, "trace": None,
                            "path": str(path), "alignment": None,
                            "error": f"Could not read file: {exc}"})
            continue
        read["alignment"] = align_read(reference, read["seq"], read["quals"], circular)
        read["error"] = None if read["alignment"] else "No match to this sequence"
        results.append(read)
    return results


def format_alignment(rows, width=60):
    """Alignment as text blocks with reference positions, like BLAST output."""
    ref_row, mid_row, read_row, first = rows
    lines = []
    position = first
    for i in range(0, len(ref_row), width):
        ref_part = ref_row[i:i + width]
        bases = sum(1 for c in ref_part if c != "-")
        lines.append(f"Ref  {position:>7,}  {ref_part}  {position + bases - 1:,}")
        lines.append(f"     {'':>7}  {mid_row[i:i + width]}")
        lines.append(f"Read {'':>7}  {read_row[i:i + width]}")
        lines.append("")
        position += bases
    return "\n".join(lines)


# ---------------------------------------------------------------- primers
PRIMER_MIN_ANNEAL = 15       # 3' bases that must match for a binding site
PRIMER_NA_MM = 50            # Tm conditions: monovalent salt (mM)
PRIMER_CONC_NM = 250         # Tm conditions: primer concentration (nM)


def clean_primer(primer):
    """Uppercase, no spaces, U -> T. Raises ValueError for anything but A/C/G/T."""
    cleaned = normalize_seq(primer).replace("U", "T")
    if not cleaned:
        raise ValueError("The primer is empty.")
    bad = sorted(set(cleaned) - set("ACGT"))
    if bad:
        raise ValueError(f"Primers may only contain A, C, G and T (found: {', '.join(bad)}).")
    return cleaned


def primer_tm(primer):
    """Melting temperature (deg C) by the nearest-neighbour method
    (Biopython Tm_NN, 50 mM Na+, 250 nM primer)."""
    from Bio.SeqUtils import MeltingTemp

    return MeltingTemp.Tm_NN(primer, Na=PRIMER_NA_MM, dnac1=PRIMER_CONC_NM, dnac2=0)


def find_primer_sites(sequence, primer, circular=False, min_anneal=PRIMER_MIN_ANNEAL):
    """Where a primer binds, on either strand.

    A site needs at least `min_anneal` bases at the primer's 3' end (or the
    whole primer, if shorter) to match the template exactly. The match is then
    extended toward the 5' end as far as bases keep matching; anything beyond
    that is a 5' tail (e.g. an added restriction site).

    Returns a list of dicts, sorted by position:
      strand     1 = primer reads along the top strand (forward primer), -1 = reverse
      start, end 1-based, inclusive, top-strand coordinates of the annealed part
                 (start > end if it crosses the origin of a circular sequence)
      anneal     number of annealed bases;  tail = primer length - anneal
      location   Biopython location of the annealed part
    """
    primer = clean_primer(primer)
    template = normalize_seq(sequence).replace("U", "T")
    length, plen = len(template), len(primer)
    seed_len = min(min_anneal, plen)
    seed = primer[-seed_len:]
    if seed_len > length:
        return []

    sites = []
    for strand, strand_seq in ((1, template), (-1, reverse_complement(template))):
        # the primer's own sequence appears in the strand it reads along
        searchable = strand_seq + strand_seq[:plen - 1] if circular else strand_seq
        i = searchable.find(seed)
        while i != -1 and i < length:
            end = i + seed_len  # exclusive, in this strand's coordinates
            anneal = seed_len
            while anneal < plen:
                j = end - anneal - 1
                if j < 0 and not circular:
                    break
                if strand_seq[j % length] != primer[plen - anneal - 1]:
                    break
                anneal += 1
                if anneal >= length:
                    break
            a, b = end - anneal, end  # annealed part in this strand's coordinates
            if strand == -1:          # convert to top-strand coordinates
                a, b = length - b, length - a
            start1 = a % length + 1
            end1 = (b - 1) % length + 1
            sites.append({
                "strand": strand, "start": start1, "end": end1,
                "anneal": anneal, "tail": plen - anneal,
                "location": make_location(start1, end1, strand, length, circular),
            })
            i = searchable.find(seed, i + 1)
    sites.sort(key=lambda site: (site["start"], -site["strand"]))
    return sites


# ---------------------------------------------------------------- digests
# Band sizes (bp) of two common, generic DNA ladders
LADDERS = {
    "1 kb ladder": [10000, 8000, 6000, 5000, 4000, 3500, 3000, 2500, 2000,
                    1500, 1000, 750, 500, 250],
    "100 bp ladder": [1500, 1200, 1000, 900, 800, 700, 600, 500, 400, 300, 200, 100],
}


def digest(sequence, enzymes, circular=False):
    """Cut a sequence with one or more enzymes.

    Returns (cuts, fragments):
      cuts       sorted list of (position, "Enzyme[/Enzyme]"), 1-based position
                 of the first base after the cut on the top strand
      fragments  list of dicts in sequence order: start, end (1-based,
                 inclusive; start > end if the fragment crosses the origin),
                 length, left and right (enzyme at each end, None for a
                 sequence end), uncut (True only for an uncut circle)

    Fragment sizes are top-strand lengths; sticky-end overhangs are ignored.
    """
    from Bio.Restriction import RestrictionBatch

    clean = normalize_seq(sequence).replace("-", "").replace("U", "T")
    length = len(clean)
    names = sorted(set(enzymes))
    by_position = {}
    if names:
        results = RestrictionBatch(names).search(Seq(clean), linear=not circular)
        for enzyme, positions in results.items():
            for pos in positions:
                pos = (pos - 1) % length + 1  # keep positions within 1..length
                by_position.setdefault(pos, []).append(str(enzyme))
    cuts = [(pos, "/".join(sorted(found))) for pos, found in sorted(by_position.items())]
    if not circular:
        cuts = [(p, e) for p, e in cuts if 1 < p <= length]  # a cut before base 1 does nothing

    fragments = []
    if circular:
        if not cuts:
            return cuts, [{"start": 1, "end": length, "length": length,
                           "left": None, "right": None, "uncut": True}]
        for k, (pos, enzyme) in enumerate(cuts):
            next_pos, next_enzyme = cuts[(k + 1) % len(cuts)]
            size = (next_pos - pos) % length or length
            fragments.append({"start": pos, "end": (next_pos - 2) % length + 1, "length": size,
                              "left": enzyme, "right": next_enzyme, "uncut": False})
    else:
        bounds = [(1, None)] + cuts + [(length + 1, None)]
        for (pos, enzyme), (next_pos, next_enzyme) in zip(bounds, bounds[1:]):
            if next_pos > pos:
                fragments.append({"start": pos, "end": next_pos - 1, "length": next_pos - pos,
                                  "left": enzyme, "right": next_enzyme, "uncut": False})
    return cuts, fragments


# ---------------------------------------------------------------- restriction
def find_restriction_sites(sequence, circular=True, commercial_only=True):
    """Return {enzyme name: [cut positions]} for every enzyme that cuts.

    Positions are 1-based: the first base after the cut on the top strand.
    """
    from Bio.Restriction import AllEnzymes, Analysis, CommOnly, RestrictionBatch

    clean = normalize_seq(sequence).replace("-", "").replace("U", "T")
    batch = RestrictionBatch(CommOnly if commercial_only else AllEnzymes)
    analysis = Analysis(batch, Seq(clean), linear=not circular)
    return {
        str(enzyme): sorted(positions)
        for enzyme, positions in analysis.full().items()
        if positions
    }


def map_sites(sites, min_site_size=6):
    """Pick the sites worth drawing on a map: enzymes that cut exactly once and
    recognise at least `min_site_size` bases. Enzymes cutting at the same spot
    (isoschizomers) are merged into one label. Returns {label: [position]}."""
    from Bio.Restriction import AllEnzymes

    by_position = {}
    for name, positions in sites.items():
        enzyme = AllEnzymes.get(name)
        if len(positions) == 1 and enzyme.size >= min_site_size:
            by_position.setdefault(positions[0], []).append(name)
    return {", ".join(sorted(names)): [pos] for pos, names in by_position.items()}


def feature_span(feature):
    """(start, end) to display, 1-based inclusive, in genome order.
    For a feature crossing the origin start > end, e.g. (2801, 80)."""
    parts = feature_parts(feature)
    if feature_strand(feature) == -1:
        parts = parts[::-1]  # reverse-strand parts are stored in biological order
    return parts[0][0] + 1, parts[-1][1]


# ---------------------------------------------------------------- NCBI
def fetch_genbank(accession, email):
    """Download one nucleotide record from NCBI as a GenBank SeqRecord."""
    Entrez.email = email  # NCBI requires an e-mail address with each request
    with Entrez.efetch(
        db="nuccore", id=accession, rettype="gbwithparts", retmode="text"
    ) as handle:
        return SeqIO.read(handle, "genbank")
