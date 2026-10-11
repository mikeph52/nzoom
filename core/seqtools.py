# N-Zoom seqtools.py by mikeph52, Init 11/10/2026
# Here there's all the seq features of nzoom.

from pathlib import Path
from Bio import SeqIO, Entrez
from Bio.Seq import Seq

DNA_CHARS = set("ACGTUNRYKMSWBDHV-")
PROTEIN_CHARS = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-")
MAX_COLORED_LENGTH = 50_000  # above this the sequence view shows plain text
GENBANK_SUFFIXES = {".gb", ".gbk", ".genbank", ".gbff"}

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
}
DEFAULT_FEATURE_COLOR = "#90a4ae"
SKIP_FEATURE_TYPES = {"source"}  # covers the whole sequence, not useful to draw and makes the app laggyu

def normalize_seq(seq):
    # Remove whitespace and force uppercase.
    return "".join(seq.split()).upper()

def detect_seq_type(seq):
    # Return 'DNA', 'Protein' or 'Unknown'. RNA counts as DNA.
    chars = set(normalize_seq(seq))
    if not chars:
        return "Unknown"
    if chars <= DNA_CHARS:
        return "DNA"
    if chars <= PROTEIN_CHARS:
        return "Protein"
    return "Unknown"

def record_sequence(record):
    # The record's sequence as a clean string
    try:
        return normalize_seq(str(record.seq))
    except Exception:  # e.g. GenBank CONTIG records without sequence data i guess
        return ""

def guess_format(path):
    return "genbank" if Path(path).suffix.lower() in GENBANK_SUFFIXES else "fasta"

def load_records(path):
    # Read every record from a FASTA or GenBank file and returns records and format
    fmt = guess_format(path)
    return list(SeqIO.parse(path, fmt)), fmt

def save_records(records, path, fmt):
    # Write records as gb or fasta
    if fmt == "genbank":
        for record in records:
            # gb output needs a molecule type unlike FASTA
            if "molecule_type" not in record.annotations:
                is_protein = detect_seq_type(record_sequence(record)) == "Protein"
                record.annotations["molecule_type"] = "protein" if is_protein else "DNA"
    SeqIO.write(records, path, fmt)


# plasmid features
def is_circular(record):
    return str(record.annotations.get("topology", "")).lower() == "circular"

def topology(record):
    return str(record.annotations.get("topology", "")).lower() or "unknown"

def display_features(record):
    # Features worth showing
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
    # 1, -1, or 0 for unknown / mixed.
    strand = feature.location.strand if feature.location is not None else None
    return strand if strand in (1, -1) else 0

# restriction sites finder
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

# fetch from NCBI
def fetch_genbank(accession, email):
    Entrez.email = email  # NCBI requires an e-mail address with each request we will see how to handle that
    with Entrez.efetch(
        db="nuccore", id=accession, rettype="gbwithparts", retmode="text"
    ) as handle:
        return SeqIO.read(handle, "genbank")