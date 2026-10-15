"""
abif_writer.py — write minimal ABIF (.ab1) chromatogram files.

N-Zoom only *reads* .ab1 files. This writer exists to make simulated
Sanger reads for the demo files and the tests, so they don't depend on
real sequencing data. It writes just the tags N-Zoom and Biopython use.
"""
import math
import random
import struct

TRACE_ORDER = "GATC"  # channel order stored in FWO_1, as on ABI machines


def _entry(tag, number, elem_type, elem_size, count, data):
    return (tag.encode(), number, elem_type, elem_size, count, data)


def write_ab1(path, name, seq, qualities, traces, peaks):
    """traces: {"A": [...], "C": [...], "G": [...], "T": [...]} (ints 0..32767);
    peaks: position of each base call in the traces."""
    entries = [
        _entry("PBAS", 2, 2, 1, len(seq), seq.encode()),
        _entry("PCON", 2, 2, 1, len(qualities), bytes(qualities)),
        _entry("PLOC", 2, 4, 2, len(peaks), struct.pack(f">{len(peaks)}h", *peaks)),
        _entry("FWO_", 1, 2, 1, 4, TRACE_ORDER.encode()),
        _entry("SMPL", 1, 18, 1, len(name) + 1, bytes([len(name)]) + name.encode()),
    ]
    for i, base in enumerate(TRACE_ORDER):
        values = traces[base]
        entries.append(_entry("DATA", 9 + i, 4, 2, len(values),
                              struct.pack(f">{len(values)}h", *values)))

    header_size = 128
    blob = bytearray()
    directory = bytearray()
    for tag, number, elem_type, elem_size, count, data in entries:
        if len(data) <= 4:  # small data lives inside the directory entry itself
            offset_field = data.ljust(4, b"\0")
        else:
            offset_field = struct.pack(">I", header_size + len(blob))
            blob += data
        directory += struct.pack(">4sI2H2I", tag, number, elem_type, elem_size, count, len(data))
        directory += offset_field + struct.pack(">I", 0)
    directory_offset = header_size + len(blob)

    header = b"ABIF" + struct.pack(">H4sI2H3I", 101, b"tdir", 1, 1023, 28,
                                   len(entries), len(directory), directory_offset)
    header += struct.pack(">I", 0)
    header = header.ljust(header_size, b"\0")
    with open(path, "wb") as handle:
        handle.write(header + bytes(blob) + bytes(directory))


def simulate_read(template, name, seed=0, spacing=12, bad_start=25, bad_end=60):
    """Make a realistic-looking Sanger read of `template` (the bases the
    sequencer reads, 5'->3'): Gaussian peaks, noise, low-quality ends with
    occasional miscalls. Returns (seq, qualities, traces, peaks)."""
    rng = random.Random(seed)
    n = len(template)
    length = n * spacing + spacing * 2
    traces = {base: [0.0] * length for base in "ACGT"}
    seq, quals, peaks = [], [], []
    for i, base in enumerate(template):
        # quality: poor at both ends, good in the middle
        edge = min(i / bad_start if bad_start else 1, (n - 1 - i) / bad_end if bad_end else 1, 1.0)
        quality = int(8 + 52 * edge ** 1.5 + rng.uniform(-4, 4))
        quality = max(4, min(60, quality))
        called = base
        if quality < 15 and rng.random() < 0.35:
            called = rng.choice("ACGTN".replace(base, ""))
        peak = spacing + i * spacing + rng.randint(-1, 1)
        height = rng.uniform(700, 1300) * (0.35 + 0.65 * edge)
        noise = (1.0 - edge) * 0.5 + 0.05
        for channel in "ACGT":
            h = height if channel == base else height * noise * rng.random()
            width = spacing * 0.33
            for x in range(max(0, peak - 3 * spacing // 2), min(length, peak + 3 * spacing // 2)):
                traces[channel][x] += h * math.exp(-((x - peak) ** 2) / (2 * width ** 2))
        seq.append(called)
        quals.append(quality)
        peaks.append(peak)
    ints = {b: [min(32767, int(v + rng.uniform(0, 15))) for v in values] for b, values in traces.items()}
    return "".join(seq), quals, ints, peaks
