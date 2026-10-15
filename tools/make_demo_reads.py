"""
make_demo_reads.py — create the simulated Sanger reads in examples/.

    python tools/make_demo_reads.py

They're reads of a pDEMO "clone" with one point mutation at position 1,000
(inside GFP), so the read-alignment tool has something to find. The script
prints which change it made.
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from core import seqtools as st  # noqa: E402
from tools.abif_writer import simulate_read, write_ab1  # noqa: E402

MUTATION = 1000  # 1-based position of the planted point mutation


def main():
    record = st.load_records(ROOT / "examples" / "pDEMO.gb")[0][0]
    plasmid = st.record_sequence(record)
    old = plasmid[MUTATION - 1]
    new = {"A": "G", "G": "A", "C": "T", "T": "C"}[old]
    clone = plasmid[:MUTATION - 1] + new + plasmid[MUTATION:]

    reads = {
        # name: (template the sequencer reads 5'->3', seed)
        "GFP_fwd": (clone[640:1460], 1),                         # forward primer upstream of GFP
        "GFP_rev": (st.reverse_complement(clone[620:1480]), 2),  # reverse primer downstream
        "ori_fwd": (clone[3040:3760], 3),                        # a clean read elsewhere
    }
    for name, (template, seed) in reads.items():
        seq, quals, traces, peaks = simulate_read(template, name, seed=seed)
        path = ROOT / "examples" / f"pDEMO_{name}.ab1"
        write_ab1(path, name, seq, quals, traces, peaks)
        print(f"wrote {path.name}: {len(seq)} bases")
    print(f"planted mutation: {old}{MUTATION}{new}")


if __name__ == "__main__":
    main()
