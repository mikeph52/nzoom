# N-Zoom

**An open-source DNA and protein sequence viewer, written in Python.**

N-Zoom is a lightweight, free alternative to commercial plasmid tools. Open a FASTA or GenBank file, or pull a sequence straight from NCBI. You can then explore it on an annotated circular map.

> **Status:** early alpha (v0.7.5). A first stable release is planned for 2027.

---

## Features

- **Open FASTA and GenBank files**, including files that hold several records
- **Circular plasmid maps** with colored, strand-aware feature arcs and auto-placed labels; click an arc to select that feature
- **Feature table** that is sortable and linked to both the map and the sequence view
- **Colored sequence view** for DNA and protein
- **Restriction site analysis** over all commercially available enzymes, with unique cutters marked on the map
- **Fetch sequences from NCBI** by accession number
- **Translate DNA, reverse complement, and GC content and composition stats**
- **Save as GenBank or FASTA**

## Requirements

- Python 3.10 or newer, with pip
- [PySide6](https://pypi.org/project/PySide6/) 6.7 or newer
- [Biopython](https://biopython.org/) 1.83 or newer

## Installation

Clone the repository, then install the dependencies:

```bash
git clone <repository-url>
cd N-Zoom
pip install -r requirements.txt
```

Using a virtual environment is recommended:

```bash
python -m venv .venv
source .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

## Quick start

```bash
python main.py
```

Then try one of these:

- **File → Open…** and pick any `.gb` or `.fasta` file.
- Click **Fetch from NCBI** and enter `L09137` to load the pUC19 plasmid. Then click **Find restriction sites**.

## Documentation

The full user manual is in **[here](docs/documentation.md)**. It covers every view and tool, color legends, file format notes, known limitations and troubleshooting.

## Project structure

```
N-Zoom/
├── core
│   ├── __init__.py
│   ├── __pycache__
│   │   ├── __init__.cpython-314.pyc
│   │   ├── seqtools.cpython-314.pyc
│   │   └── viewers.cpython-314.pyc
│   ├── seqtools.py
│   └── viewers.py
├── docs
│   └── documentation.md
├── main.py
├── README.md
├── requirements.txt
└── sample_data
    ├── L09137.2.gb
    ├── protein.faa
    ├── rna.fna
    └── ti_plasmid.fasta
```
## Contributing

Bug reports and suggestions are welcome; please open an issue. N-Zoom is in early development, so expect things to change.

---

Created by **mikeph52** _2026 GNU Public License v3_.