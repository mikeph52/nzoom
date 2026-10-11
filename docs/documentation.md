# N-Zoom User Manual

**Version 0.7.5** · by mikeph52

N-Zoom is a desktop app for viewing and analysing DNA and protein sequences. It opens FASTA and GenBank files, and draws plasmid maps with annotated features. It also finds restriction sites, translates DNA, and downloads sequences straight from NCBI.

---

## Contents

1. [Installation](#1-installation)
2. [Starting N-Zoom](#2-starting-n-zoom)
3. [The main window](#3-the-main-window)
4. [Opening sequences](#4-opening-sequences)
5. [Fetching sequences from NCBI](#5-fetching-sequences-from-ncbi)
6. [The Map view](#6-the-map-view)
7. [The Sequence view](#7-the-sequence-view)
8. [The feature table](#8-the-feature-table)
9. [Sequence information panel](#9-sequence-information-panel)
10. [Restriction site analysis](#10-restriction-site-analysis)
11. [Translation, reverse complement and copying](#11-translation-reverse-complement-and-copying)
12. [Saving and exporting](#12-saving-and-exporting)
13. [Keyboard shortcuts](#13-keyboard-shortcuts)
14. [Color reference](#14-color-reference)
15. [File format notes](#15-file-format-notes)
16. [Known limitations](#16-known-limitations)
17. [Troubleshooting](#17-troubleshooting)

---

## 1. Installation

N-Zoom needs **Python 3.9 or newer** and two libraries:

| Library | What it's used for |
|---|---|
| Biopython | Reading and writing sequence files, translation, restriction enzymes, NCBI access |
| PySide6 | The graphical interface |

Install them from the N-Zoom folder:

```
pip install -r requirements.txt
```

The N-Zoom folder should look like this:

```
N-Zoom/
├── main.py
├── requirements.txt
└── core/
    ├── __init__.py
    ├── seqtools.py
    └── viewers.py
```

## 2. Starting N-Zoom

From inside the `N-Zoom` folder, run:

```
python main.py
```

Always start the app through `main.py`. The files inside `core/` are parts of the program and can't be run on their own.

---

## 3. The main window

```
┌──────────────────────────────────────────────────────────────────────┐
│ File  Tools  About                                                   │
│ [Open] [Fetch from NCBI]  path/to/current/file.gb                    │
├────────────┬──────────────────────────────────────┬──────────────────┤
│            │  Sequence title                      │ Sequence         │
│  Record    │  ┌ Map ┬ Sequence ┐                  │ information      │
│  list      │  │                │                  │                  │
│            │  │  (viewer)      │                  │ Actions          │
│            │  ├────────────────┤                  │  Copy sequence   │
│            │  │ Feature table  │                  │  Reverse compl.  │
│            │  └────────────────┘                  │  Translate DNA   │
│            │                                      │  Find restr. ... │
├────────────┴──────────────────────────────────────┴──────────────────┤
│ Status bar                                                           │
└──────────────────────────────────────────────────────────────────────┘
```

| Area | Purpose |
|---|---|
| **Record list** (left) | Every sequence currently loaded. Click one to display it. |
| **Viewer tabs** (centre) | **Map** shows the plasmid map; **Sequence** shows the bases or amino acids. |
| **Feature table** (centre, bottom) | Every annotated feature of the current record. |
| **Information panel** (right) | Type, length, topology, GC content, feature count and composition. |
| **Actions** (right) | One-click tools for the current record. |
| **Status bar** (bottom) | Short messages: what was loaded, what was found, which feature is selected. |

You can drag the borders between the panels to resize them.

---

## 4. Opening sequences

Use **File → Open…**, the **Open** button, or **Ctrl+O**.

Supported files:

| Format | Extensions |
|---|---|
| GenBank | `.gb` `.gbk` `.genbank` `.gbff` |
| FASTA | `.fasta` `.fa` `.faa` `.fna` `.fas` |

N-Zoom picks the format from the file extension. A file with any other extension is read as FASTA.

- **Files with several records** load all of them into the record list. The first record is shown automatically.
- **Opening a file replaces** everything currently loaded. To add records to the current list without clearing it, use [Fetch from NCBI](#5-fetching-sequences-from-ncbi).
- **FASTA vs GenBank:** FASTA files contain only names and sequences. To see features on the map, you need a GenBank file.

### Sequence type detection

When you select a record, N-Zoom decides what kind of sequence it is:

| Detected type | Rule |
|---|---|
| **DNA** | Only DNA letters and IUPAC ambiguity codes: `A C G T U N R Y K M S W B D H V -` |
| **Protein** | Only amino-acid letters, `*` (stop) and `-` |
| **Unknown** | Anything else, or an empty sequence |

RNA (with `U`) is treated as DNA. Spaces and line breaks are ignored, and lowercase letters are converted to uppercase.

A short protein made only of the letters A, C, G and T would be detected as DNA. This is rare in practice.

---

## 5. Fetching sequences from NCBI

Use **Tools → Fetch from NCBI…** or the **Fetch from NCBI** button.

1. Enter a nucleotide accession number, for example `L09137` (pUC19), `NC_001416` (phage lambda) or `J01749` (pBR322).
2. **The first time only**, you'll be asked for an email address. NCBI requires one with every request, so they can contact you if something goes wrong. N-Zoom saves it on your computer and sends it only to NCBI.
3. The record is downloaded in full GenBank format, with features. It's **added to the end of the record list** and selected.

You need an internet connection. The window doesn't respond while the download runs. Large records, such as whole bacterial genomes, can take a while.

Only nucleotide records (NCBI's *nuccore* database) can be fetched. Protein accessions won't work.

---

## 6. The Map view

The **Map** tab draws the current DNA record as a circle, with position 0 at the top and positions running clockwise.

### What's on the map

| Element | Meaning |
|---|---|
| **Black circle** | The DNA backbone |
| **Short grey ticks with numbers** | Position markers. The spacing (for example every 200, 500 or 1,000 bp) adapts to the sequence length. |
| **Colored arcs outside the circle** | Features on the **forward** strand |
| **Colored arcs inside the circle** | Features on the **reverse** strand |
| **Arrowheads** | Show which way a feature reads (its 3′ end) |
| **Labels with grey leader lines** | Feature names, connected to the arc they belong to |
| **Red ticks and red labels** | Restriction sites, after you run a restriction search |
| **Centre text** | Record name and length in bp |

### Overlapping features

Features that overlap are stacked on separate rings, so they don't cover each other. Each strand can have up to four rings. Beyond that, extra features share the outermost ring.

### Features that cross the origin

On circular plasmids, a feature can start near the end of the sequence and continue past position 0. These are drawn as one continuous arc through the top of the map.

### Labels

Labels move outward automatically to avoid overlapping each other. Feature labels are placed before restriction-site labels, so feature names win when space is tight. If a label can't find free space, it's left out.

On very busy maps (more than 300 labels in total), labels are switched off completely. All features are still listed in the feature table.

### Selecting a feature on the map

**Click an arc** to select that feature. Its row is selected in the feature table, it's highlighted in the Sequence view, and its arc is drawn thicker.

### When the map is unavailable

The Map tab is disabled for protein and "Unknown" sequences. Linear DNA and FASTA records are still drawn as a circle. The Topology field in the [information panel](#9-sequence-information-panel) tells you whether the molecule is really circular.

The **source** feature, which spans the whole sequence in every GenBank file, is not drawn on the map.

---

## 7. The Sequence view

The **Sequence** tab shows the full sequence as wrapping text in a monospace font.

- **Colors:** sequences up to **50,000** characters are colored by base or amino acid (see [Color reference](#14-color-reference)). Longer sequences are shown in plain black so they load quickly.
- **Feature highlighting:** when you select a feature, its region is highlighted in yellow and the view scrolls to it. Features with several parts, such as origin-crossing features or joined exons, highlight every part.
- The view is **read-only**. You can select and copy text with the mouse, but you can't edit it.

---

## 8. The feature table

Below the viewer is a table of every annotated feature in the current record. The `source` feature is hidden, because it just covers the entire sequence.

| Column | Meaning |
|---|---|
| **Type** | GenBank feature type, for example `CDS`, `promoter`, `rep_origin`, `primer_bind` |
| **Name** | The feature's label (see below) |
| **Start** | First base, 1-based, as written in GenBank files |
| **End** | Last base, 1-based and inclusive |
| **Strand** | `→` forward, `←` reverse, `·` unknown or mixed |
| **Length** | Total length in bp. For multi-part features this is the sum of the parts. |

### How names are chosen

N-Zoom uses the first of these GenBank qualifiers that the feature has: `/label`, `/gene`, `/product`, `/locus_tag`, `/note`. If it has none of them, the feature type is shown instead.

### Origin-crossing features

For a feature that crosses position 0, **Start is larger than End**. For example, Start `3901` and End `60` on a 4,000 bp plasmid means the feature runs from 3901 to the end, then from 1 to 60.

### Using the table

- **Click a row** to select the feature on the map and highlight it in the Sequence view. The status bar shows the feature's type, name and exact location.
- **Click a column header** to sort by that column. Click it again to reverse the order.

The status bar shows the location in Biopython's internal notation, which counts from **0** and excludes the end. For example, `[2800:3000](+)` is bases 2801–3000 in normal GenBank numbering. The table always uses normal numbering.

---

## 9. Sequence information panel

| Field | Meaning |
|---|---|
| **Type** | DNA, Protein or Unknown (see [type detection](#sequence-type-detection)) |
| **Length** | Number of bases or residues |
| **Topology** | `circular` or `linear` as recorded in a GenBank file; `unknown` for FASTA |
| **GC %** | (G + C) ÷ total length × 100. DNA only. Ambiguous letters such as `N` count toward the total but not toward GC. |
| **Features** | Number of features listed in the feature table |
| **Composition** | Count of the 12 most common letters, most frequent first |

---

## 10. Restriction site analysis

Use **Tools → Find restriction sites**, or the **Find restriction sites** button. This works on DNA records only.

### What's searched

- All **commercially available** restriction enzymes in Biopython's copy of the REBASE enzyme database.
- The search honours the record's **topology**. Circular GenBank records are searched as circular, so sites spanning the origin are found. FASTA records and linear GenBank records are searched as linear.

### Results

**1. The Restriction sites window** lists every enzyme that cuts at least once:

| Column | Meaning |
|---|---|
| **Enzyme** | Enzyme name |
| **Cuts** | How many times it cuts |
| **Positions** | Cut positions (the first 12 are shown) |

The list starts with the enzymes that cut least often. Click a column header to re-sort it.

A **cut position** is the 1-based position of the first base *after* the cut on the top strand.

**2. The map** marks only the most useful sites for cloning: enzymes that **cut exactly once** and recognise **6 or more bases**. Enzymes that recognise the same site and cut at the same position (isoschizomers) are merged into one label, such as `Acc65I, KpnI`.

### Clearing sites

Use **Tools → Clear restriction sites**. Sites are also cleared automatically when you select a different record.

### Large sequences

For sequences over 2,000,000 bp, N-Zoom asks for confirmation first, because the search can take a long time.

---

## 11. Translation, reverse complement and copying

All three tools put their result on the **clipboard**, ready to paste elsewhere. None of them change the loaded sequence.

### Copy sequence

This copies the current sequence exactly as displayed: uppercase, with no spaces or line breaks.

### Reverse complement (DNA only)

This copies the reverse complement of the whole sequence. Ambiguity codes are complemented correctly, for example `R` ↔ `Y`.

### Translate DNA (DNA only)

This copies the protein translation of the whole sequence and shows a preview of up to 500 amino acids.

- Translation starts at the **first base, in frame 1**, using the **standard genetic code**.
- Gaps (`-`) are removed and `U` is read as `T` before translating.
- Leftover bases at the end that don't make a full codon are ignored.
- Stop codons appear as `*`. Translation continues past them.
- Codons with ambiguous bases translate as `X` when the amino acid can't be determined.

To translate a single gene, copy that region first and open it as its own sequence. Translating a selected feature directly isn't supported yet.

---

## 12. Saving and exporting

Use **File → Save As…** or **Ctrl+S**.

- **Everything currently in the record list is saved**, not only the selected record.
- Choose **GenBank** to keep features and annotations, or **FASTA** for sequences only.
- If you don't type an extension, `.gb` or `.fasta` is added for you.

Saving as FASTA **drops all features**. Use GenBank if you want to keep them.

If you save FASTA-loaded records as GenBank, N-Zoom fills in the molecule type (DNA or protein) automatically, because the GenBank format requires it.

---

## 13. Keyboard shortcuts

| Shortcut | Action |
|---|---|
| **Ctrl+O** | Open a file |
| **Ctrl+S** | Save As |
| **Ctrl+Q** | Exit N-Zoom |
| **Ctrl+C** | Copy selected text in the Sequence view |

On macOS, use **Cmd** instead of Ctrl.

---

## 14. Color reference

### DNA bases (Sequence view)

| Base | Color |
|---|---|
| A | Green |
| T / U | Red |
| G | Amber |
| C | Blue |
| N | Grey |
| Other ambiguity codes | Dark grey |

### Amino acids (Sequence view)

| Group | Residues | Color |
|---|---|---|
| Hydrophobic | A V L I M | Green |
| Aromatic | F W Y | Purple |
| Positive | K R H | Blue |
| Negative | D E | Red |
| Polar | S T N Q | Teal |
| Special | G P | Brown |
| Cysteine | C | Orange |
| Stop | * | Black |

### Feature types (Map)

| Feature type | Color |
|---|---|
| `CDS` | Blue |
| `gene` | Green |
| `mRNA` | Indigo |
| `promoter` | Orange |
| `terminator` | Red |
| `rep_origin` | Yellow |
| `primer_bind` | Teal |
| `protein_bind` | Brown |
| `RBS` | Purple |
| `regulatory` | Deep orange |
| `misc_feature` | Grey |
| any other type | Blue-grey |

Restriction sites are always dark red.

---

## 15. File format notes

### Writing your own GenBank features

Features in a GenBank file look like this:

```
FEATURES             Location/Qualifiers
     CDS             701..1700
                     /label="GFP"
     CDS             complement(2106..2966)
                     /label="AmpR"
     misc_feature    join(3901..4000,1..60)
                     /label="origin-spanning tag"
```

| Location | Meaning |
|---|---|
| `701..1700` | Bases 701 to 1700, forward strand |
| `complement(2106..2966)` | Same kind of range, on the reverse strand |
| `join(3901..4000,1..60)` | A feature made of several parts; used for origin-crossing features and exons |

GenBank is a column-based format: the feature type starts at column 6, and the location and qualifiers at column 22. Hand-edited files that break this layout may fail to open. It's safest to create or edit GenBank files with N-Zoom's Save As or another sequence tool.

### Circular vs linear

Topology is read from the `LOCUS` line of a GenBank file, for example `LOCUS pUC19 2686 bp DNA circular SYN`. FASTA files have no way to store topology.

---

## 16. Known limitations

These apply to version 0.5.0:

- **No editing.** Sequences and features can't be changed inside N-Zoom yet.
- **The window freezes** during NCBI downloads and long restriction searches.
- **The map is always circular.** Linear molecules are also drawn as a circle; there's no linear map view yet.
- **Translation uses frame 1 of the whole sequence** only. There's no six-frame or per-feature translation.
- **Saving GenBank can fail for long record names.** The GenBank format limits the name on the LOCUS line to 16 characters, and N-Zoom shows an error if a name is too long.
- **Very large sequences** (whole genomes) open and display, but the map becomes crowded and labels are hidden.

---

## 17. Troubleshooting

| Problem | Solution |
|---|---|
| `ModuleNotFoundError: No module named 'Bio'` or `'PySide6'` | Run `pip install -r requirements.txt` |
| `ModuleNotFoundError: No module named 'core'` | Start the app from the `N-Zoom` folder with `python main.py` |
| `ImportError: attempted relative import` | You ran a file inside `core/` directly. Run `main.py` instead. |
| The map is empty / "No circular DNA loaded" | The record isn't DNA, or it has no sequence data |
| The map shows a circle but no features | The file is FASTA, or the GenBank file has no features apart from `source` |
| "Could not read file" | The file is damaged, or its extension doesn't match its contents (for example a GenBank file named `.fasta`). Rename it to the correct extension. |
| NCBI fetch fails | Check your internet connection and the accession number. Make sure it's a nucleotide record. NCBI may also be temporarily busy; try again in a minute. |
| A wrong NCBI email was saved | Delete the saved N-Zoom settings: on Windows the registry key `HKEY_CURRENT_USER\Software\mikeph52\N-Zoom`; on Linux `~/.config/mikeph52/N-Zoom.conf`; on macOS `~/Library/Preferences/com.mikeph52.N-Zoom.plist` |
| "Could not save file" with GenBank | A record name is longer than 16 characters. Save as FASTA, or rename the record. |
| Sequence appears in black instead of colors | It's longer than 50,000 characters. This is intentional, to keep the app fast. |