# N-Zoom User Manual

**Version 0.10.1** · by mikeph52

N-Zoom is a desktop app for viewing and analysing DNA and protein sequences. It opens FASTA and GenBank files, and draws plasmid maps with annotated features. You can add and edit your own features, find open reading frames, restriction sites and primer binding sites, run virtual digests on a simulated gel, translate any region, and download sequences straight from NCBI. It has light and dark modes.

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
11. [Right-click menus](#11-right-click-menus)
12. [Translation, reverse complement and copying](#12-translation-reverse-complement-and-copying)
13. [Finding open reading frames (ORFs)](#13-finding-open-reading-frames-orfs)
14. [Adding, editing and deleting features](#14-adding-editing-and-deleting-features)
15. [Finding primers](#15-finding-primers)
16. [Virtual digest and gel](#16-virtual-digest-and-gel)
17. [Light and dark mode](#17-light-and-dark-mode)
18. [Saving and exporting](#18-saving-and-exporting)
19. [Keyboard shortcuts](#19-keyboard-shortcuts)
20. [Color reference](#20-color-reference)
21. [File format notes](#21-file-format-notes)
22. [Known limitations](#22-known-limitations)
23. [Troubleshooting](#23-troubleshooting)

---

## 1. Installation

N-Zoom needs **Python 3.10 or newer** and two libraries:

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
├── core/
│   ├── __init__.py
│   ├── dialogs.py
│   ├── digest.py
│   ├── maps.py
│   ├── primers.py
│   ├── seqtools.py
│   ├── theme.py
│   ├── viewers.py
│   └── workers.py
├── docs/
│   └── documentation.md
├── examples/
│   └── pDEMO.gb
└── tests/
```

`examples/pDEMO.gb` is a 4,000 bp demo plasmid for trying every tool. It has annotated features, two genes that translate cleanly (GFP and AmpR), and a multiple cloning site whose ten enzymes each cut only once. The genes are made of random codons, so they have the right size and shape but aren't real proteins.

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
│ File  Edit  View  Tools  About                                       │
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
│            │                                      │  Virtual digest  │
│            │                                      │  Find ORFs       │
│            │                                      │  Find primer     │
│            │                                      │  Add feature     │
├────────────┴──────────────────────────────────────┴──────────────────┤
│ Status bar                                                           │
└──────────────────────────────────────────────────────────────────────┘
```

| Area | Purpose |
|---|---|
| **Record list** (left) | Every sequence currently loaded. Click one to display it. |
| **Viewer tabs** (centre) | **Map** shows a circular or linear map; **Sequence** shows the bases or amino acids. |
| **Feature table** (centre, bottom) | Every annotated feature of the current record, plus any ORFs found. |
| **Information panel** (right) | Type, length, topology, GC content, feature count and composition. |
| **Actions** (right) | One-click tools for the current record. |
| **Status bar** (bottom) | Messages on the left (what was loaded or found, which feature is selected). On the right: a moving bar while a background task runs, and the cursor position or current selection. |

You can drag the borders between the panels to resize them.

The **window title** shows the file name. A `*` after it means there are unsaved changes.

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

You need an internet connection. The download runs in the background: a moving bar appears in the status bar, and you can keep working in the meantime. The Fetch button and menu item are disabled until it finishes. Large records, such as whole bacterial genomes, can take a while.

If the download fails, for example because the accession doesn't exist or there's no connection, an error message explains why.

Only nucleotide records (NCBI's *nuccore* database) can be fetched. Protein accessions won't work.

---

## 6. The Map view

The **Map** tab draws the current DNA record in one of two styles:

- **Circular map**, for circular DNA such as plasmids: position 0 is at the top and positions run clockwise.
- **Linear map**, for linear DNA such as PCR products, genes and FASTA sequences: a horizontal bar running from base 1 on the left to the last base on the right.

### Choosing the map style

By default N-Zoom picks the style from the record's [topology](#21-file-format-notes): circular records get the circular map, and linear, FASTA and unknown records get the linear map. To override this, use **View → Map style**:

| Option | Shows |
|---|---|
| **Automatic (from topology)** | Circular map for circular records, linear map for everything else (default) |
| **Circular** | Always the circular map |
| **Linear** | Always the linear map. Useful for looking at a plasmid region by region. |

The choice is remembered for next time. Changing the map style doesn't change the record: a circular plasmid shown as a line is still treated as circular by the ORF finder, the restriction search and so on.

### What's on the map

| Element | Circular map | Linear map |
|---|---|---|
| **The DNA** | Dark circle (backbone) | Horizontal bar with a short end mark at each end |
| **Position markers** | Ticks with numbers inside the circle | Ticks with numbers below the bar |
| **Forward-strand features** | Colored arcs outside the circle | Colored bars above the line |
| **Reverse-strand features** | Colored arcs inside the circle | Colored bars below the line |
| **Direction** | Arrowheads at each feature's 3′ end | Pointed end at each feature's 3′ end |
| **Labels** | Around the outside, with grey leader lines | Above (forward features, enzymes) and below (reverse features), with leader lines |
| **ORFs** | Thin light-green arcs | Thin light-green bars |
| **Restriction sites** | Red ticks and red labels | Red ticks across the bar, red labels above |
| **Selected region** | Yellow band on the circle | Yellow band over the bar |
| **Name and length** | In the centre | Top left |

The tick spacing (for example every 200, 500 or 1,000 bp) adapts to the sequence length. The linear map fits the whole sequence into the window's width; there's no zoom yet.

### Overlapping features

Features that overlap are stacked on separate rings (circular) or rows (linear), so they don't cover each other. Each strand can have up to four. Beyond that, extra features share the outermost one.

### Features that cross the origin

On circular plasmids, a feature can start near the end of the sequence and continue past position 0. The circular map draws it as one continuous arc through the top. The linear map draws it as two pieces, one at each end of the bar.

### Labels

Labels move outward automatically to avoid overlapping each other. Feature labels are placed before restriction-site labels, so feature names win when space is tight. If a label can't find free space, it's left out.

On very busy maps (more than 300 labels in total), labels are switched off completely. All features are still listed in the feature table.

### Selecting a feature on the map

**Click a feature** (an arc or a bar) to select it. Its row is selected in the feature table, it's highlighted in the Sequence view, and it's drawn thicker (circular) or outlined (linear).

**Right-click a feature** for its [menu](#11-right-click-menus): translate it, copy it, edit it, and so on.

### When the map is unavailable

The Map tab is disabled for protein and "Unknown" sequences.

The **source** feature, which spans the whole sequence in every GenBank file, is not drawn on the map.

---

## 7. The Sequence view

The **Sequence** tab shows the full sequence as wrapping text in a monospace font.

- **Colors:** sequences up to **50,000** characters are colored by base or amino acid (see [Color reference](#20-color-reference)). Longer sequences are shown in plain black so they load quickly.
- **Feature highlighting:** when you select a feature, its region is highlighted in yellow and the view scrolls to it. Features with several parts, such as origin-crossing features or joined exons, highlight every part.
- The view is **read-only**. You can select and copy text, but you can't edit it.

### Selecting a region

Select bases with the mouse, or with **Shift + arrow keys**. **Ctrl+A** selects everything.

The right side of the status bar shows:

| What you do | Status bar shows |
|---|---|
| Click somewhere (no selection) | `Position 42 of 4,000`: the base just after the cursor |
| Select a region of DNA | `Selected 601–657 (57 bp) · GC 57.9%` |
| Select a region of protein | `Selected 10–45 (36 aa)` |

Positions are 1-based and include both ends, matching the feature table and GenBank files.

For DNA, the selected region is also shaded **yellow on the Map**, so you can see where it lies on the plasmid.

**Ctrl+C** copies the selected bases. **Right-click** the selection for more: copy the reverse complement, translate it, or turn it into a feature (see [Right-click menus](#11-right-click-menus)).

---

## 8. The feature table

Below the viewer is a table of every annotated feature in the current record. The `source` feature is hidden, because it just covers the entire sequence.

After you run the [ORF finder](#13-finding-open-reading-frames-orfs), ORFs are listed here too, in *grey italics* with type `ORF`, so you can tell them apart from real annotations.

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
- **Double-click a row** to edit that feature. For an ORF, double-clicking turns it into a CDS feature.
- **Right-click a row** for the feature's [menu](#11-right-click-menus).
- Press **Delete** to delete the selected feature (after confirming).

The status bar shows the location in Biopython's internal notation, which counts from **0** and excludes the end. For example, `[2800:3000](+)` is bases 2801–3000 in normal GenBank numbering. The table always uses normal numbering.

---

## 9. Sequence information panel

| Field | Meaning |
|---|---|
| **Type** | DNA, Protein or Unknown (see [type detection](#sequence-type-detection)) |
| **Length** | Number of bases or residues |
| **Topology** | `circular` or `linear` as recorded in a GenBank file; `unknown` for FASTA |
| **GC %** | (G + C) ÷ total length × 100. DNA only. Ambiguous letters such as `N` count toward the total but not toward GC. |
| **Features** | Number of annotated features, plus the number of ORFs found, if any. For example `10 (+ 2 ORFs)`. |
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

### Running in the background

The search runs in the background, with a moving bar in the status bar. The button and menu item are disabled until it finishes. The results list opens in its own window, which can stay open while you keep using N-Zoom.

If you switch to another record before the search finishes, the results are discarded rather than drawn on the wrong map. The status bar tells you to run the search again.

### Large sequences

For sequences over 2,000,000 bp, N-Zoom asks for confirmation first, because the search can take a long time.

---

## 11. Right-click menus

Most tools that act on *part* of a sequence are in right-click menus.

### On a selection in the Sequence view

Select some bases, then right-click the selection:

| Menu item | What it does |
|---|---|
| **Copy** | Copies the selected bases |
| **Copy reverse complement** | Copies the reverse complement of the selection (DNA only) |
| **Translate selection** | Translates the selection, reading forward from its first base (DNA only) |
| **Translate reverse complement** | Translates the other strand of the selection (DNA only) |
| **Create feature from selection…** | Opens the [feature dialog](#14-adding-editing-and-deleting-features) with Start and End already filled in |
| **Select all** | Selects the whole sequence |

Right-clicking with nothing selected offers **Add feature…** and **Select all**.

### On a feature or ORF

Right-click a feature on the map (an arc or a bar), or its row in the feature table. The first, greyed-out line of the menu shows which feature you clicked.

| Menu item | What it does |
|---|---|
| **Select in sequence** | Switches to the Sequence tab and selects the feature's bases, so the status bar shows its length and GC%. A feature made of several pieces is highlighted instead. |
| **Copy sequence** | Copies the feature's sequence **as read on its own strand**: reverse-strand features are reverse-complemented, and pieces of joined features are put together in order. |
| **Translate** | Translates the feature (see [below](#translating-a-feature)). DNA only. |
| **Edit feature…** | Opens the [feature dialog](#14-adding-editing-and-deleting-features) |
| **Delete feature** | Deletes the feature, after asking you to confirm |
| **Add as CDS feature…** | ORFs only: turns the ORF into a real CDS feature |

---

## 12. Translation, reverse complement and copying

None of these tools change the loaded sequence.

### Copy sequence and Reverse complement (Actions panel)

**Copy sequence** copies the whole current sequence exactly as displayed: uppercase, with no spaces or line breaks.

**Reverse complement** copies the reverse complement of the whole sequence (DNA only). Ambiguity codes are complemented correctly, for example `R` ↔ `Y`.

To copy only part of a sequence, use the [right-click menus](#11-right-click-menus).

### The Translation window

Translations open in a **Translation** window showing the protein, with a **Copy** button. You can keep several of these windows open and carry on working.

The line above the protein summarises what was translated, for example:

`GFP (CDS, 701–1,417, forward strand) · 717 bp → 238 aa + stop`

It also warns you when something looks wrong:

| Note | Meaning |
|---|---|
| `last 2 base(s) ignored` | The length isn't a multiple of 3, so the incomplete codon at the end was skipped |
| `3 internal stop codon(s)` | There are stop codons (`*`) before the end: the reading frame or coordinates may be wrong |
| `genetic code 11` | A genetic code other than the standard one was used (see below) |
| `matches / differs from the /translation stored in the file` | GenBank CDS features from NCBI often store the expected protein; N-Zoom compares its result with it |

### What can be translated

| How | What is translated |
|---|---|
| **Translate DNA** button, or **Tools → Translate whole sequence** | The whole sequence, from base 1 |
| Right-click a selection → **Translate selection** | The selected bases, from the first selected base |
| Right-click a selection → **Translate reverse complement** | The reverse complement of the selection |
| Right-click a feature → **Translate** | The feature, as described below |

### Translating a feature

Feature translation follows the GenBank annotations:

- **Strand:** reverse-strand features are read on the reverse strand.
- **Pieces:** features made of several pieces (exons, or features crossing the origin) are joined before translating.
- **`/codon_start`:** if the feature has this qualifier (1, 2 or 3), translation starts at that base.
- **`/transl_table`:** if present, that genetic code is used, for example `11` for bacteria. Otherwise the standard code (1) is used.

### Translation rules

- Gaps (`-`) are removed and `U` is read as `T` before translating.
- Stop codons appear as `*`. Translation continues past them.
- Codons with ambiguous bases translate as `X` when the amino acid can't be determined.

---

## 13. Finding open reading frames (ORFs)

An **open reading frame** is a stretch of DNA that could code for a protein: it starts with a start codon (`ATG`) and runs, without interruption, to a stop codon (`TAA`, `TAG` or `TGA`).

Use **Tools → Find ORFs…** or the **Find ORFs** button (DNA only).

1. Enter the **minimum length** in amino acids. The default is 100. Shorter limits find many more ORFs, most of them occurring by chance. N-Zoom remembers your choice for next time.
2. The search runs in the background, like the restriction search.
3. ORFs appear on the **map** as thin light-green arcs, and in the **feature table** in grey italics, named `ORF1 (238 aa)`, `ORF2 (…)` and so on, in order of position.

### How ORFs are found

- All **six reading frames** are searched: three on each strand.
- Each ORF runs from the **first ATG after a stop codon** to the **next stop codon**. ATGs further inside an ORF don't produce extra, shorter ORFs.
- The stop codon **is included** in the ORF's coordinates, as in GenBank CDS features. The length in the name doesn't count it.
- On **circular** records, ORFs can run across the origin. On linear and FASTA records they can't.
- The **standard genetic code** is used.

### ORFs that are already annotated

If an ORF has exactly the same location and strand as a feature that's already in the file, it isn't shown again. The status bar tells you how many were hidden, for example `(2 already annotated, not shown.)`.

### Keeping an ORF

ORFs are **suggestions**: they're not part of your file and aren't saved. To keep one, right-click it and choose **Add as CDS feature…**, or double-click its row. The feature dialog opens with type `CDS`, the ORF's name, coordinates and strand filled in. Change the name, then click OK. The ORF is replaced by the new feature.

### Clearing ORFs

Use **Tools → Clear ORFs**. ORFs are also cleared when you select a different record.

---

## 14. Adding, editing and deleting features

All changes are made to the records in memory. Use [Save](#18-saving-and-exporting) to write them to a file.

### Adding a feature

Any of these opens the **Add feature** dialog:

- **Edit → Add feature…** (Ctrl+E) or the **Add feature** button. If bases are selected in the Sequence view, Start and End are filled in from the selection.
- Right-click a selection → **Create feature from selection…**
- Right-click an ORF → **Add as CDS feature…**

### The feature dialog

| Field | Notes |
|---|---|
| **Name** | Saved as the `/label` qualifier, the name shown everywhere in N-Zoom |
| **Type** | Pick a common GenBank feature type from the list, or type any other |
| **Start** / **End** | 1-based and inclusive, as in the feature table and GenBank files |
| **Strand** | Forward, Reverse, or None |
| **Note** | Optional; saved as the `/note` qualifier |

Below the fields, the dialog shows the feature's length.

**Features crossing the origin:** on a **circular** record, set Start larger than End. For example, Start 3,901 and End 60 on a 4,000 bp plasmid. The dialog shows `crosses the origin`. On a linear record this isn't possible: the dialog shows a red warning and OK is disabled.

### Editing a feature

**Double-click** the feature's row, right-click it → **Edit feature…**, or select it and use **Edit → Edit selected feature…**.

Features made of **several separate pieces**, such as genes with exons, can't be described by a single Start and End. For these, only the name, type and note can be changed; the location fields are disabled. Origin-crossing features count as a single piece and can be edited normally.

If a feature's name came from `/gene` or `/product` rather than `/label`, editing the name adds a `/label`. The original qualifiers are kept.

### Deleting a feature

Right-click it → **Delete feature**, select it and press **Delete**, or use **Edit → Delete selected feature**. N-Zoom asks you to confirm first.

### Undo

There's no undo yet. To throw away all changes since your last save, open the file again and choose **Discard** when asked about unsaved changes.

---

## 15. Finding primers

Use **Tools → Find primer…** (Ctrl+P) or the **Find primer** button (DNA only).

1. Type a **name**, for example `GFP-fwd`.
2. Paste the **primer sequence**, written 5′ → 3′ as you would order it. Spaces and line breaks are ignored, and lowercase is fine. Only A, C, G and T are allowed.
3. Results appear as you type: the primer's length, GC content and Tm, and a table of every **binding site**.
4. Select the sites you want and click **Add selected sites as features**. Each site is added as a `primer_bind` feature, with an arrow showing the primer's direction. Its note records the primer sequence and Tm.

### How binding sites are found

A primer binds where its **3′ end** matches the template exactly. The 3′ end is what the polymerase extends, so it must pair perfectly.

- By default the last **15 bases** at the 3′ end must match (the **Must match** setting; you can change it from 8 to 60). Shorter primers must match completely.
- From there, the match is extended toward the 5′ end for as long as bases keep matching. Bases beyond that are a **5′ tail**, such as an added restriction site or tag, which doesn't anneal in the first PCR cycle.
- **Both strands** are searched. A **forward** (→) primer has the same sequence as the top strand; a **reverse** (←) primer matches the bottom strand.
- On circular records, primers can bind across the origin.

| Column | Meaning |
|---|---|
| **Strand** | → forward or ← reverse |
| **Start / End** | The annealed part, 1-based, in top-strand numbering |
| **Annealed** | Number of primer bases that pair with the template |
| **5′ tail** | Primer bases that don't pair |
| **Tm annealed** | Melting temperature of the annealed part only. This is the Tm that matters for the first PCR cycles. |

If a primer binds **more than once**, N-Zoom warns you: PCR may give extra products.

### About the Tm

Tm is calculated with Biopython's nearest-neighbour method (`Tm_NN`), assuming **50 mM Na⁺** and **250 nM primer**, without Mg²⁺. Different calculators and buffers give values a few degrees apart, so use it to compare primers rather than as an exact annealing temperature. Your polymerase supplier's calculator is the best guide for the final value.

### Limitations

The search looks for exact 3′ matches only. Primers with a mismatch near the 3′ end, which may still bind weakly, aren't reported. There's no check for hairpins or primer dimers, and no primer design yet.

---

## 16. Virtual digest and gel

Use **Tools → Virtual digest…** (Ctrl+D) or the **Virtual digest** button (DNA only). N-Zoom first finds every enzyme that cuts the sequence, which runs in the background, then opens the digest window. If you've just run **Find restriction sites** on the same record, it opens straight away.

### Making a digest

1. **Choose enzymes** in the list on the left. Each shows how many times it cuts, for example `EcoRI (1 cut)`.
   - Type in **Filter enzymes** to find one by name.
   - The drop-down limits the list to enzymes that cut once, 1–2 times, 1–3 times, or any number of times.
2. **Check** one enzyme for a single digest, or several for a double or triple digest.
3. Click **Add lane**. The digest appears as a new lane on the gel, and the checks are cleared, ready for the next lane.

You can add up to 10 lanes. To remove one, select it in the **Lanes** list and click **Remove lane**.

### The gel

The gel shows a **ladder** (lane **L**) and one lane per digest, numbered as in the Lanes list. Choose the ladder with the **Ladder** drop-down:

| Ladder | Bands (bp) |
|---|---|
| **1 kb ladder** | 10,000 · 8,000 · 6,000 · 5,000 · 4,000 · 3,500 · 3,000 · 2,500 · 2,000 · 1,500 · 1,000 · 750 · 500 · 250 |
| **100 bp ladder** | 1,500 · 1,200 · 1,000 · 900 · 800 · 700 · 600 · 500 · 400 · 300 · 200 · 100 |

These are generic ladders; commercial ones differ slightly.

- **Band position:** small fragments run further down. Distance is proportional to the logarithm of the size, the usual approximation for agarose gels.
- **Band brightness** follows the amount of DNA. With equal numbers of molecules, longer fragments carry more DNA, so they look brighter, as with ethidium or similar stains. Fragments of the same size add up into one brighter band.
- Fragments of very similar size **merge** into one band, as they would on a real gel.

**Click a lane** on the gel to select it. **Save gel image…** saves the gel as a PNG picture.

### The fragment table

For the selected lane, the table lists each fragment:

| Column | Meaning |
|---|---|
| **#** | Fragment number, in sequence order |
| **Size (bp)** | Fragment length |
| **Start / End** | First and last base, 1-based. On circular DNA, a fragment that crosses the origin has Start larger than End. |
| **Left end / Right end** | The enzyme that made each end, or "sequence start" / "sequence end" for the ends of linear DNA |

### What the simulation assumes

- **Sizes** are measured on the top strand; sticky-end overhangs of a few bases are ignored.
- **Circular DNA** cut once gives a single linear fragment of full length. Uncut circular DNA is drawn at its full size and marked **uncut circle**. On a real gel, uncut plasmids run as several bands (supercoiled, nicked) at positions that don't match their size.
- **Topology** follows the record: circular records are digested as circles, everything else as linear DNA.
- Real gels also depend on agarose percentage, voltage and run time, and very small fragments may run off or be too faint to see.

---

## 17. Light and dark mode

Choose **View → Theme**:

| Option | What it does |
|---|---|
| **Same as system** | Follows your operating system's light or dark setting, and switches automatically when it changes (default) |
| **Light** | Always light |
| **Dark** | Always dark |

The choice is remembered for next time. The whole window switches, including the maps and the sequence colors. Base and amino-acid colors use lighter shades in dark mode so they stay readable.

The **gel** in the virtual digest always looks the same (light bands on a dark gel), as real gel photos do.

N-Zoom uses the same widget style (Qt's "Fusion") on every platform, because some native styles ignore dark colors.

---

## 18. Saving and exporting

| Command | What it does |
|---|---|
| **File → Save** (Ctrl+S) | Saves to the file you opened, in the same format |
| **File → Save As…** (Ctrl+Shift+S) | Saves to a new file. Choose **GenBank** or **FASTA**. |

- **Everything currently in the record list is saved**, not only the selected record. This includes records fetched from NCBI.
- If you don't type an extension, `.gb` or `.fasta` is added for you.

### Features and FASTA files

FASTA files can't store features. So:

- **Save** on a file that was opened as FASTA but now has features asks you where to save a **GenBank** file instead, rather than silently dropping the features.
- **Save As** with FASTA chosen warns you that all annotations will be left out, and asks you to confirm.

### Unsaved changes

The window title shows `*` when there are unsaved changes: added, edited or deleted features, or records fetched from NCBI. If you open another file or close N-Zoom with unsaved changes, you're asked whether to **Save**, **Discard** them, or **Cancel**.

### Long names in GenBank files

The GenBank format allows at most **16 characters, without spaces**, for a sequence name. Names from FASTA files are often longer. N-Zoom shortens them **in the saved file only** and tells you which names were changed. The records in N-Zoom keep their full names.

If you save FASTA-loaded records as GenBank, N-Zoom also fills in the molecule type (DNA or protein) automatically, because the GenBank format requires it.

---

## 19. Keyboard shortcuts

| Shortcut | Action |
|---|---|
| **Ctrl+O** | Open a file |
| **Ctrl+S** | Save |
| **Ctrl+Shift+S** | Save As |
| **Ctrl+E** | Add feature |
| **Ctrl+D** | Virtual digest |
| **Ctrl+P** | Find primer |
| **Delete** | Delete the selected feature (when the feature table has focus) |
| **Ctrl+Q** | Exit N-Zoom |
| **Ctrl+A** | Select the whole sequence (in the Sequence view) |
| **Ctrl+C** | Copy selected text in the Sequence view |
| **Shift + arrow keys** | Extend the selection in the Sequence view |

On macOS, use **Cmd** instead of Ctrl.

---

## 20. Color reference

The colors below are the same in light and dark mode; dark mode uses lighter shades so they stay readable.

### DNA bases (Sequence view)

| Base | Color |
|---|---|
| A | Green |
| T / U | Red |
| G | Amber |
| C | Blue |
| N | Grey |
| Other ambiguity codes | Dark grey (light grey in dark mode) |

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
| `primer_bind` (including primers you add) | Teal |
| `protein_bind` | Brown |
| `RBS` | Purple |
| `regulatory` | Deep orange |
| `misc_feature` | Grey |
| ORFs (from the ORF finder) | Light green, thinner arcs |
| any other type | Blue-grey |

Restriction sites are always dark red. The selected region is shaded yellow.

---

## 21. File format notes

### How features look in a GenBank file

Features in a GenBank file look like this:

```
FEATURES             Location/Qualifiers
     CDS             701..1417
                     /label="GFP"
     CDS             complement(2106..2966)
                     /label="AmpR"
     misc_feature    join(3901..4000,1..60)
                     /label="origin-spanning tag"
```

| Location | Meaning |
|---|---|
| `701..1417` | Bases 701 to 1417, forward strand |
| `complement(2106..2966)` | Same kind of range, on the reverse strand |
| `join(3901..4000,1..60)` | A feature made of several parts; used for origin-crossing features and exons |

GenBank is a column-based format: the feature type starts at column 6, and the location and qualifiers at column 22. Hand-edited files that break this layout may fail to open. The easiest way to add features is with N-Zoom's [feature dialog](#14-adding-editing-and-deleting-features).

### Circular vs linear

Topology is read from the `LOCUS` line of a GenBank file, for example `LOCUS pUC19 2686 bp DNA circular SYN`. FASTA files have no way to store topology, so FASTA records are treated as linear: they're shown on the linear map, and features, ORFs, primers and digest fragments can't cross the origin.

---

## 22. Known limitations

These apply to version 0.8.0:

- **The sequence itself can't be edited.** Features can be added, edited and deleted, but bases can't be inserted, removed or changed.
- **No undo** for feature changes. Discarding unsaved changes is the only way back.
- **Background tasks can't be cancelled.** If you close N-Zoom while a download or search is running, it waits for the task to finish before closing.
- **No zoom on the maps.** The whole sequence always fits the window, so dense regions of large sequences are hard to read on the map; use the Sequence view or the feature table.
- **Primer search finds exact 3′ matches only.** There's no primer design, and no check for hairpins or primer dimers.
- **The gel is an approximation.** Band positions use a simple log model, and uncut plasmids are not drawn the way they really run.
- **The ORF finder only uses ATG** as a start codon, and the standard genetic code.
- **Topology can't be changed.** A FASTA record can't be marked as circular.
- **Very large sequences** (whole genomes) open and display, but the map becomes crowded and labels are hidden.

---

## 23. Troubleshooting

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
| The ORF finder finds nothing | Lower the minimum length. Also, ORFs that are already annotated aren't shown again; check the status bar message. |
| OK is greyed out in the feature dialog | Start is after End on a linear record. Swap them, or see [Circular vs linear](#circular-vs-linear). |
| A translation has many `*` in the middle | The coordinates, strand or `/codon_start` of the feature are probably wrong |
| Sequence appears in black instead of colors | It's longer than 50,000 characters. This is intentional, to keep the app fast. |
| A plasmid is shown as a straight line | The record isn't marked circular (for example, it came from a FASTA file), or **View → Map style** is set to Linear |
| A primer you expect isn't found | Check it's written 5′ → 3′, and that its 3′ end matches the template exactly. Lower **Must match** to allow shorter annealing. |
| The enzyme you want isn't in the digest list | It doesn't cut this sequence, or the drop-down hides enzymes that cut more often. Choose **All enzymes that cut**. |
| Colors look wrong after switching theme | Choose the theme again in **View → Theme**; if that doesn't help, restart N-Zoom |
