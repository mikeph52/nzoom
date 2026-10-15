# N-Zoom by mikeph52, Init 12/9/2026
"""
main.py — the N-Zoom main window. Run with:  python main.py
"""
import sys
from collections import Counter
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QThreadPool
from PySide6.QtGui import (
    QAction, QActionGroup, QFont, QGuiApplication, QKeySequence, QShortcut, QUndoStack,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QMainWindow,
    QMenu,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from core import seqtools as st
from core import theme
from core.dialogs import FeatureDialog, RestrictionSitesDialog, TextResultDialog, number_item
from core.digest import DigestDialog
from core.maps import CircularDNAViewer, LinearDNAViewer, MapGroup, MapPanel, export_image
from core.primers import PrimerDialog
from core.reads import ReadsDialog
from core.undo import AddFeatures, DeleteFeature, EditFeature, snapshot
from core.viewers import SequenceViewer
from core.workers import Task

VERSION = "0.13.1"
MAX_FIND_HITS = 5000  # Find stops counting after this many matches
MAP_MODES = [("auto", "Automatic (from topology)"), ("circular", "Circular"), ("linear", "Linear")]
THEME_NAMES = {"system": "Same as system", "light": "Light", "dark": "Dark"}

OPEN_FILTER = (
    "Sequence files (*.gb *.gbk *.genbank *.gbff *.fasta *.fa *.faa *.fna *.fas);;"
    "GenBank files (*.gb *.gbk *.genbank *.gbff);;"
    "FASTA files (*.fasta *.fa *.faa *.fna *.fas);;"
    "All files (*)"
)
GENBANK_FILTER = "GenBank (*.gb)"
FASTA_FILTER = "FASTA (*.fasta)"
STRAND_TEXT = {1: "→", -1: "←", 0: "·"}
STRAND_WORDS = {1: "forward strand", -1: "reverse strand", 0: "no strand"}


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.resize(1300, 820)

        self.settings = QSettings("mikeph52", "N-Zoom")
        self.records = []
        self.current_record = None
        self.current_sequence = ""
        self.current_features = []   # annotated features of the current record
        self.current_orfs = []       # ORFs found by the ORF finder (not saved)
        self.current_reads = []      # aligned sequencing reads, shown as items (not saved)
        self.reads_by_record = {}    # id(record) -> read results, kept while N-Zoom runs
        self.reads_dialog = None
        self.items = []              # features + ORFs + reads: what the map and table show

        # file state
        self.current_path = None     # where Save writes to
        self.current_format = None   # "genbank" or "fasta"
        # unsaved changes: feature edits live on the undo stack (undoing back to
        # the last save makes the file "unmodified" again); records fetched from
        # NCBI aren't undoable, so they set a separate flag
        self.undo_stack = QUndoStack(self)
        self.undo_stack.cleanChanged.connect(self.on_clean_changed)
        self.unsaved_records = False
        self.find_hits = []
        self.find_index = -1

        # background work (NCBI downloads, restriction searches, ORF finding)
        self.pool = QThreadPool(self)
        self.tasks = {}              # kind -> running Task, e.g. {"ncbi": task}
        self.sites_dialog = None
        self.sites_cache = (None, None)  # (record, all enzyme sites) for the digest tool

        # appearance
        self.theme_mode = self.settings.value("theme", "system")
        if self.theme_mode not in theme.THEME_MODES:
            self.theme_mode = "system"
        self.map_mode = self.settings.value("map_mode", "auto")
        if self.map_mode not in dict(MAP_MODES):
            self.map_mode = "auto"
        theme.apply_theme(QApplication.instance(), self.theme_mode)
        hints = QGuiApplication.styleHints()
        if hasattr(hints, "colorSchemeChanged"):  # follow the system while in "system" mode
            hints.colorSchemeChanged.connect(self.on_system_color_scheme)

        self.build_ui()
        self.build_menu()
        self.update_title()

    # ================================================================ UI
    def build_menu(self):
        file_menu = self.menuBar().addMenu("File")
        self.add_action(file_menu, "Open...", self.open_file, "Ctrl+O")
        self.add_action(file_menu, "Save", self.save, "Ctrl+S")
        self.add_action(file_menu, "Save As...", self.save_as, "Ctrl+Shift+S")
        file_menu.addSeparator()
        self.add_action(file_menu, "Export map image...", self.export_map, "Ctrl+Shift+E")
        file_menu.addSeparator()
        self.add_action(file_menu, "Exit", self.close, "Ctrl+Q")

        edit_menu = self.menuBar().addMenu("Edit")
        undo_action = self.undo_stack.createUndoAction(self, "Undo")
        undo_action.setShortcut(QKeySequence.StandardKey.Undo)
        redo_action = self.undo_stack.createRedoAction(self, "Redo")
        redo_action.setShortcut(QKeySequence.StandardKey.Redo)
        edit_menu.addAction(undo_action)
        edit_menu.addAction(redo_action)
        edit_menu.addSeparator()
        find_action = self.add_action(edit_menu, "Find sequence...", self.show_find)
        find_action.setShortcut(QKeySequence.StandardKey.Find)
        next_action = self.add_action(edit_menu, "Find next", self.find_next)
        next_action.setShortcut(QKeySequence.StandardKey.FindNext)
        previous_action = self.add_action(edit_menu, "Find previous", self.find_previous)
        previous_action.setShortcut(QKeySequence.StandardKey.FindPrevious)
        edit_menu.addSeparator()
        self.add_action(edit_menu, "Add feature...", self.add_feature, "Ctrl+E")
        self.add_action(edit_menu, "Edit selected feature...", self.edit_selected_feature)
        self.add_action(edit_menu, "Delete selected feature", self.delete_selected_feature)

        view_menu = self.menuBar().addMenu("View")
        map_menu = view_menu.addMenu("Map style")
        self.map_mode_group = QActionGroup(self)
        for mode, text in MAP_MODES:
            action = self.add_action(map_menu, text, lambda _=False, m=mode: self.set_map_mode(m))
            action.setCheckable(True)
            action.setChecked(mode == self.map_mode)
            self.map_mode_group.addAction(action)
        theme_menu = view_menu.addMenu("Theme")
        self.theme_group = QActionGroup(self)
        for mode in theme.THEME_MODES:
            action = self.add_action(theme_menu, THEME_NAMES[mode],
                                     lambda _=False, m=mode: self.set_theme_mode(m))
            action.setCheckable(True)
            action.setChecked(mode == self.theme_mode)
            self.theme_group.addAction(action)

        tools_menu = self.menuBar().addMenu("Tools")
        self.ncbi_action = self.add_action(tools_menu, "Fetch from NCBI...", self.fetch_from_ncbi)
        tools_menu.addSeparator()
        self.add_action(tools_menu, "Translate whole sequence", self.show_translation)
        tools_menu.addSeparator()
        self.rs_action = self.add_action(tools_menu, "Find restriction sites", self.find_restriction_sites)
        self.add_action(tools_menu, "Clear restriction sites", lambda: self.maps.set_sites({}))
        self.digest_action = self.add_action(tools_menu, "Virtual digest...", self.open_digest, "Ctrl+D")
        tools_menu.addSeparator()
        self.orf_action = self.add_action(tools_menu, "Find ORFs...", self.find_orfs)
        self.add_action(tools_menu, "Clear ORFs", self.clear_orfs)
        tools_menu.addSeparator()
        self.add_action(tools_menu, "Find primer...", self.open_primer_dialog, "Ctrl+P")
        tools_menu.addSeparator()
        self.reads_action = self.add_action(tools_menu, "Align sequencing reads...",
                                            self.align_reads, "Ctrl+L")
        self.add_action(tools_menu, "Show aligned reads", self.show_reads_dialog)
        self.add_action(tools_menu, "Clear reads", self.clear_reads)

        about_menu = self.menuBar().addMenu("About")
        self.add_action(about_menu, "About", self.show_about)

    def add_action(self, menu, text, slot, shortcut=None):
        action = QAction(text, self)
        if shortcut:
            action.setShortcut(shortcut)
        action.triggered.connect(slot)
        menu.addAction(action)
        return action

    def build_ui(self):
        root = QWidget()
        root_layout = QVBoxLayout(root)

        # top bar
        top_bar = QHBoxLayout()
        open_button = QPushButton("Open")
        open_button.clicked.connect(self.open_file)
        top_bar.addWidget(open_button)
        self.ncbi_button = QPushButton("Fetch from NCBI")
        self.ncbi_button.clicked.connect(self.fetch_from_ncbi)
        top_bar.addWidget(self.ncbi_button)
        self.file_label = QLabel("No file loaded")
        self.file_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        top_bar.addWidget(self.file_label, 1)
        root_layout.addLayout(top_bar)

        splitter = QSplitter(Qt.Orientation.Horizontal)

        # left: record list
        self.sequence_list = QListWidget()
        self.sequence_list.setMinimumWidth(230)
        self.sequence_list.currentRowChanged.connect(self.select_record)
        splitter.addWidget(self.sequence_list)

        # center: title, viewer tabs, feature table
        center = QWidget()
        center_layout = QVBoxLayout(center)
        self.sequence_title = QLabel("Sequence")
        title_font = self.sequence_title.font()
        title_font.setBold(True)
        title_font.setPointSize(13)
        self.sequence_title.setFont(title_font)
        center_layout.addWidget(self.sequence_title)
        self.find_bar = self.build_find_bar()
        center_layout.addWidget(self.find_bar)

        center_split = QSplitter(Qt.Orientation.Vertical)

        self.view_tabs = QTabWidget()
        # both map styles are kept up to date; the stack shows one of them
        self.circular_viewer = CircularDNAViewer()
        self.linear_viewer = LinearDNAViewer()
        self.maps = MapGroup(self.circular_viewer, self.linear_viewer)
        # each map sits in a panel with zoom buttons (the linear one also gets a scrollbar)
        self.circular_panel = MapPanel(self.circular_viewer, scrollbar=False)
        self.linear_panel = MapPanel(self.linear_viewer, scrollbar=True)
        self.map_stack = QStackedWidget()
        for view in (self.circular_viewer, self.linear_viewer):
            view.featureClicked.connect(self.select_feature_row)
            view.featureContextMenu.connect(self.map_context_menu)
        self.map_stack.addWidget(self.circular_panel)
        self.map_stack.addWidget(self.linear_panel)
        self.view_tabs.addTab(self.map_stack, "Map")
        self.viewer = SequenceViewer()
        self.viewer.selectionInfo.connect(self.on_sequence_selection)
        self.viewer.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.viewer.customContextMenuRequested.connect(self.sequence_context_menu)
        self.view_tabs.addTab(self.viewer, "Sequence")
        center_split.addWidget(self.view_tabs)

        self.feature_table = QTableWidget(0, 6)
        self.feature_table.setHorizontalHeaderLabels(
            ["Type", "Name", "Start", "End", "Strand", "Length"]
        )
        self.feature_table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.feature_table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.feature_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.feature_table.verticalHeader().setVisible(False)
        # start sorted by position (Qt would default to descending by the first column)
        self.feature_table.horizontalHeader().setSortIndicator(2, Qt.SortOrder.AscendingOrder)
        self.feature_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.feature_table.itemSelectionChanged.connect(self.on_feature_selected)
        self.feature_table.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.feature_table.customContextMenuRequested.connect(self.table_context_menu)
        self.feature_table.cellDoubleClicked.connect(
            lambda row, _col: self.edit_item(self.item_index_at_row(row))
        )
        delete_shortcut = QShortcut(QKeySequence(QKeySequence.StandardKey.Delete), self.feature_table)
        delete_shortcut.setContext(Qt.ShortcutContext.WidgetShortcut)
        delete_shortcut.activated.connect(self.delete_selected_feature)
        center_split.addWidget(self.feature_table)
        center_split.setSizes([520, 200])

        center_layout.addWidget(center_split, 1)
        splitter.addWidget(center)

        # right: info + actions
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        info_box = QGroupBox("Sequence information")
        form = QFormLayout(info_box)
        self.type_value = QLabel("-")
        self.length_value = QLabel("-")
        self.topology_value = QLabel("-")
        self.gc_value = QLabel("-")
        self.features_value = QLabel("-")
        self.composition_value = QLabel("-")
        self.composition_value.setWordWrap(True)
        form.addRow("Type:", self.type_value)
        form.addRow("Length:", self.length_value)
        form.addRow("Topology:", self.topology_value)
        form.addRow("GC %:", self.gc_value)
        form.addRow("Features:", self.features_value)
        form.addRow("Composition:", self.composition_value)
        right_layout.addWidget(info_box)

        actions_box = QGroupBox("Actions")
        actions_layout = QVBoxLayout(actions_box)

        def add_button(text, slot):
            button = QPushButton(text)
            button.clicked.connect(slot)
            actions_layout.addWidget(button)
            return button

        add_button("Copy sequence", self.copy_seq)
        add_button("Reverse complement", self.show_reverse_complement)
        add_button("Translate DNA", self.show_translation)
        self.rs_button = add_button("Find restriction sites", self.find_restriction_sites)
        self.digest_button = add_button("Virtual digest", self.open_digest)
        self.orf_button = add_button("Find ORFs", self.find_orfs)
        add_button("Find primer", self.open_primer_dialog)
        self.reads_button = add_button("Align reads", self.align_reads)
        add_button("Add feature", self.add_feature)
        right_layout.addWidget(actions_box)
        right_layout.addStretch()

        splitter.addWidget(right_panel)
        splitter.setSizes([240, 800, 260])
        root_layout.addWidget(splitter, 1)

        self.setCentralWidget(root)

        # status bar: temporary messages on the left, permanent widgets on the right
        self.busy_bar = QProgressBar()
        self.busy_bar.setRange(0, 0)  # 0..0 = animated "working" bar
        self.busy_bar.setMaximumWidth(120)
        self.busy_bar.setTextVisible(False)
        self.busy_bar.hide()
        self.statusBar().addPermanentWidget(self.busy_bar)
        self.position_label = QLabel("")
        self.position_label.setMinimumWidth(260)
        self.position_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self.statusBar().addPermanentWidget(self.position_label)
        self.statusBar().showMessage("Ready")

    # ================================================================ file state
    def update_title(self):
        name = Path(self.current_path).name if self.current_path else "Untitled"
        if not self.records:
            name = "No file"
        self.setWindowTitle(f"N-Zoom {VERSION} — {name}{' *' if self.modified else ''}")

    @property
    def modified(self):
        return self.unsaved_records or not self.undo_stack.isClean()

    @modified.setter
    def modified(self, value):
        self.unsaved_records = bool(value)
        if not value:
            self.undo_stack.setClean()
        self.update_title()

    def set_modified(self, modified=True):
        self.modified = modified

    def on_clean_changed(self, _clean):
        self.update_title()

    def confirm_discard_changes(self):
        """True if it's OK to throw away the current records (saving first if asked)."""
        if not self.modified:
            return True
        answer = QMessageBox.question(
            self, "Unsaved changes",
            "You have unsaved changes. Save them first?",
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel,
            QMessageBox.StandardButton.Save,
        )
        if answer == QMessageBox.StandardButton.Save:
            return self.save()
        return answer == QMessageBox.StandardButton.Discard

    def has_annotations(self):
        return any(st.display_features(record) for record in self.records)

    # ================================================================ loading
    def open_file(self):
        if not self.confirm_discard_changes():
            return
        path, _ = QFileDialog.getOpenFileName(self, "Open sequence file", "", OPEN_FILTER)
        if not path:
            return
        try:
            records, fmt = st.load_records(path)
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Could not read file:\n{exc}")
            return
        if not records:
            QMessageBox.warning(self, "No sequences", f"No {fmt.upper()} records were found.")
            return

        self.current_path, self.current_format = path, fmt
        self.undo_stack.clear()  # old commands refer to the records being closed
        self.reads_by_record.clear()  # so do the aligned reads
        self.show_records(records)
        self.file_label.setText(str(Path(path)))
        self.set_modified(False)
        self.statusBar().showMessage(f"Loaded {len(records)} {fmt} record(s)")

    def fetch_from_ncbi(self):
        if "ncbi" in self.tasks:
            return  # a download is already running
        accession, ok = QInputDialog.getText(
            self, "Fetch from NCBI", "Nucleotide accession (e.g. L09137 for pUC19):"
        )
        accession = accession.strip()
        if not ok or not accession:
            return

        email = self.settings.value("ncbi_email", "")
        if not email:
            email, ok = QInputDialog.getText(
                self, "NCBI e-mail",
                "NCBI asks for an e-mail address with every request.\n"
                "It's saved locally and only sent to NCBI:",
            )
            email = email.strip()
            if not ok or "@" not in email:
                return
            self.settings.setValue("ncbi_email", email)

        self.start_task(
            "ncbi", f"Fetching {accession} from NCBI...",
            st.fetch_genbank, accession, email,
            on_done=self.add_fetched_record,
            on_fail=lambda message: QMessageBox.critical(
                self, "NCBI error", f"Could not fetch {accession}:\n{message}"
            ),
        )

    def add_fetched_record(self, record):
        # added to whatever is already open; not saved anywhere yet
        self.show_records(self.records + [record], select=len(self.records))
        self.set_modified()
        self.statusBar().showMessage(f"Fetched {record.id} from NCBI", 5000)

    def show_records(self, records, select=0):
        self.records = records
        self.current_record = None
        self.sequence_list.blockSignals(True)
        self.sequence_list.clear()
        for record in records:
            label = record.id
            if record.description and record.description != record.id:
                label += f" — {record.description}"
            self.sequence_list.addItem(label)
        self.sequence_list.blockSignals(False)
        self.sequence_list.setCurrentRow(select)
        self.select_record(select)
        self.update_title()

    # ================================================================ saving
    def save(self):
        """Save to the current file. Returns True if the records were saved."""
        if not self.records:
            QMessageBox.information(self, "Nothing to save", "Open or fetch a sequence first.")
            return False
        # FASTA can't hold features: if there are any, ask where to save as GenBank
        if self.current_path and (self.current_format == "genbank" or not self.has_annotations()):
            return self.write_file(self.current_path, self.current_format)
        return self.save_as()

    def save_as(self):
        if not self.records:
            QMessageBox.information(self, "Nothing to save", "Open or fetch a sequence first.")
            return False
        suggested = GENBANK_FILTER if self.has_annotations() or self.current_format != "fasta" else FASTA_FILTER
        start_dir = str(Path(self.current_path).with_suffix("")) if self.current_path else ""
        path, chosen = QFileDialog.getSaveFileName(
            self, "Save sequences", start_dir, f"{GENBANK_FILTER};;{FASTA_FILTER}", suggested
        )
        if not path:
            return False
        fmt = "fasta" if chosen == FASTA_FILTER else "genbank"
        if not Path(path).suffix:
            path += ".gb" if fmt == "genbank" else ".fasta"
        if fmt == "fasta" and self.has_annotations():
            answer = QMessageBox.question(
                self, "Features will be lost",
                "FASTA files can't store features, so all annotations will be left out "
                "of this file. Save as FASTA anyway?",
            )
            if answer != QMessageBox.StandardButton.Yes:
                return False
        return self.write_file(path, fmt)

    def write_file(self, path, fmt):
        try:
            renamed = st.save_records(self.records, path, fmt)
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Could not save file:\n{exc}")
            return False
        self.current_path, self.current_format = path, fmt
        self.file_label.setText(str(Path(path)))
        self.set_modified(False)
        self.statusBar().showMessage(f"Saved {len(self.records)} record(s) to {path}", 5000)
        if renamed:
            names = "\n".join(f"{old}  →  {new}" for old, new in renamed)
            QMessageBox.information(
                self, "Names shortened",
                f"GenBank files allow at most {st.GENBANK_NAME_LIMIT} characters, without "
                f"spaces, for a sequence name. These names were shortened in the file:\n\n{names}",
            )
        return True

    # ================================================================ display
    def select_record(self, index):
        if index < 0 or index >= len(self.records):
            return

        record = self.records[index]
        self.current_record = record
        sequence = st.record_sequence(record)
        seq_type = st.detect_seq_type(sequence)
        self.current_sequence = sequence
        self.current_features = st.display_features(record)
        self.current_orfs = []
        self.current_reads = self.read_items(self.reads_by_record.get(id(record), []))
        self.items = self.current_features + self.current_reads
        if self.reads_dialog is not None:
            self.reads_dialog.close()  # it shows another record's reads

        self.sequence_title.setText(record.description or record.id)
        self.viewer.set_sequence(sequence, seq_type)

        # map only makes sense for DNA
        if seq_type == "DNA":
            self.maps.set_record(len(sequence), record.id, self.items, sequence=sequence)
            self.view_tabs.setTabEnabled(0, True)
        else:
            self.maps.clear_seq()
            self.view_tabs.setTabEnabled(0, False)
            self.view_tabs.setCurrentIndex(1)

        self.update_map_style()
        self.fill_feature_table()

        # info panel
        self.type_value.setText(seq_type)
        self.length_value.setText(f"{len(sequence):,}")
        self.topology_value.setText(st.topology(record))
        self.update_feature_count()

        counts = Counter(sequence)
        if seq_type == "DNA" and sequence:
            self.gc_value.setText(f"{st.gc_percent(sequence):.2f}%")
        else:
            self.gc_value.setText("N/A")
        composition = ", ".join(f"{base}: {count:,}" for base, count in counts.most_common(12))
        self.composition_value.setText(composition or "-")
        if self.find_bar.isVisible():
            self.run_find()

    def update_feature_count(self):
        text = str(len(self.current_features))
        for things, word in ((self.current_orfs, "ORF"), (self.current_reads, "read")):
            if things:
                text += f"  (+ {len(things)} {word}{'' if len(things) == 1 else 's'})"
        self.features_value.setText(text)

    def refresh_features(self, select=None):
        """Redraw map and table after features or ORFs changed.
        `select`: a feature/ORF object to select afterwards."""
        if self.current_record is None:
            return
        self.current_features = st.display_features(self.current_record)
        self.items = self.current_features + self.current_orfs + self.current_reads
        if self.maps.length:
            self.maps.set_features(self.items)
        self.fill_feature_table()
        self.update_feature_count()
        if select is not None:
            for index, item in enumerate(self.items):
                if item is select:
                    self.select_feature_row(index)
                    break

    # self.items = annotated features, then ORFs, then reads
    def is_suggestion(self, index):
        """ORFs and reads: shown, but not part of the file."""
        return index >= len(self.current_features)

    def is_orf(self, index):
        return len(self.current_features) <= index < len(self.current_features) + len(self.current_orfs)

    def is_read(self, index):
        return index >= len(self.current_features) + len(self.current_orfs)

    def fill_feature_table(self):
        table = self.feature_table
        table.blockSignals(True)
        table.setSortingEnabled(False)
        table.clearContents()
        table.setRowCount(len(self.items))

        orf_font = QFont(table.font())
        orf_font.setItalic(True)
        for row, feature in enumerate(self.items):
            start, end = st.feature_span(feature)  # 1-based like GenBank
            length = sum(e - s for s, e in st.feature_parts(feature))
            cells = [
                QTableWidgetItem(feature.type),
                QTableWidgetItem(st.feature_label(feature)),
                number_item(start),
                number_item(end),
                QTableWidgetItem(STRAND_TEXT[st.feature_strand(feature)]),
                number_item(length),
            ]
            cells[0].setData(Qt.ItemDataRole.UserRole, row)  # index into self.items
            for column, cell in enumerate(cells):
                if self.is_suggestion(row):  # ORFs/reads: grey italic, not part of the file
                    cell.setForeground(Qt.GlobalColor.darkGray)
                    cell.setFont(orf_font)
                table.setItem(row, column, cell)

        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        table.blockSignals(False)
        self.on_feature_selected()

    def item_index_at_row(self, row):
        cell = self.feature_table.item(row, 0)
        return cell.data(Qt.ItemDataRole.UserRole) if cell is not None else -1

    def selected_item_index(self):
        rows = self.feature_table.selectionModel().selectedRows()
        return self.item_index_at_row(rows[0].row()) if rows else -1

    def on_feature_selected(self):
        index = self.selected_item_index()
        if index < 0:
            self.viewer.clear_highlight()
            self.maps.set_selected(-1)
            return
        feature = self.items[index]
        self.viewer.highlight_parts(st.feature_parts(feature))
        self.maps.set_selected(index)
        self.statusBar().showMessage(
            f"{feature.type}: {st.feature_label(feature)}  ({feature.location})", 8000
        )

    def select_feature_row(self, index):
        """Select the table row showing self.items[index] (e.g. after a map click)."""
        for row in range(self.feature_table.rowCount()):
            if self.item_index_at_row(row) == index:
                self.feature_table.selectRow(row)
                self.feature_table.scrollToItem(self.feature_table.item(row, 0))
                return

    def on_sequence_selection(self, start, end):
        """Show the cursor position or selected range in the status bar."""
        sequence = self.current_sequence
        if not sequence:
            self.position_label.setText("")
            self.maps.set_highlight()
            return

        is_dna = st.detect_seq_type(sequence) == "DNA"
        unit = "bp" if is_dna else "aa"
        if start == end:
            position = min(start + 1, len(sequence))
            self.position_label.setText(f"Position {position:,} of {len(sequence):,}")
            self.maps.set_highlight()
            return

        text = f"Selected {start + 1:,}–{end:,}  ({end - start:,} {unit})"
        if is_dna:
            text += f"  ·  GC {st.gc_percent(sequence[start:end]):.1f}%"
            self.maps.set_highlight(start, end)
        self.position_label.setText(text)

    # ================================================================ right-click menus
    def is_dna(self):
        return bool(self.current_sequence) and st.detect_seq_type(self.current_sequence) == "DNA"

    def sequence_context_menu(self, pos):
        self.build_selection_menu().exec(self.viewer.mapToGlobal(pos))

    def table_context_menu(self, pos):
        row = self.feature_table.rowAt(pos.y())
        if row < 0:
            return
        self.feature_table.selectRow(row)
        menu = self.build_item_menu(self.item_index_at_row(row))
        menu.exec(self.feature_table.viewport().mapToGlobal(pos))

    def map_context_menu(self, index, global_pos):
        self.select_feature_row(index)
        self.build_item_menu(index).exec(global_pos)

    def build_selection_menu(self):
        """Menu for the Sequence view: acts on the selected bases."""
        menu = QMenu(self)
        selection = self.viewer.selected_range()
        if selection:
            start, end = selection
            region = self.current_sequence[start:end]
            where = f"{start + 1:,}–{end:,}"
            menu.addAction("Copy", lambda: self.copy_text(region, f"Copied {len(region):,} bases"))
            if self.is_dna():
                menu.addAction("Copy reverse complement", lambda: self.copy_text(
                    st.reverse_complement(region), "Reverse complement copied"))
                menu.addSeparator()
                menu.addAction("Translate selection", lambda: self.show_protein(
                    region, f"Selection {where}, forward strand"))
                menu.addAction("Translate reverse complement", lambda: self.show_protein(
                    st.reverse_complement(region), f"Selection {where}, reverse strand"))
            menu.addSeparator()
            menu.addAction("Create feature from selection...",
                           lambda: self.add_feature(start=start + 1, end=end))
            menu.addSeparator()
        else:
            menu.addAction("Add feature...", self.add_feature)
        menu.addAction("Select all", self.viewer.selectAll)
        return menu

    def build_item_menu(self, index):
        """Menu for one feature, ORF or read (feature table or map)."""
        menu = QMenu(self)
        if index < 0 or index >= len(self.items):
            return menu
        item = self.items[index]
        orf = self.is_orf(index)

        title = menu.addAction(f"{item.type}: {st.feature_label(item)}")
        title.setEnabled(False)
        menu.addSeparator()
        if self.is_read(index):
            name = item.qualifiers["read_name"][0]
            menu.addAction("Show alignment and chromatogram", lambda: self.show_reads_dialog(name))
            menu.addAction("Select in sequence", lambda: self.select_item_in_sequence(index))
            menu.addSeparator()
            menu.addAction("Remove this read", lambda: self.remove_read(name))
            return menu
        menu.addAction("Select in sequence", lambda: self.select_item_in_sequence(index))
        menu.addAction("Copy sequence", lambda: self.copy_text(
            st.feature_sequence(self.current_sequence, item),
            f"Copied sequence of {st.feature_label(item)}"))
        if self.is_dna():
            menu.addAction("Translate", lambda: self.translate_item(index))
        menu.addSeparator()
        if orf:
            menu.addAction("Add as CDS feature...", lambda: self.add_orf_as_feature(index))
        else:
            menu.addAction("Edit feature...", lambda: self.edit_item(index))
            menu.addAction("Delete feature", lambda: self.delete_item(index))
        return menu

    # ================================================================ region actions
    def copy_text(self, text, message):
        QApplication.clipboard().setText(text)
        self.statusBar().showMessage(message, 3000)

    def show_protein(self, dna, description, table=1, codon_start=1, expected=None):
        """Translate `dna` and show the protein in a result window."""
        protein = st.translate_dna(dna, table=table, codon_start=codon_start)
        if not protein:
            QMessageBox.warning(self, "Too short", "There are fewer than 3 bases to translate.")
            return None
        ends_with_stop = protein.endswith("*")
        amino_acids = len(protein) - ends_with_stop
        details = [description, f"{len(dna):,} bp → {amino_acids:,} aa"
                   + (" + stop" if ends_with_stop else "")]
        if table != 1:
            details.append(f"genetic code {table}")
        if codon_start != 1:
            details.append(f"starting at base {codon_start}")
        leftover = (len(dna) - codon_start + 1) % 3
        if leftover:
            details.append(f"last {leftover} base(s) ignored")
        if "*" in protein[:-1]:
            details.append(f"{protein[:-1].count('*')} internal stop codon(s)")
        if expected:
            same = protein.rstrip("*") == expected.rstrip("*")
            details.append("matches the /translation stored in the file" if same
                           else "differs from the /translation stored in the file")
        TextResultDialog(self, "Translation", "  ·  ".join(details), protein).show()
        return protein

    def translate_item(self, index):
        item = self.items[index]
        table, codon_start = st.feature_translation_settings(item)
        start, end = st.feature_span(item)
        strand = STRAND_WORDS[st.feature_strand(item)]
        self.show_protein(
            st.feature_sequence(self.current_sequence, item),
            f"{st.feature_label(item)} ({item.type}, {start:,}–{end:,}, {strand})",
            table=table, codon_start=codon_start,
            expected=item.qualifiers.get("translation", [None])[0],
        )

    def select_item_in_sequence(self, index):
        item = self.items[index]
        parts = st.feature_parts(item)
        self.view_tabs.setCurrentWidget(self.viewer)
        if len(parts) == 1:
            self.viewer.select_range(*parts[0])  # a real selection: readout + menu work on it
        else:
            self.viewer.highlight_parts(parts)   # pieces can't be one selection; highlight them

    # ================================================================ feature editing
    def add_feature(self, start=None, end=None, strand=1, feature_type="misc_feature", name=""):
        if self.current_record is None or not self.current_sequence:
            QMessageBox.information(self, "No sequence", "Open or fetch a sequence first.")
            return None
        if start is None:  # default to the current selection, if any
            selection = self.viewer.selected_range()
            start, end = (selection[0] + 1, selection[1]) if selection else (1, 1)
        dialog = FeatureDialog(
            self, len(self.current_sequence), st.is_circular(self.current_record),
            start=start, end=end, strand=strand, feature_type=feature_type, name=name,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return None
        values = dialog.values()
        feature = st.new_feature(values["type"], values["name"], values["location"], values["note"])
        self.undo_stack.push(AddFeatures(
            self.current_record, [feature], self.on_features_changed,
            f"Add {feature.type} \u201c{st.feature_label(feature)}\u201d"))
        return feature

    def add_orf_as_feature(self, index):
        orf = self.items[index]
        start, end = st.feature_span(orf)
        number = st.feature_label(orf).split()[0]  # "ORF3"
        feature = self.add_feature(start=start, end=end, strand=st.feature_strand(orf),
                                   feature_type="CDS", name=number)
        if feature is not None and orf in self.current_orfs:
            # it's a real feature now, so drop the ORF suggestion to avoid a duplicate
            self.current_orfs = [o for o in self.current_orfs if o is not orf]
            self.refresh_features(select=feature)

    def edit_item(self, index):
        if index < 0 or index >= len(self.items):
            return
        if self.is_orf(index):
            self.add_orf_as_feature(index)
            return
        if self.is_read(index):
            self.show_reads_dialog(self.items[index].qualifiers["read_name"][0])
            return
        feature = self.items[index]
        dialog = FeatureDialog(
            self, len(self.current_sequence), st.is_circular(self.current_record), feature=feature
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        values = dialog.values()
        before = snapshot(feature)
        st.update_feature(feature, values["type"], values["name"],
                          location=values["location"], note=values["note"])
        self.undo_stack.push(EditFeature(self.current_record, feature, before, snapshot(feature),
                                         self.on_features_changed))

    def delete_item(self, index):
        if index < 0 or index >= len(self.items) or self.is_suggestion(index):
            return
        feature = self.items[index]
        answer = QMessageBox.question(
            self, "Delete feature",
            f"Delete the {feature.type} “{st.feature_label(feature)}”?",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.undo_stack.push(DeleteFeature(self.current_record, feature, self.on_features_changed))

    def on_features_changed(self, record, feature, message):
        """Called by every undo command after it changes a record's features."""
        if record is not self.current_record:
            # undoing a change to another record: show that record first
            for index, candidate in enumerate(self.records):
                if candidate is record:
                    self.sequence_list.setCurrentRow(index)
                    break
            else:
                return  # the record is no longer open
        self.refresh_features(select=feature)
        self.update_title()
        self.statusBar().showMessage(message, 5000)

    def edit_selected_feature(self):
        index = self.selected_item_index()
        if index < 0:
            QMessageBox.information(self, "No feature selected", "Select a feature in the table first.")
            return
        self.edit_item(index)

    def delete_selected_feature(self):
        index = self.selected_item_index()
        if index >= 0:
            self.delete_item(index)

    # ================================================================ background tasks
    def start_task(self, kind, message, fn, *args, on_done=None, on_fail=None, **kwargs):
        """Run fn(*args, **kwargs) on a background thread.

        `kind` names the job ("ncbi", "restriction", "orfs"); only one job of each
        kind runs at a time. on_done(result) / on_fail(message) run in the GUI thread.
        """
        if kind in self.tasks:
            self.statusBar().showMessage("Still working on the previous request...", 3000)
            return
        task = Task(fn, *args, on_done=on_done, on_fail=on_fail, **kwargs)
        task.kind = kind
        # connect to methods of this window, so Qt runs them in the GUI thread
        task.signals.finished.connect(self.on_task_finished)
        task.signals.failed.connect(self.on_task_failed)
        self.tasks[kind] = task
        self.update_busy_state()
        self.statusBar().showMessage(message)
        self.pool.start(task)

    def on_task_finished(self, payload):
        task, result = payload
        self.end_task(task)
        if task.on_done:
            task.on_done(result)

    def on_task_failed(self, payload):
        task, message = payload
        self.end_task(task)
        self.statusBar().showMessage("Failed: " + message, 8000)
        if task.on_fail:
            task.on_fail(message)

    def end_task(self, task):
        self.tasks.pop(task.kind, None)
        self.statusBar().clearMessage()
        self.update_busy_state()

    def update_busy_state(self):
        self.busy_bar.setVisible(bool(self.tasks))
        for kind, widgets in (
            ("ncbi", (self.ncbi_action, self.ncbi_button)),
            ("restriction", (self.rs_action, self.rs_button)),
            ("orfs", (self.orf_action, self.orf_button)),
            ("digest", (self.digest_action, self.digest_button)),
            ("reads", (self.reads_action, self.reads_button)),
        ):
            for widget in widgets:
                widget.setEnabled(kind not in self.tasks)

    def closeEvent(self, event):
        if not self.confirm_discard_changes():
            event.ignore()
            return
        if self.tasks:
            # can't interrupt a download midway; wait for it so we exit cleanly
            self.statusBar().showMessage("Finishing background work before closing...")
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            self.pool.waitForDone()
            QApplication.restoreOverrideCursor()
        super().closeEvent(event)

    def stale(self, record, what):
        """True (and a status message) if the user switched records while `what` ran."""
        if record is self.current_record:
            return False
        self.statusBar().showMessage(
            f"{what} for {record.id} finished, but another record is selected now. "
            "Run it again to see the results.", 8000
        )
        return True

    # ================================================================ tools
    def require_dna(self, what):
        if self.current_record is None:
            return False
        if st.detect_seq_type(self.current_sequence) != "DNA":
            QMessageBox.information(self, "DNA required", f"{what} is only available for DNA sequences.")
            return False
        return True

    def find_restriction_sites(self):
        if not self.require_dna("Restriction analysis"):
            return
        if len(self.current_sequence) > 2_000_000:
            answer = QMessageBox.question(
                self, "Large sequence",
                "This sequence is very long, so the search may take a while. Continue?",
            )
            if answer != QMessageBox.StandardButton.Yes:
                return

        record = self.current_record
        circular = st.is_circular(record)
        self.start_task(
            "restriction", f"Searching restriction sites in {record.id}...",
            st.find_restriction_sites, self.current_sequence, circular=circular,
            on_done=lambda sites: self.show_restriction_results(record, circular, sites),
            on_fail=lambda message: QMessageBox.critical(
                self, "Error", f"Restriction analysis failed:\n{message}"
            ),
        )

    def show_restriction_results(self, record, circular, sites):
        if self.stale(record, "Restriction search"):
            return
        if circular == st.is_circular(record):
            self.sites_cache = (record, sites)
        on_map = st.map_sites(sites)
        self.maps.set_sites(on_map)
        self.view_tabs.setCurrentIndex(0)
        self.statusBar().showMessage(
            f"{len(sites)} enzymes cut; {len(on_map)} unique 6+ bp sites shown on map", 8000
        )
        # non-modal: the main window stays usable while the list is open
        if self.sites_dialog is not None:
            self.sites_dialog.close()
        self.sites_dialog = RestrictionSitesDialog(self, sites, circular)
        self.sites_dialog.destroyed.connect(self.forget_sites_dialog)
        self.sites_dialog.show()

    def forget_sites_dialog(self):
        self.sites_dialog = None

    def find_orfs(self):
        if not self.require_dna("The ORF finder"):
            return
        min_aa, ok = QInputDialog.getInt(
            self, "Find ORFs",
            "Minimum ORF length (amino acids, not counting the stop codon):",
            int(self.settings.value("orf_min_aa", st.DEFAULT_MIN_ORF_AA)), 10, 100_000,
        )
        if not ok:
            return
        self.settings.setValue("orf_min_aa", min_aa)

        record = self.current_record
        circular = st.is_circular(record)
        self.start_task(
            "orfs", f"Finding ORFs in {record.id}...",
            st.find_orfs, self.current_sequence, min_aa=min_aa, circular=circular,
            on_done=lambda orfs: self.show_orfs(record, orfs, min_aa),
            on_fail=lambda message: QMessageBox.critical(
                self, "Error", f"ORF search failed:\n{message}"
            ),
        )

    def show_orfs(self, record, orfs, min_aa):
        if self.stale(record, "ORF search"):
            return
        # don't repeat ORFs that are already annotated with exactly the same location
        annotated = {
            (tuple(sorted(st.feature_parts(f))), st.feature_strand(f))
            for f in self.current_features
        }
        new = [o for o in orfs
               if (tuple(sorted(st.feature_parts(o))), st.feature_strand(o)) not in annotated]
        for number, orf in enumerate(new, 1):
            orf.qualifiers["label"] = [f"ORF{number} ({orf.qualifiers['aa'][0]} aa)"]
        hidden = len(orfs) - len(new)

        self.current_orfs = new
        self.refresh_features()
        if new:
            self.view_tabs.setCurrentIndex(0)
            message = (f"Found {len(new)} new ORF(s) of at least {min_aa} aa. "
                       "Right-click one to add it as a feature.")
        else:
            message = f"No new ORFs of at least {min_aa} aa found."
        if hidden:
            message += f" ({hidden} already annotated, not shown.)"
        self.statusBar().showMessage(message, 8000)

    def clear_orfs(self):
        self.current_orfs = []
        self.refresh_features()

    # ================================================================ find
    def build_find_bar(self):
        bar = QWidget()
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(QLabel("Find:"))
        self.find_edit = QLineEdit()
        self.find_edit.setPlaceholderText("Sequence to find, e.g. GGATCC (IUPAC codes like N, R, Y allowed)")
        self.find_edit.textChanged.connect(self.run_find)
        self.find_edit.returnPressed.connect(self.find_next)
        layout.addWidget(self.find_edit, 1)
        self.find_both = QCheckBox("Both strands")
        self.find_both.setChecked(True)
        self.find_both.toggled.connect(self.run_find)
        layout.addWidget(self.find_both)
        self.find_status = QLabel("")
        self.find_status.setMinimumWidth(170)
        layout.addWidget(self.find_status)
        for text, tip, slot in (("\u25b2", "Previous match (Shift+Enter)", self.find_previous),
                                ("\u25bc", "Next match (Enter)", self.find_next),
                                ("\u2715", "Close (Esc)", self.hide_find)):
            button = QToolButton()
            button.setText(text)
            button.setToolTip(tip)
            button.clicked.connect(slot)
            layout.addWidget(button)
        close = QShortcut(QKeySequence(Qt.Key.Key_Escape), bar)
        close.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        close.activated.connect(self.hide_find)
        previous = QShortcut(QKeySequence("Shift+Return"), self.find_edit)
        previous.setContext(Qt.ShortcutContext.WidgetShortcut)
        previous.activated.connect(self.find_previous)
        bar.hide()
        return bar

    def show_find(self):
        if not self.current_sequence:
            return
        # start from the selected bases, if a short stretch is selected
        selection = self.viewer.selected_range()
        if selection and 0 < selection[1] - selection[0] <= 100:
            self.find_edit.setText(self.current_sequence[selection[0]:selection[1]])
        self.find_bar.show()
        self.find_edit.setFocus()
        self.find_edit.selectAll()
        self.run_find()

    def hide_find(self):
        self.find_bar.hide()
        self.find_hits, self.find_index = [], -1
        self.viewer.set_matches([])
        self.maps.set_marks([])
        self.viewer.setFocus()

    def run_find(self):
        if not self.find_bar.isVisible() or not self.current_sequence:
            return
        text = self.find_edit.text()
        dna = self.is_dna()
        self.find_both.setEnabled(dna)
        self.find_hits, self.find_index = [], -1
        error = None
        if text.strip():
            try:
                self.find_hits = st.find_motif(
                    self.current_sequence, text, circular=st.is_circular(self.current_record),
                    both_strands=self.find_both.isChecked(), dna=dna, max_hits=MAX_FIND_HITS)
            except ValueError as exc:
                error = str(exc)
        parts = [(int(p.start), int(p.end)) for hit in self.find_hits for p in hit["location"].parts]
        self.viewer.set_matches(parts[:2000])  # highlighting is slow beyond this
        self.maps.set_marks(parts)
        if error:
            self.find_status.setText(f"<span style='color:#e53935'>{error}</span>")
        elif not text.strip():
            self.find_status.setText("")
        elif not self.find_hits:
            self.find_status.setText("No matches")
        else:
            self.go_to_match(0)

    def go_to_match(self, index):
        hits = self.find_hits
        if not hits:
            return
        self.find_index = index % len(hits)
        hit = hits[self.find_index]
        parts = [(int(p.start), int(p.end)) for p in hit["location"].parts]
        self.viewer.select_range(*max(parts, key=lambda p: p[1] - p[0]))
        count = f"{len(hits):,}{'+' if len(hits) >= MAX_FIND_HITS else ''}"
        strand = {1: "", -1: ", reverse strand", 0: ""}[hit["strand"]]
        self.find_status.setText(f"{self.find_index + 1:,} of {count}{strand}")

    def find_next(self):
        if not self.find_bar.isVisible():
            self.show_find()
        elif self.find_hits:
            self.go_to_match(self.find_index + 1)

    def find_previous(self):
        if self.find_hits:
            self.go_to_match(self.find_index - 1)

    # ================================================================ export
    def current_map_view(self):
        if self.map_stack.currentWidget() is self.circular_panel:
            return self.circular_viewer
        return self.linear_viewer

    def export_map(self):
        if not self.maps.length:
            QMessageBox.information(self, "No map", "Open a DNA sequence to export its map.")
            return
        name = f"{self.current_record.id}-map.png"
        path, chosen = QFileDialog.getSaveFileName(
            self, "Export map image", name, "PNG image (*.png);;SVG image (*.svg)")
        if not path:
            return
        fmt = "svg" if path.lower().endswith(".svg") or (
            chosen.startswith("SVG") and not path.lower().endswith(".png")) else "png"
        if not path.lower().endswith("." + fmt):
            path += "." + fmt
        try:
            export_image(self.current_map_view(), path, fmt, title=self.current_record.id)
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Could not export the map:\n{exc}")
            return
        self.statusBar().showMessage(f"Map saved to {path}", 5000)

    # ================================================================ primers
    def open_primer_dialog(self):
        if not self.require_dna("Primer search"):
            return
        PrimerDialog(self, self.current_sequence, st.is_circular(self.current_record),
                     on_add=self.add_primer_features).exec()

    def add_primer_features(self, features):
        count = len(features)
        self.undo_stack.push(AddFeatures(
            self.current_record, features, self.on_features_changed,
            f"Add {count} primer site{'s' if count != 1 else ''} "
            f"\u201c{st.feature_label(features[0])}\u201d"))

    # ================================================================ sequencing reads
    READ_FILTER = ("Sequencing reads (*.ab1 *.abi *.fastq *.fq *.fasta *.fa *.fas *.seq *.txt);;"
                   "All files (*)")

    def align_reads(self):
        if not self.require_dna("Read alignment"):
            return
        paths, _ = QFileDialog.getOpenFileNames(self, "Align sequencing reads", "", self.READ_FILTER)
        if not paths:
            return
        record = self.current_record
        self.start_task(
            "reads", f"Aligning {len(paths)} read{'s' if len(paths) != 1 else ''} to {record.id}...",
            st.align_reads, self.current_sequence, paths, st.is_circular(record),
            on_done=lambda results: self.add_reads(record, results),
            on_fail=lambda message: QMessageBox.critical(
                self, "Error", f"Read alignment failed:\n{message}"),
        )

    def add_reads(self, record, results):
        # a read aligned again replaces its earlier result
        kept = [r for r in self.reads_by_record.get(id(record), [])
                if r["path"] not in {new["path"] for new in results}]
        self.reads_by_record[id(record)] = kept + results
        aligned = sum(1 for r in results if r["alignment"])
        self.statusBar().showMessage(
            f"Aligned {aligned} of {len(results)} read{'s' if len(results) != 1 else ''}", 8000)
        if record is not self.current_record:
            return  # kept for when that record is shown again
        self.current_reads = self.read_items(self.reads_by_record[id(record)])
        self.refresh_features()
        self.show_reads_dialog(results[0]["name"] if results else None)

    @staticmethod
    def read_items(results):
        """Map/table items for the reads that aligned."""
        items = []
        for result in results:
            alignment = result.get("alignment")
            if alignment:
                label = f"{result['name']} ({alignment['identity']:.1f}%)"
                items.append(st.SeqFeature(alignment["location"], type="read", qualifiers={
                    "label": [label], "read_name": [result["name"]]}))
        return items

    def show_reads_dialog(self, select_name=None):
        if self.current_record is None:
            return
        results = self.reads_by_record.get(id(self.current_record), [])
        if not results:
            QMessageBox.information(self, "No reads", "No reads have been aligned to this "
                                    "sequence yet. Use Tools \u2192 Align sequencing reads.")
            return
        if self.reads_dialog is None:
            self.reads_dialog = ReadsDialog(self, self.current_record.id, results,
                                            on_jump=self.jump_to, on_add=self.align_reads)
            self.reads_dialog.destroyed.connect(self.forget_reads_dialog)
        else:
            self.reads_dialog.set_results(results)
        if select_name:
            self.reads_dialog.select_read(select_name)
        self.reads_dialog.show()
        self.reads_dialog.raise_()

    def forget_reads_dialog(self):
        self.reads_dialog = None

    def jump_to(self, start, end):
        """Select a region in the Sequence view (the map follows)."""
        self.viewer.select_range(start, end)

    def remove_read(self, name):
        if self.current_record is None:
            return
        key = id(self.current_record)
        self.reads_by_record[key] = [r for r in self.reads_by_record.get(key, []) if r["name"] != name]
        self.current_reads = self.read_items(self.reads_by_record[key])
        self.refresh_features()
        if self.reads_dialog is not None:
            if self.reads_by_record[key]:
                self.reads_dialog.set_results(self.reads_by_record[key])
            else:
                self.reads_dialog.close()

    def clear_reads(self):
        if self.current_record is None:
            return
        self.reads_by_record.pop(id(self.current_record), None)
        self.current_reads = []
        self.refresh_features()
        if self.reads_dialog is not None:
            self.reads_dialog.close()

    # ================================================================ digest
    def open_digest(self):
        if not self.require_dna("The virtual digest"):
            return
        record = self.current_record
        circular = st.is_circular(record)
        cached_record, cached_sites = self.sites_cache
        if cached_record is record:
            self.show_digest(record, circular, cached_sites)
            return
        # the enzyme list needs every site in the sequence: find them first
        self.start_task(
            "digest", f"Preparing the digest for {record.id}...",
            st.find_restriction_sites, self.current_sequence, circular=circular,
            on_done=lambda sites: self.show_digest(record, circular, sites),
            on_fail=lambda message: QMessageBox.critical(
                self, "Error", f"Restriction analysis failed:\n{message}"
            ),
        )

    def show_digest(self, record, circular, sites):
        if self.stale(record, "Preparing the digest"):
            return
        self.sites_cache = (record, sites)
        DigestDialog(self, record.id, self.current_sequence, circular, sites).exec()

    # ================================================================ appearance
    def update_map_style(self):
        """Show the circular or linear map, following the user's choice or the topology."""
        record = self.current_record
        if self.map_mode == "circular":
            circular = True
        elif self.map_mode == "linear":
            circular = False
        else:
            circular = record is None or st.is_circular(record)
        self.map_stack.setCurrentWidget(self.circular_panel if circular else self.linear_panel)

    def set_map_mode(self, mode):
        self.map_mode = mode
        self.settings.setValue("map_mode", mode)
        self.update_map_style()

    def set_theme_mode(self, mode):
        self.theme_mode = mode
        self.settings.setValue("theme", mode)
        self.apply_theme()

    def on_system_color_scheme(self, *_):
        if self.theme_mode == "system":
            self.apply_theme()

    def apply_theme(self):
        theme.apply_theme(QApplication.instance(), self.theme_mode)
        self.viewer.refresh_theme()
        self.maps.update()
        self.fill_feature_table()

    def copy_seq(self):
        if self.current_record is None:
            return
        self.copy_text(self.current_sequence, "Sequence copied to clipboard")

    def show_reverse_complement(self):
        if not self.require_dna("Reverse complement"):
            return
        self.copy_text(st.reverse_complement(self.current_sequence),
                       "Reverse complement copied to clipboard")

    def show_translation(self):
        if not self.require_dna("Translation"):
            return
        self.show_protein(self.current_sequence,
                          f"{self.current_record.id}, whole sequence, frame +1")

    def show_about(self):
        QMessageBox.about(
            self, "About N-Zoom",
            f"N-Zoom {VERSION}\nA lightweight sequence viewer by mikeph52.",
        )


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("N-Zoom")
    window = MainWindow()
    window.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
