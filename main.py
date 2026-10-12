# N-Zoom by mikeph52, Init 12/9/2026
"""
main.py — the N-Zoom main window. Run with:  python main.py
"""
import sys
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq
from PySide6.QtCore import QSettings, Qt, QThreadPool
from PySide6.QtGui import QAction
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QListWidget,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from core import seqtools as st
from core.viewers import CircularDNAViewer, SequenceViewer
from core.workers import Task

VERSION = "0.8.1"

OPEN_FILTER = (
    "Sequence files (*.gb *.gbk *.genbank *.gbff *.fasta *.fa *.faa *.fna *.fas);;"
    "GenBank files (*.gb *.gbk *.genbank *.gbff);;"
    "FASTA files (*.fasta *.fa *.faa *.fna *.fas);;"
    "All files (*)"
)
SAVE_FILTER = "GenBank (*.gb);;FASTA (*.fasta)"


def number_item(value):
    """Table cell that sorts numerically."""
    item = QTableWidgetItem()
    item.setData(Qt.ItemDataRole.DisplayRole, value)
    return item


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle(f"N-Zoom {VERSION}")
        self.resize(1300, 820)

        self.settings = QSettings("mikeph52", "N-Zoom")
        self.records = []
        self.current_record = None
        self.current_sequence = ""
        self.current_features = []

        # background work (NCBI downloads, restriction searches)
        self.pool = QThreadPool(self)
        self.tasks = {}           # kind -> running Task, e.g. {"ncbi": task}
        self.sites_dialog = None

        self.build_ui()
        self.build_menu()

    # ================================================================ UI
    def build_menu(self):
        file_menu = self.menuBar().addMenu("File")

        open_action = QAction("Open...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        save_action = QAction("Save As...", self)
        save_action.setShortcut("Ctrl+S")
        save_action.triggered.connect(self.save_as)
        file_menu.addAction(save_action)
        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+Q")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        tools_menu = self.menuBar().addMenu("Tools")

        self.ncbi_action = QAction("Fetch from NCBI...", self)
        self.ncbi_action.triggered.connect(self.fetch_from_ncbi)
        tools_menu.addAction(self.ncbi_action)
        tools_menu.addSeparator()

        translate_action = QAction("Translate sequence", self)
        translate_action.triggered.connect(self.show_translation)
        tools_menu.addAction(translate_action)

        self.rs_action = QAction("Find restriction sites", self)
        self.rs_action.triggered.connect(self.find_restriction_sites)
        tools_menu.addAction(self.rs_action)

        clear_rs_action = QAction("Clear restriction sites", self)
        clear_rs_action.triggered.connect(lambda: self.circular_viewer.set_sites({}))
        tools_menu.addAction(clear_rs_action)

        about_menu = self.menuBar().addMenu("About")
        about_action = QAction("About", self)
        about_action.triggered.connect(self.show_about)
        about_menu.addAction(about_action)

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

        center_split = QSplitter(Qt.Orientation.Vertical)

        self.view_tabs = QTabWidget()
        self.circular_viewer = CircularDNAViewer()
        self.circular_viewer.featureClicked.connect(self.select_feature_row)
        self.view_tabs.addTab(self.circular_viewer, "Map")
        self.viewer = SequenceViewer()
        self.viewer.selectionInfo.connect(self.on_sequence_selection)
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
        self.feature_table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self.feature_table.itemSelectionChanged.connect(self.on_feature_selected)
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
        for text, slot in (
            ("Copy sequence", self.copy_seq),
            ("Reverse complement", self.show_reverse_complement),
            ("Translate DNA", self.show_translation),
        ):
            button = QPushButton(text)
            button.clicked.connect(slot)
            actions_layout.addWidget(button)
        self.rs_button = QPushButton("Find restriction sites")
        self.rs_button.clicked.connect(self.find_restriction_sites)
        actions_layout.addWidget(self.rs_button)
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

    # ================================================================ loading
    def open_file(self):
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

        self.show_records(records)
        self.file_label.setText(str(Path(path)))
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
        # added to whatever is already open
        self.show_records(self.records + [record], select=len(self.records))
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

    def save_as(self):
        if not self.records:
            QMessageBox.information(self, "Nothing to save", "Open or fetch a sequence first.")
            return
        path, chosen = QFileDialog.getSaveFileName(self, "Save sequences", "", SAVE_FILTER)
        if not path:
            return
        fmt = "genbank" if chosen.startswith("GenBank") else "fasta"
        if not Path(path).suffix:
            path += ".gb" if fmt == "genbank" else ".fasta"
        try:
            st.save_records(self.records, path, fmt)
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Could not save file:\n{exc}")
            return
        self.statusBar().showMessage(f"Saved {len(self.records)} record(s) to {path}", 5000)

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

        self.sequence_title.setText(record.description or record.id)
        self.viewer.set_sequence(sequence, seq_type)

        # map only makes sense for DNA
        if seq_type == "DNA":
            self.circular_viewer.set_record(len(sequence), record.id, self.current_features)
            self.view_tabs.setTabEnabled(0, True)
        else:
            self.circular_viewer.clear_seq()
            self.view_tabs.setTabEnabled(0, False)
            self.view_tabs.setCurrentIndex(1)

        self.fill_feature_table()

        # info panel
        self.type_value.setText(seq_type)
        self.length_value.setText(f"{len(sequence):,}")
        self.topology_value.setText(st.topology(record))
        self.features_value.setText(str(len(self.current_features)))

        counts = Counter(sequence)
        if seq_type == "DNA" and sequence:
            gc = (counts.get("G", 0) + counts.get("C", 0)) / len(sequence) * 100
            self.gc_value.setText(f"{gc:.2f}%")
        else:
            self.gc_value.setText("N/A")
        composition = ", ".join(f"{base}: {count:,}" for base, count in counts.most_common(12))
        self.composition_value.setText(composition or "-")

    def fill_feature_table(self):
        table = self.feature_table
        table.blockSignals(True)
        table.setSortingEnabled(False)
        table.setRowCount(len(self.current_features))

        strand_text = {1: "→", -1: "←", 0: "·"}
        for row, feature in enumerate(self.current_features):
            parts = st.feature_parts(feature)
            start, end = st.feature_span(feature)  # 1-based like GenBank
            length = sum(e - s for s, e in parts)

            type_item = QTableWidgetItem(feature.type)
            type_item.setData(Qt.ItemDataRole.UserRole, row)  # index into current_features
            table.setItem(row, 0, type_item)
            table.setItem(row, 1, QTableWidgetItem(st.feature_label(feature)))
            table.setItem(row, 2, number_item(start))
            table.setItem(row, 3, number_item(end))
            table.setItem(row, 4, QTableWidgetItem(strand_text[st.feature_strand(feature)]))
            table.setItem(row, 5, number_item(length))

        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.horizontalHeader().setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        table.blockSignals(False)

    def on_feature_selected(self):
        rows = self.feature_table.selectionModel().selectedRows()
        if not rows:
            self.viewer.setExtraSelections([])
            self.circular_viewer.set_selected(-1)
            return
        index = self.feature_table.item(rows[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        feature = self.current_features[index]
        self.viewer.highlight_parts(st.feature_parts(feature))
        self.circular_viewer.set_selected(index)
        self.statusBar().showMessage(
            f"{feature.type}: {st.feature_label(feature)}  ({feature.location})", 8000
        )

    def select_feature_row(self, index):
        """Called when a feature is clicked on the map."""
        for row in range(self.feature_table.rowCount()):
            if self.feature_table.item(row, 0).data(Qt.ItemDataRole.UserRole) == index:
                self.feature_table.selectRow(row)
                self.feature_table.scrollToItem(self.feature_table.item(row, 0))
                return

    def on_sequence_selection(self, start, end):
        """Show the cursor position or selected range in the status bar."""
        sequence = self.current_sequence
        if not sequence:
            self.position_label.setText("")
            self.circular_viewer.set_highlight()
            return

        is_dna = st.detect_seq_type(sequence) == "DNA"
        unit = "bp" if is_dna else "aa"
        if start == end:
            position = min(start + 1, len(sequence))
            self.position_label.setText(f"Position {position:,} of {len(sequence):,}")
            self.circular_viewer.set_highlight()
            return

        text = f"Selected {start + 1:,}\u2013{end:,}  ({end - start:,} {unit})"
        if is_dna:
            region = sequence[start:end]
            gc = (region.count("G") + region.count("C")) / len(region) * 100
            text += f"  \u00b7  GC {gc:.1f}%"
            self.circular_viewer.set_highlight(start, end)
        self.position_label.setText(text)

    # ================================================================ background tasks
    def start_task(self, kind, message, fn, *args, on_done=None, on_fail=None, **kwargs):
        """Run fn(*args, **kwargs) on a background thread.

        `kind` names the job ("ncbi", "restriction"); only one job of each kind
        runs at a time. on_done(result) / on_fail(message) run in the GUI thread.
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
        ):
            for widget in widgets:
                widget.setEnabled(kind not in self.tasks)

    def closeEvent(self, event):
        if self.tasks:
            # can't interrupt a download midway; wait for it so we exit cleanly
            self.statusBar().showMessage("Finishing background work before closing...")
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            self.pool.waitForDone()
            QApplication.restoreOverrideCursor()
        super().closeEvent(event)

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
        if record is not self.current_record:
            # the user switched records while the search was running
            self.statusBar().showMessage(
                f"Restriction search for {record.id} finished, but another record "
                "is selected now. Run it again to see the results.", 8000
            )
            return

        on_map = st.map_sites(sites)
        self.circular_viewer.set_sites(on_map)
        self.view_tabs.setCurrentIndex(0)
        self.statusBar().showMessage(
            f"{len(sites)} enzymes cut; {len(on_map)} unique 6+ bp sites shown on map", 8000
        )
        self.show_sites_dialog(sites, circular)

    def show_sites_dialog(self, sites, circular):
        dialog = QDialog(self)
        dialog.setWindowTitle("Restriction sites")
        dialog.resize(520, 560)
        layout = QVBoxLayout(dialog)
        layout.addWidget(QLabel(
            f"Commercially available enzymes, {'circular' if circular else 'linear'} sequence. "
            "Unique cutters with 6+ bp sites are drawn on the map."
        ))

        table = QTableWidget(len(sites), 3)
        table.setHorizontalHeaderLabels(["Enzyme", "Cuts", "Positions"])
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.verticalHeader().setVisible(False)
        ordered = sorted(sites.items(), key=lambda item: (len(item[1]), item[0]))
        for row, (name, positions) in enumerate(ordered):
            shown = ", ".join(f"{p:,}" for p in positions[:12])
            if len(positions) > 12:
                shown += ", ..."
            table.setItem(row, 0, QTableWidgetItem(name))
            table.setItem(row, 1, number_item(len(positions)))
            table.setItem(row, 2, QTableWidgetItem(shown))
        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(table)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(dialog.reject)
        layout.addWidget(buttons)

        # non-modal: the main window stays usable while the list is open
        if self.sites_dialog is not None:
            self.sites_dialog.close()
        self.sites_dialog = dialog
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        dialog.destroyed.connect(self.forget_sites_dialog)
        dialog.show()

    def forget_sites_dialog(self):
        self.sites_dialog = None

    def copy_seq(self):
        if self.current_record is None:
            return
        QApplication.clipboard().setText(self.current_sequence)
        self.statusBar().showMessage("Sequence copied to clipboard", 3000)

    def show_reverse_complement(self):
        if not self.require_dna("Reverse complement"):
            return
        result = str(Seq(self.current_sequence).reverse_complement())
        QApplication.clipboard().setText(result)
        QMessageBox.information(self, "Reverse complement", "Reverse complement copied to clipboard.")

    def show_translation(self):
        if not self.require_dna("Translation"):
            return
        clean = self.current_sequence.replace("-", "").replace("U", "T")
        usable_length = len(clean) - (len(clean) % 3)
        if usable_length == 0:
            QMessageBox.warning(self, "Too short", "The sequence is too short to translate.")
            return

        protein = str(Seq(clean[:usable_length]).translate())
        QApplication.clipboard().setText(protein)
        preview = protein[:500] + ("\n..." if len(protein) > 500 else "")
        QMessageBox.information(
            self, "Translation", f"Protein translation copied to clipboard.\n\n{preview}"
        )

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