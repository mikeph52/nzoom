# N-Zoom by mikeph52, Init 12/9/2026
import sys
import subprocess
from pathlib import Path
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor, QTextCharFormat, QFont
from PySide6.QtWidgets import (
    QApplication,
    QMainWindow,
    QFileDialog,
    QListWidget,
    QTextEdit,
    QWidget,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSplitter,
    QMessageBox,
    QGroupBox,
    QFormLayout,
)

DNA_CHARS = set("ACGTUNRYKMSWBDHV-")
PROTEIN_CHARS = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-")
VERSION = set("0.4.1") # for build version

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

def normalize_seq(seq):
    return "".join(seq.split()).upper() # remove any whitespace and force CASES.

# detect if file contains dna or protein,!IMPORTANT! RNA returns as DNA.
def detect_seq_type(seq):
    seq = normalize_seq(seq)
    chars = set(seq)

    if chars <= DNA_CHARS:
        return "DNA"
    if chars <= PROTEIN_CHARS:
        return "Protein"
    else:
        return "Uknown"

# Sequence Viewer function
class SequenceViewer(QTextEdit):
    def __int__(self):
        super().__init__()
        self.setReadOnly(True) # Maybe change letter for editing mode.
        self.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        font = QFont("Monospace") # font settings
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)

    def set_sequence(self, sequence, seq_type):
        self.clear()
        colors = DNA_COLORS if seq_type == "DNA" else PROTEIN_COLORS
        cursor = self.textCursor()
        chunk_size = 10000

        for offset in range(0, len(sequence), chunk_size):
            chunk = sequence[offset:offset + chunk_size]
            for char in chunk:
                fmt = QTextCharFormat()
                fmt.setForeground(QColor(colors.get(char, "#424242")))
                cursor.insertText(char, fmt)
            if offset + chunk_size < len(sequence):
                cursor.insertText("\n")

        self.moveCursor(self.textCursor().MoveOperation.Start)
# Main GUI Window
class MainWindow(QMainWindow):
    def __int__(self):
        super().__init__()

        # window elements
        self.setWindowTitle("NZoom Alpha Test")
        self.resize(1200, 750)

        self.records = []
        self.current_record = None

        self.build_ui()
        self.build_menu()

    def build_menu(self):
        file_menu = self.menuBar().addMenu("&File")

        open_action = QAction("&Open FASTA...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_fasta)
        file_menu.addAction(open_action)
        file_menu.addSeparator()

        exit_action = QAction("E&xit", self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

    def build_ui(self):
        root = QWidget()
        root_layout = QVBoxLayout(root)
        top_bar = QHBoxLayout()
        open_button = QPushButton("Open FASTA")
        open_button.clicked.connect(self.open_fasta)
        top_bar.addWidget(open_button)

        self.file_label = QLabel("No file loaded")
        self.file_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        top_bar.addWidget(self.file_label, 1)
        root_layout.addLayout(top_bar)
        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.sequence_list = QListWidget()
        self.sequence_list.setMinimumWidth(230)
        self.sequence_list.currentRowChanged.connect(self.select_record)
        splitter.addWidget(self.sequence_list)

        center = QWidget()
        center_layout = QVBoxLayout(center)

        self.sequence_title = QLabel("Sequence")
        title_font = self.sequence_title.font()
        title_font.setBold(True)
        title_font.setPointSize(13)
        self.sequence_title.setFont(title_font)
        center_layout.addWidget(self.sequence_title)
        self.viewer = SequenceViewer()
        center_layout.addWidget(self.viewer, 1)
        splitter.addWidget(center)
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        info_box = QGroupBox("Sequence information")
        form = QFormLayout(info_box)

        self.type_value = QLabel("-")
        self.length_value = QLabel("-")
        self.gc_value = QLabel("-")
        self.composition_value = QLabel("-")
        self.composition_value.setWordWrap(True)
        form.addRow("Type:", self.type_value)
        form.addRow("Length:", self.length_value)
        form.addRow("GC %:", self.gc_value)
        form.addRow("Composition:", self.composition_value)

        right_layout.addWidget(info_box)
        actions_box = QGroupBox("Actions")
        actions_layout = QVBoxLayout(actions_box)

        copy_button = QPushButton("Copy sequence")
        copy_button.clicked.connect(self.copy_sequence)
        actions_layout.addWidget(copy_button)

        reverse_button = QPushButton("Reverse complement")
        reverse_button.clicked.connect(self.show_reverse_complement)
        actions_layout.addWidget(reverse_button)

        translate_button = QPushButton("Translate DNA")
        translate_button.clicked.connect(self.show_translation)
        actions_layout.addWidget(translate_button)

        right_layout.addWidget(actions_box)
        right_layout.addStretch()

        splitter.addWidget(right_panel)
        splitter.setSizes([240, 700, 260])

        root_layout.addWidget(splitter, 1)

        self.setCentralWidget(root)
        self.statusBar().showMessage("Ready")

    def open_fasta(self):
        path, _ = QFileDialog.getOpenFileName(self,"Open FASTA file","","FASTA files (*.fasta *.fa *.faa *.fna *.fas);;All files (*)")

        if not path:
            return
        try:
            records = list(SeqIO.parse(path, "fasta"))
        except Exception as exc:
            QMessageBox.critical(self, "Error", f"Could not read FASTA file:\n{exc}")
            return
        if not records:
            QMessageBox.warning(self, "No sequences", "No FASTA records were found.")
            return

        self.records = records
        self.current_record = None
        self.sequence_list.clear()

        for record in records:
            label = record.id
            if record.description and record.description != record.id:
                label += f" — {record.description}"
            self.sequence_list.addItem(label)

        self.file_label.setText(str(Path(path)))
        self.statusBar().showMessage(f"Loaded {len(records)} sequence(s)")
        self.sequence_list.setCurrentRow(0)





























def main():
    app = QApplication

if __name__ == "__main__":
    main()