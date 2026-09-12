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
        root = 









def main():
    app = QApplication

if __name__ == "__main__":
    main()