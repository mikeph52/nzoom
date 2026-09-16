# N-Zoom by mikeph52, Init 12/9/2026
import sys
import math
from pathlib import Path
from collections import Counter

from Bio import SeqIO
from Bio.Seq import Seq
from PySide6.QtCore import Qt
from PySide6.QtGui import QAction, QColor, QTextCharFormat, QFont, QPainter, QPen
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
    QTableWidget,
    QSizePolicy,
)

DNA_CHARS = set("ACGTUNRYKMSWBDHV-")
PROTEIN_CHARS = set("ABCDEFGHIKLMNPQRSTVWXYZJUO*-")
VERSION = "0.4.1" # for build version

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
        return "Unknown"

# Sequence Viewer function
class SequenceViewer(QTextEdit):
    def __init__(self):
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

# CircularDNAViwer is a class for implementing circular DNA viwer in the app. Experimental feature.
class CircularDNAViewer(QWidget):
    def __init__(self):
        super().__init__()

        self.sequence = ""
        self.title = "Circular DNA"
        self.setMinimumSize(400, 400)

    def set_seq(self, sequence, title="Circular DNA"):
        self.sequence = normalize_seq(sequence)
        self.title = title
        self.update()

    def clear_seq(self):
        self.sequence = ""
        self.update()

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        # set background
        painter.fillRect(self.rect(), QColor("#ffffff"))

        if not self.sequence:
            painter.setPen(QColor("#777777"))
            painter.drawText(
                self.rect(),
                Qt.AlignmentFlag.AlignCenter,
                "No circular DNA loaded",
            )
            return

        width = self.width()
        height = self.height()

        center_x = width / 2
        center_y = height / 2

        radius = min(width, height) * 0.32
        # draw dna
        pen = QPen(QColor("#333333"))
        pen.setWidth(5)
        painter.setPen(pen)

        painter.drawEllipse(
            int(center_x - radius),
            int(center_y - radius),
            int(radius * 2),
            int(radius * 2),
        )
        # tickmarks
        sequence_length = len(self.sequence)

        if sequence_length <= 1000:
            tick_interval = 100
        elif sequence_length <= 5000:
            tick_interval = 500
        elif sequence_length <= 10000:
            tick_interval = 1000
        elif sequence_length <= 50000:
            tick_interval = 5000
        else:
            tick_interval = 10000

        for position in range(0, sequence_length, tick_interval):

            fraction = position / sequence_length
            angle = fraction * 360 - 90
            angle_rad = math.radians(angle)
            inner_radius = radius - 8
            outer_radius = radius + 8

            x1 = center_x + math.cos(angle_rad) * inner_radius
            y1 = center_y + math.sin(angle_rad) * inner_radius
            x2 = center_x + math.cos(angle_rad) * outer_radius
            y2 = center_y + math.sin(angle_rad) * outer_radius

            tick_pen = QPen(QColor("#555555"))
            tick_pen.setWidth(2)

            painter.setPen(tick_pen)
            painter.drawLine(
                int(x1),
                int(y1),
                int(x2),
                int(y2),
            )
            #tick labels
            label_radius = radius + 30
            label_x = center_x + math.cos(angle_rad) * label_radius
            label_y = center_y + math.sin(angle_rad) * label_radius

            painter.setPen(QColor("#333333"))
            painter.drawText(
                int(label_x - 30),
                int(label_y - 10),
                60,
                20,
                Qt.AlignmentFlag.AlignCenter,
                f"{position:,}",
            )
            # center title
            title_font = QFont()
            title_font.setPointSize(12)
            title_font.setBold(True)

            painter.setFont(title_font)
            painter.setPen(QColor("#222222"))
            painter.drawText(
                int(center_x - radius),
                int(center_y - 35),
                int(radius * 2),
                30,
                Qt.AlignmentFlag.AlignCenter,
                self.title,
            )
            # seq length
            length_font = QFont()
            length_font.setPointSize(11)

            painter.setFont(length_font)
            painter.setPen(QColor("#555555"))

            painter.drawText(
                int(center_x - radius),
                int(center_y + 5),
                int(radius * 2),
                30,
                Qt.AlignmentFlag.AlignCenter,
                f"{sequence_length:,} bp",
            )
# Main GUI Window
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        # window elements
        self.setWindowTitle("NZoom Alpha Test")
        self.resize(1200, 750)

        self.records = []
        self.current_record = None

        self.build_ui()
        self.build_menu()

    def build_menu(self):
        # file submenu
        file_menu = self.menuBar().addMenu("File")

        open_action = QAction("Open FASTA...", self)
        open_action.setShortcut("Ctrl+O")
        open_action.triggered.connect(self.open_fasta)
        file_menu.addAction(open_action)
        file_menu.addSeparator()

        exit_action = QAction("Exit", self)
        exit_action.setShortcut("Ctrl+C")
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

        # tools submenu
        tools_menu = self.menuBar().addMenu("Tools")

        ncbi_action = QAction("Fetch data from NCBI", self)
        tools_menu.addAction(ncbi_action)
        tools_menu.addSeparator()

        translate_action = QAction("Translate sequence", self)
        tools_menu.addAction(translate_action)
        tools_menu.addSeparator()

        rs_action = QAction("Find restriction sites", self)
        tools_menu.addAction(rs_action)

        # about submenu
        about_menu = self.menuBar().addMenu("About")
        about_action = QAction("About", self)
        about_menu.addAction(about_action)
        about_menu.addSeparator()

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
        copy_button.clicked.connect(self.copy_seq)
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

    def select_record(self, index: int):
        if index < 0 or index >= len(self.records):
            return

        record = self.records[index]
        self.current_record = record

        sequence = normalize_seq(str(record.seq))
        seq_type = detect_seq_type(sequence)

        self.sequence_title.setText(record.description or record.id)
        self.viewer.set_sequence(sequence, seq_type)
        self.type_value.setText(seq_type)
        self.length_value.setText(f"{len(sequence):,}")

        counts = Counter(sequence)
        total = len(sequence)

        if seq_type == "DNA" and total:
            gc = ((counts.get("G", 0) + counts.get("C", 0)) / total) * 100
            self.gc_value.setText(f"{gc:.2f}%")
        else:
            self.gc_value.setText("N/A")

        common = counts.most_common(12)
        composition = ", ".join(f"{base}: {count:,}" for base, count in common)
        self.composition_value.setText(composition or "-")

    def copy_seq(self):
        if self.current_record is None:
            return

        QApplication.clipboard().setText(normalize_seq(str(self.current_record.seq)))
        self.statusBar().showMessage("Sequence copied to clipboard", 3000)

    def show_reverse_complement(self):
        if self.current_record is None:
            return

        sequence = normalize_seq(str(self.current_record.seq))
        if detect_seq_type(sequence) != "DNA":
            QMessageBox.information(self,"DNA required","Reverse complement is only available for DNA sequences.")
            return

        result = str(Seq(sequence).reverse_complement())
        dialog = QMessageBox(self)
        dialog.setWindowTitle("Reverse complement")
        dialog.setText("Reverse complement copied to clipboard.")
        QApplication.clipboard().setText(result)
        dialog.exec()

    def show_translation(self):
        if self.current_record is None:
            return

        sequence = normalize_seq(str(self.current_record.seq))

        if detect_seq_type(sequence) != "DNA":
            QMessageBox.information(
                self,
                "DNA required",
                "Translation is only available for DNA sequences.",
            )
            return
        clean = sequence.replace("-", "").replace("U", "T")
        usable_length = len(clean) - (len(clean) % 3)
        if usable_length == 0:
            QMessageBox.warning(self, "Too short", "The sequence is too short to translate.")
            return

        protein = str(Seq(clean[:usable_length]).translate())
        QApplication.clipboard().setText(protein)

        preview = protein[:500]
        if len(protein) > 500:
            preview += "\n..."

        QMessageBox.information(
            self,
            "Translation",
            f"Protein translation copied to clipboard.\n\n{preview}",
        )

def main():
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()