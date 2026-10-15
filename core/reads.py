"""
reads.py — the "Sequencing reads" window: Sanger reads aligned to the
current sequence, their differences, the alignment text and chromatograms.

TraceView    draws an .ab1 chromatogram (four colored traces + base calls)
ReadsDialog  table of reads, table of differences, alignment and trace tabs
"""
from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHeaderView,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from . import seqtools as st
from . import theme
from .dialogs import monospace_font, number_item

TRACE_COLORS_LIGHT = {"A": "#2e7d32", "C": "#1565c0", "G": "#212121", "T": "#c62828", "N": "#9e9e9e"}
TRACE_COLORS_DARK = {"A": "#66bb6a", "C": "#42a5f5", "G": "#e0e0e0", "T": "#ef5350", "N": "#9e9e9e"}


def status_of(result):
    """Short verdict for a read: matches / differences / no match."""
    alignment = result.get("alignment")
    if alignment is None:
        return "no match"
    real = [d for d in alignment["differences"] if not d["low_quality"]]
    if not real:
        return "matches"
    return f"{len(real)} difference{'s' if len(real) != 1 else ''}"


# =================================================================== chromatogram
class TraceView(QWidget):
    """A chromatogram. Lives in a QScrollArea: it is as wide as the whole trace."""
    X_SCALE = 2.0  # pixels per trace point

    def __init__(self):
        super().__init__()
        self.read = None
        self.marked = -1  # base index to highlight
        self.setMinimumHeight(190)

    def set_read(self, read):
        self.read = read
        self.marked = -1
        width = 200
        if read and read.get("trace"):
            width = int(len(next(iter(read["trace"]["channels"].values()))) * self.X_SCALE) + 40
        self.setFixedWidth(width)
        self.update()

    def peak_x(self, index):
        return 20 + self.read["trace"]["peaks"][index] * self.X_SCALE

    def mark(self, index):
        self.marked = index
        self.update()

    def sizeHint(self):
        return QSize(self.width(), 220)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        colors = theme.current()
        painter.fillRect(event.rect(), QColor(colors.background))
        read = self.read
        if not read or not read.get("trace"):
            painter.setPen(QColor(colors.muted))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "No chromatogram: this read has no trace data (only .ab1 files do).")
            return
        palette = TRACE_COLORS_DARK if colors.dark else TRACE_COLORS_LIGHT
        trace = read["trace"]
        channels, peaks = trace["channels"], trace["peaks"]
        top, base_line = 34, self.height() - 26
        tallest = max(max(values) for values in channels.values()) or 1
        visible = event.rect()
        first = max(0, int((visible.left() - 20) / self.X_SCALE) - 2)
        last = int((visible.right() - 20) / self.X_SCALE) + 2

        # trimmed (unused) ends in grey
        trim_start, trim_end = read.get("trim", (0, len(read["seq"])))
        shade = QColor(colors.muted)
        shade.setAlpha(40)
        if trim_start > 0:
            painter.fillRect(QRectF(0, 0, self.peak_x(trim_start) - 6, self.height()), shade)
        if trim_end < len(peaks):
            x = self.peak_x(trim_end - 1) + 6
            painter.fillRect(QRectF(x, 0, self.width() - x, self.height()), shade)

        # highlighted base
        if 0 <= self.marked < len(peaks):
            x = self.peak_x(self.marked)
            painter.fillRect(QRectF(x - 7, 0, 14, self.height()), QColor(*colors.highlight))

        # the four traces
        for base, values in channels.items():
            path = QPainterPath()
            end = min(last, len(values))
            for i in range(first, end):
                point = QPointF(20 + i * self.X_SCALE,
                                base_line - values[i] / tallest * (base_line - top))
                if i == first:
                    path.moveTo(point)
                else:
                    path.lineTo(point)
            painter.setPen(QPen(QColor(palette.get(base, palette["N"])), 1.2))
            painter.drawPath(path)

        # base calls, quality bars and base numbers
        font = QFont("Monospace")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(9)
        painter.setFont(font)
        metrics = QFontMetrics(font)
        quals = read.get("quals") or []
        for i, peak in enumerate(peaks):
            x = self.peak_x(i)
            if x < visible.left() - 10 or x > visible.right() + 10:
                continue
            base = read["seq"][i] if i < len(read["seq"]) else "N"
            if i < len(quals):  # grey bar: taller = better quality
                h = min(quals[i], 60) / 60 * 18
                bar = QColor(colors.muted)
                bar.setAlpha(90)
                painter.fillRect(QRectF(x - 3, top - 4 - h, 6, h), bar)
            painter.setPen(QColor(palette.get(base, palette["N"])))
            painter.drawText(QRectF(x - 6, 2, 12, metrics.height()), Qt.AlignmentFlag.AlignCenter, base)
            if (i + 1) % 10 == 0:
                painter.setPen(QColor(colors.muted))
                painter.drawText(QRectF(x - 20, self.height() - 20, 40, 16),
                                 Qt.AlignmentFlag.AlignCenter, str(i + 1))
        painter.end()


class TraceScrollArea(QScrollArea):
    """Scrolls the chromatogram sideways; its height follows the window."""

    def resizeEvent(self, event):
        super().resizeEvent(event)
        if self.widget() is not None:
            self.widget().setFixedHeight(max(190, self.viewport().height()))


# =================================================================== dialog
class ReadsDialog(QDialog):
    """Non-modal window listing aligned reads.

    on_jump(start, end): select a region (0-based, end exclusive) in the main window.
    on_add(): ask the main window to align more reads.
    """
    READ_COLUMNS = ["Read", "Strand", "Covers", "Aligned", "Identity", "Result"]
    DIFF_COLUMNS = ["Position", "Type", "Reference", "Read", "Quality", "Note"]

    def __init__(self, parent, title, results, on_jump, on_add):
        super().__init__(parent)
        self.setWindowTitle(f"Sequencing reads — {title}")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(1050, 720)
        self.results = []
        self.on_jump = on_jump

        layout = QVBoxLayout(self)
        self.summary = QLabel()
        self.summary.setWordWrap(True)
        layout.addWidget(self.summary)

        splitter = QSplitter(Qt.Orientation.Vertical)
        self.read_table = self._table(self.READ_COLUMNS)
        self.read_table.itemSelectionChanged.connect(self.show_selected_read)
        splitter.addWidget(self.read_table)

        lower = QSplitter(Qt.Orientation.Horizontal)
        self.diff_table = self._table(self.DIFF_COLUMNS)
        self.diff_table.itemSelectionChanged.connect(self.show_selected_difference)
        lower.addWidget(self.diff_table)

        self.tabs = QTabWidget()
        self.trace = TraceView()
        self.trace_scroll = TraceScrollArea()
        self.trace_scroll.setWidget(self.trace)
        self.trace_scroll.setWidgetResizable(False)
        self.trace_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.tabs.addTab(self.trace_scroll, "Chromatogram")
        self.alignment_text = QPlainTextEdit()
        self.alignment_text.setReadOnly(True)
        self.alignment_text.setFont(monospace_font(10))
        self.alignment_text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self.tabs.addTab(self.alignment_text, "Alignment")
        lower.addWidget(self.tabs)
        lower.setSizes([420, 630])
        splitter.addWidget(lower)
        splitter.setSizes([200, 480])
        layout.addWidget(splitter, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        add_button = QPushButton("Add reads...")
        add_button.clicked.connect(on_add)
        buttons.addButton(add_button, QDialogButtonBox.ButtonRole.ActionRole)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

        self.set_results(results)

    @staticmethod
    def _table(columns):
        table = QTableWidget(0, len(columns))
        table.setHorizontalHeaderLabels(columns)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        table.horizontalHeader().setStretchLastSection(True)
        return table

    # ------------------------------------------------------------------ reads
    def set_results(self, results, select=0):
        self.results = list(results)
        table = self.read_table
        table.blockSignals(True)
        table.setRowCount(len(self.results))
        for row, result in enumerate(self.results):
            alignment = result.get("alignment")
            status = status_of(result)
            if alignment:
                strand = "→ forward" if alignment["strand"] == 1 else "← reverse"
                covers = f"{alignment['start']:,}–{alignment['end']:,}"
                cells = [result["name"], strand, covers, number_item(alignment["columns"]),
                         f"{alignment['identity']:.1f}%", status]
            else:
                cells = [result["name"], "", "", "", "", result.get("error") or status]
            for column, cell in enumerate(cells):
                item = cell if isinstance(cell, QTableWidgetItem) else QTableWidgetItem(cell)
                if column == 5:
                    color = {"matches": "#2e7d32", "no match": "#c62828"}.get(status, "#ef6c00")
                    item.setForeground(QColor(color))
                table.setItem(row, column, item)
        table.blockSignals(False)
        good = sum(1 for r in self.results if status_of(r) == "matches")
        aligned = sum(1 for r in self.results if r.get("alignment"))
        self.summary.setText(
            f"{len(self.results)} read{'s' if len(self.results) != 1 else ''}: {aligned} aligned, "
            f"{good} without high-quality differences. Differences at bases with quality "
            f"below {st.LOW_QUALITY} are marked “low quality”: they are probably "
            "sequencing errors, not real changes. Click a difference to see it in the "
            "chromatogram and the main window.")
        if self.results:
            table.selectRow(min(select, len(self.results) - 1))
        else:
            self.show_selected_read()

    def selected_result(self):
        rows = self.read_table.selectionModel().selectedRows()
        return self.results[rows[0].row()] if rows else None

    def select_read(self, name):
        for row, result in enumerate(self.results):
            if result["name"] == name:
                self.read_table.selectRow(row)
                return

    def show_selected_read(self):
        result = self.selected_result()
        alignment = result.get("alignment") if result else None
        self.trace.set_read(None)
        self.diff_table.setRowCount(0)
        if result is None:
            self.alignment_text.setPlainText("")
            return
        read_view = dict(result)
        read_view["trim"] = alignment["trim"] if alignment else (0, len(result["seq"]))
        self.trace.set_read(read_view)
        if alignment is None:
            self.alignment_text.setPlainText(
                result.get("error") or "This read doesn't match the sequence.")
            return
        self.alignment_text.setPlainText(st.format_alignment(alignment["rows"]))

        differences = alignment["differences"]
        table = self.diff_table
        table.blockSignals(True)
        table.setRowCount(len(differences))
        for row, d in enumerate(differences):
            quality = "" if d["quality"] is None else str(d["quality"])
            note = "low quality" if d["low_quality"] else ("" if d["quality"] is None else "check")
            position = f"{d['ref_pos']:,}" + ("½" if d["kind"] == "insertion" else "")
            cells = [position, d["kind"], d["ref"], d["read"], quality, note]
            for column, text in enumerate(cells):
                item = QTableWidgetItem(text)
                if d["low_quality"]:
                    item.setForeground(QColor("#9e9e9e"))
                table.setItem(row, column, item)
        table.blockSignals(False)
        # scroll the chromatogram to where the read starts being used
        if result.get("trace"):
            self.trace_scroll.horizontalScrollBar().setValue(
                max(0, int(self.trace.peak_x(alignment["trim"][0])) - 40))

    def show_selected_difference(self):
        result = self.selected_result()
        rows = self.diff_table.selectionModel().selectedRows()
        if not result or not result.get("alignment") or not rows:
            return
        d = result["alignment"]["differences"][rows[0].row()]
        if result.get("trace"):
            self.trace.mark(d["read_index"])
            self.tabs.setCurrentWidget(self.trace_scroll)
            x = int(self.trace.peak_x(d["read_index"]))
            self.trace_scroll.horizontalScrollBar().setValue(
                max(0, x - self.trace_scroll.viewport().width() // 2))
        start = d["ref_pos"] - 1
        length = len(d["ref"]) if d["kind"] == "deletion" else 1
        self.on_jump(start, start + length)
