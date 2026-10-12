"""
digest.py — virtual restriction digest with a simulated agarose gel.

GelView       draws a ladder and one lane per digest
DigestDialog  pick enzymes, add lanes, see fragment sizes, save the gel image

Band positions use the usual approximation that DNA moves a distance
proportional to -log(size). Real gels also depend on agarose %, voltage and
DNA topology (uncut supercoiled plasmids run faster than their size suggests).
"""
import math

from PySide6.QtCore import QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from . import seqtools as st
from .dialogs import number_item

MAX_LANES = 10
CUT_FILTERS = [("Enzymes that cut once", 1), ("Cut 1–2 times", 2),
               ("Cut 1–3 times", 3), ("All enzymes that cut", 10**9)]


# =================================================================== gel
class GelView(QWidget):
    """Simulated agarose gel. Gels look the same in light and dark mode."""
    laneClicked = Signal(int)  # sample lane index (0 = first lane after the ladder)

    BACKGROUND = QColor("#14181c")
    GEL = QColor("#1f2c35")
    WELL = QColor("#0b0f12")
    BAND = (225, 245, 255)
    TEXT = QColor("#c9d6df")

    def __init__(self):
        super().__init__()
        self.ladder_name = next(iter(st.LADDERS))
        self.lanes = []      # list of (label, [fragment sizes])
        self.selected = -1
        self._lane_rects = []
        self.setMinimumSize(320, 380)

    def set_lanes(self, lanes, ladder_name, selected=-1):
        self.lanes = list(lanes)
        self.ladder_name = ladder_name
        self.selected = selected
        self.update()

    def size_range(self):
        sizes = list(st.LADDERS[self.ladder_name])
        for _, fragment_sizes in self.lanes:
            sizes.extend(fragment_sizes)
        return max(10, min(sizes) * 0.8), max(sizes) * 1.25

    def y_for(self, size, top, bottom):
        low, high = self.size_range()
        fraction = (math.log10(high) - math.log10(size)) / (math.log10(high) - math.log10(low))
        return top + min(max(fraction, 0.0), 1.0) * (bottom - top)

    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), self.BACKGROUND)

        font = QFont()
        font.setPointSize(8)
        painter.setFont(font)
        metrics = QFontMetrics(font)

        axis_w = metrics.horizontalAdvance("10,000 bp") + 12
        gel = QRectF(axis_w, 28, self.width() - axis_w - 10, self.height() - 40)
        painter.fillRect(gel, self.GEL)
        count = len(self.lanes) + 1  # + ladder
        lane_w = gel.width() / max(count, 1)
        band_w = min(lane_w * 0.72, 70)
        top, bottom = gel.top() + 22, gel.bottom() - 10

        all_lanes = [("L", st.LADDERS[self.ladder_name])] + self.lanes
        max_mass = max((size for _, sizes in all_lanes for size in sizes), default=1)
        self._lane_rects = []
        for index, (label, sizes) in enumerate(all_lanes):
            cx = gel.left() + lane_w * (index + 0.5)
            lane_rect = QRectF(cx - lane_w / 2, gel.top(), lane_w, gel.height())
            if index > 0:
                self._lane_rects.append(lane_rect)
            if index - 1 == self.selected and index > 0:
                painter.fillRect(lane_rect, QColor(255, 255, 255, 18))

            # well and lane number
            painter.fillRect(QRectF(cx - band_w / 2, gel.top() + 6, band_w, 6), self.WELL)
            painter.setPen(self.TEXT)
            painter.drawText(QRectF(cx - lane_w / 2, 6, lane_w, 18),
                             Qt.AlignmentFlag.AlignCenter, label if index == 0 else str(index))

            # bands: brightness follows DNA mass (equal molar amounts -> longer is brighter);
            # identical sizes add up
            masses = {}
            for size in sizes:
                masses[size] = masses.get(size, 0) + size
            for size, mass in masses.items():
                y = self.y_for(size, top, bottom)
                alpha = int(70 + 185 * math.sqrt(min(mass / max_mass, 1.0)))
                painter.fillRect(QRectF(cx - band_w / 2, y - 2, band_w, 4), QColor(*self.BAND, alpha))

        # ladder sizes on the left
        painter.setPen(self.TEXT)
        last_y = -1e9
        for size in sorted(st.LADDERS[self.ladder_name], reverse=True):
            y = self.y_for(size, top, bottom)
            if y - last_y < metrics.height():
                continue
            painter.drawText(QRectF(0, y - metrics.height() / 2, axis_w - 6, metrics.height()),
                             Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                             f"{size:,}")
            last_y = y
        painter.end()

    def mousePressEvent(self, event):
        for index, rect in enumerate(self._lane_rects):
            if rect.contains(event.position()):
                self.laneClicked.emit(index)
                return
        super().mousePressEvent(event)


# =================================================================== dialog
class DigestDialog(QDialog):
    """`sites`: {enzyme: [positions]} for every enzyme that cuts the sequence."""

    def __init__(self, parent, name, sequence, circular, sites):
        super().__init__(parent)
        self.setWindowTitle(f"Virtual digest — {name}")
        self.resize(1000, 640)
        self.sequence = sequence
        self.circular = circular
        self.sites = sites
        self.lanes = []  # list of {"enzymes": [...], "fragments": [...]}

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            f"{name}: {len(sequence):,} bp, {'circular' if circular else 'linear'}. "
            "Check one or more enzymes, then click Add lane."
        ))
        splitter = QSplitter(Qt.Orientation.Horizontal)

        # ---- left: enzyme picker
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText("Filter enzymes")
        self.filter_edit.textChanged.connect(self.fill_enzymes)
        left_layout.addWidget(self.filter_edit)
        self.cut_filter = QComboBox()
        for text, _ in CUT_FILTERS:
            self.cut_filter.addItem(text)
        self.cut_filter.setCurrentIndex(1)
        self.cut_filter.currentIndexChanged.connect(self.fill_enzymes)
        left_layout.addWidget(self.cut_filter)
        self.enzyme_list = QListWidget()
        self.enzyme_list.itemChanged.connect(self.update_buttons)
        left_layout.addWidget(self.enzyme_list, 1)
        self.add_lane_button = QPushButton("Add lane")
        self.add_lane_button.clicked.connect(self.add_lane_from_checked)
        left_layout.addWidget(self.add_lane_button)
        splitter.addWidget(left)

        # ---- middle: gel
        middle = QWidget()
        middle_layout = QVBoxLayout(middle)
        middle_layout.setContentsMargins(0, 0, 0, 0)
        ladder_row = QHBoxLayout()
        ladder_row.addWidget(QLabel("Ladder:"))
        self.ladder_combo = QComboBox()
        self.ladder_combo.addItems(list(st.LADDERS))
        self.ladder_combo.currentIndexChanged.connect(lambda _index: self.refresh())
        ladder_row.addWidget(self.ladder_combo, 1)
        middle_layout.addLayout(ladder_row)
        self.gel = GelView()
        self.gel.laneClicked.connect(self.select_lane)
        middle_layout.addWidget(self.gel, 1)
        splitter.addWidget(middle)

        # ---- right: lanes and fragments
        right = QWidget()
        right_layout = QVBoxLayout(right)
        right_layout.setContentsMargins(0, 0, 0, 0)
        right_layout.addWidget(QLabel("Lanes"))
        self.lane_list = QListWidget()
        self.lane_list.setMaximumHeight(140)
        self.lane_list.currentRowChanged.connect(self.show_fragments)
        right_layout.addWidget(self.lane_list)
        self.remove_lane_button = QPushButton("Remove lane")
        self.remove_lane_button.clicked.connect(self.remove_lane)
        right_layout.addWidget(self.remove_lane_button)
        self.fragment_label = QLabel("Fragments")
        right_layout.addWidget(self.fragment_label)
        self.fragment_table = QTableWidget(0, 6)
        self.fragment_table.setHorizontalHeaderLabels(
            ["#", "Size (bp)", "Start", "End", "Left end", "Right end"])
        self.fragment_table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.fragment_table.verticalHeader().setVisible(False)
        self.fragment_table.horizontalHeader().setSortIndicator(0, Qt.SortOrder.AscendingOrder)
        self.fragment_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        right_layout.addWidget(self.fragment_table, 1)
        splitter.addWidget(right)
        splitter.setSizes([230, 330, 440])
        layout.addWidget(splitter, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        save_button = QPushButton("Save gel image...")
        save_button.clicked.connect(self.save_image)
        buttons.addButton(save_button, QDialogButtonBox.ButtonRole.ActionRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self.fill_enzymes()
        self.refresh()

    # ------------------------------------------------------------------ enzymes
    def checked_enzymes(self):
        return [self.enzyme_list.item(i).data(Qt.ItemDataRole.UserRole)
                for i in range(self.enzyme_list.count())
                if self.enzyme_list.item(i).checkState() == Qt.CheckState.Checked]

    def fill_enzymes(self):
        keep = set(self.checked_enzymes())
        text = self.filter_edit.text().strip().lower()
        max_cuts = CUT_FILTERS[self.cut_filter.currentIndex()][1]
        self.enzyme_list.blockSignals(True)
        self.enzyme_list.clear()
        for name, positions in sorted(self.sites.items(), key=lambda kv: (len(kv[1]), kv[0])):
            if len(positions) > max_cuts or (text and text not in name.lower()):
                continue
            cuts = len(positions)
            item = QListWidgetItem(f"{name}   ({cuts} cut{'s' if cuts != 1 else ''})")
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if name in keep else Qt.CheckState.Unchecked)
            self.enzyme_list.addItem(item)
        self.enzyme_list.blockSignals(False)
        self.update_buttons()

    def set_checked(self, names):
        """Check exactly these enzymes (used by tests and keyboard users alike)."""
        wanted = set(names)
        for i in range(self.enzyme_list.count()):
            item = self.enzyme_list.item(i)
            state = Qt.CheckState.Checked if item.data(Qt.ItemDataRole.UserRole) in wanted else Qt.CheckState.Unchecked
            item.setCheckState(state)

    # ------------------------------------------------------------------ lanes
    def add_lane(self, enzymes):
        if len(self.lanes) >= MAX_LANES:
            return False
        _, fragments = st.digest(self.sequence, enzymes, self.circular)
        self.lanes.append({"enzymes": sorted(enzymes), "fragments": fragments})
        self.refresh(select=len(self.lanes) - 1)
        return True

    def add_lane_from_checked(self):
        enzymes = self.checked_enzymes()
        if enzymes and self.add_lane(enzymes):
            self.set_checked([])

    def remove_lane(self):
        row = self.lane_list.currentRow()
        if 0 <= row < len(self.lanes):
            del self.lanes[row]
            self.refresh(select=min(row, len(self.lanes) - 1))

    def select_lane(self, index):
        self.lane_list.setCurrentRow(index)

    def refresh(self, select=None):
        if select is None:
            select = self.lane_list.currentRow()
        self.lane_list.blockSignals(True)
        self.lane_list.clear()
        for number, lane in enumerate(self.lanes, 1):
            sizes = ", ".join(f"{f['length']:,}" for f in sorted(
                lane["fragments"], key=lambda f: -f["length"]))
            self.lane_list.addItem(f"{number}: {' + '.join(lane['enzymes'])}   ({sizes} bp)")
        self.lane_list.blockSignals(False)
        gel_lanes = [(str(n), [f["length"] for f in lane["fragments"]])
                     for n, lane in enumerate(self.lanes, 1)]
        if self.lanes:
            select = min(max(select, 0), len(self.lanes) - 1)
            self.lane_list.setCurrentRow(select)
        self.gel.set_lanes(gel_lanes, self.ladder_combo.currentText(),
                           select if self.lanes else -1)
        self.show_fragments(select if self.lanes else -1)
        self.update_buttons()

    def show_fragments(self, row):
        self.gel.selected = row
        self.gel.update()
        table = self.fragment_table
        table.setSortingEnabled(False)
        if not 0 <= row < len(self.lanes):
            table.setRowCount(0)
            self.fragment_label.setText("Fragments")
            return
        lane = self.lanes[row]
        fragments = lane["fragments"]
        self.fragment_label.setText(
            f"Fragments in lane {row + 1} ({' + '.join(lane['enzymes'])}): {len(fragments)}")
        table.setRowCount(len(fragments))
        for i, fragment in enumerate(fragments):
            if fragment["uncut"]:
                left = right = "uncut circle"
            else:
                left = fragment["left"] or "sequence start"
                right = fragment["right"] or "sequence end"
            cells = [number_item(i + 1), number_item(fragment["length"]),
                     number_item(fragment["start"]), number_item(fragment["end"]),
                     QTableWidgetItem(left), QTableWidgetItem(right)]
            for column, cell in enumerate(cells):
                table.setItem(i, column, cell)
        table.setSortingEnabled(True)

    def update_buttons(self):
        self.add_lane_button.setEnabled(bool(self.checked_enzymes()) and len(self.lanes) < MAX_LANES)
        self.add_lane_button.setText("Add lane" if len(self.lanes) < MAX_LANES
                                     else f"Add lane (max {MAX_LANES})")
        self.remove_lane_button.setEnabled(self.lane_list.currentRow() >= 0)

    def save_image(self):
        path, _ = QFileDialog.getSaveFileName(self, "Save gel image", "gel.png", "PNG image (*.png)")
        if path:
            if not path.lower().endswith(".png"):
                path += ".png"
            self.gel.grab().save(path)
