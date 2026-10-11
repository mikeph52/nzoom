"""
N-Zoom seqtools.py by mikeph52, Init 11/10/2026
All view functions:
- SequenceViewer      colored text view of the sequence, can highlight regions
- CircularDNAViewer   plasmid map with features, ticks and restriction sites
"""
import math
from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontMetrics,
    QPainter,
    QPen,
    QPolygonF,
    QTextCharFormat,
    QTextCursor,
)
from PySide6.QtWidgets import QTextEdit, QWidget
from .seqtools import (
    DNA_COLORS,
    MAX_COLORED_LENGTH,
    PROTEIN_COLORS,
    feature_color,
    feature_label,
    feature_parts,
    feature_strand,
    normalize_seq,
) # dot for loading from core

# text view
class SequenceViewer(QTextEdit):
    def __init__(self):
        super().__init__()
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)

        font = QFont("Monospace")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)

    def set_sequence(self, sequence, seq_type):
        self.clear()
        self.setExtraSelections([])
        sequence = normalize_seq(sequence)

        if len(sequence) > MAX_COLORED_LENGTH:
            self.setPlainText(sequence)  # big sequences: plain text, fast
        else:
            colors = DNA_COLORS if seq_type == "DNA" else PROTEIN_COLORS
            cursor = self.textCursor()
            cursor.beginEditBlock()
            # insert runs of identical letters with one format call each
            i, n = 0, len(sequence)
            while i < n:
                char = sequence[i]
                j = i + 1
                while j < n and sequence[j] == char:
                    j += 1
                fmt = QTextCharFormat()
                fmt.setForeground(QColor(colors.get(char, "#424242")))
                cursor.insertText(sequence[i:j], fmt)
                i = j
            cursor.endEditBlock()

        self._move_to(0)

    def highlight_parts(self, parts, color="#fff59d"):
        # Give each (start, end) region a background color and scroll to it.
        last = max(0, self.document().characterCount() - 1)
        selections = []
        for start, end in parts:
            cursor = QTextCursor(self.document())
            cursor.setPosition(min(max(start, 0), last))
            cursor.setPosition(min(max(end, 0), last), QTextCursor.MoveMode.KeepAnchor)
            fmt = QTextCharFormat()
            fmt.setBackground(QColor(color))
            selection = QTextEdit.ExtraSelection()
            selection.cursor = cursor
            selection.format = fmt
            selections.append(selection)
        self.setExtraSelections(selections)
        if parts:
            self._move_to(parts[0][0])

    def _move_to(self, position):
        last = max(0, self.document().characterCount() - 1)
        cursor = QTextCursor(self.document())
        cursor.setPosition(min(max(position, 0), last))
        self.setTextCursor(cursor)
        self.ensureCursorVisible()

# map view
MAX_LANES = 4        # feature rings per strand before features start to stack
LANE_GAP = 11        # pixels between rings
FEATURE_WIDTH = 8    # thickness of a feature arc
MAX_LABELS = 300     # skip labels on very busy maps (e.g. whole genomes)

def _overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]

def nice_interval(length, target=10):
    # Tick spacing of 1, 2 or 5 x 10^n giving about `target` ticks.
    raw = length / target
    if raw <= 1:
        return 1
    magnitude = 10 ** math.floor(math.log10(raw))
    for step in (1, 2, 5, 10):
        if step * magnitude >= raw:
            return int(step * magnitude)
    return int(10 * magnitude)

class CircularDNAViewer(QWidget):
    featureClicked = Signal(int)  # index into the features list

    def __init__(self):
        super().__init__()
        self.length = 0
        self.title = "Circular DNA"
        self.features = []
        self.sites = {}
        self.selected = -1
        self._lanes = {}
        self._hit_regions = []
        self._cx = self._cy = 0.0
        self.setMinimumSize(450, 450)

    # public API
    def set_record(self, length, title, features=None):
        self.length = length
        self.title = title
        self.features = list(features or [])
        self.sites = {}
        self.selected = -1
        self._lanes = self._assign_lanes()
        self.update()

    def set_sites(self, sites):
        # sites: {enzyme name: [positions]}
        self.sites = dict(sites)
        self.update()

    def set_selected(self, index):
        self.selected = index
        self.update()

    def clear_seq(self):
        self.set_record(0, "", [])

    # geometry
    def _angle(self, position):
        # Radians; 0 bp is at 12 o'clock, positions run clockwise.
        return position / self.length * 2 * math.pi - math.pi / 2

    def _point(self, position, r):
        a = self._angle(position)
        return QPointF(self._cx + math.cos(a) * r, self._cy + math.sin(a) * r)

    def _is_forward(self, index):
        return feature_strand(self.features[index]) != -1

    def _assign_lanes(self):
        # Put each feature on a ring so overlapping features don't collide.
        lanes = {}
        for forward in (True, False):
            group = [
                (i, f) for i, f in enumerate(self.features)
                if (feature_strand(f) != -1) == forward
            ]
            # longest first, so big features sit on the inner lanes
            group.sort(key=lambda item: -sum(e - s for s, e in feature_parts(item[1])))
            lane_contents = []
            for i, feature in group:
                parts = feature_parts(feature)
                for lane, contents in enumerate(lane_contents):
                    if not any(_overlaps(p, q) for p in parts for q in contents):
                        contents.extend(parts)
                        lanes[i] = lane
                        break
                else:
                    if len(lane_contents) < MAX_LANES:
                        lane_contents.append(list(parts))
                        lanes[i] = len(lane_contents) - 1
                    else:
                        lane_contents[-1].extend(parts)
                        lanes[i] = MAX_LANES - 1
        return lanes

    def _lane_count(self, forward):
        used = [lane for i, lane in self._lanes.items() if self._is_forward(i) == forward]
        return max(used) + 1 if used else 0

    def _feature_mid(self, feature):
        parts = feature_parts(feature)
        start, end = max(parts, key=lambda p: p[1] - p[0])
        return (start + end) / 2

    # painting
    def paintEvent(self, event):
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#ffffff"))

        if self.length <= 0:
            painter.setPen(QColor("#777777"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No circular DNA loaded")
            return

        self._cx = self.width() / 2
        self._cy = self.height() / 2
        radius = min(self.width(), self.height()) * 0.27
        fwd_lanes = self._lane_count(True)
        rev_lanes = self._lane_count(False)
        inner_r = radius - 12 - rev_lanes * LANE_GAP

        # backbone
        painter.setPen(QPen(QColor("#333333"), 3))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QPointF(self._cx, self._cy), radius, radius)

        self._draw_ticks(painter, radius, inner_r - 14)

        # features: forward strand outside the backbone, reverse inside
        labels = []
        self._hit_regions = []
        for i, feature in enumerate(self.features):
            if i not in self._lanes:
                continue
            offset = 12 + self._lanes[i] * LANE_GAP
            r = radius + offset if self._is_forward(i) else radius - offset
            self._draw_feature(painter, i, feature, r)
            self._hit_regions.append((i, r, feature_parts(feature)))
            labels.append((
                0,  # priority: features get placed before enzymes
                self._feature_mid(feature),
                feature_label(feature),
                QColor(feature_color(feature)).darker(140),
                r,
            ))

        # restriction sites
        site_pen = QPen(QColor("#c62828"), 1.5)
        for name, positions in self.sites.items():
            for pos in positions:
                painter.setPen(site_pen)
                painter.drawLine(self._point(pos, radius - 5), self._point(pos, radius + 6))
                labels.append((1, pos, name, QColor("#c62828"), radius + 6))

        if len(labels) <= MAX_LABELS:
            self._draw_labels(painter, labels, radius + 12 + fwd_lanes * LANE_GAP + 14)

        self._draw_center_text(painter, inner_r)

    def _draw_ticks(self, painter, radius, label_r):
        interval = nice_interval(self.length)
        font = QFont()
        font.setPointSize(8)
        painter.setFont(font)
        for pos in range(0, self.length, interval):
            painter.setPen(QPen(QColor("#555555"), 1))
            painter.drawLine(self._point(pos, radius - 4), self._point(pos, radius + 4))
            if self.length - pos < interval / 2 and pos != 0:
                continue  # would crowd the "0" label
            p = self._point(pos, label_r)
            painter.setPen(QColor("#666666"))
            painter.drawText(
                QRectF(p.x() - 30, p.y() - 7, 60, 14),
                Qt.AlignmentFlag.AlignCenter,
                f"{pos:,}",
            )

    def _draw_feature(self, painter, index, feature, r):
        color = QColor(feature_color(feature))
        width = FEATURE_WIDTH + (4 if index == self.selected else 0)
        pen = QPen(color, width)
        pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        painter.setPen(pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)

        rect = QRectF(self._cx - r, self._cy - r, 2 * r, 2 * r)
        parts = feature_parts(feature)
        for start, end in parts:
            # Qt angles: 1/16 degree, 0 = 3 o'clock, positive = counter-clockwise
            start_deg = 90 - start / self.length * 360
            span_deg = -(end - start) / self.length * 360
            painter.drawArc(rect, round(start_deg * 16), round(span_deg * 16))

        # arrowhead at the 3' end
        strand = feature_strand(feature)
        if strand and parts:
            start, end = parts[-1]
            theta = self._angle(end if strand == 1 else start)
            tx, ty = -math.sin(theta) * strand, math.cos(theta) * strand  # along the strand
            nx, ny = math.cos(theta), math.sin(theta)                     # radial
            px, py = self._cx + nx * r, self._cy + ny * r
            half = width / 2 + 3
            arrow = QPolygonF([
                QPointF(px + tx * 7, py + ty * 7),
                QPointF(px + nx * half, py + ny * half),
                QPointF(px - nx * half, py - ny * half),
            ])
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(color)
            painter.drawPolygon(arrow)
            painter.setBrush(Qt.BrushStyle.NoBrush)

    def _draw_labels(self, painter, labels, base_r):
        font = QFont()
        font.setPointSize(8)
        painter.setFont(font)
        metrics = QFontMetrics(font)
        placed = []

        for _, pos, text, color, anchor_r in sorted(labels, key=lambda item: (item[0], item[1])):
            theta = self._angle(pos)
            cos_t, sin_t = math.cos(theta), math.sin(theta)
            w, h = metrics.horizontalAdvance(text) + 2, metrics.height()
            right_side = cos_t >= 0

            # push the label outwards until it no longer overlaps another one
            r = base_r
            for _ in range(6):
                x, y = self._cx + cos_t * r, self._cy + sin_t * r
                rect = QRectF(x if right_side else x - w, y - h / 2, w, h)
                if not any(rect.intersects(other) for other in placed):
                    break
                r += h
            else:
                continue  # no free spot: skip this label rather than overlap
            placed.append(rect)

            painter.setPen(QPen(QColor("#bdbdbd"), 1))
            painter.drawLine(self._point(pos, anchor_r), QPointF(x, y))
            painter.setPen(color)
            align = Qt.AlignmentFlag.AlignLeft if right_side else Qt.AlignmentFlag.AlignRight
            painter.drawText(rect, align | Qt.AlignmentFlag.AlignVCenter, text)

    def _draw_center_text(self, painter, inner_r):
        box_w = max(inner_r * 2, 140)
        left = self._cx - box_w / 2

        title_font = QFont()
        title_font.setPointSize(12)
        title_font.setBold(True)
        painter.setFont(title_font)
        painter.setPen(QColor("#222222"))
        painter.drawText(QRectF(left, self._cy - 30, box_w, 26),
                         Qt.AlignmentFlag.AlignCenter, self.title)

        length_font = QFont()
        length_font.setPointSize(10)
        painter.setFont(length_font)
        painter.setPen(QColor("#555555"))
        painter.drawText(QRectF(left, self._cy, box_w, 22),
                         Qt.AlignmentFlag.AlignCenter, f"{self.length:,} bp")

    # enaBLE clicking
    def mousePressEvent(self, event):
        if self.length > 0:
            p = event.position()
            dx, dy = p.x() - self._cx, p.y() - self._cy
            dist = math.hypot(dx, dy)
            theta = math.atan2(dy, dx)
            pos = ((theta + math.pi / 2) % (2 * math.pi)) / (2 * math.pi) * self.length
            for index, r, parts in self._hit_regions:
                if abs(dist - r) <= FEATURE_WIDTH / 2 + 3 and any(s <= pos <= e for s, e in parts):
                    self.featureClicked.emit(index)
                    return
        super().mousePressEvent(event)