"""
maps.py — the Map views.

MapView            shared base: features, lanes, sites, selection, clicking
CircularDNAViewer  plasmid map (a circle)
LinearDNAViewer    the same information on a horizontal bar
MapGroup           sends every update to several map views at once

Positions are 0-based and end-exclusive everywhere in this file.
"""
import math

from PySide6.QtCore import QPoint, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import QWidget

from . import theme
from .seqtools import feature_color, feature_label, feature_parts, feature_strand

MAX_LANES = 4        # feature rows per strand before features start to stack
MAX_LABELS = 300     # skip labels on very busy maps (e.g. whole genomes)


def _overlaps(a, b):
    return a[0] < b[1] and b[0] < a[1]


def nice_interval(length, target=10):
    """Tick spacing of 1, 2 or 5 x 10^n giving about `target` ticks."""
    raw = length / target
    if raw <= 1:
        return 1
    magnitude = 10 ** math.floor(math.log10(raw))
    for step in (1, 2, 5, 10):
        if step * magnitude >= raw:
            return int(step * magnitude)
    return int(10 * magnitude)


def small_font(size=8, bold=False):
    font = QFont()
    font.setPointSize(size)
    font.setBold(bold)
    return font


# =================================================================== base class
class MapView(QWidget):
    featureClicked = Signal(int)              # index into the features list
    featureContextMenu = Signal(int, QPoint)  # index, global position for the menu

    def __init__(self):
        super().__init__()
        self.length = 0
        self.title = ""
        self.features = []
        self.sites = {}
        self.selected = -1
        self.highlight = None  # (start, end) of the text selection, or None
        self._lanes = {}

    # ---------------------------------------------------------- public API
    def set_record(self, length, title, features=None):
        self.length = length
        self.title = title
        self.features = list(features or [])
        self.sites = {}
        self.selected = -1
        self.highlight = None
        self._lanes = self._assign_lanes()
        self.update()

    def set_features(self, features):
        """Replace the features (after an edit) but keep sites and highlight."""
        self.features = list(features or [])
        self.selected = -1
        self._lanes = self._assign_lanes()
        self.update()

    def set_sites(self, sites):
        """sites: {enzyme name: [positions]}"""
        self.sites = dict(sites)
        self.update()

    def set_selected(self, index):
        self.selected = index
        self.update()

    def set_highlight(self, start=None, end=None):
        """Shade the region start..end; no arguments clears it."""
        new = (start, end) if start is not None and end is not None and end > start else None
        if new != self.highlight:
            self.highlight = new
            self.update()

    def clear_seq(self):
        self.set_record(0, "", [])

    # ---------------------------------------------------------- lanes
    def _is_forward(self, index):
        return feature_strand(self.features[index]) != -1

    def _assign_lanes(self):
        """Put each feature on a lane (ring/row) so overlapping features don't collide."""
        lanes = {}
        for forward in (True, False):
            group = [
                (i, f) for i, f in enumerate(self.features)
                if (feature_strand(f) != -1) == forward
            ]
            # longest first, so big features sit next to the backbone
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

    @staticmethod
    def _feature_mid(feature):
        parts = feature_parts(feature)
        start, end = max(parts, key=lambda p: p[1] - p[0])
        return (start + end) / 2

    # ---------------------------------------------------------- shared painting
    def _begin_paint(self, empty_text):
        """Painter with the background filled; None if there's nothing to draw."""
        colors = theme.current()
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor(colors.background))
        if self.length <= 0:
            painter.setPen(QColor(colors.muted))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, empty_text)
            painter.end()
            return None
        return painter

    # ---------------------------------------------------------- clicking
    def feature_at(self, point):
        """Index of the feature drawn under a widget point, or -1."""
        raise NotImplementedError

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            index = self.feature_at(event.position())
            if index >= 0:
                self.featureClicked.emit(index)
                return
        super().mousePressEvent(event)

    def contextMenuEvent(self, event):
        index = self.feature_at(QPointF(event.pos()))
        if index >= 0:
            self.featureContextMenu.emit(index, event.globalPos())


# =================================================================== circular
RING_GAP = 11        # pixels between rings
ARC_WIDTH = 8        # thickness of a feature arc


class CircularDNAViewer(MapView):
    def __init__(self):
        super().__init__()
        self._hit_regions = []
        self._cx = self._cy = 0.0
        self.setMinimumSize(450, 450)

    def _angle(self, position):
        """Radians; 0 bp is at 12 o'clock, positions run clockwise."""
        return position / self.length * 2 * math.pi - math.pi / 2

    def _point(self, position, r):
        a = self._angle(position)
        return QPointF(self._cx + math.cos(a) * r, self._cy + math.sin(a) * r)

    def paintEvent(self, event):
        painter = self._begin_paint("No DNA loaded")
        if painter is None:
            return
        colors = theme.current()

        self._cx = self.width() / 2
        self._cy = self.height() / 2
        radius = min(self.width(), self.height()) * 0.27
        fwd_lanes = self._lane_count(True)
        rev_lanes = self._lane_count(False)
        inner_r = radius - 12 - rev_lanes * RING_GAP

        # backbone
        painter.setPen(QPen(QColor(colors.backbone), 3))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawEllipse(QPointF(self._cx, self._cy), radius, radius)

        # selected region from the Sequence tab, as a translucent band
        if self.highlight:
            start, end = self.highlight
            band = QPen(QColor(*colors.highlight), 18)
            band.setCapStyle(Qt.PenCapStyle.FlatCap)
            painter.setPen(band)
            rect = QRectF(self._cx - radius, self._cy - radius, 2 * radius, 2 * radius)
            span = max((end - start) / self.length * 360, 0.5)  # stay visible when tiny
            painter.drawArc(rect, round((90 - start / self.length * 360) * 16), round(-span * 16))

        self._draw_ticks(painter, radius, inner_r - 14)

        # features: forward strand outside the backbone, reverse inside
        labels = []
        self._hit_regions = []
        for i, feature in enumerate(self.features):
            if i not in self._lanes:
                continue
            offset = 12 + self._lanes[i] * RING_GAP
            r = radius + offset if self._is_forward(i) else radius - offset
            self._draw_feature(painter, i, feature, r)
            self._hit_regions.append((i, r, feature_parts(feature)))
            labels.append((0, self._feature_mid(feature), feature_label(feature),
                           colors.label_color(feature_color(feature)), r))

        # restriction sites
        site_pen = QPen(QColor(colors.site), 1.5)
        for name, positions in self.sites.items():
            for pos in positions:
                painter.setPen(site_pen)
                painter.drawLine(self._point(pos, radius - 5), self._point(pos, radius + 6))
                labels.append((1, pos, name, QColor(colors.site), radius + 6))

        if len(labels) <= MAX_LABELS:
            self._draw_labels(painter, labels, radius + 12 + fwd_lanes * RING_GAP + 14)

        self._draw_center_text(painter, inner_r)
        painter.end()

    def _draw_ticks(self, painter, radius, label_r):
        colors = theme.current()
        interval = nice_interval(self.length)
        painter.setFont(small_font())
        for pos in range(0, self.length, interval):
            painter.setPen(QPen(QColor(colors.tick), 1))
            painter.drawLine(self._point(pos, radius - 4), self._point(pos, radius + 4))
            if self.length - pos < interval / 2 and pos != 0:
                continue  # would crowd the "0" label
            p = self._point(pos, label_r)
            painter.setPen(QColor(colors.muted))
            painter.drawText(QRectF(p.x() - 30, p.y() - 7, 60, 14),
                             Qt.AlignmentFlag.AlignCenter, f"{pos:,}")

    def _draw_feature(self, painter, index, feature, r):
        color = QColor(feature_color(feature))
        base = ARC_WIDTH - 3 if feature.type == "ORF" else ARC_WIDTH
        width = base + (4 if index == self.selected else 0)
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
        colors = theme.current()
        font = small_font()
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

            painter.setPen(QPen(QColor(colors.leader), 1))
            painter.drawLine(self._point(pos, anchor_r), QPointF(x, y))
            painter.setPen(color)
            align = Qt.AlignmentFlag.AlignLeft if right_side else Qt.AlignmentFlag.AlignRight
            painter.drawText(rect, align | Qt.AlignmentFlag.AlignVCenter, text)

    def _draw_center_text(self, painter, inner_r):
        colors = theme.current()
        box_w = max(inner_r * 2, 140)
        left = self._cx - box_w / 2
        painter.setFont(small_font(12, bold=True))
        painter.setPen(QColor(colors.text))
        painter.drawText(QRectF(left, self._cy - 30, box_w, 26),
                         Qt.AlignmentFlag.AlignCenter, self.title)
        painter.setFont(small_font(10))
        painter.setPen(QColor(colors.muted))
        painter.drawText(QRectF(left, self._cy, box_w, 22),
                         Qt.AlignmentFlag.AlignCenter, f"{self.length:,} bp")

    def feature_at(self, point):
        if self.length <= 0:
            return -1
        dx, dy = point.x() - self._cx, point.y() - self._cy
        dist = math.hypot(dx, dy)
        theta = math.atan2(dy, dx)
        pos = ((theta + math.pi / 2) % (2 * math.pi)) / (2 * math.pi) * self.length
        for index, r, parts in self._hit_regions:
            if abs(dist - r) <= ARC_WIDTH / 2 + 3 and any(s <= pos <= e for s, e in parts):
                return index
        return -1


# =================================================================== linear
SIDE_MARGIN = 50     # pixels left/right of the bar
ROW_GAP = 14         # pixels between feature rows
BAR_HEIGHT = 9       # height of a feature bar (ORFs are thinner)


class LinearDNAViewer(MapView):
    def __init__(self):
        super().__init__()
        self._hit_rects = []  # (index, QRectF)
        self.setMinimumSize(450, 300)

    def _x(self, position):
        usable = max(1, self.width() - 2 * SIDE_MARGIN)
        return SIDE_MARGIN + position / self.length * usable

    def paintEvent(self, event):
        painter = self._begin_paint("No DNA loaded")
        if painter is None:
            return
        colors = theme.current()
        fwd_rows = self._lane_count(True)
        rev_rows = self._lane_count(False)

        # vertical layout: labels / forward rows / bar / reverse rows / ticks / labels
        bar_y = self.height() / 2 + (fwd_rows - rev_rows) * ROW_GAP / 2
        x0, x1 = self._x(0), self._x(self.length)

        self._draw_title(painter)

        # selected region from the Sequence tab
        if self.highlight:
            start, end = self.highlight
            left, right = self._x(start), self._x(end)
            painter.fillRect(QRectF(left, bar_y - 9, max(right - left, 2), 18),
                             QColor(*colors.highlight))

        # backbone with end caps
        painter.setPen(QPen(QColor(colors.backbone), 3))
        painter.drawLine(QPointF(x0, bar_y), QPointF(x1, bar_y))
        painter.setPen(QPen(QColor(colors.backbone), 2))
        for x in (x0, x1):
            painter.drawLine(QPointF(x, bar_y - 6), QPointF(x, bar_y + 6))

        tick_label_y = bar_y + 12 + rev_rows * ROW_GAP + 4
        self._draw_ticks(painter, bar_y, tick_label_y)

        # features: forward strand above the bar, reverse below
        top_labels, bottom_labels = [], []
        self._hit_rects = []
        for i, feature in enumerate(self.features):
            if i not in self._lanes:
                continue
            forward = self._is_forward(i)
            offset = 14 + self._lanes[i] * ROW_GAP
            y = bar_y - offset if forward else bar_y + offset
            self._draw_feature(painter, i, feature, y)
            label = (0, self._x(self._feature_mid(feature)), feature_label(feature),
                     colors.label_color(feature_color(feature)), y)
            (top_labels if forward else bottom_labels).append(label)

        # restriction sites: ticks across the bar, labels on top
        painter.setPen(QPen(QColor(colors.site), 1.5))
        for name, positions in self.sites.items():
            for pos in positions:
                x = self._x(pos)
                painter.drawLine(QPointF(x, bar_y - 7), QPointF(x, bar_y + 7))
                top_labels.append((1, x, name, QColor(colors.site), bar_y - 7))

        if len(top_labels) + len(bottom_labels) <= MAX_LABELS:
            top_base = bar_y - 14 - max(fwd_rows, 1) * ROW_GAP - 4
            self._draw_labels(painter, top_labels, top_base, upward=True)
            self._draw_labels(painter, bottom_labels, tick_label_y + 22, upward=False)
        painter.end()

    def _draw_title(self, painter):
        colors = theme.current()
        painter.setFont(small_font(12, bold=True))
        painter.setPen(QColor(colors.text))
        metrics = QFontMetrics(painter.font())
        painter.drawText(QPointF(12, 10 + metrics.ascent()), self.title)
        offset = metrics.horizontalAdvance(self.title) + 10
        painter.setFont(small_font(10))
        painter.setPen(QColor(colors.muted))
        painter.drawText(QPointF(12 + offset, 10 + metrics.ascent()), f"{self.length:,} bp")

    def _draw_ticks(self, painter, bar_y, label_y):
        colors = theme.current()
        interval = nice_interval(self.length)
        painter.setFont(small_font())
        metrics = QFontMetrics(painter.font())
        positions = list(range(0, self.length, interval)) + [self.length]
        last_right = -1e9
        for pos in positions:
            x = self._x(pos)
            painter.setPen(QPen(QColor(colors.tick), 1))
            painter.drawLine(QPointF(x, bar_y + 3), QPointF(x, bar_y + 7))
            text = f"{max(pos, 1):,}"  # show the first base as 1
            w = metrics.horizontalAdvance(text)
            if x - w / 2 < last_right + 6:
                continue  # skip labels that would touch the previous one
            painter.setPen(QColor(colors.muted))
            painter.drawText(QRectF(x - w / 2 - 2, label_y, w + 4, metrics.height()),
                             Qt.AlignmentFlag.AlignCenter, text)
            last_right = x + w / 2

    def _draw_feature(self, painter, index, feature, y):
        color = QColor(feature_color(feature))
        height = BAR_HEIGHT - 3 if feature.type == "ORF" else BAR_HEIGHT
        half = height / 2
        strand = feature_strand(feature)
        parts = feature_parts(feature)

        painter.setBrush(color)
        if index == self.selected:
            painter.setPen(QPen(QColor(theme.current().text), 2))
        else:
            painter.setPen(Qt.PenStyle.NoPen)
        for number, (start, end) in enumerate(parts):
            left, right = self._x(start), self._x(end)
            right = max(right, left + 2)
            arrow = min(8.0, (right - left) * 0.6)
            is_last = number == len(parts) - 1  # the 3' end is in the last part
            if strand == 1 and is_last:
                shape = [QPointF(left, y - half), QPointF(right - arrow, y - half),
                         QPointF(right, y), QPointF(right - arrow, y + half),
                         QPointF(left, y + half)]
            elif strand == -1 and is_last:
                shape = [QPointF(right, y - half), QPointF(left + arrow, y - half),
                         QPointF(left, y), QPointF(left + arrow, y + half),
                         QPointF(right, y + half)]
            else:
                shape = [QPointF(left, y - half), QPointF(right, y - half),
                         QPointF(right, y + half), QPointF(left, y + half)]
            painter.drawPolygon(QPolygonF(shape))
            self._hit_rects.append((index, QRectF(left, y - half, right - left, height)))
        painter.setBrush(Qt.BrushStyle.NoBrush)

    def _draw_labels(self, painter, labels, base_y, upward):
        """Place labels in rows above (upward) or below the features, avoiding overlaps."""
        colors = theme.current()
        font = small_font()
        painter.setFont(font)
        metrics = QFontMetrics(font)
        h = metrics.height()
        placed = []
        for _, x, text, color, anchor_y in sorted(labels, key=lambda item: (item[0], item[1])):
            w = metrics.horizontalAdvance(text) + 4
            left = min(max(x - w / 2, 2), self.width() - w - 2)
            for row in range(8):
                top = base_y - h - row * h if upward else base_y + row * h
                rect = QRectF(left, top, w, h)
                if rect.top() < 0 or rect.bottom() > self.height():
                    rect = None
                    break
                if not any(rect.intersects(other) for other in placed):
                    break
            else:
                rect = None
            if rect is None:
                continue  # no free spot: skip this label rather than overlap
            placed.append(rect)
            painter.setPen(QPen(QColor(colors.leader), 1))
            edge = rect.bottom() if upward else rect.top()
            painter.drawLine(QPointF(x, anchor_y), QPointF(min(max(x, rect.left()), rect.right()), edge))
            painter.setPen(color)
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, text)

    def feature_at(self, point):
        for index, rect in reversed(self._hit_rects):  # topmost drawn last
            if rect.adjusted(-2, -3, 2, 3).contains(point):
                return index
        return -1


# =================================================================== group
class MapGroup:
    """Forwards every call (set_record, set_sites, ...) to all map views,
    so switching between circular and linear needs no refresh. Reading an
    attribute (e.g. .sites) returns the first view's value."""

    def __init__(self, *views):
        self.views = views

    def __getattr__(self, name):
        first = getattr(self.views[0], name)
        if not callable(first):
            return first

        def call_all(*args, **kwargs):
            result = None
            for view in self.views:
                result = getattr(view, name)(*args, **kwargs)
            return result
        return call_all

    def update(self):
        for view in self.views:
            view.update()
