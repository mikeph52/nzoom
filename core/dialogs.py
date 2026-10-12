"""
dialogs.py — pop-up windows for N-Zoom.

FeatureDialog           add or edit a feature
TextResultDialog        show a sequence/translation with a Copy button
RestrictionSitesDialog  table of enzymes and cut positions
"""
from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QTextOption
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from . import seqtools as st

STRAND_CHOICES = [("Forward  →", 1), ("Reverse  ←", -1), ("None", 0)]


def number_item(value):
    """Table cell that sorts numerically."""
    item = QTableWidgetItem()
    item.setData(Qt.ItemDataRole.DisplayRole, value)
    return item


def monospace_font(size=11):
    font = QFont("Monospace")
    font.setStyleHint(QFont.StyleHint.Monospace)
    font.setPointSize(size)
    return font


# =================================================================== features
class FeatureDialog(QDialog):
    """Add a new feature, or edit an existing one.

    Coordinates are typed 1-based and inclusive, like in GenBank files.
    After exec() returns Accepted, read the result with values().
    """

    def __init__(self, parent, length, circular, feature=None,
                 start=1, end=1, strand=1, feature_type="misc_feature", name="", note=""):
        super().__init__(parent)
        self.length = length
        self.circular = circular
        self.setWindowTitle("Edit feature" if feature is not None else "Add feature")
        self.setMinimumWidth(420)

        # a feature made of several pieces (e.g. exons) can't be described by a
        # single start/end, so its location is shown but not editable here
        self.location_editable = feature is None or st.has_simple_location(feature, length)
        if feature is not None:
            feature_type = feature.type
            name = feature.qualifiers.get("label", [""])[0] or st.feature_label(feature)
            note = st.feature_note(feature)
            start, end = st.feature_span(feature)
            strand = st.feature_strand(feature)

        layout = QVBoxLayout(self)
        form = QFormLayout()

        self.name_edit = QLineEdit(name)
        self.name_edit.setPlaceholderText("e.g. AmpR, lac promoter")
        form.addRow("Name:", self.name_edit)

        self.type_combo = QComboBox()
        self.type_combo.setEditable(True)  # any GenBank feature key is allowed
        types = list(st.COMMON_FEATURE_TYPES)
        if feature_type not in types:
            types.insert(0, feature_type)
        self.type_combo.addItems(types)
        self.type_combo.setCurrentText(feature_type)
        form.addRow("Type:", self.type_combo)

        self.start_spin = QSpinBox()
        self.end_spin = QSpinBox()
        for spin, value in ((self.start_spin, start), (self.end_spin, end)):
            spin.setRange(1, max(1, length))
            spin.setGroupSeparatorShown(True)
            spin.setValue(min(max(1, value), max(1, length)))
            spin.valueChanged.connect(self.update_hint)
        form.addRow("Start:", self.start_spin)
        form.addRow("End:", self.end_spin)

        self.strand_combo = QComboBox()
        for text, value in STRAND_CHOICES:
            self.strand_combo.addItem(text, value)
        self.strand_combo.setCurrentIndex([v for _, v in STRAND_CHOICES].index(strand))
        form.addRow("Strand:", self.strand_combo)

        self.note_edit = QLineEdit(note)
        self.note_edit.setPlaceholderText("optional")
        form.addRow("Note:", self.note_edit)
        layout.addLayout(form)

        self.hint = QLabel()
        self.hint.setWordWrap(True)
        layout.addWidget(self.hint)

        if not self.location_editable:
            for widget in (self.start_spin, self.end_spin, self.strand_combo):
                widget.setEnabled(False)
            parts = len(st.feature_parts(feature))
            self.hint.setText(
                f"This feature is made of {parts} separate pieces ({feature.location}). "
                "Its name, type and note can be changed here, but not its location."
            )

        self.buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self.buttons.accepted.connect(self.accept)
        self.buttons.rejected.connect(self.reject)
        layout.addWidget(self.buttons)

        self.update_hint()
        self.name_edit.setFocus()

    def update_hint(self):
        if not self.location_editable:
            return
        start, end = self.start_spin.value(), self.end_spin.value()
        ok = True
        if start <= end:
            text = f"Length: {end - start + 1:,} bp"
        elif self.circular:
            text = f"Length: {self.length - start + 1 + end:,} bp — crosses the origin"
        else:
            text = ("Start is after end. Only features on circular sequences can "
                    "cross the origin.")
            ok = False
        self.hint.setText(text)
        self.hint.setStyleSheet("" if ok else "color: #c62828;")
        self.buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)

    def values(self):
        """dict with type, name, note and location (None if the location wasn't editable)."""
        location = None
        if self.location_editable:
            location = st.make_location(
                self.start_spin.value(), self.end_spin.value(),
                self.strand_combo.currentData(), self.length, self.circular,
            )
        return {
            "type": self.type_combo.currentText().strip() or "misc_feature",
            "name": self.name_edit.text().strip(),
            "note": self.note_edit.text().strip(),
            "location": location,
        }


# =================================================================== text
class TextResultDialog(QDialog):
    """Non-modal window showing a sequence or protein, with a Copy button."""

    def __init__(self, parent, title, summary, text):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(560, 380)
        self.text = text

        layout = QVBoxLayout(self)
        label = QLabel(summary)
        label.setWordWrap(True)
        label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(label)

        self.text_view = QPlainTextEdit(text)
        self.text_view.setReadOnly(True)
        self.text_view.setFont(monospace_font())
        self.text_view.setWordWrapMode(QTextOption.WrapMode.WrapAnywhere)
        layout.addWidget(self.text_view, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        copy_button = QPushButton("Copy")
        copy_button.clicked.connect(self.copy)
        buttons.addButton(copy_button, QDialogButtonBox.ButtonRole.ActionRole)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)

    def copy(self):
        QApplication.clipboard().setText(self.text)
        self.setWindowTitle(self.windowTitle().removesuffix(" (copied)") + " (copied)")


# =================================================================== enzymes
class RestrictionSitesDialog(QDialog):
    """Non-modal table of every enzyme that cuts the sequence."""

    def __init__(self, parent, sites, circular):
        super().__init__(parent)
        self.setWindowTitle("Restriction sites")
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        self.resize(520, 560)

        layout = QVBoxLayout(self)
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
        table.horizontalHeader().setSortIndicator(1, Qt.SortOrder.AscendingOrder)  # fewest cuts first
        table.setSortingEnabled(True)
        table.resizeColumnsToContents()
        table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(table)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        layout.addWidget(buttons)
