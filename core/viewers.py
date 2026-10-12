"""
N-Zoom seqtools.py by mikeph52, Init 11/10/2026
All view functions:
- SequenceViewer      colored text view of the sequence, can highlight regions
- CircularDNAViewer   plasmid map with features, ticks and restriction sites
"""
from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QColor, QFont, QTextCharFormat, QTextCursor
from PySide6.QtWidgets import QTextEdit

from . import theme
from .seqtools import MAX_COLORED_LENGTH, normalize_seq


class SequenceViewer(QTextEdit):
    # (start, end) in 0-based, end-exclusive sequence coordinates.
    # start == end means nothing is selected and the cursor is at `start`.
    selectionInfo = Signal(int, int)

    def __init__(self):
        super().__init__()
        self.setReadOnly(True)
        self.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        # allow Shift+arrow keys as well as the mouse for selecting
        self.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
            | Qt.TextInteractionFlag.TextSelectableByKeyboard
        )
        self.selectionChanged.connect(self._emit_selection)
        self.cursorPositionChanged.connect(self._emit_selection)
        self._sequence = ""
        self._seq_type = "Unknown"
        self._highlighted_parts = []

        font = QFont("Monospace")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)

    def set_sequence(self, sequence, seq_type):
        self._sequence = normalize_seq(sequence)
        self._seq_type = seq_type
        self._highlighted_parts = []
        self._render()
        self._move_to(0)

    def refresh_theme(self):
        """Re-color the text for the current theme, keeping cursor and highlight."""
        cursor = self.textCursor()
        anchor, position = cursor.anchor(), cursor.position()
        parts = self._highlighted_parts
        self._render()
        if parts:
            self.highlight_parts(parts)
        restored = QTextCursor(self.document())
        restored.setPosition(anchor)
        restored.setPosition(position, QTextCursor.MoveMode.KeepAnchor)
        self.setTextCursor(restored)

    def _render(self):
        self.clear()
        self.setExtraSelections([])
        sequence = self._sequence
        colors = theme.current()
        if len(sequence) > MAX_COLORED_LENGTH:
            self.setPlainText(sequence)  # big sequences: plain text, fast
            return
        palette = colors.dna_colors if self._seq_type == "DNA" else colors.protein_colors
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
            fmt.setForeground(QColor(palette.get(char, colors.other_letter)))
            cursor.insertText(sequence[i:j], fmt)
            i = j
        cursor.endEditBlock()

    def highlight_parts(self, parts, color=None):
        """Give each (start, end) region a background color and scroll to it."""
        self._highlighted_parts = list(parts)
        color = color or theme.current().feature_highlight
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

    def clear_highlight(self):
        self._highlighted_parts = []
        self.setExtraSelections([])

    def select_range(self, start, end):
        """Select start..end (0-based, end exclusive) and scroll to it."""
        last = max(0, self.document().characterCount() - 1)
        cursor = QTextCursor(self.document())
        cursor.setPosition(min(max(start, 0), last))
        cursor.setPosition(min(max(end, 0), last), QTextCursor.MoveMode.KeepAnchor)
        self.setTextCursor(cursor)
        self.ensureCursorVisible()

    def selected_range(self):
        """(start, end) of the current selection, or None."""
        cursor = self.textCursor()
        if cursor.hasSelection():
            return cursor.selectionStart(), cursor.selectionEnd()
        return None

    def _emit_selection(self):
        # the text is the bare sequence (no newlines), so document positions
        # are exactly sequence positions
        cursor = self.textCursor()
        self.selectionInfo.emit(cursor.selectionStart(), cursor.selectionEnd())

    def _move_to(self, position):
        last = max(0, self.document().characterCount() - 1)
        cursor = QTextCursor(self.document())
        cursor.setPosition(min(max(position, 0), last))
        self.setTextCursor(cursor)
        self.ensureCursorVisible()
