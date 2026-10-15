"""
theme.py — light and dark color schemes.

Standard widgets (buttons, tables, menus) get their colors from a Qt palette.
The widgets N-Zoom draws itself (maps, sequence colors) read the active
Theme from `current()` every time they paint, so a theme switch only needs a
repaint.

    apply_theme(app, "dark")   # "system", "light" or "dark"
    current().background       # colors for custom drawing
"""
from dataclasses import dataclass, field

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QGuiApplication, QPalette

from .seqtools import DNA_COLORS, PROTEIN_COLORS

THEME_MODES = ("system", "light", "dark")


@dataclass(frozen=True)
class Theme:
    name: str
    dark: bool
    background: str          # map canvas
    text: str                # titles on maps
    muted: str               # tick labels, secondary text
    backbone: str            # the DNA line/circle
    tick: str
    leader: str              # lines from labels to features
    site: str                # restriction sites
    highlight: tuple         # RGBA band showing the selected region on maps
    feature_highlight: str   # background of a selected feature in the Sequence view
    other_letter: str        # sequence letters without a specific color
    match: str               # Find matches marked on the maps
    match_highlight: str     # background of Find matches in the Sequence view
    dna_colors: dict = field(default_factory=dict)
    protein_colors: dict = field(default_factory=dict)

    def label_color(self, feature_hex):
        """Feature label color that stays readable on this background."""
        color = QColor(feature_hex)
        return color.lighter(140) if self.dark else color.darker(140)


LIGHT = Theme(
    name="light", dark=False,
    background="#ffffff", text="#222222", muted="#666666", backbone="#333333",
    tick="#555555", leader="#bdbdbd", site="#c62828",
    highlight=(255, 193, 7, 160), feature_highlight="#fff59d", other_letter="#424242",
    match="#8e24aa", match_highlight="#ffcc80",
    dna_colors=dict(DNA_COLORS), protein_colors=dict(PROTEIN_COLORS),
)

DARK = Theme(
    name="dark", dark=True,
    background="#1e1f22", text="#e8e8e8", muted="#a0a0a0", backbone="#c8c8c8",
    tick="#9a9a9a", leader="#55585e", site="#ef5350",
    highlight=(255, 202, 40, 110), feature_highlight="#5c4f12", other_letter="#bdbdbd",
    match="#ce93d8", match_highlight="#6d4318",
    dna_colors={"A": "#66bb6a", "T": "#ef5350", "U": "#ef5350", "G": "#ffca28",
                "C": "#42a5f5", "N": "#9e9e9e"},
    protein_colors={
        **{aa: "#81c784" for aa in "AVLIM"}, **{aa: "#ce93d8" for aa in "FWY"},
        **{aa: "#64b5f6" for aa in "KRH"}, **{aa: "#e57373" for aa in "DE"},
        **{aa: "#4dd0e1" for aa in "STNQ"}, **{aa: "#bcaaa4" for aa in "GP"},
        "C": "#ffb74d", "*": "#f5f5f5",
    },
)

_current = LIGHT


def current():
    """The theme custom-drawn widgets should paint with."""
    return _current


def system_prefers_dark():
    hints = QGuiApplication.styleHints()
    if hasattr(hints, "colorScheme"):  # Qt 6.5+
        return hints.colorScheme() == Qt.ColorScheme.Dark
    return False


def dark_palette():
    window, base, text = QColor("#2b2d30"), QColor("#1e1f22"), QColor("#e8e8e8")
    disabled = QColor("#7a7a7a")
    palette = QPalette()
    role = QPalette.ColorRole
    palette.setColor(role.Window, window)
    palette.setColor(role.WindowText, text)
    palette.setColor(role.Base, base)
    palette.setColor(role.AlternateBase, window)
    palette.setColor(role.ToolTipBase, window)
    palette.setColor(role.ToolTipText, text)
    palette.setColor(role.Text, text)
    palette.setColor(role.Button, QColor("#3a3d41"))
    palette.setColor(role.ButtonText, text)
    palette.setColor(role.BrightText, QColor("#ff6b6b"))
    palette.setColor(role.Link, QColor("#64b5f6"))
    palette.setColor(role.Highlight, QColor("#2f65ca"))
    palette.setColor(role.HighlightedText, QColor("#ffffff"))
    palette.setColor(role.PlaceholderText, QColor("#8a8a8a"))
    palette.setColor(role.Light, QColor("#4a4d52"))
    palette.setColor(role.Midlight, QColor("#3f4246"))
    palette.setColor(role.Mid, QColor("#2f3134"))
    palette.setColor(role.Dark, QColor("#18191b"))
    palette.setColor(role.Shadow, QColor("#000000"))
    for r in (role.WindowText, role.Text, role.ButtonText):
        palette.setColor(QPalette.ColorGroup.Disabled, r, disabled)
    return palette


def apply_theme(app, mode):
    """Switch the whole app to `mode` ("system", "light" or "dark").
    Returns the Theme now in use."""
    global _current
    dark = mode == "dark" or (mode == "system" and system_prefers_dark())
    # Fusion draws the same on every platform and fully respects the palette
    # (some native styles ignore dark palettes)
    app.setStyle("Fusion")
    app.setPalette(dark_palette() if dark else app.style().standardPalette())
    _current = DARK if dark else LIGHT
    return _current
