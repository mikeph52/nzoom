"""
primers.py — the "Find primer" dialog: paste a primer, see where it binds
on either strand with its Tm, and add the sites to the map as features.
"""
from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
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
from .dialogs import monospace_font, number_item


class PrimerDialog(QDialog):
    """`on_add(features)` is called with new primer_bind features when the
    user clicks "Add as features"."""

    COLUMNS = ["Strand", "Start", "End", "Annealed", "5′ tail", "Tm annealed (°C)"]

    def __init__(self, parent, sequence, circular, on_add):
        super().__init__(parent)
        self.setWindowTitle("Find primer")
        self.resize(640, 520)
        self.sequence = sequence
        self.circular = circular
        self.on_add = on_add
        self.sites = []
        self.primer = ""

        layout = QVBoxLayout(self)
        form = QFormLayout()
        self.name_edit = QLineEdit()
        self.name_edit.setPlaceholderText("e.g. GFP-fwd")
        form.addRow("Name:", self.name_edit)

        self.primer_edit = QPlainTextEdit()
        self.primer_edit.setFont(monospace_font())
        self.primer_edit.setPlaceholderText("Primer sequence, 5′ → 3′")
        self.primer_edit.setFixedHeight(56)
        self.primer_edit.textChanged.connect(self.search)
        form.addRow("Sequence:", self.primer_edit)

        self.anneal_spin = QSpinBox()
        self.anneal_spin.setRange(8, 60)
        self.anneal_spin.setValue(st.PRIMER_MIN_ANNEAL)
        self.anneal_spin.setSuffix(" bases at the 3′ end")
        self.anneal_spin.valueChanged.connect(self.search)
        form.addRow("Must match:", self.anneal_spin)
        layout.addLayout(form)

        self.info = QLabel(" ")
        self.info.setWordWrap(True)
        self.info.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        layout.addWidget(self.info)

        self.table = QTableWidget(0, len(self.COLUMNS))
        self.table.setHorizontalHeaderLabels(self.COLUMNS)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.verticalHeader().setVisible(False)
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self.table.itemSelectionChanged.connect(self.update_buttons)
        layout.addWidget(self.table, 1)

        note = QLabel(
            f"Tm: nearest-neighbour method, {st.PRIMER_NA_MM} mM Na⁺, "
            f"{st.PRIMER_CONC_NM} nM primer. Treat it as an estimate."
        )
        note.setStyleSheet("color: gray;")
        note.setWordWrap(True)
        layout.addWidget(note)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        self.add_button = QPushButton("Add selected sites as features")
        self.add_button.clicked.connect(self.add_selected)
        buttons.addButton(self.add_button, QDialogButtonBox.ButtonRole.ActionRole)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self.update_buttons()

    # ------------------------------------------------------------------
    def search(self):
        text = self.primer_edit.toPlainText()
        self.sites = []
        self.primer = ""
        self.table.setRowCount(0)
        if not text.strip():
            self.info.setText(" ")
            self.update_buttons()
            return
        try:
            primer = st.clean_primer(text)
        except ValueError as error:
            self.info.setText(f"<span style='color:#c62828'>{error}</span>")
            self.update_buttons()
            return
        self.primer = primer

        details = [f"{len(primer)} nt", f"GC {st.gc_percent(primer):.0f}%"]
        if len(primer) >= 8:
            details.append(f"Tm {st.primer_tm(primer):.1f} °C (whole primer)")
        if len(primer) < 10:
            details.append("too short to search")
            self.info.setText("  ·  ".join(details))
            self.update_buttons()
            return

        self.sites = st.find_primer_sites(self.sequence, primer, self.circular,
                                          min_anneal=self.anneal_spin.value())
        count = len(self.sites)
        details.append("no binding sites" if not count else
                       f"{count} binding site{'s' if count != 1 else ''}")
        if count > 1:
            details.append("binds more than once: PCR may give extra products")
        self.info.setText("  ·  ".join(details))

        self.table.setRowCount(count)
        for row, site in enumerate(self.sites):
            annealed = primer[-site["anneal"]:]
            cells = [
                QTableWidgetItem("→ forward" if site["strand"] == 1 else "← reverse"),
                number_item(site["start"]),
                number_item(site["end"]),
                number_item(site["anneal"]),
                number_item(site["tail"]),
                QTableWidgetItem(f"{st.primer_tm(annealed):.1f}"),
            ]
            for column, cell in enumerate(cells):
                self.table.setItem(row, column, cell)
        self.table.selectAll()  # most primers bind once: ready to add
        self.update_buttons()

    def update_buttons(self):
        self.add_button.setEnabled(bool(self.table.selectionModel().selectedRows()))

    def selected_sites(self):
        rows = sorted({index.row() for index in self.table.selectionModel().selectedRows()})
        return [self.sites[row] for row in rows]

    def make_features(self, sites):
        name = self.name_edit.text().strip() or "primer"
        tm = st.primer_tm(self.primer)
        features = []
        for site in sites:
            note = f"5'-{self.primer}-3', Tm {tm:.1f} C"
            if site["tail"]:
                note += f", {site['tail']} nt 5' tail"
            features.append(st.new_feature("primer_bind", name, site["location"], note))
        return features

    def add_selected(self):
        sites = self.selected_sites()
        if sites:
            self.on_add(self.make_features(sites))
            self.accept()
