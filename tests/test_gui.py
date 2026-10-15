"""
Tests for the main window. Slower than the seqtools tests, so keep them to
whole workflows: open -> do something -> check what the user would see.

Pop-up dialogs would block a test, so tests replace them using pytest's
`monkeypatch` (the `window` fixture in conftest.py already handles message
boxes).
"""
import time

import pytest
from PySide6.QtWidgets import QDialog, QFileDialog, QInputDialog

from core import seqtools as st
from core.dialogs import FeatureDialog, TextResultDialog
from conftest import DEMO


def open_demo(window, monkeypatch, tmp_path):
    """Open a copy of pDEMO.gb (so saving doesn't change the real example)."""
    path = tmp_path / "pDEMO.gb"
    path.write_bytes(DEMO.read_bytes())
    monkeypatch.setattr(QFileDialog, "getOpenFileName", staticmethod(lambda *a, **k: (str(path), "")))
    window.open_file()
    return path


def fill_feature_dialog(monkeypatch, **fields):
    """Make the next FeatureDialog fill in `fields` and press OK."""
    seen = {}

    def fake_exec(dialog):
        seen.update(start=dialog.start_spin.value(), end=dialog.end_spin.value(),
                    type=dialog.type_combo.currentText(), name=dialog.name_edit.text())
        if "name" in fields:
            dialog.name_edit.setText(fields["name"])
        if "start" in fields:
            dialog.start_spin.setValue(fields["start"])
        if "end" in fields:
            dialog.end_spin.setValue(fields["end"])
        return QDialog.DialogCode.Accepted

    monkeypatch.setattr(FeatureDialog, "exec", fake_exec)
    return seen


def wait_for_tasks(window, qapp, timeout=30):
    deadline = time.time() + timeout
    while window.tasks and time.time() < deadline:
        qapp.processEvents()
        time.sleep(0.01)
    assert not window.tasks, "background task didn't finish"
    qapp.processEvents()


def labels(window):
    return [st.feature_label(f) for f in window.current_features]


# ------------------------------------------------------------------ opening
def test_open_file(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    assert window.feature_table.rowCount() == 10
    assert window.type_value.text() == "DNA"
    assert not window.modified


# ------------------------------------------------------------------ selection
def test_selection_readout_and_translation(window, qapp, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    window.viewer.select_range(600, 657)  # the multiple cloning site
    assert window.position_label.text().startswith("Selected 601–657  (57 bp)")
    assert window.circular_viewer.highlight == (600, 657)

    menu = window.build_selection_menu()
    [translate] = [a for a in menu.actions() if a.text() == "Translate selection"]
    translate.trigger()
    qapp.processEvents()
    result = window.findChildren(TextResultDialog)[-1]
    assert result.text == st.translate_dna(window.current_sequence[600:657])


# ------------------------------------------------------------------ features
def test_add_feature_from_selection_and_save(window, monkeypatch, tmp_path):
    path = open_demo(window, monkeypatch, tmp_path)
    window.viewer.select_range(100, 250)
    seen = fill_feature_dialog(monkeypatch, name="my feature")

    window.add_feature()

    assert (seen["start"], seen["end"]) == (101, 250)  # pre-filled from the selection
    assert "my feature" in labels(window)
    assert window.modified and window.windowTitle().endswith("*")

    assert window.save()
    assert not window.modified
    back = st.load_records(path)[0][0]
    assert "my feature" in [st.feature_label(f) for f in st.display_features(back)]


def test_edit_feature_across_origin(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    index = labels(window).index("AmpR")
    fill_feature_dialog(monkeypatch, name="AmpR moved", start=3901, end=60)

    window.edit_item(index)

    feature = window.current_features[labels(window).index("AmpR moved")]
    assert st.wraps_origin(feature, 4000)
    assert st.feature_strand(feature) == -1  # strand kept


def test_delete_feature(window, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    open_demo(window, monkeypatch, tmp_path)
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    window.delete_item(labels(window).index("T7 terminator"))
    assert "T7 terminator" not in labels(window)


def test_linear_sequence_blocks_origin_crossing(qapp):
    dialog = FeatureDialog(None, 4000, circular=False, start=3900, end=50)
    assert not dialog.buttons.button(dialog.buttons.StandardButton.Ok).isEnabled()


# ------------------------------------------------------------------ ORFs
def test_orf_finder_hides_annotated_genes(window, qapp, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    monkeypatch.setattr(QInputDialog, "getInt", staticmethod(lambda *a, **k: (100, True)))
    window.find_orfs()
    wait_for_tasks(window, qapp)
    # pDEMO's only long ORFs are GFP and AmpR, which are already annotated
    assert window.current_orfs == []
    assert "2 already annotated" in window.statusBar().currentMessage()


def test_orf_to_cds(window, qapp, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    monkeypatch.setattr(QInputDialog, "getInt", staticmethod(lambda *a, **k: (40, True)))
    window.find_orfs()
    wait_for_tasks(window, qapp)
    assert window.current_orfs

    first_orf = len(window.current_features)  # ORFs come after features in window.items
    seen = fill_feature_dialog(monkeypatch)
    window.add_orf_as_feature(first_orf)

    assert seen["type"] == "CDS" and seen["name"] == "ORF1"
    assert "ORF1" in labels(window)


# ------------------------------------------------------------------ restriction sites
def test_restriction_sites_on_map(window, qapp, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    window.find_restriction_sites()
    wait_for_tasks(window, qapp)
    on_map = " ".join(window.circular_viewer.sites)
    for enzyme in ("EcoRI", "BamHI", "HindIII"):
        assert enzyme in on_map


# ------------------------------------------------------------------ background tasks
def test_ncbi_fetch_runs_in_background(window, qapp, monkeypatch, tmp_path):
    """Uses a fake download, so the test works offline and never hits NCBI."""
    open_demo(window, monkeypatch, tmp_path)

    def fake_fetch(accession, email):
        time.sleep(0.5)  # pretend to be slow
        record = st.load_records(DEMO)[0][0]
        record.id = accession
        return record

    monkeypatch.setattr(st, "fetch_genbank", fake_fetch)
    answers = iter([("TEST123", True)])
    monkeypatch.setattr(QInputDialog, "getText", staticmethod(lambda *a, **k: next(answers)))
    window.settings.setValue("ncbi_email", "test@example.com")

    window.fetch_from_ncbi()
    assert not window.ncbi_button.isEnabled()  # busy
    wait_for_tasks(window, qapp)

    assert [r.id for r in window.records] == ["pDEMO", "TEST123"]
    assert window.ncbi_button.isEnabled()


# ------------------------------------------------------------------ map styles
def linear_copy_of_demo():
    record = st.load_records(DEMO)[0][0]
    record.annotations["topology"] = "linear"
    record.id = "linear-demo"
    return record


def test_map_style_follows_topology(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    assert window.map_stack.currentWidget() is window.circular_viewer
    window.show_records(window.records + [linear_copy_of_demo()], select=1)
    assert window.map_stack.currentWidget() is window.linear_viewer

    window.set_map_mode("circular")  # the user can override it
    assert window.map_stack.currentWidget() is window.circular_viewer
    window.set_map_mode("auto")
    assert window.map_stack.currentWidget() is window.linear_viewer


def test_both_maps_get_the_same_data(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    window.viewer.select_range(600, 657)
    for view in (window.circular_viewer, window.linear_viewer):
        assert view.highlight == (600, 657)
        assert len(view.features) == len(window.items)


def test_click_feature_on_linear_map(window, qapp, monkeypatch, tmp_path):
    window.show_records([linear_copy_of_demo()])
    view = window.linear_viewer
    view.grab()  # paint once: clickable areas are worked out while painting
    gfp = labels(window).index("GFP")
    [rect] = [r for i, r in view._hit_rects if i == gfp]
    assert view.feature_at(rect.center()) == gfp
    view.featureClicked.emit(view.feature_at(rect.center()))
    assert window.selected_item_index() == gfp


# ------------------------------------------------------------------ themes
def test_dark_and_light_theme(window, qapp, monkeypatch, tmp_path):
    from core import theme
    open_demo(window, monkeypatch, tmp_path)
    window.set_theme_mode("dark")
    qapp.processEvents()  # widgets pick up a new palette through events
    assert theme.current().dark
    assert window.palette().window().color().lightness() < 100
    window.circular_viewer.grab()  # paints without errors
    window.linear_viewer.grab()
    window.set_theme_mode("light")
    qapp.processEvents()
    assert not theme.current().dark
    assert window.palette().window().color().lightness() > 150
    assert window.settings.value("theme") == "light"


# ------------------------------------------------------------------ primers
def test_primer_dialog_adds_feature(window, monkeypatch, tmp_path):
    from core.primers import PrimerDialog
    open_demo(window, monkeypatch, tmp_path)
    seq = window.current_sequence
    dialog = PrimerDialog(window, seq, True, on_add=window.add_primer_features)
    dialog.name_edit.setText("GFP-fwd")
    dialog.primer_edit.setPlainText(seq[700:722])

    assert len(dialog.sites) == 1
    assert "1 binding site" in dialog.info.text()
    dialog.add_selected()

    feature = window.current_features[labels(window).index("GFP-fwd")]
    assert feature.type == "primer_bind"
    assert st.feature_span(feature) == (701, 722)
    assert window.modified


def test_primer_dialog_rejects_bad_letters(window, monkeypatch, tmp_path):
    from core.primers import PrimerDialog
    open_demo(window, monkeypatch, tmp_path)
    dialog = PrimerDialog(window, window.current_sequence, True, on_add=lambda f: None)
    dialog.primer_edit.setPlainText("ACGTXXACGTACGT")
    assert dialog.sites == [] and "only contain" in dialog.info.text()
    assert not dialog.add_button.isEnabled()


# ------------------------------------------------------------------ digest
def test_digest_dialog(window, qapp, monkeypatch, tmp_path):
    from core.digest import DigestDialog
    open_demo(window, monkeypatch, tmp_path)
    opened = []
    monkeypatch.setattr(DigestDialog, "exec", lambda self: opened.append(self) or 0)

    window.open_digest()
    wait_for_tasks(window, qapp)
    [dialog] = opened

    dialog.set_checked(["EcoRI", "HindIII"])
    dialog.add_lane_from_checked()
    assert dialog.lanes[0]["enzymes"] == ["EcoRI", "HindIII"]
    assert sorted(f["length"] for f in dialog.lanes[0]["fragments"]) == [51, 3949]
    assert dialog.fragment_table.rowCount() == 2
    assert dialog.gel.lanes == [("1", [51, 3949])]
    assert dialog.checked_enzymes() == []  # checks are cleared for the next lane

    dialog.add_lane(["BamHI"])
    dialog.lane_list.setCurrentRow(0)
    dialog.remove_lane()
    assert [lane["enzymes"] for lane in dialog.lanes] == [["BamHI"]]

    image = tmp_path / "gel.png"
    monkeypatch.setattr(QFileDialog, "getSaveFileName", staticmethod(lambda *a, **k: (str(image), "")))
    dialog.save_image()
    assert image.exists() and image.stat().st_size > 1000


def test_digest_reuses_restriction_search(window, qapp, monkeypatch, tmp_path):
    from core.digest import DigestDialog
    open_demo(window, monkeypatch, tmp_path)
    window.find_restriction_sites()
    wait_for_tasks(window, qapp)
    monkeypatch.setattr(DigestDialog, "exec", lambda self: 0)
    window.open_digest()
    assert "digest" not in window.tasks  # opened straight away from the cached search
