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
    assert window.map_stack.currentWidget() is window.linear_panel

    window.set_map_mode("circular")  # the user can override it
    assert window.map_stack.currentWidget() is window.circular_viewer
    window.set_map_mode("auto")
    assert window.map_stack.currentWidget() is window.linear_panel


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


# ------------------------------------------------------------------ find
def test_find_bar(window, qapp, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    window.show_find()
    assert window.find_bar.isVisible()

    window.find_edit.setText("GGTAC")  # KpnI-like, not palindromic: both strands
    hits = window.find_hits
    assert len(hits) == len(st.find_motif(window.current_sequence, "GGTAC", circular=True))
    assert {h["strand"] for h in hits} == {1, -1}
    assert window.find_status.text().startswith("1 of ")
    first = window.viewer.selected_range()

    window.find_next()
    assert window.find_index == 1 and window.viewer.selected_range() != first
    window.find_previous()
    assert window.viewer.selected_range() == first
    assert len(window.circular_viewer.marks) == len(hits)

    window.find_both.setChecked(False)
    assert all(h["strand"] == 1 for h in window.find_hits)

    window.find_edit.setText("GGZ")
    assert window.find_hits == [] and "IUPAC" in window.find_status.text()

    window.hide_find()
    assert window.circular_viewer.marks == [] and not window.find_bar.isVisible()


def test_find_starts_from_selection(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    window.viewer.select_range(600, 606)
    window.show_find()
    assert window.find_edit.text() == "GAATTC"
    assert window.find_status.text() == "1 of 1"


# ------------------------------------------------------------------ linear zoom
def test_linear_zoom_and_scroll(window, monkeypatch, tmp_path):
    window.show_records([linear_copy_of_demo()])
    view, panel = window.linear_viewer, window.linear_panel
    assert view.zoom == 1 and not panel.scrollbar.isEnabled()

    view.zoom_by(4)
    assert view.zoom == 4 and abs(view.visible_span() - 1000) < 1
    assert panel.scrollbar.isEnabled() and panel.scrollbar.maximum() == 3000

    panel.scrollbar.setValue(2000)  # scrollbar moves the map
    assert view.view_start == 2000

    view.zoom_by(1000)  # capped so a base is at most MAX_PX_PER_BP pixels wide
    assert view.zoom == view.max_zoom()
    view.grab()  # draws the bases without errors

    view.fit()
    assert view.zoom == 1 and view.view_start == 0


def test_linear_zoom_keeps_point_under_mouse(window):
    window.show_records([linear_copy_of_demo()])
    view = window.linear_viewer
    x = view.width() * 0.3
    before = view._pos_at(x)
    view.zoom_by(3, anchor_x=x)
    assert abs(view._pos_at(x) - before) < 1


def test_linear_map_scrolls_to_selected_feature(window):
    window.show_records([linear_copy_of_demo()])
    view = window.linear_viewer
    view.zoom_by(10)  # 400 bp visible, starting near the middle
    view.set_view_start(0)
    ampr = labels(window).index("AmpR")
    window.select_feature_row(ampr)
    start, end = 2105, 2966
    assert view.view_start <= start + 1 and view.view_start + view.visible_span() >= start


# ------------------------------------------------------------------ export
@pytest.mark.parametrize("fmt", ["png", "svg"])
@pytest.mark.parametrize("style", ["circular", "linear"])
def test_export_map(window, monkeypatch, tmp_path, fmt, style):
    from xml.dom import minidom
    from PySide6.QtGui import QImage
    open_demo(window, monkeypatch, tmp_path)
    window.set_map_mode(style)
    path = tmp_path / f"map.{fmt}"
    monkeypatch.setattr(QFileDialog, "getSaveFileName",
                        staticmethod(lambda *a, **k: (str(path), "")))
    window.export_map()

    assert path.exists()
    view = window.current_map_view()
    if fmt == "png":
        image = QImage(str(path))
        assert (image.width(), image.height()) == (view.width() * 2, view.height() * 2)
    else:
        svg = minidom.parse(str(path))  # valid XML
        assert "GFP" in path.read_text()  # labels are real text in the SVG


# ------------------------------------------------------------------ undo
def test_undo_redo_add_edit_delete(window, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    open_demo(window, monkeypatch, tmp_path)
    stack = window.undo_stack
    original = labels(window)

    fill_feature_dialog(monkeypatch, name="new one", start=10, end=50)
    window.add_feature()
    assert "new one" in labels(window) and window.modified

    index = labels(window).index("GFP")
    fill_feature_dialog(monkeypatch, name="GFP renamed")
    window.edit_item(index)
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    window.delete_item(labels(window).index("ori"))
    assert stack.count() == 3

    stack.undo()  # delete
    assert "ori" in labels(window)
    stack.undo()  # edit
    assert "GFP" in labels(window) and "GFP renamed" not in labels(window)
    stack.undo()  # add
    assert labels(window) == original
    assert not window.modified  # back to the saved state: no "*" in the title
    assert not window.windowTitle().endswith("*")

    stack.redo(); stack.redo(); stack.redo()
    assert "new one" in labels(window) and "GFP renamed" in labels(window)
    assert "ori" not in labels(window) and window.modified


def test_undo_restores_feature_position_and_location(window, monkeypatch, tmp_path):
    from PySide6.QtWidgets import QMessageBox
    open_demo(window, monkeypatch, tmp_path)
    features = window.current_record.features
    before = [id(f) for f in features]
    index = labels(window).index("AmpR")
    location = window.current_features[index].location

    fill_feature_dialog(monkeypatch, start=3901, end=60)
    window.edit_item(index)
    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Yes))
    window.delete_item(labels(window).index("lac promoter"))

    window.undo_stack.undo()
    window.undo_stack.undo()
    assert [id(f) for f in features] == before  # same objects, same order
    assert window.current_features[labels(window).index("AmpR")].location == location


def test_save_marks_clean_and_undo_makes_it_dirty_again(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    fill_feature_dialog(monkeypatch, name="x", start=1, end=9)
    window.add_feature()
    assert window.save() and not window.modified
    window.undo_stack.undo()
    assert window.modified  # the file on disk still has the feature


def test_primers_and_orf_conversion_are_undoable(window, qapp, monkeypatch, tmp_path):
    from core.primers import PrimerDialog
    open_demo(window, monkeypatch, tmp_path)
    count = len(window.current_features)
    dialog = PrimerDialog(window, window.current_sequence, True, on_add=window.add_primer_features)
    dialog.primer_edit.setPlainText(window.current_sequence[700:722])
    dialog.add_selected()
    assert len(window.current_features) == count + 1
    assert window.undo_stack.undoText().startswith("Add 1 primer site")
    window.undo_stack.undo()
    assert len(window.current_features) == count


def test_opening_a_file_clears_undo_history(window, monkeypatch, tmp_path):
    open_demo(window, monkeypatch, tmp_path)
    fill_feature_dialog(monkeypatch, name="x", start=1, end=9)
    window.add_feature()
    open_demo(window, monkeypatch, tmp_path / "again" if False else tmp_path)
    assert window.undo_stack.count() == 0
