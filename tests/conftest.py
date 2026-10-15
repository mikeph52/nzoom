"""
Shared setup for all tests. pytest loads this file automatically.

- Makes `import core...` and `import main` work from inside tests/.
- Runs Qt without opening real windows ("offscreen"), so tests also work
  on servers and in CI.
- Provides fixtures: small reusable pieces a test asks for by name.
"""
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

DEMO = ROOT / "examples" / "pDEMO.gb"


@pytest.fixture(scope="session")
def qapp():
    """One QApplication for the whole test run (Qt allows only one)."""
    from PySide6.QtWidgets import QApplication
    return QApplication.instance() or QApplication([])


@pytest.fixture
def demo_record():
    """A fresh copy of the demo plasmid for each test."""
    from core import seqtools as st
    return st.load_records(DEMO)[0][0]


@pytest.fixture
def window(qapp, monkeypatch, tmp_path):
    """A main window that never blocks on a pop-up.

    Questions ("Save changes?") answer Discard, and information boxes are
    collected in window.infos instead of being shown. Settings go to a
    temporary file, so tests never touch your real N-Zoom settings.
    """
    from PySide6.QtCore import QSettings
    from PySide6.QtWidgets import QMessageBox
    import main

    monkeypatch.setattr(QMessageBox, "question",
                        staticmethod(lambda *a, **k: QMessageBox.StandardButton.Discard))
    infos = []
    monkeypatch.setattr(QMessageBox, "information", staticmethod(lambda *a, **k: infos.append(a[2])))
    monkeypatch.setattr(QMessageBox, "critical", staticmethod(lambda *a, **k: infos.append(a[2])))

    w = main.MainWindow()
    w.settings = QSettings(str(tmp_path / "settings.ini"), QSettings.Format.IniFormat)
    w.infos = infos
    w.resize(1300, 820)
    w.show()
    yield w
    w.modified = False  # don't ask about saving on close
    w.close()
