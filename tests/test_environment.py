"""
Checks that the installed libraries work correctly, separately from N-Zoom's
own code. If one of these fails, the problem is in a library version, and the
fix is usually in requirements.txt.
"""
import gc
import sys

import pytest


@pytest.mark.skipif(sys.version_info >= (3, 12),
                    reason="None can't be freed on Python 3.12+, so the leak is harmless there")
def test_pyside_does_not_leak_none(qapp):
    """PySide6 6.12.0 on Python 3.10/3.11 drops a reference to None on every
    method call that returns nothing, which soon crashes Python with
    'Fatal Python error: none_dealloc'. See requirements.txt."""
    from PySide6.QtWidgets import QLabel

    label = QLabel()
    gc.collect()
    before = sys.getrefcount(None)
    for _ in range(1000):
        label.setText("x")
    gc.collect()
    lost = before - sys.getrefcount(None)
    assert lost < 100, (
        f"This PySide6 version loses references to None ({lost} in 1,000 calls) "
        "and will crash on this Python version. Install an older PySide6, as "
        "listed in requirements.txt."
    )