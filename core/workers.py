"""
workers.py — run slow functions (network, heavy analysis) off the GUI thread,
so the window keeps responding.

Usage from a QObject (e.g. the main window):

    task = Task(slow_function, arg1, arg2,
                on_done=lambda result: ...,      # called in the GUI thread
                on_fail=lambda message: ...)     # called in the GUI thread
    task.signals.finished.connect(self.on_task_finished)
    task.signals.failed.connect(self.on_task_failed)
    pool.start(task)

`on_task_finished` / `on_task_failed` must be methods of a QObject (the main
window). Qt then runs them in the GUI thread, and they call task.on_done /
task.on_fail there. Never touch widgets from inside `slow_function` itself:
it runs on the background thread.
"""
from PySide6.QtCore import QObject, QRunnable, Signal


class TaskSignals(QObject):
    finished = Signal(object)  # payload: (task, result)
    failed = Signal(object)    # payload: (task, error message)


class Task(QRunnable):
    def __init__(self, fn, *args, on_done=None, on_fail=None, **kwargs):
        super().__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.on_done = on_done
        self.on_fail = on_fail
        self.signals = TaskSignals()
        self.setAutoDelete(False)  # Python keeps the object alive, not Qt

    def run(self):  # runs on a background thread
        try:
            result = self.fn(*self.args, **self.kwargs)
        except Exception as exc:  # report every error back to the GUI
            self.signals.failed.emit((self, f"{type(exc).__name__}: {exc}"))
        else:
            self.signals.finished.emit((self, result))