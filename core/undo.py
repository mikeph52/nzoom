"""
undo.py — undoable feature changes for Qt's QUndoStack.

Each command changes a record's features in redo() and changes them back in
undo(). QUndoStack.push() calls redo() straight away, so making a change and
recording it for Undo are the same step.

After every redo/undo the command calls on_change(record, feature, message),
which the main window uses to redraw. `feature` is the feature to select
afterwards, or None.
"""
import copy

from PySide6.QtGui import QUndoCommand

from . import seqtools as st


def snapshot(feature):
    """Everything an edit can change, copied so later edits don't alter it."""
    return feature.type, copy.deepcopy(feature.qualifiers), feature.location


def restore(feature, state):
    feature.type, qualifiers, feature.location = state
    feature.qualifiers = copy.deepcopy(qualifiers)


def remove_by_identity(features, feature):
    """Remove `feature` itself (two features can be equal without being the same)."""
    for position, candidate in enumerate(features):
        if candidate is feature:
            del features[position]
            return position
    return None


class FeatureCommand(QUndoCommand):
    def __init__(self, text, record, on_change):
        super().__init__(text)
        self.record = record
        self.on_change = on_change


class AddFeatures(FeatureCommand):
    def __init__(self, record, features, on_change, text):
        super().__init__(text, record, on_change)
        self.features = list(features)

    def redo(self):
        self.record.features.extend(self.features)
        self.on_change(self.record, self.features[0], self.text())

    def undo(self):
        for feature in self.features:
            remove_by_identity(self.record.features, feature)
        self.on_change(self.record, None, f"Undid: {self.text()}")

class DeleteFeature(FeatureCommand):
    def __init__(self, record, feature, on_change):
        super().__init__(f"Delete {feature.type} “{st.feature_label(feature)}”",
                         record, on_change)
        self.feature = feature
        self.position = None

    def redo(self):
        self.position = remove_by_identity(self.record.features, self.feature)
        self.on_change(self.record, None, self.text())

    def undo(self):
        position = len(self.record.features) if self.position is None else self.position
        self.record.features.insert(position, self.feature)  # back in its old place
        self.on_change(self.record, self.feature, f"Undid: {self.text()}")

class EditFeature(FeatureCommand):
    """Records a change already made to `feature` (before -> after states)."""

    def __init__(self, record, feature, before, after, on_change):
        super().__init__(f"Edit {feature.type} “{st.feature_label(feature)}”",
                         record, on_change)
        self.feature = feature
        self.before = before
        self.after = after

    def redo(self):
        restore(self.feature, self.after)
        self.on_change(self.record, self.feature, self.text())

    def undo(self):
        restore(self.feature, self.before)
        self.on_change(self.record, self.feature, f"Undid: {self.text()}")
