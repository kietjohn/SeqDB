from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s


class TestListView(QtGui.QListWidget):
    def __init__(self, type, parent=None):
        super(TestListView, self).__init__(parent)
        self.setAcceptDrops(True)
        self.setIconSize(QtCore.QSize(72, 72))
        self.file_list = []

    def dragEnterEvent(self, event):
        if event.mimeData().hasUrls:
            event.accept()
        else:
            event.ignore()

    def dragMoveEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
        else:
            event.ignore()

    def dropEvent(self, event):
        if event.mimeData().hasUrls:
            event.setDropAction(QtCore.Qt.CopyAction)
            event.accept()
            for url in event.mimeData().urls():
                self.file_list.append(str(url.toLocalFile()))
            self.emit(QtCore.SIGNAL("dropped"), self.file_list)
        else:
            event.ignore()


class ManualListView(QtGui.QListWidget):
    def __init__(self, type, parent=None):
        super(ManualListView, self).__init__(parent)
        self.setAcceptDrops(False)
        self.setIconSize(QtCore.QSize(72, 72))
        self.file_list = []