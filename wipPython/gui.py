import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QGridLayout
from PyQt5.QtWidgets import QLineEdit
from PyQt5.QtWidgets import QPushButton
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QHBoxLayout
from PyQt5.QtWidgets import QStackedWidget
from PyQt5.QtWidgets import QAction
from functools import partial

__version__ = 0.1
__author__ = 'Jose Mauricio Rodas Rodriguez'

class UIWindow(QWidget):
    def __init__(self, parent=None):
        super(UIWindow, self).__init__(parent)
        # mainwindow.setWindowIcon(QtGui.QIcon('PhotoIcon.png'))
        self.ToolsBTN = QPushButton('text', self)
        self.ToolsBTN.move(50, 350)


class UIToolTab(QWidget):
    def __init__(self, parent=None):
        super(UIToolTab, self).__init__(parent)
        self.CPSBTN = QPushButton("text2", self)
        self.CPSBTN.move(100, 350)


class MainWindow(QMainWindow):
    def __init__(self, parent=None):
        super(MainWindow, self).__init__(parent)
        self.setGeometry(50, 50, 400, 450)
        self.setFixedSize(400, 450)
        self.startUIWindow()

    def startUIToolTab(self):
        self.ToolTab = UIToolTab(self)
        self.setWindowTitle("UIToolTab")
        self.setCentralWidget(self.ToolTab)
        self.ToolTab.CPSBTN.clicked.connect(self.startUIWindow)
        self.show()

    def startUIWindow(self):
        self.Window = UIWindow(self)
        self.setWindowTitle("UIWindow")
        self.setCentralWidget(self.Window)
        self.Window.ToolsBTN.clicked.connect(self.startUIToolTab)
        self.show()


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    sys.exit(app.exec_())