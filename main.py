from rdkit import Chem
from rdkit.Chem import rdchem, Draw, AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize

from PyQt5.QtGui import QPixmap, QPalette, QColor, QImage, QPainter
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QFileDialog, QHBoxLayout, QLabel, QDialog, QTabWidget, QComboBox, QAction
from PyQt5.QtSvg import QSvgWidget, QSvgRenderer
from PyQt5.QtCore import Qt, QByteArray

from muropeptide import Muropeptide
from murotab import MuroTab

class muroview(QMainWindow):
    def __init__(self):
        super().__init__()

        # set window title
        self.setWindowTitle("muroview")
        self.resize(800, 600)

        self.currentMuro = Muropeptide()

        # create central widget
        self.central_widget = QWidget(self)

        # create 2D canvas
        self.canvas = QLabel(self.central_widget)
        self.canvas.setAlignment(Qt.AlignTop)
        self.canvas.setFixedSize(600, 600)

        # create layout
        self.display_layout = QVBoxLayout()
        self.display_layout.addWidget(self.canvas)

        # set white background
        pal = QPalette()
        pal.setColor(QPalette.Background, QColor(255, 255, 255))
        self.canvas.setAutoFillBackground(True)
        self.canvas.setPalette(pal)

        self.list_tab = QWidget(self.central_widget)
        self.list_layout = QVBoxLayout(self.list_tab)
        self.list_layout.setAlignment(Qt.AlignTop)
        self.list_layout.setSpacing(0)

        self.murotab = MuroTab(self.currentMuro, self.central_widget, self.canvas)

        self.tab_widget = QTabWidget(self.central_widget)
        self.tab_widget.addTab(self.murotab, "Muropeptide")
        #self.tab_widget.addTab(self.list_tab, "List")

        # create main layout
        self.main_layout = QHBoxLayout(self.central_widget)
        self.main_layout.addLayout(self.display_layout)
        self.main_layout.addWidget(self.tab_widget)

        # set central widget
        self.setCentralWidget(self.central_widget)

        # create menu bar
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("File")

        # create submenu for opening files
        save_as_action = QAction("Save as...", self)
        save_as_action.triggered.connect(self.save_as)
        file_menu.addAction(save_as_action)

    def save_as(self):
        print("TODO")


if __name__ == "__main__":
    app = QApplication([])
    viewer = muroview()
    viewer.show()
    app.exec_()
