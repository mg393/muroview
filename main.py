from rdkit import Chem
from rdkit.Chem import rdchem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Descriptors

from PyQt5.QtGui import QPixmap, QPalette, QColor, QImage, QPainter
from PyQt5.QtWidgets import QApplication, QMainWindow, QMenu, QVBoxLayout, QSizePolicy, QMessageBox, QWidget, QPushButton, QFileDialog, QHBoxLayout, QLabel, QDialog, QTabWidget, QComboBox
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

        self.murotab = MuroTab(self.currentMuro, self.central_widget)

        self.tab_widget = QTabWidget(self.central_widget)
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
        save_menu = QMenu("Save as SMILES", self)
        file_menu.addMenu(save_menu)

    def draw_molecule2D(self, muro):
        AllChem.Compute2DCoords(muro)
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
        drawer.ClearDrawing()
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, muro)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        svg_data = QByteArray(svg.encode('utf-8'))
        svg_renderer = QSvgRenderer(svg_data)
        img = QImage(self.canvas.size(), QImage.Format_ARGB32)
        img.fill(Qt.white)
        painter = QPainter(img)
        svg_renderer.render(painter)
        painter.end()
        pixmap = QPixmap.fromImage(img)
        self.canvas.setPixmap(pixmap)
        
    def load_molecule(self, muro):
        # clean up and uncharge molecule
        mol = rdMolStandardize.Cleanup(muro)
        uncharger = rdMolStandardize.Uncharger()
        mol = uncharger.uncharge(mol)

        self.murotab = MuroTab(muro, self.central_widget)
        self.murotab.molecularWeightLabel.setText("<b>Average mass:</b> " + str(round(self.murotab.muropeptide.averageMass, 4)))
        self.murotab.monoisotopicMwLabel.setText("<b>Monoisotopic mass:</b> " + str(round(self.murotab.muropeptide.monoisotopicMass, 4)))
        self.murotab.m1_pve_Label.setText("[M+H]<sup>+</sup>: " + str(round(self.murotab.muropeptide.m1_pve, 4)))
        self.murotab.m2_pve_Label.setText("[M+2H]<sup>2+</sup>: " + str(round(self.murotab.muropeptide.m2_pve, 4)))
        self.murotab.m3_pve_Label.setText("[M+3H]<sup>3+</sup>: " + str(round(self.murotab.muropeptide.m3_pve, 4)))
        self.murotab.mNa_pve_Label.setText("[M+Na]<sup>+</sup>: " + str(round(self.murotab.muropeptide.mNa_pve, 4)))
        self.murotab.mK_pve_label.setText("[M+K]<sup>+</sup>: " + str(round(self.murotab.muropeptide.mK_pve, 4)))
        self.murotab.m1_nve_Label.setText("[M-H]<sup>-</sup>: " + str(round(-1*self.murotab.muropeptide.m1_nve, 4)))
        self.murotab.m2_nve_Label.setText("[M-2H]<sup>2-</sup>: " + str(round(-1*self.murotab.muropeptide.m2_nve, 4)))
        self.murotab.m3_nve_Label.setText("[M-3H]<sup>3-</sup>: " + str(round(-1*self.murotab.muropeptide.m3_nve, 4)))

        self.tab_widget.addTab(self.murotab, "Muropeptide")
        self.draw_molecule2D(self.murotab.muropeptide)

if __name__ == "__main__":
    app = QApplication([])
    viewer = muroview()
    viewer.show()
    viewer.load_molecule(muro=Muropeptide())
    app.exec_()
