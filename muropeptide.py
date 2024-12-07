from rdkit.Chem import rdchem, rdmolfiles, Descriptors, rdMolDescriptors
from rdkit.Chem import rdchem, Draw, AllChem, Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

from pyteomics import mass

from PyQt5.QtSvg import QSvgWidget, QSvgRenderer
from PyQt5.QtGui import QPixmap, QPalette, QColor, QImage, QPainter
from PyQt5.QtCore import Qt, QByteArray

class GlcNAc():
    def __init__(self):
        self.NDeAc = False
        self.OAc = False

    def setNAc(self, deacetylated):
        self.NDeAc = deacetylated

    def setOAc(self, oacetylated):
        self.OAc = oacetylated

    def smiles(self):
        if self.NDeAc == True:
            return "N[C@@H]1[C@H]([C@@H]([C@H](OC1O(*))CO)O)O"
        else: 
            return "CC(=O)N[C@@H]1[C@H]([C@@H]([C@H](OC1O(*))CO)O)O"

class MurNAc():
    def __init__(self):
        self.NDeAc = False
        self.OAc = False
        self.anhydro = False
        self.reduced = True

    def setNAc(self, deacetylated):
        self.NDeAc = deacetylated

    def setOAc(self, oacetylated):
        self.OAc = oacetylated

    def anhydro(self, anhydrom):
        self.anhydro = anhydrom

    def reducec(self, reducedm):
        self.reduced = reducedm

    def smiles(self):
        if self.NDeAc:
            print("yeah boi")
        elif self.OAc:
            print("yeah boi")
        elif self.NDeAc and self.OAc:
            print("yeah boi")
        else: 
            return "O=C(O)[C@H](O[C@H]1[C@H](*)[C@H](OC(O)[C@@H]1NC(=O)C)CO)C"

class AminoAcid():
    def __init__(self):
        print("TODO")


class Muropeptide(rdchem.Mol):
    def __init__(self):
        self.glcnac = GlcNAc()
        self.murnac = MurNAc()
        self.aa1 = AminoAcid()
        self.calc()


    def calc(self):
        print(self.glcnac.smiles())
        combinedSmiles = self.glcnac.smiles() + "." + self.murnac.smiles()
        muro = rdmolfiles.MolFromSmiles(combinedSmiles.replace('(*)', '9'))
        print(rdmolfiles.MolToSmiles(muro))
        rdchem.Mol.__init__(self, muro)

        self.monoisotopicMass = Descriptors.ExactMolWt(muro)
        self.averageMass = Descriptors.MolWt(muro)
        self.formula = rdMolDescriptors.CalcMolFormula(muro, False)
        self.m1_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 1)
        self.m2_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 2)
        self.m3_pve = mass.calculate_mass(formula = self.formula, ion_type="M", charge = 3)
        self.mNa_pve = mass.calculate_mass(formula = self.formula, charge_carrier="Na+", ion_type="M", charge = 1)
        self.mK_pve = mass.calculate_mass(formula = self.formula, charge_carrier="K+", ion_type="M", charge = 1)
        self.m1_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -1)
        self.m2_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -2)
        self.m3_nve = -mass.calculate_mass(formula = self.formula, ion_type="M", charge = -3)

    def makeSVG(self):
        AllChem.Compute2DCoords(self)
        drawer = rdMolDraw2D.MolDraw2DSVG(600, 600)
        drawer.ClearDrawing()
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, self)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        return QByteArray(svg.encode('utf-8'))

    def draw(self, canvas):
        svg_renderer = QSvgRenderer(self.makeSVG())
        img = QImage(canvas.size(), QImage.Format_ARGB32)
        img.fill(Qt.white)
        painter = QPainter(img)
        svg_renderer.render(painter)
        painter.end()
        return QPixmap.fromImage(img)

        

        
