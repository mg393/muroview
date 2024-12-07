from rdkit.Chem import rdchem, rdmolfiles, Descriptors, rdMolDescriptors
from pyteomics import mass

class GlcNAc():
    def __init__(self):
        self.NDeAc = False
        self.OAc = False

    def setNAc(self, deacetylated):
        self.NDeAc = deacetylated

    def setOAc(self, oacetylated):
        self.OAc = oacetylated

    def smiles(self):
        if self.NDeAc:
            print("yeah boi")
        elif self.OAc:
            print("yeah boi")
        elif self.NDeAc and self.OAc:
            print("yeah boi")
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
        

class Muropeptide(rdchem.Mol):
    def __init__(self):
        self.glcnac = GlcNAc()
        self.murnac = MurNAc()

        combinedSmiles = self.glcnac.smiles() + "." + self.murnac.smiles()
        print(combinedSmiles)
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
        

        
