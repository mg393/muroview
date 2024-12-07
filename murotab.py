from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QLabel, QHBoxLayout, QVBoxLayout, QComboBox

class MuroTab(QWidget):
	def __init__(self, muro, central_widget, canvas):
		super().__init__()
		self.muropeptide = muro
		self.canvas = canvas
		self.molecularWeightLabel = QLabel()
		self.monoisotopicMwLabel = QLabel()
		self.space1Label = QLabel("")
		self.pve_charge_Label = QLabel("<u><b>Positive Charge States</b></u>")
		self.m1_pve_Label = QLabel()
		self.mNa_pve_Label = QLabel()
		self.mK_pve_label = QLabel()
		self.space2Label = QLabel("")
		self.m2_pve_Label = QLabel()
		self.m3_pve_Label = QLabel()
		self.space3Label = QLabel("")
		self.nve_charge_Label = QLabel("<u><b>Negative Charge States</b></u>")
		self.m1_nve_Label = QLabel()
		self.m2_nve_Label = QLabel()
		self.m3_nve_Label = QLabel()

		self.main_layout = QVBoxLayout(central_widget)
		self.main_layout.setAlignment(Qt.AlignTop)
		self.main_layout.setSpacing(0)

		self.GlcNAc_layout = QHBoxLayout()
		self.GlcNAc_label = QLabel("GlcNAc: ")
		self.GlcNAc_dropdown = QComboBox()
		self.GlcNAc_dropdown.addItems(["Unmodified", "N-deacetylated", "O-acetylated"])
		self.GlcNAc_layout.addWidget(self.GlcNAc_label)
		self.GlcNAc_layout.addWidget(self.GlcNAc_dropdown)

		self.MurNAc_layout = QHBoxLayout()
		self.MurNAc_label = QLabel("MurNAc: ")
		self.MurNAc_dropdown = QComboBox()
		self.MurNAc_dropdown.addItems(["Unmodified", "N-deacetylated", "O-acetylated"])
		self.MurNAc_layout.addWidget(self.MurNAc_label)
		self.MurNAc_layout.addWidget(self.MurNAc_dropdown)

		self.GlcNAc_dropdown.currentIndexChanged.connect(self.on_GlcNAc_dropdown_changed)
		self.MurNAc_dropdown.currentIndexChanged.connect(self.on_MurNAc_dropdown_changed)

		self.main_layout.addLayout(self.GlcNAc_layout)
		self.main_layout.addLayout(self.MurNAc_layout)
		self.main_layout.addWidget(self.molecularWeightLabel)
		self.main_layout.addWidget(self.monoisotopicMwLabel)
		self.main_layout.addWidget(self.space1Label)
		self.main_layout.addWidget(self.pve_charge_Label)
		self.main_layout.addWidget(self.m1_pve_Label)
		self.main_layout.addWidget(self.mNa_pve_Label)
		self.main_layout.addWidget(self.mK_pve_label)
		self.main_layout.addWidget(self.space2Label)
		self.main_layout.addWidget(self.m2_pve_Label)
		self.main_layout.addWidget(self.m3_pve_Label)
		self.main_layout.addWidget(self.space3Label)
		self.main_layout.addWidget(self.nve_charge_Label)
		self.main_layout.addWidget(self.m1_nve_Label)
		self.main_layout.addWidget(self.m2_nve_Label)
		self.main_layout.addWidget(self.m3_nve_Label)
		self.setLayout(self.main_layout)
		self.update()

	def update(self):
		self.muropeptide.calc()
		self.molecularWeightLabel.setText("<b>Average mass:</b> " + str(round(self.muropeptide.averageMass, 4)))
		self.monoisotopicMwLabel.setText("<b>Monoisotopic mass:</b> " + str(round(self.muropeptide.monoisotopicMass, 4)))
		self.m1_pve_Label.setText("[M+H]<sup>+</sup>: " + str(round(self.muropeptide.m1_pve, 4)))
		self.m2_pve_Label.setText("[M+2H]<sup>2+</sup>: " + str(round(self.muropeptide.m2_pve, 4)))
		self.m3_pve_Label.setText("[M+3H]<sup>3+</sup>: " + str(round(self.muropeptide.m3_pve, 4)))
		self.mNa_pve_Label.setText("[M+Na]<sup>+</sup>: " + str(round(self.muropeptide.mNa_pve, 4)))
		self.mK_pve_label.setText("[M+K]<sup>+</sup>: " + str(round(self.muropeptide.mK_pve, 4)))
		self.m1_nve_Label.setText("[M-H]<sup>-</sup>: " + str(round(-1*self.muropeptide.m1_nve, 4)))
		self.m2_nve_Label.setText("[M-2H]<sup>2-</sup>: " + str(round(-1*self.muropeptide.m2_nve, 4)))
		self.m3_nve_Label.setText("[M-3H]<sup>3-</sup>: " + str(round(-1*self.muropeptide.m3_nve, 4)))
		self.canvas.setPixmap(self.muropeptide.draw(self.canvas))


	def on_GlcNAc_dropdown_changed(self, index):
		selected_item = self.GlcNAc_dropdown.itemText(index)
		if selected_item == "N-deacetylated":
			self.muropeptide.glcnac.setNAc(True)
		else: 
			self.muropeptide.glcnac.setNAc(False)
		print(self.muropeptide.glcnac.NDeAc)
		self.update()
		print(f"GlcNAc dropdown changed: Selected {selected_item}")

	def on_MurNAc_dropdown_changed(self, index):
		selected_item = self.MurNAc_dropdown.itemText(index)
		print(f"MurNAc dropdown changed: Selected {selected_item}")
