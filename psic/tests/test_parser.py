import unittest
import sys
sys.path.append("..")
import psic

class TestParser(unittest.TestCase):
	def setUp(self):
		self.structure = psic.Parser().parse("../data/test.pdb")
		self.chain = self.structure.chains[0]
		self.residue = self.structure.residues[0]
		self.atom = self.structure.atoms[0]
		# self.mol_structure = psic.Parser().parse("../data/crystal_ligand.mol2")
		# self.mol_chain = self.mol_structure.chains[0]
		# self.mol_residue = self.mol_structure.residues[0]
		# self.mol_atom = self.mol_structure.atoms[0]

	def test_structure(self):
		self.assertEqual(self.structure.name, "test")
		self.assertEqual(len(self.structure.atoms), 4783)
		self.assertEqual(len(self.structure.residues), 760)
		self.assertEqual(len(self.structure.chains), 4)

	def test_chain(self):
		self.assertEqual(self.structure.name, "test")
		self.assertEqual(len(self.chain.atoms), 2256)
		self.assertEqual(self.chain.atoms[0].name, "N")
		self.assertEqual(self.chain.atoms[-1].name, "O")
		self.assertEqual(len(self.chain.residues), 367)
		self.assertEqual(self.chain.residues[-1].name, "HOH")
		self.assertEqual(self.chain.residues[0].name, "ALA")
		self.assertEqual(self.chain.structure.name, "test")
		

	def test_residue(self):
		self.assertEqual(self.residue.structure.name, "test")
		self.assertEqual(self.residue.chain.name, "A")
		self.assertEqual(self.residue.number, "523")
		self.assertEqual(self.residue.name, "ALA")
		self.assertEqual(len(self.residue.atoms), 5)

	def test_atom(self):
		self.assertEqual(self.atom.identifier, "ATOM")
		self.assertEqual(self.atom.name, "N")
		self.assertEqual(self.atom.residue_name, "ALA")
		self.assertEqual(self.atom.chain_name, "A")
		self.assertEqual(self.atom.residue_number, "523")
		self.assertEqual(self.atom.x, "65.336")
		self.assertEqual(self.atom.y, "11.629")
		self.assertEqual(self.atom.z, "-15.337")
		self.assertEqual(self.atom.occupancy, "1.00")
		self.assertEqual(self.atom.temperature_factor, "126.58")
		self.assertEqual(self.atom.segment_id, "")
		self.assertEqual(self.atom.element, "N")
		self.assertEqual(self.atom.structure.name, "test")
		self.assertEqual(self.atom.chain.name, "A")
		self.assertEqual(self.atom.residue.name, "ALA")

	# def test_mol_structure(self):
	# 	self.assertEqual(self.mol_structure.name, "crystal_ligand")
	# 	self.assertEqual(len(self.mol_structure.atoms), 44)
	# 	self.assertEqual(len(self.mol_structure.residues), 1)
	# 	self.assertEqual(len(self.mol_structure.chains), 1)

	# def test_mol_chain(self):
	# 	self.assertEqual(len(self.mol_chain.atoms), 44)
	# 	self.assertEqual(self.mol_chain.atoms[0].name, "C1")
	# 	self.assertEqual(self.mol_chain.atoms[-1].name, "H15")
	# 	self.assertEqual(len(self.mol_chain.residues), 1)
	# 	self.assertEqual(self.mol_chain.residues[-1].name, "UNK")
	# 	self.assertEqual(self.mol_chain.residues[0].name, "UNK")
	# 	self.assertEqual(self.mol_chain.structure.name, "crystal_ligand")
		

	# def test_mol_residue(self):
	# 	self.assertEqual(self.mol_residue.structure.name, "crystal_ligand")
	# 	self.assertEqual(self.mol_residue.chain.name, "")
	# 	self.assertEqual(self.mol_residue.number, "1")
	# 	self.assertEqual(self.mol_residue.name, "UNK")
	# 	self.assertEqual(len(self.mol_residue.atoms), 44)

	# def test_atom(self):
	# 	self.assertEqual(self.mol_atom.identifier, "ATOM")
	# 	self.assertEqual(self.mol_atom.name, "C1")
	# 	self.assertEqual(self.mol_atom.residue_name, "UNK")
	# 	self.assertEqual(self.mol_atom.chain_name, "")
	# 	self.assertEqual(self.mol_atom.residue_number, "1")
	# 	self.assertEqual(self.mol_atom.x, "9.2220")
	# 	self.assertEqual(self.mol_atom.y, "11.4610")
	# 	self.assertEqual(self.mol_atom.z, "1.9550")
	# 	self.assertEqual(self.mol_atom.occupancy, "")






if __name__ == "__main__":
	unittest.main()