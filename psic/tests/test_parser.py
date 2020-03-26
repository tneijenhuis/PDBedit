import unittest
import sys
sys.path.append("..")
import psic

class TestParser(unittest.TestCase):
	def setUp(self):
		self.structure = psic.Parser("../data/test.pdb").parse()
		self.chain = self.structure.chains[0]
		self.residue = self.structure.residues[0]
		self.atom = self.structure.atoms[0]

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


if __name__ == "__main__":
	unittest.main()