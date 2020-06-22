import numpy as np

class Structure:
	
	def __init__(self, name):
		self.name = name
		self.chains =  np.empty([0])
		self.residues =  np.empty([0])
		self.atoms = np.empty([0])

class Chain:

	def __init__(self, name, structure):
		self.residues = np.empty([0])
		self.atoms = np.empty([0])
		self.name = name
		self.structure = structure
	
class Residue:

	def __init__(self, name, structure, number=0, chain=None):
		self.atoms = np.empty([0])
		self.name = name
		self.number = number
		self.structure = structure
		self.chain = chain


class Atom:
	"""Contains all information of each atom from the PDB file"""

	def  __init__(self, identifier, name, residue_name, chain_name, residue_number, x, y, z, 
			occupancy, temperature_factor, segment_id, element, structure = None, chain = None, 
			residue = None):

		self.identifier = identifier
		self.name = name
		self.residue_name = residue_name
		self.chain_name = chain_name
		self.residue_number = residue_number
		self.x = x
		self.y = y
		self.z = z
		self.occupancy = occupancy
		self.temperature_factor = temperature_factor
		self.segment_id = segment_id
		self.element = element
		self.structure = structure
		self.chain = chain
		self.residue = residue

class Voxel:

	"""Represents a cubic box in space
	----------------------------------
	arguments:
	size = the length of the sides
	x, y, z coordinates
	"""

	def __init__(self, size, x, y, z):
		self.content = []
		self.empty = True
		self.size = size
		self.x = x
		self.y = y
		self.z = z

	def add_content(self, content):
		self.content.append(content)
		self.empty = False

class Pocket:
	def __init__(self):
		self.atoms = np.empty([0])


class Parser:

	def parse(self, file, local=True, identifier="ATOM", residue_name="UNK", chain_name=""):

		name = file.split("/")[-1].split(".")[0]

		if file[-4:] == ".pdb":
			if local:
				with open(file) as pdb:
					structure = self._read_pdb(pdb, name)
			else:
				import urllib3
				
				http = urllib3.PoolManager()
				page = http.request('GET', "https://files.rcsb.org/view/{}".format(file))
				f = page.data.decode("utf-8")
				pdb = f.split("\n")
				structure = self._read_pdb(pdb, name)

		elif file[-5:] == ".mol2":
			with open(file) as mol2:
				structure = self._read_mol2(mol2, name, identifier, residue_name, chain_name)

		else:
			raise ValueError("File extention not recognized by parser")

		return structure

	def _read_pdb(self, pdb, name):
		"""Function takes a PDB formatted file and returns a Structure object which contains all atom information of the file"""

		current_structure = Structure(name)
		self.current_chain = ""
		self.current_residue = ""
		self.residue_numbers = {}
		self.chain_names = []

		for line in pdb:
		   
			identifier = line[0:6].strip() #identifier 

			if identifier == "ATOM" or identifier == "HETATM": 

		 
				number = line[6:12].strip()     #atomnmbr
				name = line[12:17].strip()      #atom
				residue_name = line[17:20].strip()   #residue
				chain_name = line[21].strip()        #chain
				residue_number = line[22:26].strip()   #resnmbr
				col7 = line[26].strip()         #?
				x = float(line[30:38].strip())         #X
				y = float(line[38:46].strip())         #Y      
				z = float(line[46:54].strip())         #z
				occupancy = line[54:60].strip()     #occupancy
				temperature_factor = line[60:66].strip()     #Temperature factor
				segid = line[72:76].strip()     #Segment identifer
				element = line[76:78].strip()   #element symbol

			   
				atom = Atom(identifier, name, residue_name, chain_name, residue_number, x, y, z, occupancy,
				temperature_factor, segid, element)

				self._add_atom(atom, current_structure)

			   
		return current_structure

	def _read_mol2(self, mol2, name, identifier, residue_name, chain_name):
		
		current_structure = Structure(name)
		self.current_chain = ""
		self.current_residue = ""
		self.residue_numbers = {}
		self.chain_names = []

		reading = False
		
		for line in mol2:

			if line.strip() == "@<TRIPOS>ATOM":
				reading = True
			elif len(line) > 0:  
				if line.strip().startswith("@"):
					reading = False

			if reading and line.strip() != "@<TRIPOS>ATOM":
				line = line.split()
				identifier = identifier
				number = line[0]     #atomnmbr
				name = line[1]     #atom
				residue_name = residue_name   #residue
				chain_name = chain_name        #chain
				residue_number = line[6] #resnmbr
				x = line[2]         #X
				y = line[3]         #Y      
				z = line[4]         #z
				occupancy = ""     #occupancy
				temperature_factor = ""     #Temperature factor
				segid = chain_name     #Segment identifer
				element = line[5].split(".")[0]   #element symbol

				atom = Atom(identifier, name, residue_name, chain_name, residue_number, x, y, z, occupancy,
				temperature_factor, segid, element)

				self._add_atom(atom, current_structure)


		return current_structure
		


	def _add_atom(self, atom, structure):

		addarray = np.array([atom])
		structure.atoms = np.concatenate((structure.atoms, addarray))
		atom.structure = structure

		# Find the chain of the atom
		if atom.chain_name not in self.chain_names:
			self.current_chain = self._make_new_chain(structure, atom)
			
		elif atom.chain_name != self.current_chain.name:

			self.current_chain = structure.chains[self.chain_names.index(atom.chain_name)]

		atom.chain = self.current_chain

		self._add_atom_to_chain(atom, self.current_chain)

		# Find the residue of the atom
		
		if atom.residue_number not in self.residue_numbers[self.current_chain.name]:
			self.current_residue = self._make_new_residue(structure, atom)
		
		atom.residue = self.current_residue
		self._add_atom_to_residue(atom, self.current_residue)

	def _make_new_chain(self, structure, atom):
		self.current_chain = Chain(atom.chain_name, atom.structure)
		self._add_chain_to_structure(self.current_chain, structure)
		self.chain_names.append(atom.chain_name)
		self.residue_numbers[self.current_chain.name] = []    
		return self.current_chain

	def _add_atom_to_chain(self, atom, chain):
		addarray = np.array([atom])
		chain.atoms = np.concatenate([chain.atoms, addarray])

	def _add_chain_to_structure(self, chain, structure):
		addarray = np.array([chain])
		structure.chains = np.concatenate((structure.chains, addarray))

	def _make_new_residue(self, structure, atom):
		self.current_residue = Residue(atom.residue_name, atom.structure, atom.residue_number, atom.chain)
		self._add_residue_to_structure(self.current_residue, structure)
		self.residue_numbers[atom.chain_name].append(atom.residue_number)
		self._add_residue_to_chain(self.current_residue, self.current_chain)
		return self.current_residue
	
	def _add_residue_to_structure(self, residue, structure):
		addarray = np.array([residue])
		structure.residues = np.concatenate((structure.residues, addarray))

	def _add_atom_to_residue(self, atom, residue):
		addarray = np.array([atom])
		residue.atoms = np.concatenate([residue.atoms, addarray])

	def _add_residue_to_chain(self, residue, chain):
		addarray = np.array([residue])
		chain.residues = np.concatenate([chain.residues, addarray])
	


class Builder(Parser):

	def build_dummy_atom(self, x, y, z, chain_name = "A"):
		identifier = "ATOM"
		name = "X"
		residue_name = "DUM"
		chain_name = chain_name
		residue_number = 0 
		occupancy = ""
		segid = ""
		temperature_factor = ""
		element = "X"

		atom = Atom(identifier, name, residue_name, chain_name, residue_number, x, y, z, occupancy,
				temperature_factor, segid, element)

		return atom

		
	def add_dummy_to_pocket(self, dummy, pocket):
		addarray = np.array([dummy])
		pocket.atoms = np.concatenate((pocket.atoms, addarray))

	def plaster_structure(self):
		pass

		
def make_dummy(x, y, z):

	"""Creates a dummy atom as a Atom object"""

	atom = Atom("DUM", "X")
	atom.set_ident("ATOM")
	atom.set_residue("DUM", 0)
	atom.set_chain(" ")
	atom.set_pos(x, y, z)
	atom.set_tempf("")
	atom.set_occupancy("")
	atom.set_segid("")

	return atom
		
class Writer():

	def write_pdb(structure, filename):
		pass

def write_pdb(structure, filename):

	"""
	takes a Structure object to generate a PDB formatted file using the Atoms information
	-----------------------------------
	arguments:
	structure = Structure object
	filename = PATH + name of the to be generated file
	"""

	import data

	if ".pdb" not in filename:
		print("can only make files witb pdb extention")
		
	else:
		with open(filename, 'w') as f:
			atoms = structure.atoms
			atom_nmbr = 0
			for atom in atoms:
				
				if atom_nmbr > 0:
					if col4.strip() in data.protein_aa and atom.residue_name in data.protein_aa and col5.strip() != atom.chain:
						
						atom_nmbr += 1
						
						col2 = str(atom_nmbr)
						for i in range((5 - len(col2))):
							col2 = " " + col2
						
						f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))
					elif col4.strip() in data.protein_aa and atom.residue_name not in data.protein_aa:

						atom_nmbr += 1
						
						col2 = str(atom_nmbr)
						for i in range((5 - len(col2))):
							col2 = " " + col2
						
						f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))

				atom_nmbr +=1
				
				col1 = atom.identifier
				for i in range((6 - len(col1))):
					col1 = col1 + " "
				
				col2 = str(atom_nmbr)
				for i in range((5 - len(col2))):
					col2 = " " + col2
				
				col3 = atom.name   
				for i in range((5 - len(col3))):
					if i < 2:
						col3 = " " + col3
					else:
						col3 = col3 + " "
				
				col4 = atom.residue_name
				for i in range((5 - len(col4))):
					if i < 1:
						col4 = " " + col4
					else:
						col4 = col4 + " "
				
				col5 = atom.chain_name

				col6 = str(atom.residue_name)
				for i in range((4 - len(col6))):
					col6 = " " + col6
				
				col7 = str(atom.x)
				if len(col7) > 11:
					col7 = col7[:7]
				for i in range((11 - len(col7))):
					col7 = " " + col7

				col8 = str(atom.y)
				if len(col8) > 8:
					col8 = col8[:7]
				for i in range(8 - len(col8)):
					col8 = " " + col8
				
				col9 = str(atom.z)
				if len(col9) > 8:
					col9 = col9[:7]
				for i in range(8 - len(col9)):
					col9 = " " + col9 
				
				col10 = atom.occupancy
				for i in range(6-len(col10)):
					col10 = " " + col10
				
				col11 = atom.temperature_factor
				for i in range(6-len(col11)):
					col11 = " " + col11
				
				col12 = atom.segment_id
				for i in range(4-len(atom.segment_id)):
					col12 = " " + col12

				col13 = atom.element
				for i in range(2-len(col13)):
					col13 = " " + col13

				f.write("{}{}{}{}{}{} {}{}{}{}{}      {}{}\n".format(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13))
			f.write("END")


def make_atom_array(atoms):
	""" Takes a list of atom objects and returns a list of its positions
	"""
	if len(atoms) == 0:
		print("An empty list cannot be converted")

	else:
		for i, atom in enumerate(atoms):
			if i == 0:
				array = np.array([[atom, atom.x, atom.y, atom.z]])
			else:
				addarray = np.array([[atom, atom.x, atom.y, atom.z]])
				array = np.concatenate((array, addarray))

		return array


def make_3d_grid(array, voxel_size = 5):

	"""
	Takes an array to generate a 3d array containing Voxel object [z[y[[x[Voxel]]]]
	--------------------------------
	arguments:
	array = a 2d array containing [object, x, y, z]

	optional:
	voxel_size (default 5) = the cubic size of each voxel
	"""

	# print("Separating atoms into voxels of size {} A^3".format(voxel_size))
	
	axis1 = 0
	axis2 = 1
	axis3 = 2

	maxs = np.amax(array[:,1:], axis = 0)

	mins = np.amin(array[:,1:], axis=0)

	centrs = (mins + maxs) / 2

	axis1_start = centrs[axis1]
	axis2_start = centrs[axis2]
	axis3_start = centrs[axis3]
	axis1_count = 0
	axis2_count = 0
	axis3_count = 0 

	# initializing search space
	while axis1_start > mins[axis1]:
		axis1_start -= voxel_size

	while axis2_start > mins[axis2]:
		axis2_start -= voxel_size

	while axis3_start > mins[axis3]:
		axis3_start -= voxel_size

	# array dimention calculation
	axis3_curr = axis3_start
	while axis3_curr < maxs[axis3]:
		axis3_count += 1
		axis3_curr += voxel_size
	axis2_curr = axis2_start
	while axis2_curr < maxs[axis2]:
		axis2_count += 1
		axis2_curr += voxel_size
	axis1_curr = axis1_start
	while axis1_curr < maxs[axis1]:
		axis1_count += 1
		axis1_curr += voxel_size
	   
	# generatin the 3d grid

	axis3_curr = axis3_start

	axis3_array = np.empty([0, axis2_count, axis1_count])

	while axis3_curr < maxs[axis3]:

		allanal = 0

		axis3_slice = np.empty([0, 4])
		for atom in array:
			if atom[axis3+1] >= axis3_curr and atom[axis3+1] < (axis3_curr + voxel_size):
				addarray = np.array([atom])
				axis3_slice = np.concatenate((axis3_slice, addarray))
				allanal +=1


		axis2_array = np.empty([0, axis1_count])
		axis2_curr = axis2_start

		while axis2_curr < maxs[axis2]:

			axis2_line = np.empty([0, 4])
			for atom in axis3_slice:
				if atom[axis2+1] >= axis2_curr and atom[axis2+1] < axis2_curr + voxel_size:
					addarray = np.array([atom])
					axis2_line = np.concatenate((axis2_line, addarray))            

			axis1_array = np.empty([0])
			axis1_curr = axis1_start

			while axis1_curr < maxs[axis1]:

				x = (axis1_curr * 2 + voxel_size) / 2
				y = (axis2_curr * 2 + voxel_size) / 2
				z = (axis3_curr * 2 + voxel_size) / 2

				voxel = Voxel(voxel_size, x, y, z)

				for atom in axis2_line:
					if atom[axis1 +1] >= axis1_curr and atom[axis1+1] < axis1_curr + voxel_size:
						voxel.add_content(atom[0])

				#print(voxel_array.shape)
				add_array = np.array([voxel])
				axis1_array = np.concatenate((axis1_array, add_array))

				axis1_curr += voxel_size

			#print(axis1_array.shape)
			add_array = np.array([axis1_array])
			axis2_array = np.concatenate((axis2_array, add_array))


			axis2_curr += voxel_size

		add_array = np.array([axis2_array])
		axis3_array = np.concatenate((axis3_array, add_array))


		axis3_curr += voxel_size
	return axis3_array


def surf_from_grid(gird):
	"""
	Calculates surface residues
	"""

	surfs = np.empty([0])

	# this fuction will check the content of a specific voxel and compairs it with its neighbours.
	# if a neibour is empty the voxel is positioned at the surface

	print("Calculating surface")    

	for current_slice_numb in range(grid.shape[0]):
		print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
		for current_line_numb in range(grid.shape[1]):
			for current_voxel_numb in range(grid.shape[2]):
				current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

				if (current_slice_numb == 0 or current_slice_numb == grid.shape[0] - 1 or current_line_numb == 0 or 
					current_line_numb == grid.shape[1] -1 or current_voxel_numb == 0 or current_voxel_numb == grid.shape[2] -1 ):
					if current_voxel.empty == False:
						for atom in current_voxel.content:   
							if atom.Residue not in surfs:

								addarray = np.array([atom.Residue])
								surfs = np.concatenate((surfs, addarray)) 

				
				else:
					neighbour1 = grid[current_slice_numb - 1, current_line_numb, current_voxel_numb]
					neighbour2 = grid[current_slice_numb + 1, current_line_numb, current_voxel_numb]
					neighbour3 = grid[current_slice_numb, current_line_numb - 1, current_voxel_numb]
					neighbour4 = grid[current_slice_numb, current_line_numb + 1, current_voxel_numb]
					neighbour5 = grid[current_slice_numb, current_line_numb, current_voxel_numb - 1]
					neighbour6 = grid[current_slice_numb, current_line_numb, current_voxel_numb + 1]  

					if (neighbour1.empty or neighbour2.empty or neighbour3.empty or neighbour4.empty or neighbour5.empty or neighbour6.
					empty) and current_voxel.empty == False: 
						for atom in current_voxel.content:   
							if atom.Residue not in surfs:

								addarray = np.array([atom.Residue])
								surfs = np.concatenate((surfs, addarray)) 

	return surfs

def find_interface(grid):
	"""
	Finds the interface between proteins
	"""
	interface_voxels = np.empty([0])
	print("searching for interface")


	for current_slice_numb in range(grid.shape[0]):

		print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
		for current_line_numb in range(grid.shape[1]):
			for current_voxel_numb in range(grid.shape[2]):
				current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

				chain_names = []
				for atom in current_voxel.content:
					if atom.chain not in chain_names:
						chain_names.append(atom.chain)

				if len(chain_names) > 1:
					addarray = np.array([current_voxel])
					interface_voxels = np.concatenate((interface_voxels, addarray))

	return interface_voxels


######################################################
# Everything under this line is under development 

def fill_void(grid, neg_array, max_distance):
	"""
	finds pockets in a protein structure
	"""
	import math

	empty_voxels = np.empty([0])
	# print("Searching for pockets")


	for current_slice_numb in range(grid.shape[0]):

		# print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
		for current_line_numb in range(grid.shape[1]):
			for current_voxel_numb in range(grid.shape[2]):
				current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

				
				if (current_slice_numb == 0 or current_slice_numb == grid.shape[0] - 1 or 
						current_line_numb == 0 or current_line_numb == grid.shape[1] - 1 or
						current_voxel_numb == 0 or current_voxel_numb == grid.shape[2] - 1):
					pass

				else: 

					neighbour1 = grid[current_slice_numb - 1, current_line_numb, current_voxel_numb]
					neighbour2 = grid[current_slice_numb + 1, current_line_numb, current_voxel_numb]
					neighbour3 = grid[current_slice_numb, current_line_numb - 1, current_voxel_numb]
					neighbour4 = grid[current_slice_numb, current_line_numb + 1, current_voxel_numb]
					neighbour5 = grid[current_slice_numb, current_line_numb, current_voxel_numb - 1]
					neighbour6 = grid[current_slice_numb, current_line_numb, current_voxel_numb + 1] 


					if current_voxel.empty and (neighbour1.empty and neighbour2.empty and
							neighbour3.empty and neighbour4.empty and neighbour5.empty and
							neighbour6.empty): 

						in_struct = True

						for neg in neg_array:
							distance = math.sqrt(math.pow(neg[1] - current_voxel.x, 2) + math.pow(neg[2] - current_voxel.y, 2) + math.pow(neg[3] - current_voxel.z, 2))
							if distance <= max_distance:
								in_struct = False

						if in_struct:

							addarray = np.array([current_voxel])
							empty_voxels = np.concatenate((empty_voxels, addarray))
 
	return empty_voxels

def plaster_structure(grid):
	"""
	plasters the extirior of the protein structure with dummy atoms
	"""

	empty_voxels = np.empty([0])
	#print("plastering")
  
	for current_slice_numb in range(grid.shape[0]):
		# print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
		for current_line_numb in range(grid.shape[1]):
			for current_voxel_numb in range(grid.shape[2]):
				current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

				if current_voxel.empty:
					
					if current_slice_numb == 0:
						neighbour1 = Voxel(0, 0, 0, 0)
					else:                    
						neighbour1 = grid[current_slice_numb - 1, current_line_numb, current_voxel_numb]
					
					if current_slice_numb == grid.shape[0] - 1:
						neighbour2 = Voxel(0, 0, 0, 0)
					else:
						neighbour2 = grid[current_slice_numb + 1, current_line_numb, current_voxel_numb]
					
					if current_line_numb == 0:
						neighbour3 = Voxel(0, 0, 0, 0)
					else:
						neighbour3 = grid[current_slice_numb, current_line_numb - 1, current_voxel_numb]
					
					if current_line_numb == grid.shape[1] - 1:
						neighbour4 = Voxel(0, 0, 0, 0)
					else:
						neighbour4 = grid[current_slice_numb, current_line_numb + 1, current_voxel_numb]
					
					if current_voxel_numb == 0:
						neighbour5 = Voxel(0, 0, 0, 0)
					else:
						neighbour5 = grid[current_slice_numb, current_line_numb, current_voxel_numb - 1]
					
					if current_voxel_numb == grid.shape[2] - 1:
						neighbour6 = Voxel(0, 0, 0, 0)
					else:
						neighbour6 = grid[current_slice_numb, current_line_numb, current_voxel_numb + 1] 
					
					neighbours = [neighbour1, neighbour2, neighbour3, neighbour4, neighbour5, neighbour6]

					empty_neighbours = 0
					for neighbour in neighbours:
						if neighbour.empty:
							empty_neighbours += 1
					
					if empty_neighbours > 4:
					
						addarray = np.array([current_voxel])
						empty_voxels = np.concatenate((empty_voxels, addarray))
 
	plaster = Pocket()
	
	for voxel in empty_voxels:
		dum = Builder().build_dummy_atom(voxel.x, voxel.y, voxel.z)
		Builder().add_dummy_to_pocket(dum, plaster)
	
	return plaster

def cluster_voxels(voxels):
	import math
	#print("Clustering pockets")
	# I use a normal list for pockets as this does not require similar axis lengths
	pockets = []
	neighbours = np.empty([0])
	while len(voxels) > 0:
		
		if len(neighbours) == 0:
			current_dummy = voxels[0]
			current_pocket = np.array([current_dummy])
			voxels = np.delete(voxels, np.where(voxels == current_dummy))

			for dummy in voxels:
				distance = math.sqrt(math.pow(dummy.x - current_dummy.x, 2) + math.pow(dummy.y - current_dummy.y, 2) + math.pow(dummy.z - current_dummy.z, 2))
				if distance == dummy.size*2:
					addarray = np.array([dummy])
					current_pocket = np.concatenate((current_pocket, addarray))
					neighbours = np.concatenate((neighbours, addarray))
					voxels = np.delete(voxels, np.where(voxels == dummy))

		else:
			for neighbour in neighbours:
				current_dummy = neighbour
				for dummy in voxels:
					distance = math.sqrt(math.pow(dummy.x - current_dummy.x, 2) + math.pow(dummy.y - current_dummy.y, 2) + math.pow(dummy.z - current_dummy.z, 2))
					if distance == dummy.size:
						addarray = np.array([dummy])
						current_pocket = np.concatenate((current_pocket, addarray))
						neighbours = np.concatenate((neighbours, addarray))
						voxels = np.delete(voxels, np.where(voxels == dummy))

				neighbours = np.delete(neighbours, np.where(neighbours == current_dummy))
		
		if len(neighbours) == 0:
			
			if len(current_pocket) > 5:
				pockets.append(current_pocket)

			if len(voxels) == 1:
				voxels = []
			else:
				voxels = np.delete(voxels, np.where(voxels == current_pocket))

	if len(current_pocket) > 5:
		pockets.append(current_pocket)
	
	pockets = sorted(pockets, key = lambda x : len(x), reverse = True )
	return pockets


def find_pockets(structure):
	"""This function calls all functions for the pocket detection"""

	array = make_atom_array(structure.atoms)
	first_grid = make_3d_grid(array, 4)
	plaster = plaster_structure(first_grid)

	plaster_array = make_atom_array(plaster.atoms)
	grid = make_3d_grid(array, 2)

	pocket_potential = fill_void(grid, plaster_array, 4)

	pockets = cluster_voxels(pocket_potential)

	return pockets[:10]

def compare_pocket_to_ligand(pocket, ligand):

	import math

	distances = []
	for atom in ligand:
		distance = 100000
		for dum in pocket:
			current_distance = math.sqrt(math.pow(dum.x - atom.x, 2) + math.pow(dum.y - atom.y, 2) + math.pow(dum.z - atom.z, 2))
			if distance > current_distance:
				distance = current_distance
		distances.append(distance)
	mean = (sum(distances) / len(distances))
	return mean

################################################################3##
	
# struct= Parser().parse("data/receptor.pdb")


# pockets = find_pockets(struct)


# best_pocket = []
# for pocket in pockets:
#     if len(best_pocket) < len(pocket):
#         best_pocket = pocket

# best_pocket_struct = Pocket()
# for voxel in best_pocket:
#     dum = Builder().build_dummy_atom(voxel.x, voxel.y, voxel.z)
#     Builder().add_dummy_to_pocket(dum, best_pocket_struct)


# all_pocket_struct = Structure("all")
# for pocket in pockets:
#     for voxel in pocket:
# 	    dum = Builder().build_dummy_atom(voxel.x, voxel.y, voxel.z)
# 	    Builder().add_dummy_to_pocket(dum, all_pocket_struct)


# write_pdb(all_pocket_struct, "pockets.pdb")
# write_pdb(best_pocket_struct, "best.pdb")

# ligand = open_mol2("data/crystal_ligand.mol2", "ligand")

# write_pdb(ligand, "ligand.pdb")

#################
#corona#
#################

# struct = open_web_pdb("6Y84")

# no_hoh = Structure("no_hoh")
# for atom in struct.Atoms:
#     if atom.residue_name != "HOH":
#         no_hoh.add_atom(atom)


# pockets = find_pockets(no_hoh)

# all_pocket_struct = Structure("all")
# for pocket in pockets:
#     for voxel in pocket:
#         dum = make_dummy(voxel.x, voxel.y, voxel.z)
#         all_pocket_struct.add_atom(dum)

# write_pdb(all_pocket_struct, "pockets.pdb")
# write_pdb(struct, "6Y84.pdb")
