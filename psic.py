import data
import numpy as np

class Current_open():
    structures = []

    def add_struct(self, struct):
        self.structures.append(structatom)

class Structure:
    
    def __init__(self, name):
        self.name = name
        self.Chains =  np.empty([0])
        self.Residues =  np.empty([0])
        self.Atoms = np.empty([0])
        self.chainnames =  []

    def add_chain(self, chain):
        addarray = np.array([chain])
        self.Chains = np.concatenate((self.Chains, addarray))
    
    def add_residue(self, residue):
        addarray = np.array([residue])
        self.Residues = np.concatenate((self.Residues, addarray))
    
    def add_atom(self, atom):
        addarray = np.array([atom])
        self.Atoms = np.concatenate((self.Atoms, addarray))
        atom.set_struct(self)

    def add_chainname(self, chainname):
        self.chainnames.append(chainname)    

   
class Chain:

    def __init__(self):
        self.Residues = []
        self.Atoms = []
        self.name = ""
    
    def add_atom(self, atom):
        self.Atoms.append(atom)

    def add_residue(self, residue):
        self.Residues.append(residue)

    def set_name(self, name):
        self.name = name
        
        return self.name

    
class Residue:

    def __init__(self):
        self.Atoms = []
        self.name = ""
        self.nmbr = 0

    def add_atom(self, atom):
        self.Atoms.append(atom)

    def set_nmbr(self, nmbr):
        self.nmbr = nmbr

        return self.nmbr

    def set_name(self, name):
        self.name = name

        return self.name



class Atom:

    def  __init__(self, name, element,):
        self.name = name
        self.element = element

    def set_ident(self, ident):
        self.ident = ident

    def set_chain(self, chain):
        self.chain = chain
        
    def set_residue(self, resname, resnmbr):
        self.resname = resname
        self.resnmbr = resnmbr

    def set_pos(self, xcor, ycor, zcor):
        self.xcor = float(xcor)
        self.ycor = float(ycor)
        self.zcor = float(zcor)
        

    def set_tempf(self, tempf):
        self.tempf = tempf
        
    def set_segid(self,segid):
        self.segid = segid

    def set_occupancy(self, occupancy):
        self.occupancy = occupancy

    def set_struct(self, structure):
        self.Structure = structure

    def set_Residue(self, residue):
        self.Residue = residue
        

class Voxel():

    def __init__(self, size, xcor, ycor, zcor):
        self.content = []
        self.empty = True
        self.size = size
        self.xcor = xcor
        self.ycor = ycor
        self.zcor = zcor

    def add_content(self, content):
        self.content.append(content)
        self.empty = False

def open_pdb(pdb, name):
    """
    Function takes a PDB formatted file and returns a Structure object which contains all atom information of the file
    ---------------
    arguments :
    pdb = iteratable PDB file
    name = name of the generated structure
    """

    build_res = Residue()
    build_chain = Chain()
    build_structure = Structure(name)

    for line in pdb:
       
        identifier = line[0:6].strip() #identifier 
        
        if identifier == "ATOM" or identifier == "HETATM":      
       
            atomnmbr = line[6:12].strip()       #atomnmbr
            atomname = line[12:17].strip()      #atom
            resname = line[17:20].strip()      #residue
            chain = line[21].strip()         #chain
            resnmbr = line[22:26].strip()      #resnmbr
            col7 = line[26].strip()         #?
            xcor = line[30:38].strip()      #X
            ycor = line[38:46].strip()      #Y      
            zcor = line[46:54].strip()     #z
            occup = line[54:60].strip()     #occupancy
            tempf = line[60:66].strip()     #Temperature factor
            segid = line[72:76].strip()     #Segment identifer
            element = line[76:78].strip()     #element symbol
            
            atom = Atom(atomname, element)
            atom.set_struct(build_structure)
            atom.set_Residue(build_res)
            atom.set_ident(identifier)
            atom.set_residue(resname, resnmbr)
            atom.set_chain(chain)
            atom.set_pos(xcor, ycor, zcor)
            atom.set_tempf(tempf)
            atom.set_occupancy(occup)
            atom.set_segid(segid)


            #building residues from atoms
            if len(build_res.Atoms) == 0  or build_res.Atoms[0].resnmbr == atom.resnmbr and build_chain.Atoms[0].chain == atom.chain:
                build_res.add_atom(atom)
                build_chain.add_atom(atom)
                build_structure.add_atom(atom)
                
            #building chains from residues
            elif len(build_chain.Residues) == 0 or build_chain.Atoms[0].chain == atom.chain:
                build_res.set_nmbr(build_res.Atoms[0].resnmbr)
                build_res.set_name(build_res.Atoms[0].resname)
                build_chain.add_residue(build_res)
                build_structure.add_residue(build_res)
                build_res = Residue()
                build_res.add_atom(atom)
                build_chain.add_atom(atom)
                build_structure.add_atom(atom)

            elif atom.chain in build_structure.chainnames:
                build_res.set_nmbr(build_res.Atoms[0].resnmbr)
                build_res.set_name(build_res.Atoms[0].resname)
                build_chain.add_residue(build_res)
                build_chain.set_name(build_res.Atoms[0].chain)
                
                if build_chain.name not in build_structure.chainnames:
                    build_structure.add_chain(build_chain)
                    build_structure.add_chainname(build_chain.name)
                
                build_chain = build_structure.Chains[build_structure.chainnames.index(atom.chain)]
                build_res = Residue()
                build_res.add_atom(atom)
                build_chain.add_atom(atom)
                build_structure.add_atom(atom)

            else:
                build_res.set_nmbr(build_res.Atoms[0].resnmbr)
                build_res.set_name(build_res.Atoms[0].resname)
                build_chain.add_residue(build_res)
                build_chain.set_name(build_res.Atoms[0].chain)
                #build_structure.add_chain(build_chain)
                
                if build_chain.name not in build_structure.chainnames:
                    build_structure.add_chain(build_chain)
                    build_structure.add_chainname(build_chain.name)

                build_chain = Chain()
                build_res = Residue()
                build_res.add_atom(atom)
                build_chain.add_atom(atom)
                build_structure.add_atom(atom)


    build_res.set_nmbr(build_res.Atoms[0].resnmbr)
    build_res.set_name(build_res.Atoms[0].resname)
    build_chain.add_residue(build_res)
    build_structure.add_residue(build_res)
    build_chain.set_name(build_res.Atoms[0].chain)

    if atom.chain not in build_structure.chainnames:
        build_structure.add_chain(build_chain)
        build_structure.add_chainname(build_chain.name)

    Current_open.structures.append(build_structure)
    
    return build_structure



def open_local_pdb(file):
    """
    This function opens a local PDB file and generates a Structure object
    ---------
    arguments:
    file = path to be parsed PDB file
    """

    name = file[-8:-4]
    
    with open(file) as f:

        structure = open_pdb(f, name)

    return structure


def open_web_pdb(pdb):
    """
    Parses a PDB structure from the web and returns a Structure object
    ------------------------
    arguments:
    pdb = name of the pdb structure
    """

    import urllib3
    http = urllib3.PoolManager()

    if ".pdb" not in pdb:
        pdb += ".pdb"

    name = pdb[:-4]
    page = http.request('GET', "https://files.rcsb.org/view/{}".format(pdb))
    
    f = page.data.decode("utf-8")
    file =f.split("\n")

    structure = open_pdb(file, name)

    return structure
    
def make_dummy(xcor, ycor, zcor):

    atom = Atom("DUM", "X")
    atom.set_ident("ATOM")
    atom.set_residue("DUM", 0)
    atom.set_chain(" ")
    atom.set_pos(xcor, ycor, zcor)
    atom.set_tempf("")
    atom.set_occupancy("")
    atom.set_segid("")

    return atom
        
def write_pdb(structure, filename):

    """
    takes a Structure object to generate a PDB formatted file using the Atoms information
    -----------------------------------
    arguments:
    structure = Structure object
    filename = PATH + name of the to be generated file
    """

    if ".pdb" not in filename:
        print("can only make files witb pdb extention")
        
    else:
        with open(filename, 'w') as f:
            atoms = structure.Atoms
            atom_nmbr = 0
            for atom in atoms:
                
                if atom_nmbr > 0:
                    if col4.strip() in data.protein_aa and atom.resname in data.protein_aa and col5.strip() != atom.chain:
                        
                        atom_nmbr += 1
                        
                        col2 = str(atom_nmbr)
                        for i in range((5 - len(col2))):
                            col2 = " " + col2
                        
                        f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))
                    elif col4.strip() in data.protein_aa and atom.resname not in data.protein_aa:

                        atom_nmbr += 1
                        
                        col2 = str(atom_nmbr)
                        for i in range((5 - len(col2))):
                            col2 = " " + col2
                        
                        f.write("TER   {}     {}{}{}\n".format(col2, col4, col5, col6))

                atom_nmbr +=1
                
                col1 = atom.ident
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
                
                col4 = atom.resname
                for i in range((5 - len(col4))):
                    if i < 1:
                        col4 = " " + col4
                    else:
                        col4 = col4 + " "
                
                col5 = atom.chain

                col6 = str(atom.resnmbr)
                for i in range((4 - len(col6))):
                    col6 = " " + col6
                
                col7 = str(atom.xcor)
                for i in range((11 - len(col7))):
                    col7 = " " + col7

                col8 = str(atom.ycor)
                for i in range(8 - len(col8)):
                    col8 = " " + col8
                
                col9 = str(atom.zcor)
                for i in range(8 - len(col9)):
                    col9 = " " + col9
                
                col10 = atom.occupancy
                for i in range(6-len(col10)):
                    col10 = " " + col10
                
                col11 = atom.tempf
                for i in range(6-len(col11)):
                    col11 = " " + col11
                
                col12 = atom.segid
                for i in range(4-len(atom.segid)):
                    col12 = " " + col12

                col13 = atom.element
                for i in range(2-len(col13)):
                    col13 = " " + col13

                f.write("{}{}{}{}{}{} {}{}{}{}{}      {}{}\n".format(col1, col2, col3, col4, col5, col6, col7, col8, col9, col10, col11, col12, col13))
            f.write("END")


def make_atom_array(atoms):
    """ Takes a list of atom objects and returns a list of its positions
    """

    for i, atom in enumerate(atoms):
        if i == 0:
            array = np.array([[atom, atom.xcor, atom.ycor, atom.zcor]])
        else:
            addarray = np.array([[atom, atom.xcor, atom.ycor, atom.zcor]])
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

    print("Separating atoms into voxels of size {} A^3".format(voxel_size))
    
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


######################################################
# Everything under this line is under development 

def find_pockets(gird, neg_array):
    """
    finds pockets in a protein structure
    """
    import math

    empty_voxels = np.empty([0])
    print("Searching for pockets")


    for current_slice_numb in range(grid.shape[0]):

        print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
        for current_line_numb in range(grid.shape[1]):
            for current_voxel_numb in range(grid.shape[2]):
                current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

                if (current_slice_numb == 0 or current_slice_numb == grid.shape[0] - 1 or current_line_numb == 0 or 
                    current_line_numb == grid.shape[1] -1 or current_voxel_numb == 0 or current_voxel_numb == grid.shape[2] -1 ):
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
                            distance = math.sqrt(math.pow(neg[1] - current_voxel.xcor, 2) + math.pow(neg[2] - current_voxel.ycor, 2) + math.pow(neg[3] - current_voxel.zcor, 2))
                            if distance <= 5:
                                in_struct = False

                        if in_struct:

                            addarray = np.array([current_voxel])
                            empty_voxels = np.concatenate((empty_voxels, addarray))
 
    return empty_voxels

def plaster_structure(gird):
    """
    finds pockets in a protein structure
    """

    empty_voxels = np.empty([0])
    print("plastering")


    for current_slice_numb in range(grid.shape[0]):

        print("{}%".format((float(current_slice_numb)/float(grid.shape[0])*100)))
        for current_line_numb in range(grid.shape[1]):
            for current_voxel_numb in range(grid.shape[2]):
                current_voxel = grid[current_slice_numb, current_line_numb, current_voxel_numb]

                if (current_slice_numb == 0 or current_slice_numb == grid.shape[0] - 1 or current_line_numb == 0 or current_voxel.empty or
                    current_line_numb == grid.shape[1] -1 or current_voxel_numb == 0 or current_voxel_numb == grid.shape[2] -1 ):

                        addarray = np.array([current_voxel])
                        empty_voxels = np.concatenate((empty_voxels, addarray))
 
    return empty_voxels


struct = open_local_pdb("total.pdb")

alphas = np.empty([0])
for atom in struct.Atoms:
    #print(atom.ident)

    if atom.ident != "HETATM":
        addarray = np.array([atom])
        alphas = np.concatenate((alphas, addarray))



arr = make_atom_array(alphas)
grid = make_3d_grid(arr, 6)


pocket = plaster_structure(grid)
not_poc = Structure("dummy")
for voxel in pocket:

    x = voxel.xcor
    y = voxel.ycor 
    z = voxel.zcor 

    dum = make_dummy(x, y, z)
    not_poc.add_atom(dum)


out_arr = make_atom_array(not_poc.Atoms)
arr = make_atom_array(alphas)
grid = make_3d_grid(arr, 2)

pocket = find_pockets(grid, out_arr)
dummy = Structure("dummy")
for voxel in pocket:

    dum = make_dummy(voxel.xcor, voxel.ycor, voxel.zcor)
    dummy.add_atom(dum)

write_pdb(dummy, "dummy.pdb")
write_pdb(not_poc, "plastered.pdb")
# surfs = surf_from_grid(grid)

# new_struct = Structure("Name")

# for res in surfs:
#     for atom in res.Atoms:
#         new_struct.add_atom(atom)

# print(len(surfs))
# write_pdb(struct, "total.pdb")
# write_pdb(new_struct, "surf.pdb")


