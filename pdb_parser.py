# this script will be used to parse the pdb, 
# here PDB files will be read and documented
import data

class Current_open():
    structures = []

    def add_struct(self, struct):
        self.structures.append(structatom)



class Structure:
    
    def __init__(self, name):
        self.name = name
        self.Chains = []
        self.Residues = []
        self.Atoms = []
        self.chainnames = []

        

    def add_chain(self, chain):
        self.Chains.append(chain)
    
    def add_residue(self, residue):
        self.Residues.append(residue)
    
    def add_atom(self, atom):
        self.Atoms.append(atom)

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
        self.xcor = xcor
        self.ycor = ycor
        self.zcor = zcor
        

    def set_tempf(self, tempf):
        self.tempf = tempf
        
    def set_segid(self,segid):
        self.segid = segid

    def set_occupancy(self, occupancy):
        self.occupancy = occupancy

    def set_struct(self, structure):
        self.Structure = structure
        

def open_pdb(pdb, name):

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
    name = file[-8:-4]
    
    with open(file) as f:

        structure = open_pdb(f, name)

    return structure


def open_web_pdb(pdb):
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
    

        
def write_pdb(structure, filename):

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
                
                col7 = atom.xcor
                for i in range((11 - len(col7))):
                    col7 = " " + col7

                col8 = atom.ycor
                for i in range(8 - len(col8)):
                    col8 = " " + col8
                
                col9 = atom.zcor
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

    import numpy as np
    for i, atom in enumerate(atoms):
        if i == 0:
            array = np.array([[atom, atom.xcor, atom.ycor, atom.zcor]])
        else:
            addarray = np.array([[atom, atom.xcor, atom.ycor, atom.zcor]])
            array = np.concatenate((array, addarray))

    return(array)

struct = open_local_pdb("test/3J95.pdb")

make_atom_array(struct.Atoms)


