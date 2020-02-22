# this script will be used to parse the pdb, 
# here PDB files will be read and documented

class Structure():
    pass

class Chain():
    pass

class Redidue():
    pass

class Atom(number):
    
    Atom.number = number

def open_pdb(file):
    with open(file) as f:
        for line in f:
            identifier = line[0:6].strip() #identifier 
            if identifier == "ATOM" or "HETATM":      
                col2 = line[6:12].strip()       #atomnmbr

                atom = Atom(col2)
                col3 = line[12:17].strip()      #atom
                col4 = line[17:20].strip()      #residue
                col5 = line[21].strip()         #chain
                col6 = line[22:26].strip()      #resnmbr
                col7 = line[26].strip()         #?
                col8 = line[30:38].strip()      #X
                col9 = line[38:46].strip()      #Y      
                col10 = line[46:54].strip()     #z
                col11 = line[54:60].strip()     #occupancy
                col12 = line[60:66].strip()     #Temperature factor
                col13 = line[72:76].strip()     #Segment identifer
                col14 = line[76:78].strip()     #element symbol
                print(line)


pdb_path = open_pdb("test/part.pdb")