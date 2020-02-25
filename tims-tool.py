import pdb_parser as pp 
import numpy as np

struct = pp.open_local_pdb('test/3J95.pdb')

chains = struct.chains

for i in range(len(chains)):
    for residue in chains[i].residues:
        for atom in residue.atoms:
            refat = np.array([float(atom.xcor), float(atom.ycor), float(atom.zcor)])
            for j in range(i , len(chains)):
                for altres in chains[j].residues:
                    for altatm in residue.atoms:
                        alta = np.array([float(altatm.xcor), float(altatm.ycor), float(altatm.zcor)])


