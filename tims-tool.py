import pdb_parser as pp 

struct = pp.open_local_pdb('test/3J95.pdb')

chains = struct.chains

chain = chains[-2]

for residue in chain.residues:
	print(residue.name)