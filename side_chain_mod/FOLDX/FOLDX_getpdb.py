import sys

def get_pdb():

	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)

	from Bio.PDB.Polypeptide import one_to_three as ott
	from Bio.PDB.PDBIO import PDBIO
	from Bio.PDB.PDBIO import Select

	#pathmmcif = "/Users/tarun/Documents/mmCIF"
	pathmmcif = "/Volumes/BIOINFO/mmCIF"

	pdb = sys.argv[1]
	C = sys.argv[2]	# CHAIN

	fol = pdb[1:3]

	pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
	tar = gzip.open("{}".format(pdbfile),"rb")
	out = open("pdbprocess.cif","wb")
	out.write(tar.read())
	tar.close()
	out.close()

	structure_id = "{}".format(pdb)
	filename = "pdbprocess.cif"
	structure = parser.get_structure(structure_id,filename)

	model = structure[0]

	chain = model["{}".format(C)]
	c1 = chain.get_list()		# LIST ALL THE RESIDUES
		
	io = PDBIO()
	io.set_structure(chain)
	io.save("WT.pdb")

get_pdb()

