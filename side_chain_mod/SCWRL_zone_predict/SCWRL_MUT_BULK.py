import sys

def side_chain_rem():

	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.PDBParser import PDBParser
	parser1 = PDBParser(PERMISSIVE=0,QUIET=True)
	from Bio.PDB.Polypeptide import one_to_three as ott
	from Bio.PDB.PDBIO import PDBIO
	from Bio.PDB.PDBIO import Select

	#pathmmcif = "/Users/tarun/Documents/mmCIF"
	pathmmcif = "/Volumes/BIOINFO/mmCIF"

	pdb = sys.argv[1]
	C = sys.argv[2]	# CHAIN
	pos = int(sys.argv[3])
	wtres = ott(sys.argv[4])
	mutres = ott(sys.argv[5])
	zone = sys.argv[6]	# ZONE

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
		
	class atom_type(Select):
		def accept_atom(self, r1):
			r2 = r1.get_id()
			p1 = r1.get_parent()
			p = p1.get_id()
			if p[0] == " " and "{}".format(p[1]) in zone:
				if r2 == "CA" or r2 == "N" or r2 == "C" or r2 == "O":
					return 1
				else:
					return 0
			elif p[0] == " ":
				return 1
			elif p[0] != " ":
				return 0

	io = PDBIO()
	io.set_structure(chain)
	io.save("WT.pdb", atom_type())

	print(wtres,mutres,pos)

	# MUTATING THE RESIDUE

	f = open("WT.pdb","r")
	ft = f.readlines()
	f.close()

	g = open("MUT.pdb","w")

	k = 0
	list1 = []
	while k < len(ft):
		ft1 = ft[k].split()
		k1 = 0
		while k1 < len(ft1):
			if len(ft1) > 2:
				if wtres in ft1[k1]:
					if "{}".format(pos) in ft1[(k1+2)]:
						list1.append(k)
			k1 = k1 + 1
		k = k + 1
	
	k = 0
	while k < len(ft):
		t1 = ft[k]
		if k in list1:
			t1 = ft[k].replace("{}".format(wtres),"{}".format(mutres))
		g.write("{}".format(t1))

		k = k + 1

	g.close()
		
side_chain_rem()

