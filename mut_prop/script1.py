def RMSD(var1,var2,var3):

	import math
	
	k = 0
	l = 0.0
	while k < len(var1):
		x1 = var1[k][0]
		y1 = var1[k][1]
		z1 = var1[k][2]

		x2 = var2[k][0]
		y2 = var2[k][1]
		z2 = var2[k][2]

		l = l + (((x1-x2)**2) + ((y1-y2)**2) + ((z1-z2)**2))

		k = k + 1

	l = l / var3
	RMSD = math.sqrt(l)
	RMSD = round(RMSD,2)

	return(RMSD)

def func1():
	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	#from Bio.PDB.PDBParser import PDBParser
	#parser = PDBParser(PERMISSIVE=0,QUIET=True)

	from Bio.PDB.PDBIO import PDBIO

	pathmmcif = "/Users/tarun/Documents/mmCIF"
	#pathmmcif = "/data/pdb/divided/mmCIF"
	#pathmmcif = "/Volumes/BIOINFO/mmCIF"
	#pathmmcif = "/Volumes/RCSB_DATA/pdb"

	#count = 0	
	#if count == 0:
	try:
		pdb1 = "{}".format(sys.argv[2])
		fol = pdb1[1:3]		
		c1 = "{}".format(sys.argv[3])
		pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb1)
		#pdbfile = "{}/{}/pdb{}.ent.gz".format(pathmmcif,fol,pdb1)
		tar = gzip.open("{}".format(pdbfile),"rb")
		out = open("pdbprocess.cif","wb")
		#out = open("pdbprocess.pdb","wb")
		out.write(tar.read())
		tar.close()
		out.close()
		structure_id = "{}".format(pdb1)
		filename = "pdbprocess.cif"
		#filename = "pdbprocess.pdb"
		structure = parser.get_structure(structure_id,filename)
		model = structure[0]
		chain = model["{}".format(c1)]

		io = PDBIO()
		io.set_structure(chain)
		io.save("chain1.pdb")

		pdb2 = "{}".format(sys.argv[4])
		fol = pdb2[1:3]		
		c2 = "{}".format(sys.argv[5])

		pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb2)
		#pdbfile = "{}/{}/pdb{}.ent.gz".format(pathmmcif,fol,pdb2)
		tar = gzip.open("{}".format(pdbfile),"rb")
		out = open("pdbprocess.cif","wb")
		#out = open("pdbprocess.pdb","wb")
		out.write(tar.read())
		tar.close()
		out.close()
		structure_id = "{}".format(pdb2)
		filename = "pdbprocess.cif"
		#filename = "pdbprocess.pdb"
		structure = parser.get_structure(structure_id,filename)
		model = structure[0]
		chain = model["{}".format(c2)]

		io = PDBIO()
		io.set_structure(chain)
		io.save("chain2.pdb")
	except:
		print("FILE NOT FOUND")

def func2(arg1):

	from Bio.PDB.PDBParser import PDBParser
	parser = PDBParser(QUIET=True)

	# CONVERTING RASMOL FILES INTO PDB FORMAT

	f = open("output_all","r")
	ft = f.readlines()
	f.close()
	
	g = open("file1.pdb","w")		# ONLY C-ALPHA

	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		if ft1[0] == "ATOM" or ft1[0].strip("\n") == "TER" or ft1[0] == "REMARK":
			g.write("{}".format(ft[k]))
		k = k + 1
	g.close()

	f = open("output_all_atm","r")
	ft = f.readlines()
	f.close()

	g = open("file2.pdb","w")		# ALL ATOMS

	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		if ft1[0] == "ATOM" or ft1[0].strip("\n") == "TER" or ft1[0] == "REMARK":
			g.write("{}".format(ft[k]))
		k = k + 1
	g.close()

	# READING THE GLOBAL RMSD VALUE AND THE COVERAGE VALUE

	f = open("temp","r")
	ft = f.readlines()
	f.close()

	k = 0
	ft1 = ft[k].split()
	try:
		t1 = ft1[0]
	except:
		t1 = "NA"
	while t1 != "Length":
		k = k + 1
		ft1 = ft[k].split()
		try:
			t1 = ft1[0]
		except:
			t1 = "NA"

	ft1 = ft[k].split()
	lenwt = ft1[3]
	k = k + 1
	ft1 = ft[k].split()
	lenmut = ft1[3]

	ft1 = ft[k].split()
	try:
		t1 = ft1[0]
	except:
		t1 = "NA"
	while t1 != "Aligned":
		k = k + 1
		ft1 = ft[k].split()
		try:
			t1 = ft1[0]
		except:
			t1 = "NA"

	ft1 = ft[k].split(",")
	t1 = ft1[1].split()
	grmsd = t1[1]
	t1 = ft1[0].split()
	coverage = t1[2]

	g = open("results","w")

	g.write("{} {} {} {}\n".format(grmsd,lenwt,lenmut,coverage))

	if arg1 == 1:

		reswt1 = []
		reswt2 = []
		resmut1 = []
		resmut2 = []
		cdwt1 = []
		cdmut1 = []
		cdwt2 =[]
		cdmut2 = []

		list1 = sys.argv[2].split(",")
		wt = dict()
		k = 1
		while k < len(list1):
			wt[int(list1[k])] = "in"
			k = k + 1
	
		# LOCAL RMSD BASED ON C-ALPHA CARBONS AND SIDE CHAINS

		pdb = "file2"
		structure_id = "{}".format(pdb)
		filename = "{}.pdb".format(pdb)
		structure = parser.get_structure(structure_id,filename)
		model = structure[0]

		# FOR WILD TYPE
		
		chain = model["A"]
		c1 = chain.get_list() 	# LIST ALL THE RESIDUES
		k1 = 0
		while k1 < len(c1):
			c2 = c1[k1].get_id()
			resid = c2[1]
			if resid in wt.keys():
				reswt1.append(resid)
				residue = chain[c2]	
				r1 = residue.get_list()		# LIST ALL THE ATOMS

				k2 = 0
				while k2 < len(r1):
					r2 = r1[k2].get_id()
					if r2 == "CA":
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdwt1.append(a1)
					if r2 != "CA" and r2 != "N" and r2 != "C" and r2 != "O" and r2[0:1] != "H":

						# ONLY SIDE CHAIN

						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdwt2.append(a1)

					k2 = k2 + 1

			k1 = k1 + 1
							
		# FOR MUTANT
		
		chain = model["B"]
		c1 = chain.get_list() 	# LIST ALL THE RESIDUES
		k1 = 0
		while k1 < len(c1):
			c2 = c1[k1].get_id()
			resid = c2[1]
			if resid in wt.keys():
				resmut1.append(resid)
				residue = chain[c2]	
				r1 = residue.get_list()		# LIST ALL THE ATOMS

				k2 = 0
				while k2 < len(r1):
					r2 = r1[k2].get_id()
					if r2 == "CA":
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdmut1.append(a1)
					if r2 != "CA" and r2 != "N" and r2 != "C" and r2 != "O" and r2[0:1] != "H":

						# ONLY SIDE CHAIN
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdmut2.append(a1)

					k2 = k2 + 1

			k1 = k1 + 1
							
		# CALCULATING LOCAL RMSD FOR C-ALPHA AND SIDE CHAINS

		if len(cdwt1) == len(cdmut1):
			lrmsd1 = RMSD(cdwt1,cdmut1,len(reswt1))
		else:
			lrmsd1 = "NA"
		if len(cdwt2) == len(cdmut2):
			lrmsd2 = RMSD(cdwt2,cdmut2,len(reswt1))
		else:
			lrmsd2 = "NA"

		g.write("{} {}".format(lrmsd1,lrmsd2))

	g.close()


import sys

if int(sys.argv[1]) == 0:
	func1()
elif int(sys.argv[1]) == 1:
	func2(0)
elif int(sys.argv[1]) == 2:
	func2(1)
