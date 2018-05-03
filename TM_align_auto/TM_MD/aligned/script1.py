def dihedral(var1,var2):

	import numpy as np

	v1 = np.array(var1)
	v2 = np.array(var2)

	coord1 = v1[1] - v1[0]
	coord2 = v1[2] - v1[1]
	coord3 = v1[3] - v1[2]

	x1 = np.cross(coord1,coord2)
	x2 = np.cross(coord2,coord3)

	x1n = x1 / np.linalg.norm(x1)
	x2n = x2 / np.linalg.norm(x2)

	xdot = np.dot(x1n,x2n)
	xdot = round(xdot,3)
	theta1 = np.degrees(np.arccos(xdot))

	coord1 = v2[1] - v2[0]
	coord2 = v2[2] - v2[1]
	coord3 = v2[3] - v2[2]

	x1 = np.cross(coord1,coord2)
	x2 = np.cross(coord2,coord3)

	x1n = x1 / np.linalg.norm(x1)
	x2n = x2 / np.linalg.norm(x2)

	xdot = np.dot(x1n,x2n)
	xdot = round(xdot,3)
	theta2 = np.degrees(np.arccos(xdot))

	theta = round((theta1 - theta2),2)

	return(theta)

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

def RMSD_SC(var1,var2,var3,var4,var5,var6,var7,d1,d2,d3):

	import math

	# VAR1 = Coord WT
	# VAR2 = Coord MUT
	# VAR3 = NUM RES
	# VAR4 = AT WT
	# VAR5 = AT MUT
	# VAR6 = RESNAME
	# VAR7 = RESID
	# d1 = dict
	# d2 =dict

	# FORMING DICT OF WT AND RES ATOMTYPE AND RESIDUE NAME

	count = 0
	wt = dict()
	mut = dict()
	k = 0
	k1 = 0
	while k < len(var4):
		t1 = var7[k1]
		if var6[k1] != "G":
			numat = len(d1[var6[k1]])
			k2 = 0
			while k2 < numat:
				t2 = var4[k]
				t3 = var5[k]
				key1 = (t1,t2)
				key2 = (t1,t3)
				if key1 not in wt.keys():
					wt[key1] = var1[k]
				else:
					wt[key1].append(var1[k])
				if key2 not in mut.keys():
					mut[key2] = var2[k]
				else:
					mut[key2].append(var2[k])
				k = k + 1
				k2 = k2 + 1
		else:
			count = count + 1
		k1 = k1 + 1				

	# ORDERING THE LISTS FOR WILDTYPE AND MUTANT

	listWT = []
	listMUT= []
	k = 0
	while k < var3:
		t1 = var7[k]
		numat = len(d1[var6[k]])
		numrot = len(d2[var6[k]])
		if numat != 0:
			k1 = 0
			count1 = 0
			while k1 < numat:
				t2 = d1[var6[k]][k1]
				if numrot > 0:
					if t2 in d2[var6[k]] and count1 == 0:

						count1 = count1 + 1

						# CHECKING FOR ALTERNATE POSITIONS

						if var6[k][0] == t2:
							tkey1 = (t1,d2[var6[k]][0])
							tkey2 = (t1,d2[var6[k]][1])
						else:
							tkey2 = (t1,d2[var6[k]][0])
							tkey1 = (t1,d2[var6[k]][1])

						if var6[k] != "Y" and var6[k] != "F":
							coordwt1 = wt[tkey1]
							coordmut1 = mut[tkey1]
							coordwt2 = wt[tkey2]
							coordmut2 = mut[tkey2]


							dis1 = ((coordwt1[0]-coordmut1[0])**2) + ((coordwt1[1]-coordmut1[1])**2) + ((coordwt1[2]-coordmut1[2])**2)

							dis2 = ((coordwt2[0]-coordmut2[0])**2) + ((coordwt2[1]-coordmut2[1])**2) + ((coordwt2[2]-coordmut2[2])**2)
							
							dis11 = ((coordwt1[0]-coordmut2[0])**2) + ((coordwt1[1]-coordmut2[1])**2) + ((coordwt1[2]-coordmut2[2])**2)

							dis22 = ((coordwt2[0]-coordmut1[0])**2) + ((coordwt2[1]-coordmut1[1])**2) + ((coordwt2[2]-coordmut1[2])**2)
							
							if (dis11+dis22) < (dis1+dis2):
								mut[tkey1] = coordmut2
								mut[tkey2] = coordmut1
						else:

							# TO ACCOUNT FOR RING ROTATION

							tkey3 = (t1,d3[var6[k]][0])
							tkey4 = (t1,d3[var6[k]][1])
							coordwt1 = wt[tkey1]
							coordmut1 = mut[tkey1]
							coordwt2 = wt[tkey2]
							coordmut2 = mut[tkey2]
							coordwt3 = wt[tkey3]
							coordmut3 = mut[tkey3]
							coordwt4 = wt[tkey4]
							coordmut4 = mut[tkey4]

							dis1 = ((coordwt1[0]-coordmut1[0])**2) + ((coordwt1[1]-coordmut1[1])**2) + ((coordwt1[2]-coordmut1[2])**2)

							dis2 = ((coordwt2[0]-coordmut2[0])**2) + ((coordwt2[1]-coordmut2[1])**2) + ((coordwt2[2]-coordmut2[2])**2)

							dis3 = ((coordwt3[0]-coordmut3[0])**2) + ((coordwt3[1]-coordmut3[1])**2) + ((coordwt3[2]-coordmut3[2])**2)

							dis4 = ((coordwt4[0]-coordmut4[0])**2) + ((coordwt4[1]-coordmut4[1])**2) + ((coordwt4[2]-coordmut4[2])**2)
							
							dis11 = ((coordwt1[0]-coordmut2[0])**2) + ((coordwt1[1]-coordmut2[1])**2) + ((coordwt1[2]-coordmut2[2])**2)

							dis22 = ((coordwt2[0]-coordmut1[0])**2) + ((coordwt2[1]-coordmut1[1])**2) + ((coordwt2[2]-coordmut1[2])**2)

							dis33 = ((coordwt3[0]-coordmut4[0])**2) + ((coordwt3[1]-coordmut4[1])**2) + ((coordwt3[2]-coordmut4[2])**2)

							dis44 = ((coordwt4[0]-coordmut3[0])**2) + ((coordwt4[1]-coordmut3[1])**2) + ((coordwt4[2]-coordmut3[2])**2)
							
							if (dis11+dis22+dis33+dis44) < (dis1+dis2+dis3+dis4):
								mut[tkey1] = coordmut2
								mut[tkey2] = coordmut1
								mut[tkey3] = coordmut4
								mut[tkey4] = coordmut3
						
				key1 = (t1,t2)
				listWT.append(wt[key1])
				listMUT.append(mut[key1])
				k1 = k1 + 1
		k = k + 1
	
	k = 0
	l = 0.0
	while k < len(listWT):

		x1 = listWT[k][0]
		y1 = listWT[k][1]
		z1 = listWT[k][2]

		x2 = listMUT[k][0]
		y2 = listMUT[k][1]
		z2 = listMUT[k][2]

		l = l + (((x1-x2)**2) + ((y1-y2)**2) + ((z1-z2)**2))

		k = k + 1

	l = l / (len(listWT))
	RMSD = math.sqrt(l)
	RMSD = round(RMSD,2)

	return(RMSD)

def RMSD_dihedral(var1,var2,var3,var4,var5,var6,var7,d1,d2,d3,d4,c1,c2,c3,c4,c5):

	import math

	# VAR1 = Coord WT
	# VAR2 = Coord MUT
	# VAR3 = NUM RES
	# VAR4 = AT WT
	# VAR5 = AT MUT
	# VAR6 = RESNAME
	# VAR7 = RESID
	# d1 = dict
	# d2 = dict rotation
	# d3 = dict rotation2
	# d4 = dict all atom types
 	# c1 = dict chi1
 	# c2 = dict chi2
 	# c3 = dict chi3
 	# c4 = dict chi4
 	# c5 = dict chi5

	# FORMING DICT OF WT AND RES ATOMTYPE AND RESIDUE NAME

	wt = dict()
	mut = dict()
	k = 0
	k1 = 0
	while k < len(var4):
		t1 = var7[k1]
		numat = len(d4[var6[k1]])
		k2 = 0
		while k2 < numat:
			t2 = var4[k]
			t3 = var5[k]
			key1 = (t1,t2)
			key2 = (t1,t3)
			if key1 not in wt.keys():
				wt[key1] = var1[k]
			else:
				wt[key1].append(var1[k])
			if key2 not in mut.keys():
				mut[key2] = var2[k]
			else:
				mut[key2].append(var2[k])
			k = k + 1
			k2 = k2 + 1

		k1 = k1 + 1			

	# ORDERING THE LISTS FOR WILDTYPE AND MUTANT

	k = 0
	while k < var3:
		t1 = var7[k]
		numat = len(d1[var6[k]])
		numrot = len(d2[var6[k]])
		if numat != 0:
			k1 = 0
			count1 = 0
			while k1 < numat:
				t2 = d1[var6[k]][k1]
				if numrot > 0:
					if t2 in d2[var6[k]] and count1 == 0:

						count1 = count1 + 1

						# CHECKING FOR ALTERNATE POSITIONS

						if var6[k][0] == t2:
							tkey1 = (t1,d2[var6[k]][0])
							tkey2 = (t1,d2[var6[k]][1])
						else:
							tkey2 = (t1,d2[var6[k]][0])
							tkey1 = (t1,d2[var6[k]][1])

						if var6[k] != "Y" and var6[k] != "F":
							coordwt1 = wt[tkey1]
							coordmut1 = mut[tkey1]
							coordwt2 = wt[tkey2]
							coordmut2 = mut[tkey2]

							dis1 = ((coordwt1[0]-coordmut1[0])**2) + ((coordwt1[1]-coordmut1[1])**2) + ((coordwt1[2]-coordmut1[2])**2)

							dis2 = ((coordwt2[0]-coordmut2[0])**2) + ((coordwt2[1]-coordmut2[1])**2) + ((coordwt2[2]-coordmut2[2])**2)
							
							dis11 = ((coordwt1[0]-coordmut2[0])**2) + ((coordwt1[1]-coordmut2[1])**2) + ((coordwt1[2]-coordmut2[2])**2)

							dis22 = ((coordwt2[0]-coordmut1[0])**2) + ((coordwt2[1]-coordmut1[1])**2) + ((coordwt2[2]-coordmut1[2])**2)
							
							if (dis11+dis22) < (dis1+dis2):
								mut[tkey1] = coordmut2
								mut[tkey2] = coordmut1
						else:

							# TO ACCOUNT FOR RING ROTATION

							tkey3 = (t1,d3[var6[k]][0])
							tkey4 = (t1,d3[var6[k]][1])
							coordwt1 = wt[tkey1]
							coordmut1 = mut[tkey1]
							coordwt2 = wt[tkey2]
							coordmut2 = mut[tkey2]
							coordwt3 = wt[tkey3]
							coordmut3 = mut[tkey3]
							coordwt4 = wt[tkey4]
							coordmut4 = mut[tkey4]

							dis1 = ((coordwt1[0]-coordmut1[0])**2) + ((coordwt1[1]-coordmut1[1])**2) + ((coordwt1[2]-coordmut1[2])**2)

							dis2 = ((coordwt2[0]-coordmut2[0])**2) + ((coordwt2[1]-coordmut2[1])**2) + ((coordwt2[2]-coordmut2[2])**2)

							dis3 = ((coordwt3[0]-coordmut3[0])**2) + ((coordwt3[1]-coordmut3[1])**2) + ((coordwt3[2]-coordmut3[2])**2)

							dis4 = ((coordwt4[0]-coordmut4[0])**2) + ((coordwt4[1]-coordmut4[1])**2) + ((coordwt4[2]-coordmut4[2])**2)
							
							dis11 = ((coordwt1[0]-coordmut2[0])**2) + ((coordwt1[1]-coordmut2[1])**2) + ((coordwt1[2]-coordmut2[2])**2)

							dis22 = ((coordwt2[0]-coordmut1[0])**2) + ((coordwt2[1]-coordmut1[1])**2) + ((coordwt2[2]-coordmut1[2])**2)

							dis33 = ((coordwt3[0]-coordmut4[0])**2) + ((coordwt3[1]-coordmut4[1])**2) + ((coordwt3[2]-coordmut4[2])**2)

							dis44 = ((coordwt4[0]-coordmut3[0])**2) + ((coordwt4[1]-coordmut3[1])**2) + ((coordwt4[2]-coordmut3[2])**2)
							
							if (dis11+dis22+dis33+dis44) < (dis1+dis2+dis3+dis4):
								mut[tkey1] = coordmut2
								mut[tkey2] = coordmut1
								mut[tkey3] = coordmut4
								mut[tkey4] = coordmut3

				k1 = k1 + 1
		k = k + 1

	# LOCAL RMSD CHI ANGLES OF THE SIDE CHAIN

	k = 0
	numc1 = 0
	numc2 = 0
	numc3 = 0
	numc4 = 0
	numc5 = 0
	vc1 = 0.0
	vc2 = 0.0
	vc3 = 0.0
	vc4 = 0.0
	vc5 = 0.0
	while k < len(var6):
		t1 = var6[k]
		t2 = var7[k]

		# CHI1
		N = c1[t1]
		list1 = []
		list2 = []
		k1 = 0
		while k1 < len(N):
			at = N[k1]
			key1 = (t2,at)
			cwt = wt[key1]
			cmut = mut[key1]
			list1.append(cwt)
			list2.append(cmut)
			k1 = k1 + 1

		if len(N) != 0:
			numc1 = numc1 + 1
			dih = dihedral(list1,list2)
			dih = dih * dih
			vc1 = vc1 + dih

		# CHI2
		N = c2[t1]
		list1 = []
		list2 = []
		k1 = 0
		while k1 < len(N):
			at = N[k1]
			key1 = (t2,at)
			cwt = wt[key1]
			cmut = mut[key1]
			list1.append(cwt)
			list2.append(cmut)
			k1 = k1 + 1

		if len(N) != 0:
			numc2 = numc2 + 1
			dih = dihedral(list1,list2)
			dih = dih * dih
			vc2 = vc2 + dih

		# CHI3
		N = c3[t1]
		list1 = []
		list2 = []
		k1 = 0
		while k1 < len(N):
			at = N[k1]
			key1 = (t2,at)
			cwt = wt[key1]
			cmut = mut[key1]
			list1.append(cwt)
			list2.append(cmut)
			k1 = k1 + 1

		if len(N) != 0:
			numc3 = numc3 + 1
			dih = dihedral(list1,list2)
			dih = dih * dih
			vc3 = vc3 + dih

		# CHI4
		N = c4[t1]
		list1 = []
		list2 = []
		k1 = 0
		while k1 < len(N):
			at = N[k1]
			key1 = (t2,at)
			cwt = wt[key1]
			cmut = mut[key1]
			list1.append(cwt)
			list2.append(cmut)
			k1 = k1 + 1

		if len(N) != 0:
			numc4 = numc4 + 1
			dih = dihedral(list1,list2)
			dih = dih * dih
			vc4 = vc4 + dih

		# CHI5
		N = c5[t1]
		list1 = []
		list2 = []
		k1 = 0
		while k1 < len(N):
			at = N[k1]
			key1 = (t2,at)
			cwt = wt[key1]
			cmut = mut[key1]
			list1.append(cwt)
			list2.append(cmut)
			k1 = k1 + 1

		if len(N) != 0:
			numc5 = numc5 + 1
			dih = dihedral(list1,list2)
			dih = dih * dih
			vc5 = vc5 + dih
	
		k = k + 1

	if numc1 != 0:
		vc1 = vc1 / numc1
		vc1 = math.sqrt(vc1)
		vc1 = round(vc1,2)
	else:
		vc1 = "NA"

	if numc2 != 0:
		vc2 = vc2 / numc2
		vc2 = math.sqrt(vc2)
		vc2 = round(vc2,2)
	else:
		vc2 = "NA"

	if numc3 != 0:
		vc3 = vc3 / numc3
		vc3 = math.sqrt(vc3)
		vc3 = round(vc3,2)
	else:
		vc3 = "NA"

	if numc4 != 0:
		vc4 = vc4 / numc4
		vc4 = math.sqrt(vc4)
		vc4 = round(vc4,2)
	else:
		vc4 = "NA"

	if numc5 != 0:
		vc5 = vc5 / numc5
		vc5 = math.sqrt(vc5)
		vc5 = round(vc5,2)
	else:
		vc5 = "NA"

	return(vc1,vc2,vc3,vc4,vc5)

def func1():
	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.PDBParser import PDBParser
	parser1 = PDBParser(PERMISSIVE=0,QUIET=True)

	from Bio.PDB.PDBIO import PDBIO

	#pathmmcif = "/Users/tarun/Documents/mmCIF"
	#pathmmcif = "/data/pdb/divided/mmCIF"
	#pathmmcif = "/Volumes/RCSB_DATA/pdb"
	pathmmcif = "/Volumes/BIOINFO/mmCIF/"
	pathmmcif1 = "/Volumes/BIOINFO/MDoutput/NO_WATER/"

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

		pdbfile = "{}/{}.pdb".format(pathmmcif1,pdb2)
		#pdbfile = "{}/{}/pdb{}.ent.gz".format(pathmmcif,fol,pdb2)
		#tar = gzip.open("{}".format(pdbfile),"rb")
		#out = open("pdbprocess.cif","wb")
		#out = open("pdbprocess.pdb","wb")
		#out.write(tar.read())
		#tar.close()
		#out.close()
		structure_id = "{}".format(pdb2)
		filename = pdbfile
		#filename = "pdbprocess.pdb"
		structure = parser1.get_structure(structure_id,filename)
		modid = "{}".format(sys.argv[6])
		model = structure[int(modid)]
		#chain = model["{}".format(c2)]

		io = PDBIO()
		io.set_structure(model)
		io.save("chain2.pdb")
	except:
		print("FILE NOT FOUND")

def func2(arg1):

	from Bio.PDB.Polypeptide import three_to_one as tto

	AA = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]

	AT = [["CB"],["CB","CG","CD","NE","CZ","NH1","NH2"],["CB","CG","OD1","ND2"],["CB","CG","OD1","OD2"],["CB","SG"],["CB","CG","CD","OE1","OE2"],["CB","CG","CD","OE1","NE2"],[],["CB","CG","ND1","CE1","NE2","CD2"],["CB","CG1","CG2","CD1"],["CB","CG","CD1","CD2"],["CB","CG","CD","CE","NZ"],["CB","CG","SD","CE"],["CB","CG","CD1","CE1","CZ","CE2","CD2"],["CB","CG","CD"],["CB","OG"],["CB","CG2","OG1"],["CB","CG","CD1","NE1","CE2","CD2","CE3","CZ3","CH2","CZ2"],["CB","CG","CD1","CE1","CZ","OH","CE2","CD2"],["CB","CG1","CG2"]]

	BA = ["CA","C","N","O"]

	IN = [[],["NH1","NH2"],[],["OD1","OD2"],[],["OE1","OE2"],[],[],[],[],["CD1","CD2"],[],[],["CD1","CD2"],[],[],[],[],["CD1","CD2"],["CG1","CG2"]]

	IN1 = [[],[],[],[],[],[],[],[],[],[],[],[],[],["CE1","CE2"],[],[],[],[],["CE1","CE2"],[]]

	chi1 = [[],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","SG"],["N","CA","CB","CG"],["N","CA","CB","CG"],[],["N","CA","CB","CG"],["N","CA","CB","CG1"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","OG"],["N","CA","CB","OG1"],["N","CA","CB","CG"],["N","CA","CB","CG"],["N","CA","CB","CG1"]]

	chi2 = [[],["CA","CB","CG","CD"],["CA","CB","CG","OD1"],["CA","CB","CG","OD1"],[],["CA","CB","CG","CD"],["CA","CB","CG","CD"],[],["CA","CB","CG","ND1"],["CA","CB","CG1","CD1"],["CA","CB","CG","CD1"],["CA","CB","CG","CD"],["CA","CB","CG","SD"],["CA","CB","CG","CD1"],["CA","CB","CG","CD"],[],[],["CA","CB","CG","CD1"],["CA","CB","CG","CD1"],[]]

	chi3 = [[],["CB","CG","CD","NE"],[],[],[],["CB","CG","CD","OE1"],["CB","CG","CD","OE1"],[],[],[],[],["CB","CG","CD","CE"],["CB","CG","SD","CE"],[],[],[],[],[],[],[]]

	chi4 = [[],["CG","CD","NE","CZ"],[],[],[],[],[],[],[],[],[],["CG","CD","CE","NZ"],[],[],[],[],[],[],[],[]]

	chi5 = [[],["CD","NE","CZ","NH1"],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]

	k = 0
	datall = dict()
	dat = dict()
	drot = dict()
	drot1 = dict()
	chia1 = dict()
	chia2 = dict()
	chia3 = dict()
	chia4 = dict()
	chia5 = dict()
	while k < len(AA):
		AT_all = BA + AT[k]
		datall[AA[k]] = AT_all
		dat[AA[k]] = AT[k]
		drot[AA[k]] = IN[k]
		drot1[AA[k]] = IN1[k]
		chia1[AA[k]] = chi1[k]
		chia2[AA[k]] = chi2[k]
		chia3[AA[k]] = chi3[k]
		chia4[AA[k]] = chi4[k]
		chia5[AA[k]] = chi5[k]
		k = k + 1

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
		atwt2 = []
		atmut2 = []
		resnamewt1 = []
		cdwt3 = []
		cdmut3 = []
		cdwt3 =[]
		cdmut3 = []
		atwt3 = []
		atmut3 = []
		
		list1 = sys.argv[2].split(",")
		wt = dict()
		k = 1
		while k < len(list1):
			wt[int(list1[k])] = "in"
			k = k + 1

		# FOR MD TRAJECTORY LOCAL ZONE

		list1 = sys.argv[3].split(",")
		wt1 = dict()
		k = 1
		while k < len(list1):
			wt1[int(list1[k])] = "in"
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
				tresname = residue.get_resname()
				resname = tto("{}".format(tresname))
				resnamewt1.append(resname)
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
						atwt2.append(r2)
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdwt2.append(a1)

					if r2[0:1] != "H":

						# ALL ATOM TYPES
						atwt3.append(r2)
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdwt3.append(a1)

					k2 = k2 + 1

			k1 = k1 + 1
							
		# FOR MUTANT
		
		chain = model["B"]
		c1 = chain.get_list() 	# LIST ALL THE RESIDUES
		k1 = 0
		while k1 < len(c1):
			c2 = c1[k1].get_id()
			resid = c2[1]
			if resid in wt1.keys():
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
						atmut2.append(r2)
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdmut2.append(a1)

					if r2[0:1] != "H":
			
						# ALL ATOM TYPES
						atmut3.append(r2)
						atom = residue["{}".format(r2)]
						a1 = atom.get_coord()
						cdmut3.append(a1)


					k2 = k2 + 1

			k1 = k1 + 1
							
		# CALCULATING LOCAL RMSD FOR C-ALPHA AND SIDE CHAINS

		if len(cdwt1) == len(cdmut1):
			lrmsd1 = RMSD(cdwt1,cdmut1,len(reswt1))
		else:
			lrmsd1 = "NA"
		if len(cdwt2) == len(cdmut2):
			lrmsd2 = RMSD_SC(cdwt2,cdmut2,len(reswt1),atwt2,atmut2,resnamewt1,reswt1,dat,drot,drot1)
		else:
			lrmsd2 = "NA"
		if len(cdwt3) == len(cdmut3):
			lrmsd3 = RMSD_dihedral(cdwt3,cdmut3,len(reswt1),atwt3,atmut3,resnamewt1,reswt1,dat,drot,drot1,datall,chia1,chia2,chia3,chia4,chia5)
		else:
			lrmsd3 = ["NA","NA","NA","NA","NA"]

		g.write("{} {} {} {} {} {} {}".format(lrmsd1,lrmsd2,lrmsd3[0],lrmsd3[1],lrmsd3[2],lrmsd3[3],lrmsd3[4]))

	g.close()


import sys

if int(sys.argv[1]) == 0:
	func1()
elif int(sys.argv[1]) == 1:
	func2(0)
elif int(sys.argv[1]) == 2:
	func2(1)
