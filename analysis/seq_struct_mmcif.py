# STEP 2: GOING FROM SEQUENCE TO STRUCTURE

def COM(var1,var2):

	# TO DETERMINE THE CENTRE OF MASS OF THE ATOMS DEFINED BY VARIABLE var1 AND VARIABLE var2

	elem = ['C','N','O','S','P','H']
	mass = [12.0107,14.0067,15.9994,32.065,30.9738,1.0079]	# CORRESPONDING MASS OF THE ABOVE ELEMENTS
	avgx = 0.0
	avgy = 0.0
	avgz = 0.0
	M = 0.0
	for i in range(0,len(var1)):
		x = var2[i][0]
		y = var2[i][1]
		z = var2[i][2]
		e = var1[i][0:1]

		k = 0
		count = 0
		while k < len(elem):
			if elem[k] == e:
				count = count + 1
				w = mass[k]
			k = k + 1 
		if count == 0:
			w = 12.0

		avgx = avgx + (w*x)
		avgy = avgy + (w*y)
		avgz = avgz + (w*z)
		M = M + w

	if M > 0:
		comx = avgx / M
		comy = avgy / M
		comz = avgz / M
		return [comx,comy,comz]
	else:
		print ("UNKNOWN ATOM TYPE FOUND :: QUITING THE COM CALCULATION")
		quit()

def neighbours():
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.Polypeptide import three_to_one as tto

	# ZONE VARIABLE
	zdis = 10.0

	# PATH FOR THE PDB/mmCIF FILES

	pathmmcif = "/bmm/data/pdbmmcif/data/structures/all/mmCIF"
	#pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"
	#pathPDB = "/bmm/home/tkhanna1/Documents/Database/First_10000/test_set"

	f = open("number_of_distinct_mutants.csv","r")
	ft = f.readlines()
	f.close()

	g = open("number_of_distinct_mutations.csv","r")
	gt = g.readlines()
	g.close()

	h = open("structure.txt","w")

	# DETERMINING THE DISTINCT MUTATIONS AND MUTANTS

	dis = open("distinct_mutants.txt","w")
	k = 0
	while k < len(ft):
		dis_mut = []
		mu=ft[k].split(',')
		mut=gt[k].split(',')

		if len(mu) > 1:			
			pdb = mu[0]
			dis_mut.append(pdb)
			k1 = 1
			while k1 < len(mut):
				dumstr2 = mut[k1].strip("\n")
				mut_res=dumstr2[1:len(dumstr2)-1]
				k2 = 1
				count = 0
				while k2 < len(dis_mut):
					if dis_mut[k2] == mut_res:
						count = count + 1
					k2 = k2 + 1
				if count == 0:
					dis_mut.append(mut_res)
				k1 = k1 + 1
			dis.write("{}".format(dis_mut))
			dis.write("\n")
		k = k + 1
				
	dis.close()

	dis = open("distinct_mutants.txt","r")
	ht = dis.readlines()
	dis.close()

	# DETERMINING THE ZONE AROUND THE MUTANT SITE TO DETERMINE ANY STRUCTURAL CHANGE TAKING PLACE DUE TO THE MUTATION

	z = open("zone.txt","w")
	temp = open("COM.txt","w")

	k=0
	while k < len(ht): # end = len(ht)
		mutant = []
		mu=ht[k].split(',')

		pdbid=mu[0].strip('[|\,|\'|]')
		pdb=pdbid[0:4]		# PDB NAME
		C=pdbid[5:6]		# CHAIN
		
		print("*** {} :: {} of {} ***" .format(pdb,k,len(ht)))

		mutant.append(pdbid)

		# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

		try:
			pdbfile = "{}/{}.cif.gz".format(pathmmcif,pdb)
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

			k1 = 0
			resid_list = []
			resname_list = []
			com_list = []
			while k1 < len(c1):
				c2 = c1[k1].get_id()
				resid = c2[1]
				if c2[0] == " ":
					residue = chain[c2]
					tresname = residue.get_resname()
					try:
						resname = tto("{}".format(tresname))
					except:
						rename = "X"
					r1 = residue.get_list() # LIST ALL THE ATOMS OF A PARTICULAR RESIDUE

					k2 = 0
					res = []	# ATOM NAMES
					cd = []		# ATOM COORDINATES
					while k2 < len(r1):
						r2 = r1[k2].get_id()
						# COM OF THE BACKBONE
						if r2 == "CA" or r2 == "N" or r2 == "C" or r2 == "O":
							res.append(r2)

							atom = residue['{}'.format(r2)]
							a1 = atom.get_coord()
							cd.append(a1)
						k2 = k2 + 1
			
					CM = COM(res,cd)
					resid_list.append(resid)
					resname_list.append(resname)
					com_list.append(CM)
			
				k1 = k1 + 1
			temp.write("{}".format(com_list))
			temp.write("\n")

			# DETERMINING THE ZONE AROUND THE MUTATED RESIDUE

			k1 = 1
			while k1 < len(mu):
				pos= mu[k1].strip(' |[|,|]|\'|\n')
				pos = int(pos)
				z.write("{}".format(pdbid))
				z.write("\n")
				zres = "{}".format(pos)
				k2 = 0
				count = 0
				while k2 < len(resid_list):
					posc = int(resid_list[k2])
					if posc == pos:
						count = count + 1
						rcpos = k2
						k2 = len(resid_list)
					k2 = k2 + 1
			
				if count != 0:
					zresname = "{}".format(resname_list[rcpos])
					v1 = com_list[rcpos]
					# COMPUTING THE DISTANCES

					k3 = 0
					while k3 < len(com_list):
						if k3 != rcpos:
							v2 = com_list[k3]
							dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
							if dis < (zdis*zdis):
								zres = zres + ",{}".format(resid_list[k3])
								zresname = zresname +  "{}".format(resname_list[k3])
						k3 = k3 + 1
			
				z.write("{}".format(zres))
				z.write("\n")
				z.write("{}".format(zresname))
				z.write("\n")
				k1 = k1 +1
		except:
			print("FILE NOT FOUND")
			z.write("NA")
			z.write("\n")
			z.write("NA")
			z.write("\n")

		k = k + 1

	temp.close()
	z.close()
							
neighbours()


