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
	from Bio.PDB.PDBParser import PDBParser
	parser = PDBParser(PERMISSIVE=0)

	# PATH FOR THE PDB/mmCIF FILES

	pathmmcif = "/bmm/data/rcsb/data/structures/all/mmCIF"
	pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"

	f = open("mutants.txt","r")
	ft = f.readlines()
	f.close()

	g = open("mutations.txt","r")
	gt = g.readlines()
	g.close()

	h = open("structure","w")

	# DETERMINING THE DISTINCT MUTATIONS AND MUTANTS

	dis = open("distinct_mutants","w")
	k = 0
	while k < len(ft):
		dis_mut = []
		mu=ft[k].split(',')
		mut=gt[k].split(',')

		if len(mu) > 1:			
			pdb = mu[0].strip("\n")
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
	quit()

	# DUMMY VARIABLE TO REMOVE REPETITION IN THE CALCULATIONS

	mutant = []

	k=1
	while k < len(f):
		mu=ft[k].split(',')
		mut=gt[k].split(',')
		dumstr = mu
	
		k1=1
		while k1 < len(dumstr):
			pdb=mu[k1][0:4]
			C=mu[k1][4:5]
			dumstr2 = mut[k1].strip("\n")
			mut_res=dumstr2[1:len(dumstr2)-1]

			# PREVENTING FROM PROCESSING THE PDB AGAIN

			k2=0
			count=0
			while k2 < len(mutant):
				if mutant[k2] == pdb:
					count = count + 1
					pos = k2
				k2=k2+1

			if count == 0:

				# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

				pdbfile = "{}/pdb{}.ent.gz".format(pathPDB,pdb)
				tar = gzip.open("{}".format(pdbfile),"rb")
				out = open("pdbprocess.pdb","wb")
				out.write(tar.read())
				tar.close()
				out.close()

				mutant.append(pdb)

				structure_id = "{}".format(pdb)
				filename = "pdb{}.ent".format(pdb)
				structure = parser.get_structure(structure_id,filename)

				model = structure[0]
				m1 = model.get_list()	# THIS GIVES THE CHAINS IN THE PDB FILE

				chain = model["{}".format(C)]

				residue = chain[mut_res]

				r1 = residue.get_resname()
				r2 = residue.is_disordered()
				r3 = residue.get_unpacked_list() # LIST ALL THE ATOMS


				
neighbours()





















