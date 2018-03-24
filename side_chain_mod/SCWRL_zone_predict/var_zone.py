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
		return []

def neighbours():
	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.Polypeptide import three_to_one as tto

	# ZONE VARIABLE
	#zdis = 10.0

	# PATH FOR THE PDB/mmCIF FILES

	pathmmcif = "/Volumes/BIOINFO/mmCIF"
	#pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"
	#pathPDB = "/bmm/home/tkhanna1/Documents/Database/First_10000/test_set"

	# DETERMINING THE ZONE AROUND THE MUTANT SITE TO DETERMINE ANY STRUCTURAL CHANGE TAKING PLACE DUE TO THE MUTATION

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	z = open("{}.txt".format(sys.argv[4]),"w")
	#temp = open("COM.txt","w")

	start = sys.argv[2]
	end = (sys.argv[3])
	if end == "END" or end == "end":
		end = len(ft)
	end = int(end)
	k=int(start)
	while k < end:
		mutant = []
		mu=ft[k].split()

		pdbid=mu[0]
		pdb=pdbid[0:4]		# PDB NAME
		C=pdbid[5:6]		# CHAIN

		pdbidmut = mu[1]
		pdbmut = pdbidmut[0:4]
		Cmut = pdbidmut[5:6]
		reswt = mu[3]
		resmut = mu[4].strip("\n")
		
		print("*** {} :: {} of {} ***" .format(pdb,k,end))

		mutant.append(pdbid)

		# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

		#ccount = 0
		#if ccount == 0:
		try:
			fol = pdb[1:3]
			pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess{}.cif".format(start),"wb")
			out.write(tar.read())
			tar.close()
			out.close()

			structure_id = "{}".format(pdb)
			filename = "pdbprocess{}.cif".format(start)
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
						resname = "X"
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
			#temp.write("{}".format(com_list))
			#temp.write("\n")

			# DETERMINING THE ZONE AROUND THE MUTATED RESIDUE

			pos= mu[2].strip("\n")
			pos = int(pos)
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
				list1 = [pos]
				for zdis in range(4,11,1):
					zres = "{}".format(pos)
					v1 = com_list[rcpos]
					# COMPUTING THE DISTANCES
					k3 = 0
					while k3 < len(com_list):
						if k3 != rcpos:
							v2 = com_list[k3]
							dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
							if dis < (zdis*zdis):
								zres = zres + ",{}".format(resid_list[k3])
						k3 = k3 + 1
					list1.append(zres)

				z.write("{} {} {} {} {} {} {} {} {} {} {} {} {} {}\n".format(pdb,C,pdbmut,Cmut,list1[0],list1[1],list1[2],list1[3],list1[4],list1[5],list1[6],list1[7],reswt,resmut))
		except:
			print("FILE NOT FOUND")
			#z.write("{}".format(pdbid))
			#z.write("\n")
			#z.write("NA")
			#z.write("\n")
			#z.write("NA")
			#z.write("\n")

		k = k + 1

	#temp.close()
	z.close()

#distinct_nei()				
neighbours()


