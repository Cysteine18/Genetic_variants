def dihedral(p):
	# THIS FUNCTION DETERMINE THE DIHEDRAL ANGLE BETWEEN 4 POINTS PRESENT IN P

	import numpy as np

	import math

	b1 = p[0] - p[1]
	normb1 = np.linalg.norm(b1)
	b1 = b1 / normb1

	b2 = p[2] - p[1]
	normb2 = np.linalg.norm(b2)
	b2 = b2 / normb2

	b3 = p[3] - p[2]
	normb3 = np.linalg.norm(b3)
	b3 = b3 / normb3

	b1_c_b2 = np.cross(b1,b2)
	b3_c_b2 = np.cross(b3,b2)

	ang = np.dot(b1_c_b2,b3_c_b2)
	ang = math.acos(ang)
	ang = round(ang,3)
	ang = (ang * 180) / 3.14
	ang = round(ang,2)

	return(ang)

def disulphide_bond():
	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.Polypeptide import three_to_one as tto

	import numpy as np

	# ZONE VARIABLE
	zdis = 2.2

	# PATH FOR THE PDB/mmCIF FILES

	#pathmmcif = "/bmm/data/pdbmmcif/data/structures/all/mmCIF"
	pathmmcif = "/Volumes/BIOINFO/mmCIF"
	#pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"
	#pathPDB = "/bmm/home/tkhanna1/Documents/Database/First_10000/test_set"

	dis = open("mut_data.txt","r")
	ht = dis.readlines()
	dis.close()

	# DETERMINING THE ZONE AROUND THE MUTANT SITE TO DETERMINE ANY STRUCTURAL CHANGE TAKING PLACE DUE TO THE MUTATION

	z = open("{}".format(sys.argv[3]),"w")
	#temp = open("COM.txt","w")

	start = sys.argv[1]
	end = sys.argv[2]
	if end == "END" or end == "end":
		end = len(ht)
	end = int(end)
	k=int(start)
	while k < end: # end = len(ht)
		mutant = []
		mu=ht[k].split()

		check1 = mu[3].strip("\n")
		check2 = mu[4].strip("\n")
		pos = mu[2].strip("\n")

		count3 = 0
		if check1 == "C":
			pdbid=mu[0].strip("\n")
			pdb=pdbid[0:4]		# PDB NAME
			C=pdbid[5:6]		# CHAIN
			count3 = count3 + 1
		elif check2 == "C":
			pdbid=mu[1].strip("\n")
			pdb=pdbid[0:4]		# PDB NAME
			C=pdbid[5:6]		# CHAIN
			count3 = count3 + 1

		if count3 != 0:
			
			print("*** {} :: {} of {} ***" .format(pdb,k,len(ht)))

			# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

			#cc = 0
			#if cc == 0:
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
				com_list1 = []
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
						if resname == "C":
							r1 = residue.get_list() # LIST ALL THE ATOMS OF A PARTICULAR RESIDUE

							k2 = 0
							res = []	# ATOM NAMES
							cd = []		# ATOM COORDINATES CB
							cd1 = []		# ATOM COORDINATES SG
							count = 0
							while k2 < len(r1):
								r2 = r1[k2].get_id()
								# COM OF THE BACKBONE
								if r2 == "CB":
									res.append(r2)
									atom = residue['{}'.format(r2)]
									a1 = atom.get_coord()
									list1 = [a1[0],a1[1],a1[2]]
									cd.append(list1)
									count = count + 1
								if r2 == "SG":
									res.append(r2)
									atom = residue['{}'.format(r2)]
									a1 = atom.get_coord()
									list1 = [a1[0],a1[1],a1[2]]
									cd1.append(list1)
									count = count + 1
							
								k2 = k2 + 1
						
							if count == 2:
								resid_list.append(resid)
								resname_list.append(resname)
								list1 = [cd[0][0],cd[0][1],cd[0][2]]
								com_list.append(list1)
								list1 = [cd1[0][0],cd1[0][1],cd1[0][2]]
								com_list1.append(list1)
					
					k1 = k1 + 1

				# DETERMINING THE ZONE AROUND THE MUTATED RESIDUE

				k1 = 0
				count = 0
				while k1 < len(resid_list):
					posc = resid_list[k1]
					count = 0
					if int(posc) == int(pos):
						count = count + 1
						rcpos = k1
						k1 = len(resid_list)
					k1 = k1 + 1
				
				if count != 0:
					zresname = "{}".format(resname_list[rcpos])
					v1 = com_list1[rcpos]
					v11 = com_list[rcpos]

					# COMPUTING THE DISTANCES

					k3 = 0
					list1 = []
					while k3 < len(com_list1):
						if k3 != rcpos:
							v2 = com_list1[k3]
							v21 = com_list[k3]
							dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
							if dis < (zdis*zdis):
								p = np.array([v1,v11,v2,v21])
								ang = dihedral(p)
								# DIHEDRAL CRITERIA
								if ang > 75.0 and ang < 105.0:
									list1.append(resid_list[k3])
						k3 = k3 + 1

					if len(list1) > 0:
						z.write("{} {} {} {} {} YES {}\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1))
						print("{} {} {} {} {} YES {}".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1))
					else:
						z.write("{} {} {} {} {} NO\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
				else:
						z.write("{} {} {} {} {} NO\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
		
			except:
				print("FILE NOT FOUND")
				z.write("{} {} {} {} {} NO\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))

		else:
			#print("NOT CYSTINE")
			z.write("{} {} {} {} {} NO\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))

		k = k + 1

	z.close()

#distinct_nei()				
disulphide_bond()


