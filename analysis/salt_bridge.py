def salt_bridge():
	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.Polypeptide import three_to_one as tto

	import numpy as np

	ALL = dict()
	ALL["D"] = "SB"
	ALL["E"] = "SB"
	ALL["K"] = "SB"
	ALL["R"] = "SB"
	ALL["H"] = "SB"

	ALLN = dict()
	ALLN["D"] = "SB"
	ALLN["E"] = "SB"

	ALLP = dict()
	ALLP["K"] = "SB"
	ALLP["R"] = "SB"
	ALLP["H"] = "SB"

	neg = ["D","E"]
	negat = ["CG","CD"]
	posi = ["K","R","H"]
	posat = ["NZ","NE",["ND1","NE2"]]

	# ZONE VARIABLE
	zdis = 4.0

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
		pos = mu[2].strip("\n")
		check1 = mu[3].strip("\n")
		check2 = mu[4].strip("\n")

		kk = 0
		while kk < 2:
			if check1 in ALL.keys() or check2 in ALL.keys():

				pdbid=mu[(kk+0)].strip("\n")
				pdb=pdbid[0:4]		# PDB NAME
				C=pdbid[5:6]		# CHAIN
				
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

							count = 0
							kn = 0
							while kn < len(neg):
								if resname == neg[kn]:
									count = count + 1
									ksel = kn
								kn = kn + 1 

							if count == 0:
								count1 = 0
								kn = 0
								while kn < len(posi):
									if resname == posi[kn]:
										count1 = count1 + 1
										ksel = kn
									kn = kn + 1 

							if count != 0:
								r1 = residue.get_list() # LIST ALL THE ATOMS OF A PARTICULAR RESIDUE

								k2 = 0
								res = []	# ATOM NAMES
								cd = []		# ATOM COORDINATES 
								countc = 0
								while k2 < len(r1):
									r2 = r1[k2].get_id()
									if r2 == negat[ksel]:

										res.append(r2)
										atom = residue['{}'.format(r2)]
										a1 = atom.get_coord()
										list1 = [a1[0],a1[1],a1[2]]
										cd.append(list1)
										countc = countc + 1
								
									k2 = k2 + 1
							
								if countc == 1:
									resid_list.append(resid)
									resname_list.append(resname)
									list1 = [cd[0][0],cd[0][1],cd[0][2]]
									com_list.append(list1)

							if count1 != 0:
								r1 = residue.get_list() # LIST ALL THE ATOMS OF A PARTICULAR RESIDUE
								k2 = 0
								res = []	# ATOM NAMES
								cd = []		# ATOM COORDINATES 
								countc = 0
								while k2 < len(r1):
									r2 = r1[k2].get_id()
									if ksel != 2:
										if r2 == posat[ksel]:
											res.append(r2)
											atom = residue['{}'.format(r2)]
											a1 = atom.get_coord()
											list1 = [a1[0],a1[1],a1[2]]
											cd.append(list1)
											countc = countc + 1
									else:
										if r2 == posat[ksel][0] or r2 == posat[ksel][1]:
											res.append(r2)
											atom = residue['{}'.format(r2)]
											a1 = atom.get_coord()
											list1 = [a1[0],a1[1],a1[2]]
											cd.append(list1)
											countc = countc + 1
									
									k2 = k2 + 1
							
								if countc == 1:
									resid_list.append(resid)
									resname_list.append(resname)
									list1 = [cd[0][0],cd[0][1],cd[0][2]]
									com_list.append(list1)
								if countc == 2:
									resid_list.append(resid)
									resname_list.append(resname)
									list1 = [cd[0][0],cd[0][1],cd[0][2],cd[1][0],cd[1][1],cd[1][2]]
									com_list.append(list1)
						
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
						hd = -1
						hd1 = -1
						zresname = "{}".format(resname_list[rcpos])
						if zresname == "H":
							v1 = com_list[rcpos]
							# COMPUTING THE DISTANCES
							list1 = []
							for x in range(0,2):
								skip = (x*3)
								k3 = 0
								while k3 < len(com_list):
									if k3 != rcpos:
										v2 = com_list[k3]
										zresnamec = "{}".format(resname_list[k3])
										if zresnamec in ALLN.keys():
											if zresnamec == "H" and int(hd) != int(resid_list[k3]):
												for y in range(0,2):
													skip2 = y*3
													dis = (v1[0+skip]-v2[0+skip2])**2 + (v1[1+skip]-v2[1+skip2])**2 + (v1[2+skip]-v2[2+skip2])**2
													if dis < (zdis*zdis):
														list1.append(resid_list[k3])
														hd = resid_list[k3]
														y = 3
											elif zresnamec != "H" and hd1 != 3:
												dis = (v1[0+skip]-v2[0])**2 + (v1[1+skip]-v2[1])**2 + (v1[2+skip]-v2[2])**2
												if dis < (zdis*zdis):
													list1.append(resid_list[k3])
													hd1 = 3
									k3 = k3 + 1
						else:
							v1 = com_list[rcpos]
							if zresname in ALLP.keys():
								k3 = 0
								list1 = []
								while k3 < len(com_list):
									if k3 != rcpos:
										v2 = com_list[k3]
										zresnamec = "{}".format(resname_list[k3])
										if zresnamec in ALLN.keys():
											if zresnamec == "H" and int(hd) != int(resid_list[k3]):
												for x in range(0,2):
													skip = (x*3)
													dis = (v1[0]-v2[0+skip])**2 + (v1[1]-v2[1+skip])**2 + (v1[2]-v2[2+skip])**2
													if dis < (zdis*zdis):
														list1.append(resid_list[k3])
														hd = resid_list[k3]
														x = 2
											elif zresnamec != "H":
												dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
												if dis < (zdis*zdis):
													list1.append(resid_list[k3])
									k3 = k3 + 1

							elif zresname in ALLN.keys():
								k3 = 0
								list1 = []
								while k3 < len(com_list):
									if k3 != rcpos:
										v2 = com_list[k3]
										zresnamec = "{}".format(resname_list[k3])
										if zresnamec in ALLP.keys():
											if zresnamec == "H" and int(hd) != int(resid_list[k3]):
												for x in range(0,2):
													skip = (x*3)
													dis = (v1[0]-v2[0+skip])**2 + (v1[1]-v2[1+skip])**2 + (v1[2]-v2[2+skip])**2
													if dis < (zdis*zdis):
														list1.append(resid_list[k3])
														hd = resid_list[k3]
														x = 2
											elif zresnamec != "H":
												dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
												if dis < (zdis*zdis):
													list1.append(resid_list[k3])
									k3 = k3 + 1

						if len(list1) > 0:
							if kk == 0:
								z.write("{} {} {} {} {} YES WT {} {}\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1,len(list1)))
								print("{} {} {} {} {} YES {} WT".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1))
							else:
								z.write("{} {} {} {} {} YES MUT {} {}\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1,len(list1)))
								print("{} {} {} {} {} YES {} MUT".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n"),list1))

						else:
							if kk == 0:
								z.write("{} {} {} {} {} NO WT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
							else:
								z.write("{} {} {} {} {} NO MUT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
					else:
						if kk == 0:
							z.write("{} {} {} {} {} NO WT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
						else:
							z.write("{} {} {} {} {} NO MUT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
			
				except:
					print("FILE NOT FOUND")
					if kk == 0:
						z.write("{} {} {} {} {} NO WT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
					else:
						z.write("{} {} {} {} {} NO MUT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))

				kk = kk + 1

			else:
				kk = 2
				#print("NOT SALT BRIDGE RESIDUES")
				z.write("{} {} {} {} {} NO WT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))
				z.write("{} {} {} {} {} NO MUT\n".format(mu[0],mu[1],mu[2],mu[3],mu[4].strip("\n")))

		k = k + 1
	
	z.close()

#distinct_nei()				
salt_bridge()


