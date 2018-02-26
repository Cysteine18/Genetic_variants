from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import gzip

pathfiles = "/Users/tarunkhanna/Documents/Bioinformatics/file_links/"

def res_filter():
	f = open("{}/PDB_classifier.txt".format(pathfiles),"r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	d1 = dict()
	d2 = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
		date = ft1[2].strip("\n") # DATE
		year = date.split("-")
		t3 = year[0]	# YEAR
		t4 = ft1[4].strip("\n")	# TYPE OF EXPERIMENT 
		if t1 not in d.keys():
			d["{}".format(t1)] = t2
			d1["{}".format(t1)] = t3
			d2["{}".format(t1)] = t4

		k = k + 1

	return(d,d1,d2)

def HETATM_KW():

	# DETERMINATION OF HETATM BASED ON KEYWORD

	keywords = ["BOUND","COMPLEXED","COMPLEX"]
	# IF THE KEYWORD IS "COMPLEX" THEN IT HAS TO BE FOLLOWED BY "IN"

	f = open("{}/entries.idx".format(pathfiles),"r")
	ft = f.readlines()
	f.close()

	#g = open("bound.txt","w")

	d = dict()
	k = 0
	while k < len(ft):
		#print("CHECKING ROW {} OF {}".format(k,end))
		ft1 = ft[k].split()
		k1 = 0
		count = 0
		while k1 < len(ft1):
			t2 = ft1[k1].strip("\n")
			k2 = 0
			while k2 < len(keywords):
				if t2 == keywords[k2]:
					if t2 == "COMPLEX":
						if ft1[(k1-2)] == "IN":
							#g.write("{} IN {}\n".format(ft1[0],t2))
							d["{}".format(ft1[0].lower())] = "BOUND"
							count = count + 1
							k1 = len(ft1)
					else:
						#g.write("{} {}\n".format(ft1[0],t2))
						d["{}".format(ft1[0].lower())] = "BOUND"
						count = count + 1
						k1 = len(ft1)
				k2 = k2 + 1
			k1 = k1 + 1
		
		if count == 0:
			if len(ft1[0]) == 4:
				d["{}".format(ft1[0].lower())] = "UNBOUND"

		k = k + 1

	#g.close()
	return(d)

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

def seqres_atom_map(mmcif_dict):

	category = "_pdbx_poly_seq_scheme"
	seq_len = len(mmcif_dict[category + ".seq_id"])
	seqres = {}
	for i in range(seq_len):
		seqres_index = mmcif_dict["_pdbx_poly_seq_scheme.seq_id"][i]
		pdb_seq_id = int(mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"][i])
		chain = mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"][i]
		key1 = (seqres_index,chain)

		seqres[key1] = pdb_seq_id
	return seqres

def mut_dat(file1):

	import sys

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close() 

	k = 0
	ft1 = ft[k].split(",")
	t1 = ft1[0].strip("\n")
	while t1 != sys.argv[1]:
		k = k + 1
		if k > len(ft):
			print("{} NOT FOUND".format(sys.argv[1]))
			quit()
		ft1 = ft[k].split(",")
		t1 = ft1[0].strip("\n")
		
	d = dict()
	k1 = 0
	wt = ft1[k1].strip("\n")
	k1 = k1 + 1
	while k1 < len(ft1):
		t2 = ft1[k1].strip("\n")
		if wt not in d.keys():
			d["{}".format(wt)] = [t2]
		else:
			d["{}".format(wt)].append(t2)
		#if mut not in d.keys():
			#d["{}".format(mut)] = 1
		k1 = k1 + 1

	print("FOUND {} SIMILAR STRUCTURE FOR {}".format(len(ft1),wt))
	return(d)

def main_func():
	

	# CRITERIA
	r = 2.5

	sd = mut_dat("{}/seq_sim.txt".format(pathfiles))

	resolution = res_filter()

	HKW = HETATM_KW()

	import sys
	import re
	import gzip
	from Bio.PDB.MMCIFParser import MMCIFParser
	parser = MMCIFParser(QUIET=True)
	from Bio.PDB.Polypeptide import three_to_one as tto

	# ZONE VARIABLE
	zdis = 10.0

	# PATH FOR THE PDB/mmCIF FILES

	#pathmmcif = "/Users/tarun/Documents/mmCIF"
	pathmmcif = "/Volumes/BIOINFO/mmCIF"

	# LOCAL RMSD OF THE POLYPEPTIDE PRESENT IN DICTIONARY sd

	# DETERMINING THE ZONE AROUND THE MUTANT SITE TO DETERMINE ANY STRUCTURAL CHANGE TAKING PLACE DUE TO THE MUTATION

	z = open("{}".format(sys.argv[3]),"w")
	#temp = open("COM.txt","w")

	count2 = 0
	for x in sd.keys():
		pdbid = "{}".format(x)
		pdb = pdbid[0:4]
		cw = pdbid[5:6]

		typ = HKW["{}".format(pdb)]
		

		# CHECKING IF THE PDB's ARE ALIGNED
		count1 = 0
		if count1 == 0:
		#try:
			fol = pdb[1:3]		
			pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess1.cif","wb")
			out.write(tar.read())
			tar.close()
			out.close()

			mmcif = MMCIF2Dict("pdbprocess1.cif")
			idmap1 = seqres_atom_map(mmcif)

			kk = 0
			while kk < len(sd["{}".format(x)]):
				pdbid1 = "{}".format(sd["{}".format(x)][kk])
				#print(pdbid1)
				pdb1 = pdbid1[0:4]
				cm = pdbid1[5:6]

				typ1 = HKW["{}".format(pdb1)]
				rescr = resolution[0]["{}".format(pdb1)]

				if rescr == "None":
					rescr = 50.0

				if typ1 == typ and float(rescr) < r:

					fol = pdb1[1:3]		
					pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb1)
					tar = gzip.open("{}".format(pdbfile),"rb")
					out = open("pdbprocess2.cif","wb")
					out.write(tar.read())
					tar.close()
					out.close()

					mmcif = MMCIF2Dict("pdbprocess2.cif")
					idmap2 = seqres_atom_map(mmcif)

					count = 0
					for i in idmap1.keys():
						if i[1] == cw:
							for m in idmap2.keys():
								if m[1] == cm and i[0] == m[0]:
									if idmap2[m] != idmap1[i]:
										count = count + 1

					if count == 0:
						if count2 == 0:
							count2 = count2 + 1
							#print("*** {} AND {} :: {} of {} ***" .format(pdbid,pdbid1,kk,len(sd["{}".format(x)])))

							# FOR WILDTYPE

							# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

							structure_id = "{}".format(pdb)
							filename = "pdbprocess1.cif"
							structure = parser.get_structure(structure_id,filename)

							model = structure[0]

							chain = model["{}".format(cw)]
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
							
							k1 = 0
							while k1 < len(resid_list):
								if resid_list[k1] == int(sys.argv[2]):
									v1 = com_list[k1]
									zres = "{}".format(resid_list[k1])
									zresname = "{}".format(resname_list[k1])
									k2 = 0
									while k2 < len(resid_list):
										if k2 != k1:
											v2 = com_list[k2]
											dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
											if dis < (zdis*zdis):
												zres = zres + ",{}".format(resid_list[k2])
												zresname = zresname +  "{}".format(resname_list[k2])
										k2 = k2 + 1
									k1 = len(resid_list)
								k1 = k1 + 1

						z.write("{} {} {} {} {}\n".format(pdb,cw,pdb1,cm,zres))
					else:
						print("WT AND MUT MISMATCH")
					print("{} of {} ;; {} {}".format(kk,len(sd["{}".format(x)]),pdbid,pdbid1))
				kk = kk + 1

		#except:
			#print("FILE NOT FOUND")

	z.close()


main_func()




			
			
