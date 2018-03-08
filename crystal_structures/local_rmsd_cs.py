from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import gzip

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

def PFAM(file1):

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close()

	# DICTIONARY OF ALL THE UNIQUE PDBS

	p = dict()
	k = 0
	while k < len(ft):
		#print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t1 = t1.lower()
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")
		t4 = ft1[3].strip("\n")
		t5 = ft1[4].strip("\n")
		key1 = "{}_{}".format(t1,t2)
		list1 = (t3,t4,t5)
		if key1 not in p.keys():
			p["{}".format(key1)] = [list1]
		else:
			p["{}".format(key1)].append(list1)
		k = k + 1

	return(p)

def cath(file1):

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close()

	# DICTIONARY OF ALL THE UNIQUE PDBS

	p = dict()
	k = 0
	while k < len(ft):
		#print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		temp = ft1[0].strip("\n")
		t1 = temp[0:4]
		t2 = temp[4:(len(temp)-2)]
		t3 = temp[len(temp):(len(temp)-2)]
		key1 = "{}_{}".format(t1,t2)
		list1 = t3
		if key1 not in p.keys():
			p["{}".format(key1)] = [list1]
		else:
			p["{}".format(key1)].append(list1)
		k = k + 1

	p1 = dict()
	for x in p.keys():
		t1 = len(p["{}".format(x)])
		p1["{}".format(x)] = t1

	return(p1)

def mut_dat(file1):

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close() 

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")
		t4 = ft1[3].strip("\n")
		wt = "{}_{}".format(t1,t2)
		mut = "{}_{}".format(t3,t4)
		if wt not in d.keys():
			d["{}".format(wt)] = [mut]
		else:
			d["{}".format(wt)].append(mut)
		#if mut not in d.keys():
			#d["{}".format(mut)] = 1
		k = k + 1

	return(d)

def main_func():
	
	pf = PFAM("pdb_pfam_mapping.txt")
	ca = cath("cath_domain.txt")
	md = mut_dat("set1_sc.txt")

	f = open("PFAM_dis.txt","w")

	g = open("CATH_dis.txt","w")

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
	#pathmmcif = "/Volumes/BIOINFO/mmCIF"
	pathmmcif = "/data/pdb/divided/mmCIF"

	d1 = dict()
	sd = dict()
	nf = 0
	for x in md:
		try:
			t1 = "{}".format(ca["{}".format(x)])
			if t1 not in d1.keys() and t1 == "1":
				d1["{}".format(t1)] = [x]
				if x not in sd.keys():
					sd["{}".format(x)] = md["{}".format(x)]
			elif t1 == "1":
				d1["{}".format(t1)].append(x)
				if x not in sd.keys():
					sd["{}".format(x)] = md["{}".format(x)]
		except:
			nf = nf + 1
			#print("{} NOT FOUND IN CATH".format(x))

	#print("{} OUT OF {} ARE NOT FOUND IN CATH FILE".format(nf,len(md)))

	k = 1
	for x in d1:
		t1 = len(d1["{}".format(x)])
		g.write("{} {} {}\n".format(k,x,t1))
		k = k + 1

	d = dict()
	nf = 0
	for x in sd:
		try:
			k1 = 0
			while k1 < len(pf["{}".format(x)]):
				#print(pf["{}".format(x))
				t1 = pf["{}".format(x)][k1][2]
				if t1 not in d.keys():
					d["{}".format(t1)] = 	[x]
				else:
					d["{}".format(t1)].append(x)
				k1 = k1 + 1
		except:
			nf = nf + 1
			#print("{} NOT FOUND IN PFAM".format(x))

	#print("{} OUT OF {} ARE NOT FOUND IN PFAM FILE".format(nf,len(md)))
	
	k = 1
	for x in d:
		t1 = len(d["{}".format(x)])
		f.write("{} {} {}\n".format(k,x,t1))
		k = k + 1

	f.close()
	g.close()

	# LOCAL RMSD OF THE POLYPEPTIDE PRESENT IN DICTIONARY sd

	# DETERMINING THE ZONE AROUND THE MUTANT SITE TO DETERMINE ANY STRUCTURAL CHANGE TAKING PLACE DUE TO THE MUTATION

	z = open("{}".format(sys.argv[3]),"w")
	#temp = open("COM.txt","w")

	start = sys.argv[1]
	end = sys.argv[2]
	if end == "END" or end == "end":
		end = len(sd)
	end = int(end)
	start =int(start)
	k = 0
	for x in sd.keys():
		print("# {} of {} #".format(k,end))
		if k >= start:
			pdbid = "{}".format(x)
			pdb = pdbid[0:4]
			cw = pdbid[5:6]

			# CHECKING IF THE PDB's ARE ALIGNED
			#count1 = 0
			#if count1 == 0:
			try:
				fol = pdb[1:3]		
				pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
				tar = gzip.open("{}".format(pdbfile),"rb")
				out = open("pdbprocess1{}.cif".format(start),"wb")
				out.write(tar.read())
				tar.close()
				out.close()

				mmcif = MMCIF2Dict("pdbprocess1{}.cif".format(start))
				idmap1 = seqres_atom_map(mmcif)

				kk = 0
				while kk < len(sd["{}".format(x)]):
					pdbid1 = "{}".format(sd["{}".format(x)][kk])
					#print(pdbid1)
					pdb1 = pdbid1[0:4]
					cm = pdbid1[5:6]

					fol = pdb1[1:3]		
					pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb1)
					tar = gzip.open("{}".format(pdbfile),"rb")
					out = open("pdbprocess2{}.cif".format(start),"wb")
					out.write(tar.read())
					tar.close()
					out.close()

					mmcif = MMCIF2Dict("pdbprocess2{}.cif".format(start))
					idmap2 = seqres_atom_map(mmcif)

					count = 0
					for i in idmap1.keys():
						if i[1] == cw:
							for m in idmap2.keys():
								if m[1] == cm and i[0] == m[0]:
									if idmap2[m] != idmap1[i]:
										count = count + 1

					count2 = 0
					if count == 0:

						if count2 == 0:
							count2 = count2 + 1
							#print("*** {} AND {} :: {} of {} ***" .format(pdbid,pdbid1,kk,len(sd["{}".format(x)])))

							# FOR WILDTYPE

							# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

							structure_id = "{}".format(pdb)
							filename = "pdbprocess1{}.cif".format(start)
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
						
							k1 = 20
							ss = int((len(resid_list) - 40) / 5)
							if ss < 6:
								ss = 6
							while k1 < (len(resid_list) -20):
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

								k1 = k1 + ss

								z.write("{} {} {} {} {}\n".format(pdb,cw,pdb1,cm,zres))

								#z.write("{}".format(pdbid))
								#z.write("\n")
								#z.write("{}".format(zres))
								#z.write("\n")
								#z.write("{}".format(zresname))
								#z.write("\n")



						# FOR MUTANT

						#structure_id = "{}".format(pdb1)
						#filename = "pdbprocess2.cif"
						#structure = parser.get_structure(structure_id,filename)

						#model = structure[0]

						#chain = model["{}".format(cm)]
						#c1 = chain.get_list()		# LIST ALL THE RESIDUES

						#k1 = 0
						#resid_list = []
						#resname_list = []
						#com_list = []
						#while k1 < len(c1):
							#c2 = c1[k1].get_id()
							#resid = c2[1]
							#if c2[0] == " ":
								#residue = chain[c2]
								#tresname = residue.get_resname()
								#try:
									#resname = tto("{}".format(tresname))
								#except:
									#resname = "X"
								#r1 = residue.get_list() # LIST ALL THE ATOMS OF A PARTICULAR RESIDUE

								#k2 = 0
								#res = []	# ATOM NAMES
								#cd = []		# ATOM COORDINATES
								#while k2 < len(r1):
									#r2 = r1[k2].get_id()
									# COM OF THE BACKBONE
									#if r2 == "CA" or r2 == "N" or r2 == "C" or r2 == "O":
										#res.append(r2)

										#atom = residue['{}'.format(r2)]
										#a1 = atom.get_coord()
										#cd.append(a1)
									#k2 = k2 + 1
							
								#CM = COM(res,cd)
								#resid_list.append(resid)
								#resname_list.append(resname)
								#com_list.append(CM)
						
							#k1 = k1 + 1
						
						#k1 = 20
						#ss = int((len(resid_list) - 40) / 5)
						#if ss < 6:
							#ss = 6
						#while k1 < (len(resid_list) -20):
							#v1 = com_list[k1]
							#zres = "{}".format(resid_list[k1])
							#zresname = "{}".format(resname_list[k1])
							#k2 = 0
							#while k2 < len(resid_list):
								#if k2 != k1:
									#v2 = com_list[k2]
									#dis = (v1[0]-v2[0])**2 + (v1[1]-v2[1])**2 + (v1[2]-v2[2])**2
									#if dis < (zdis*zdis):
										#zres = zres + ",{}".format(resid_list[k2])
										#zresname = zresname +  "{}".format(resname_list[k2])
								#k2 = k2 + 1

							#z.write("{}".format(pdbid1))
							#z.write("\n")
							#z.write("{}".format(zres))
							#z.write("\n")
							#z.write("{}".format(zresname))
							#z.write("\n")

							#k1 = k1 + ss

					else:
						print("WT AND MUT MISMATCH")
					kk = kk + 1

			except:
				print("FILE NOT FOUND")
				#z.write("{}".format(pdbid))
				#z.write("\n")
				#z.write("NA")
				#z.write("\n")
				#z.write("NA")
				#xz.write("\n")

		k = k + 1
		if k > end:
				break

	z.close()


main_func()




			
			
