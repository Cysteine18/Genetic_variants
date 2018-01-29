# ASSIGNMENT OF EACH MUTATIONS AS BURIED OR EXPOSED

def mutations():
	import sys

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("mut_data.txt","w")

	k = 0
	d = dict()
	nmut = 1
	cid = 1
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0]
		try:
			t2 = ft1[1]
		except:
			t2 = "NF"
		if t1 == "#" and t2 == "{}".format(cid):
			cid = cid + 1
			k = k + 1
			ft1 = ft[k].split()
			wt = ft1[0].strip("('|',|\n")
			k = k + 1
			ft1 = ft[k].split()
			mut = ft1[0].strip("(,|',|\n")
			while mut != "#" and k < (len(ft)-1):
				pos = ft1[(len(ft1)-1)].strip("'|')")
				if pos != "NOT_FOUND":
					pos1 = pos.split(",")
					mut_pos = pos1[0]
					d[nmut] = (wt,mut,mut_pos)
					nmut = nmut + 1
				k = k + 1
				ft1 = ft[k].split()
				mut = ft1[0].strip("(,|',|\n")
		else:
			k = k + 1

	list1 = (d,nmut)
	for x in range(1,nmut):
		t1 = d[x][0]
		t2 = d[x][1]
		t3 = d[x][2]
		g.write("{} {} {}".format(t1,t2,t3))
		g.write("\n")
	g.close()
	return(list1)

def DSSP_ass():

	path = "/Volumes/RCSB_DATA/pdb/"

	import subprocess
	import gzip
	import os
	import time
	# THIS FUNTION DIVIDES ALL MUTATIONS INTO BURIED AND EXPOSED
	# OUTPUT FORMAT :: wt (b/i/e) mut (b/i/e) if change '1' otherwise '0'

	x = mutations()
	nmut = x[1]
	mut = x[0]
	
	f = open("solvent_ass_3.csv","w")
	f.write("#wt,mut,chwt,chmut,reswt,resmut,pos,wt_acc,mut_acc")
	f.write("\n")

	for y in range(4000,nmut):
		#print(mut[y])
		pdb1 = mut[y][0][0:4]
		c1 = mut[y][0][5:(len(mut[y][0]))]
		fol1 = mut[y][0][1:3]
		pdb2 = mut[y][1][0:4]
		c2 = mut[y][1][5:(len(mut[y][1]))]
		fol2 = mut[y][1][1:3]
		pos = mut[y][2]

		print("{} of {}".format(y,nmut))

		try:

		#count = 0
		#if count == 0:

			# WILDTYPE

			pathPDB = "{}/{}".format(path,fol1)
			pdbfile = "{}/pdb{}.ent.gz".format(pathPDB,pdb1)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess.pdb","wb")
			out.write(tar.read())
			tar.close()
			out.close()


			subprocess.Popen(['mkdssp', '-i','pdbprocess.pdb', '-o','output{}.dssp'.format(pdb1)])
			time.sleep(0.5)


			#subprocess.Popen(['mkdssp', '-i','pdbprocess.pdb', '-o','output{}.dssp'.format(pdb1)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
			#os.rename("pdbprocess.pdb","{}.pdb".format(pdb1))
			#os.remove("pdbprocess.pdb")

			g = open("output{}.dssp".format(pdb1),"r")
			gt = g.readlines()
			g.close()

			#os.remove("output{}.dssp".format(pdb1))

			k = 0
			gt1 = gt[k].split()
			t1 = gt1[0]
			while t1 != "#":
				k = k + 1
				gt1 = gt[k].split()
				t1 = gt1[0]
			k = k + 1
			gt1 = gt[k].split()
			t2 = gt1[1]
			chwt = gt1[2]
			while t2 != pos or chwt != c1:
				k = k + 1
				gt1 = gt[k].split()
				t2 = gt1[1]
				chwt = gt1[2]
			reswt = gt1[3]
			gt2 = gt[k].split(",")
			gt3 = gt2[0].split()
			t3 = gt3[(len(gt3)-2)]
			wt = pdb1
			wt_acc = t3
		

			# MUTANT

			pathPDB = "{}/{}".format(path,fol2)
			pdbfile = "{}/pdb{}.ent.gz".format(pathPDB,pdb2)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess.pdb","wb")
			out.write(tar.read())
			tar.close()
			out.close()

			subprocess.Popen(['mkdssp', '-i','pdbprocess.pdb', '-o','output{}.dssp'.format(pdb2)])
			time.sleep(0.5)

			#subprocess.Popen(['mkdssp', '-i','pdbprocess.pdb', '-o','output{}.dssp'.format(pdb2)], stdout=open(os.devnull, 'wb'), stderr=open(os.devnull, 'wb'))
			#os.rename("pdbprocess.pdb","{}.pdb".format(pdb2))
			#os.remove("pdbprocess.pdb")

			g = open("output{}.dssp".format(pdb2),"r")
			gt = g.readlines()
			g.close()

			#os.remove("output{}.dssp".format(pdb2))

			k = 0
			gt1 = gt[k].split()
			t1 = gt1[0]
			while t1 != "#":
				k = k + 1
				gt1 = gt[k].split()
				t1 = gt1[0]
			k = k + 1
			gt1 = gt[k].split()
			t2 = gt1[1]
			chmut = gt1[2]
			while t2 != pos or chmut != c2:
				k = k + 1
				gt1 = gt[k].split()
				t2 = gt1[1]
				chmut = gt1[2]
			resmut = gt1[3]
			gt2 = gt[k].split(",")
			gt3 = gt2[0].split()
			t3 = gt3[(len(gt3)-2)]
			mut1 = pdb2
			mut_acc = t3
			f.write("{},{},{},{},{},{},{},{},{}".format(pdb1,c1,pdb2,c2,reswt,resmut,pos,wt_acc,mut_acc))
			f.write("\n")
				
		except:
			print("#### ERROR ####")
			f.write("{},{},{},{},ERROR,ERROR,{},ERROR,ERROR".format(pdb1,c1,pdb2,c2,pos))
			f.write("\n")

	f.close()

def mut_prop():

	# THIS CODE ROSS AND SANDER DEFINATION OF MAXIMUM ACCESSIBILITY TO CALCULATE THE BURIED AND EXPOSED RESIDUES

	list1 = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','a']
	list2 = [106,160,135,163,194,197,84,184,169,205,164,188,157,136,198,248,130,142,142,227,180,222,196,135]

	acc = dict()
	for x in range(0,len(list1)):
		acc["{}".format(list1[x])] = list2[x]


	# BINARY MODEL 'B' 'E'
	modB1 = 0.16
	BN = 0
	BY = 0
	BB = 0
	BE = 0
	
	# TERNARY MODEL 'B' 'I' 'E'
	modT1 = 0.09
	modT2 = 0.36
	TN = 0
	TY = 0
	TB = 0
	TI = 0
	TE = 0

	f = open("solvent_ass.csv","r")
	ft = f.readlines()
	f.close()

	g = open("buried_exposed.csv","w")
	g.write("#wt,mut,chwt,chmut,reswt,resmut,pos,wt_acc,mut_acc,B/E_WT,B/E_MUT,change,B/I/E_WT,B/I/E_MUT,change")
	g.write("\n")

	k = 1
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split(",")
		check = ft1[4].strip("\n")
		if check != "ERROR":
			count = 0
			t1 = ft1[0].strip("\n")
			t2 = ft1[1].strip("\n")
			t3 = ft1[2].strip("\n")
			t4 = ft1[3].strip("\n")
			t5 = ft1[4].strip("\n")
			t6 = ft1[5].strip("\n")
			t7 = ft1[6].strip("\n")
			t8 = ft1[7].strip("\n")
			t9 = ft1[8].strip("\n")

			try:
				wtra = int(t8) / acc["{}".format(t5)]	# WILDTYPE
				mutra = int(t9) / acc["{}".format(t6)]	# MUTANT
			except:
				count = count + 1

			if count == 0:
			
				# BINARY ASSIGNMENT

				if wtra < modB1:
					t10 = "B"
					BB = BB + 1
				else:
					t10 = "E"
					BE = BE + 1
			
				if mutra < modB1:
					t11 = "B"
					BB = BB + 1
				else:
					t11 = "E"
					BE = BE + 1

				if t10 == t11:
					t12 = "NO"
					BN = BN + 1
				else:
					t12 = "YES"
					BY = BY + 1

				# TERNARY ASSIIGNMENT

				if wtra < modT1:
					t13 = "B"
					TB = TB + 1
				elif wtra > modT1 and wtra < modT2:
					t13 = "I"
					TI = TI + 1
				elif wtra > modT2:
					t13 = "E"
					TE = TE + 1

				if mutra < modT1:
					t14 = "B"
					TB = TB + 1
				elif mutra > modT1 and mutra < modT2:
					t14 = "I"
					TI = TI + 1
				elif mutra > modT2:
					t14 = "E"
					TE = TE + 1

				if t13 == t14:
					t15 = "NO"
					TN = TN + 1
				else:
					t15 = "YES"
					TY = TY + 1

				g.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15))
				g.write("\n")
			else:
				g.write("{},{},{},{},{},{},{},{},{},ERROR,ERROR,ERROR,ERROR,ERROR,ERROR".format(t1,t2,t3,t4,t5,t6,t7,t8,t9))
				g.write("\n")
		else:
				t1 = ft1[0].strip("\n")
				t2 = ft1[1].strip("\n")
				t3 = ft1[2].strip("\n")
				t4 = ft1[3].strip("\n")
				g.write("{},{},{},{},ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR".format(t1,t2,t3,t4))
				g.write("\n")

		k = k + 1

	g.close()
	h = open("buried_exposed_stats","w")
	h.write("# BINARY ASSIGNMENT\n")
	h.write("BURIED :: {}\n".format(BB))
	h.write("EXPOSED :: {}\n".format(BE))
	h.write("CHANGE :: {}\n".format(BY))
	h.write("NO_CHANGE :: {}\n".format(BN))
	nf = len(ft) - int(BY) - int(BN)
	h.write("CAN'T ASSIGN :: {}\n".format(nf))
	h.write("\n")

	h.write("# TERNARY ASSIGNMENT\n")
	h.write("BURIED :: {}\n".format(TB))
	h.write("INTERMEDIATE :: {}\n".format(TI))
	h.write("EXPOSED :: {}\n".format(TE))
	h.write("CHANGE :: {}\n".format(TY))
	h.write("NO_CHANGE :: {}\n".format(TN))
	nf = len(ft) - int(TY) - int(TN)
	h.write("CAN'T ASSIGN :: {}\n".format(nf))

	h.close()
	

mut_prop()













				
			


