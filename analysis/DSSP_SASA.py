def mutations():
	import sys

	f = open("{}".format(sys.argv[1]),"r")		# CLUSTER FILE
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
				posw = ft1[(len(ft1)-2)].strip("'|')|',")
				resm = ft1[(len(ft1)-3)].strip("'|')|',")
				resw = ft1[(len(ft1)-4)].strip("'|')|',")
				if pos != "NOT_FOUND" and posw != "NOT_FOUND" and resm != "NOT_FOUND" and resw != "NOT_FOUND" and resm != "NA" and resw != "NA":
					pos1 = pos.split(",")
					mut_pos = pos1[0]
					mut_res = resm[0]
					wt_res = resw[0]
					d[nmut] = (wt,mut,mut_pos,wt_res,mut_res)
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
		t4 = d[x][3]
		t5 = d[x][4]
		g.write("{} {} {} {} {}".format(t1,t2,t3,t4,t5))
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
	
	f = open("solvent_ass.csv","w")
	f.write("#wt,mut,chwt,chmut,reswt,resmut,pos,wt_acc,mut_acc")
	f.write("\n")

	for y in range(1,nmut):
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

DSSP_ass()

