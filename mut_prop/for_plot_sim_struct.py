import sys

f = open("{}".format(sys.argv[1]),"r")		# MUTATIONS
ft = f.readlines()
f.close()

g = open("{}".format(sys.argv[2]),"r")		# ENSEMBLES
gt = g.readlines()
g.close()

h = open("ensemble.txt","w")
h1 = open("ensemble.csv","w")
h.write("#WT	MUT	GRMSD	GRMSD_ENSEMBLE	GRMSD_ENSEMBLE_SD	LRMSD	LRMSD_ENSEMBLE	LRMSD_ENSEMBLE_SD	LRMSDSC	LRMSDSC_ENSEMBLE	LRMSDSC_ENSEMBLE_SD\n")
h1.write("#WT,MUT,GRMSD,GRMSD_ENSEMBLE	,GRMSD_ENSEMBLE_SD,LRMSD,LRMSD_ENSEMBLE,LRMSD_ENSEMBLE_SD,LRMSDSC,LRMSDSC_ENSEMBLE,LRMSDSC_ENSEMBLE_SD\n")

k = 1
while k < len(ft):
	ft1 = ft[k].split()
	t1 = ft1[2]
	t2 = ft1[3]
	t3 = ft1[4]
	t4 = ft1[5]

	t1c = "{}_{}".format(t1,t2)
	t2c = "{}_{}".format(t3,t4)

	grmsd = ft1[8]
	lrmsd = ft1[9]
	lrmsdsc = ft1[10].strip("\n")
	
	if grmsd == "ERROR" or grmsd == "NA":
		grmsd = 0.0

	if lrmsd == "ERROR" or lrmsd == "NA":
		lrmsd = 0.0

	if lrmsdsc == "ERROR" or lrmsdsc == "NA":
		lrmsdsc = 0.0
	
	k1 = 1
	gt1 = gt[k1].split()
	temp1 = gt1[0]
	temp2 = gt1[1]
	count = 0
	while temp1 != t1c or temp2 != t2c:
		k1 = k1 + 1
		if k1 >= len(gt):
			break
			count = count + 1
		gt1 = gt[k1].split()
		temp1 = gt1[0]
		temp2 = gt1[1]

	if count == 0:
		nwt = int(gt1[2])
		nmut = int(gt1[3])
		print(t1c,t2c)
		if nwt > 0:
			genwt = float(gt1[4])
			genwtsd = float(gt1[5])
			try:
				lenwt = float(gt1[8])
				lenwtsd = float(gt1[9])
			except:
				lenwt = 0.0
				lenwtsd = 0.0
			try:
				lenscwt = float(gt1[12])
				lenscwtsd = float(gt1[13])
			except:
				lenscwt = 0.0
				lenscwtsd = 0.0
		else:
			genwt = 0.0
			lenwt = 0.0
			lenscwt = 0.0
		if nmut > 0:
			genmut = float(gt1[6])
			genmutsd = float(gt1[7])
			try:
				lenmut = float(gt1[10])
				lenmutsd = float(gt1[11])
			except:
				lenmut = 0.0
				lenmutsd = 0.0
			try:
				lenscmut = float(gt1[14])
				lenscmutsd = float(gt1[15])
			except:
				lenscmut = 0.0
				lenscmutsd = 0.0
		else:
			genmut = 0.0
			lenmut = 0.0
			lenscmut = 0.0
	
		if nwt != 0 or nmut != 0:
			if genwt >= genmut:
				mgen = genwt
				mgensd = genwtsd
			else:
				mgen = genmut
				mgensd = genmutsd

			if lenwt >= lenmut:
				mlen = lenwt
				mlensd = lenwtsd
			else:
				mlen = lenmut
				mlensd = lenmutsd

			if lenscwt >= lenscmut:
				mlensc = lenscwt
				mlenscsd = lenscwtsd
			else:
				mlensc = lenscmut
				mlenscsd = lenscmutsd

			h.write("{}	{}	{}	{}	{}	{}	{}	{}	{}	{}	{}\n".format(t1c,t2c,grmsd,mgen,mgensd,lrmsd,mlen,mlensd,lrmsdsc,mlensc,mlenscsd))
			h1.write("{},{},{},{},{},{},{},{},{},{},{}\n".format(t1c,t2c,grmsd,mgen,mgensd,lrmsd,mlen,mlensd,lrmsdsc,mlensc,mlenscsd))

	k = k + 1

h.close()
h1.close()









		

	
