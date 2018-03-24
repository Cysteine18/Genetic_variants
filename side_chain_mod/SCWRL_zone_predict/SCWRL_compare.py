import sys

f = open("{}".format(sys.argv[1]),"r")		# RMSD SCWRL
ft = f.readlines()
f.close()

g = open("{}".format(sys.argv[2]),"r")		# RMSD CRYSTAL STRUCTURE
gt = g.readlines()
g.close()

h = open("SCWRL_compare.csv","w")
h.write("#CRY_MUT,CHAIN,PRED_MUT,CHAIN,LRMSD_SC(0),LRMSD_SC(4),LRMSD_SC(5),LRMSD_SC(6), LRMSD_SC(7),LRMSD_SC(8),LRMSD_SC(9),LRMSD_SC(10),LRMSD_SC(crystal),LRMSD\n")

k = 1
while k < len(ft):
	print("{} of {}".format(k,len(ft)))
	ft1 = ft[k].split(",")
	wt = ft1[0]
	wtc = ft1[1]
	mut = ft1[2]
	mutc = ft1[3]
	lenwt = ft1[4]
	lenmut = ft1[5]
	allen = ft1[6]

	temp = mut.split("-")
	wt1 = temp[0][0:4]
	wtc1 = temp[0][4:(len(temp[0]))]
	mut1 = temp[1][0:4]
	mutc1 = temp[1][4:(len(temp[0]))]
	k1 = 0
	gt1 = gt[k1].split(",")
	t2 = gt1[0]
	t3 = gt1[1]
	t4 = gt1[2]
	t5 = gt1[3]

	if lenwt == allen and lenmut == allen:
		count = 0
		while t2 != wt1 or t3 != wtc1 or t4 != mut1 or t5 != mutc1:
			k1 = k1 + 1
			if k1 >= len(gt):
				count = count + 1
				break
			gt1 = gt[k1].split(",")
			t2 = gt1[0]
			t3 = gt1[1]
			t4 = gt1[2]
			t5 = gt1[3]

		if count == 0:
			lrmsdsc = gt1[9].strip("\n")
			lrmsd = gt1[8]

			if lrmsdsc != "NA":
				list1 = []
				for x in range(0,8):
					ft1 = ft[(k+x)].split(",")
					t1 = ft1[9]
					if t1 != "NA":
						list1.append(t1)
				if len(list1) == 8:
					h.write("{},{},{}-{},{},{},{},{},{},{},{},{},{},{},{}\n".format(wt,wtc,wt1,mut1,mutc,list1[0],list1[1],list1[2],list1[3],list1[4],list1[5], list1[6],list1[7],lrmsdsc,lrmsd))
	
	k = k + 8

h.close()

