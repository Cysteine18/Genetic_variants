import sys

f = open("{}".format(sys.argv[1]),"r")
ft = f.readlines()
f.close()

g = open("cluster_high_res.txt","r")
gt = g.readlines()
g.close()

h = open("human_variants_tarun.csv","w")
h.write("#WT,WT_CHAIN,MUT,MUT_CHAIN,UNIPROT,WT_RES,MUT_RES,MUT_POS,MUT_TYPE")
h.write("\n")

h1 = open("represented_PDB.txt","w")

ft1 = ft[0].split(", ")
k = 0
k1 = 0
while k < len(ft1):
	print("{} OF {}".format(k,len(ft1)))
	t1 = ft1[k].strip("['|'|']|\n")
	#print(t1)
	gt1 = gt[k1].split()
	t2 = gt1[0].strip("\n")
	try:
		t3 = gt1[1].strip("\n")
	except:
		t3 = "NF"
	while t2 != "#" or t3 != t1:
		k1 = k1 + 1
		gt1 = gt[k1].split()
		t2 = gt1[0].strip("\n")
		try:
			t3 = gt1[1].strip("\n")
		except:
			t3 = "NF"

	k1 = k1 + 1
	gt1 = gt[k1].split()
	t1 = gt1[0]
	wt = gt1[0].strip("('|',")
	wtpdb = wt[0:4]
	wtchain = wt[5:len(wt)]
	h1.write("{} {}".format(wtpdb,wtchain))
	h1.write("\n")
	uwt = gt1[2].strip("['|']|)|\n")
	k1 = k1 + 1
	gt1 = gt[k1].split()
	t1 = gt1[0]
	while t1 != "#":
		mut = gt1[0].strip("('|',")
		mutpdb = mut[0:4]
		mutchain = mut[5:len(mut)]
		umut = gt1[2].strip("['|'],|)|\n")
		count = 0
		if umut == uwt:
			check = gt1[3].strip("'|',|\n")
			if check !=  "NOT_FOUND" and check !=  "NA":
				reswt = gt1[3][1]
			else:
				reswt = "NR"
				count = count + 1
			check = gt1[4].strip("'|',|\n")
			if check !=  "NOT_FOUND" and check !=  "NA":
				resmut = gt1[4][1]
			else:
				resmut = "NR"
				count = count + 1
			check = gt1[5].strip("'|',|\n")
			if check !=  "NOT_FOUND" and check !=  "NA":
				pos = gt1[5].split(",")
				pos1 = pos[0].strip("'")
			else:
				#check = gt1[6].strip("'|',|\n")
				pos1 = "NR"
				if check !=  "NOT_FOUND" and check !=  "NA":
					pos = gt1[6].split(",")
					pos1 = pos[0].strip("')")
				else:
					pos1 = "NR"
					count = count + 1
	
			if count == 0:
				h.write("{},{},{},{},{},{},{},{}".format(wtpdb,wtchain,mutpdb,mutchain,uwt,reswt,resmut,pos1))
				h.write("\n")
			
		k1 = k1 + 1
		gt1 = gt[k1].split()
		t1 = gt1[0]

	k = k + 1

h.close()
h1.close()













		
