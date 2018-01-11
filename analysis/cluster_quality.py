# THIS CODE LOOKS AT EACH PAIR OF WILDTYPE AND MUTANT TO DETERMINE IF THEY REPRESENT THE SAME UNIPROT ID

f = open("cluster.txt","r")
ft = f.readlines()
f.close()

lenfilter = 30

g = open("wrong_clusters","w")
g1 = open("wrong_mutants","w")
g2 = open("capri_fit_data.txt","w")

k = 0
cid = 1
while k < len(ft):
	ft1 = ft[k].split()	
	if ft1[0] == "#" and ft1[1] == "{}".format(cid):
		cl = ft1[4].strip("'|,")
		print("CHECKING CLUSTER {} OF SIZE {}".format(cid,cl))
		k = k + 1
		t1 = ft[k].split(", ")
		pdbwt = t1[0].strip("('|'")
		pdbwtn = pdbwt[0:4]
		chainwt = pdbwt[5:len(pdbwt)]
		orgwt = t1[1].strip("'")
		uid = t1[2].strip("[|'|]|)|\n")
		k = k + 1
		ft1 = ft[k].split()
		count = 0
		while ft1[0] != "#" and k < (len(ft) - 1):
			t1 = ft[k].split(", ")
			t2 = t1[2].strip("[|'|]|)|\n")
			pdbmut = t1[0].strip("('|'")
			pdbmutn = pdbmut[0:4]
			chainmut = pdbmut[5:len(pdbmut)]
			orgmut = t1[1].strip("'")
			if t2 != uid and int(cl) > lenfilter:
				count = count + 1
				g1.write("{}	{}".format(pdbwt,pdbmut))
				g1.write("\n")
				g2.write("{} {} {} {}".format(pdbwtn,chainwt,pdbmutn,chainmut))
				g2.write("\n")
			k = k + 1
			t1 = ft[k].split(", ")
			uid = t1[2].strip("[|'|]|)|\n")
			ft1 = ft[k].split()		

		if count != 0 and int(cl) > lenfilter:
			g.write("{}".format(cid))
			g.write("\n")
		cid = int(cid) + 1
	else:
		k = k + 1

g.close()
g1.close()
g2.close()

