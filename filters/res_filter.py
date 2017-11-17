import sys

f = open("seq_sim.txt","r")
ft = f.readlines()
f.close()

g = open("PDB_classifier.txt","r")
gt = g.readlines()
g.close()

ne = 1
while ne < len(sys.argv):
	pdb = sys.argv[ne]
	k = 0
	ft1 = ft[k].split(",")
	t1 = ft1[0].strip("\n")

	while t1 != pdb:
		k = k + 1
		if k > len(ft):
			print("{} NOT FOUND".format(pdb))
			quit()
		ft1 = ft[k].split(",")
		t1 = ft1[0].strip("\n")	

	k1 = 0
	minr = 100.0
	while k1 < len(ft1):
		t1c = ft1[k1].strip("\n")
		t1 = t1c[0:4]
		k2 = 0

		gt1 = gt[k2].split()
		t2 = gt1[0]
		count = 0
		while t2 != t1:
			k2 = k2 + 1
			if k2 > len(gt):
				print("{} NOT FOUND IN CLASSIFIER".format(t1))
				count = count + 1
			gt1 = gt[k2].split()
			t2 = gt1[0]
			t3 = gt1[1]
		if count == 0:
			if t3 != "None":
				if float(t3) < minr:
					pdbc = t1c
	
		k1 = k1 + 1

	print("FOR PDBID {} THE PDBID {} IS THE LOWEST IN RESOLUTION".format(pdb,pdbc))

	ne = ne + 1	
