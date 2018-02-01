# THIS CODE DIVIDES THE 100 SIMILAR STRUCTURES INTO THREE SETS BASED ON VARIOUS CRYTALLOGRAPY PARAMETERS WITH THE AIM TO QUANTIFY THE CRYTALLOGRAPHY NOISE

def set1():
	
	# SETS OF PDB CHAINS FROM THE SAME PDB

	import sys

	f = open("{}".format(sys.argv[1]),"r")	# INPUT WHICH CONTAINS ALL THE PAIRS
	ft = f.readlines()
	f.close()

	g = open("{}".format(sys.argv[2]),"w")

	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")
		t4 = ft1[3].strip("\n")

		if t1 == t3:
			g.write("{} {} {} {}".format(t1,t2,t3,t4))
			g.write("\n")
		
		k = k + 1

	g.close()

def set2():
	
	# SETS OF PDB CHAINS FROM DIFFERENT PDB's BUT SAME AUTHORS

	import sys

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	h = open("author.idx","r")		# CONTAINING THE AUTHOR LIST
	ht = h.readlines()
	h.close()

	# FORMING A DICTIONARY OF AUTHORS

	d = dict()
	k = 5
	while k < (len(ht)-1):
		ht1 = ht[k].split(" ; ")
		t1 = ht1[1].split(",")
		aut = t1[0]
		pdb = ht1[0]

		d["{}".format(pdb)] = aut
		k = k + 1

	g = open("{}".format(sys.argv[2]),"w")

	k = 0
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")
		t4 = ft1[3].strip("\n")

		t1u = t1.upper()
		t3u = t3.upper()

		try:
			if d["{}".format(t1u)] == d["{}".format(t3u)] and t1u != t3u:
				g.write("{} {} {} {}".format(t1,t2,t3,t4))
				g.write("\n")
		except:
			print("DICT_ERROR")

		k = k + 1

	g.close()
	
set2()
