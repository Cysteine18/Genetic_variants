# THIS CODE DIVIDES THE 100 SIMILAR STRUCTURES INTO THREE SETS BASED ON VARIOUS CRYTALLOGRAPY PARAMETERS WITH THE AIM TO QUANTIFY THE CRYTALLOGRAPHY NOISE

def set1():
	
	# SETS OF PDB CHAINS FROM THE SAME PDB

	import sys

	f = open("{}".format(sys.argv[1]),"r")
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

set1()
