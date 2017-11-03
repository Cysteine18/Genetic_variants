f = open("number_of_distinct_mutants.txt","r")
ft = f.readlines()
f.close()

g = open("pdb_seqres.txt","r")
gt = g.readlines()
g.close()

h = open("cluster_length.txt","w")

k = 0
while k < len(ft):
	s1 = "{}".format(ft[k])
	ft1 = s1.split()
	t1 = ft1[0]
	print(t1)
	k1 = 0
	gt1 = gt[k1].split()
	t2 = gt1[0].strip(">")
	while t2 != t1:
		k1 = k1 + 2
		gt1 = gt[k1].split()
		t2 = gt1[0].strip(">")
	
	t1len = len(gt[k1+1])
	t1len = t1len - 1

	h.write("{} {} {}".format(t1,len(ft1),t1len))
	h.write("\n")

	k = k + 1

h.close()
