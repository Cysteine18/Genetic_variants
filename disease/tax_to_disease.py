# MAPPING HUMAN MUTANTS TO DISEASES

f = open("disease_causing_mutations.txt","r")
ft = f.readlines()
f.close()

g = open("unique_humans_mutants.txt","r")
gt = g.readlines()
g.close()

h = open("disease_causing_pdbIDs.txt","w")

k = 0
nmut = 0
while k < len(gt):
	t1 = gt[k].strip("\n")
	print ("{} :: {} of {}".format(t1,k,len(gt)))
	k1 = 0
	count = 0
	while k1 < len(ft):
		ft1 = ft[k1].split(",")
		ft2 = ft1[1].split(";")
		k2 = 0
		while k2 < len(ft2):
			t2 = ft2[k2]
			if t1 == t2:
				nmut = nmut + 1
				count = count + 1
				k2 = len(ft2)
				k1 = len(ft)
			k2 = k2 + 1
		k1 = k1 + 1
	if count != 0:
		if nmut > 1:
			h.write("\n")
		h.write("{}".format(t1))
	k = k + 1

print ("	### NUMBER OF DISEASE CAUSING MUTANTS = {}".format(nmut))
h.close()
