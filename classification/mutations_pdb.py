# THIS CODE CONVERTS THE MUTATIONS FROM SEQRES FILE TO CORRESPONDING PDB FILES

f = open("number_of_distinct_mutations.csv","r")
ft = f.readlines()
f.close()

g = open("distinct_mutants_pdb.txt","r")
gt = g.readlines()
g.close()

h = open("distinct_mutants_seqres.txt","r")
ht = h.readlines()
h.close()

# forming a dixtionary with pdb chain as key and mutation site in seqres and PDB as values

mut_site = dict()

k = 0
k1 = 0
oldk1 = 0
while k < len(ht):
	ht1 = ht[k].split(", ")
	t1 = ht1[0].strip("'|[|]|\n")

	oldk1 = k1
	count = 0
	gt1 = gt[k1].split(", ")
	t2 = gt1[0].strip("'|[|]|\n")
	while t2 != t1:
		k1 = k1 + 1
		if k1 >= len(gt):
			count = count + 1
			break
		gt1 = gt[k1].split(", ")
		t2 = gt1[0].strip("'|[|]|\n")

	if count != 0:
		k1 = oldk1

	k2 = 1
	while k2 < len(ht1):
		t2 = ht1[k2].strip("'|[|]|\n")
		if count == 0:
			t3 = gt1[k2].strip("'|[|]|\n")
		else:
			t3 = t2
		list1 = [(t2,t3)]
		if "{}".format(t1) in mut_site.keys():
			mut_site["{}".format(t1)].append(list1)
		else:
			mut_site["{}".format(t1)] = list1

		k2 = k2 + 1

	k = k + 1

# forming a number_of_distant_mutations.csv file based on the PDB mutation sites

m = open("number_of_distint_mutations_pdb.csv","w")

k = 0
while k < len(ft):
	ft1 = ft[k].split(",")
	t1 = ft1[0].strip("\n")
	print(t1)
	m.write("{}".format(t1))

	k1 = 1
	while k1 < len(ft1):
		t2 = ft1[k1].strip("\n")
		res1 = t2[0:1]
		res2 = t2[(len(t2)-1):len(t2)]
		pos = t2[1:(len(t2)-1)]

		value = mut_site["{}".format(t1)]

		k2 = 0
		while k2 < len(value):
			if "{}".format(value[k2][0]) == "{}".format(pos):
				t3 = value[k2][1]
			k2 = k2 + 1

		t4 = "{}{}{}".format(res1,t3,res2)

		m.write(",{}".format(t4))

		k1 = k1 + 1
	
	m.write("\n")

	k = k + 1
		
m.close()
