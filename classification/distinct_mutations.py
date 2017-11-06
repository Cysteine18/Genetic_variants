# FORMING A CSV FILE CONTAINING DISTINCT MUTATIONS

f = open("all_mutations_mod.txt","r")
ft = f.readlines()
f.close()

g = open("all_mutants.txt","r")
gt = g.readlines()
g.close()

h = open("number_of_distinct_mutants.txt","r")
ht = h.readlines()
h.close()

m = open("number_of_distinct_mutations.csv","w")
n = open("number_of_distinct_mutants.csv","w")

end = len(ht) - 2
k = 0
while k < end:
	tempt1 = "{}".format(ht[k])
	ht1 = tempt1.split()
	k1 = 0
	t1 = ht1[k1].strip("}|{|\n")
	print(t1)
	m.write("{}".format(t1))
	n.write("{}".format(t1))

	k2 = 0
	gt1 = gt[k2].split(",")
	t2 = gt1[0]
	while t2 != t1:
		k2 = k2 + 1
		gt1 = gt[k2].split(",")
		t2 = gt1[0]
	ft1 = ft[k2].split(",")
	k1 = k1 + 1	
	while k1 < len(ht1):
		t1 = ht1[k1].strip("}|{|\n")
		n.write(",{}".format(t1))
		k3 = 1
		t2 = gt1[k3].strip("\n")
		while t2 != t1:
			k3 = k3 + 1
			t2 = gt1[k3].strip("\n")
		t3 = ft1[k3].strip("\n")
		m.write(",{}".format(t3))
		k1 = k1 + 1
	m.write("\n")
	n.write("\n")
	k = k + 1

m.close()
n.close()	
