# TO FORM THE DISTINCT MUTANT FILE WITH ONLY HIGHEST RESOLUTION PDB'S

f = open("number_of_distinct_mutants.csv","r")
ft = f.readlines()
f.close()

g = open("PDB_classifier.txt","r")
gt = g.readlines()
g.close()

h = open("seq_sim.txt","r")
ht = h.readlines()
h.close()

m = open("number_of_distinct_mutants_after_res_filter.csv","w")

k = 0
while k < len(ft):
	print("{} of {}".format(k,len(ft)))
	ft1 = ft[k].split(",")
	k1 = 0
	list1 = ""
	while k1 < len(ft1):
		t1 = ft1[k1].strip("\n")
		k2 = 0
		ht1 = ht[k2].split(",")
		t2 = ht1[0].strip("\n")
		while t1 != t2:
			k2 = k2 + 1
			if k2 >= len(ht):
				print("ERROR")
				quit()
			ht1 = ht[k2].split(",")
			t2 = ht1[0].strip("\n")
	
		pdb = t1
		res = 100.0
		k3 = 0
		while k3 < len(ht1):
			t2 = ht1[k3].strip("\n")

			k4 = 0
			gt1 = gt[k4].split()
			t3 = gt1[0].strip("\n")
			count = 0
			while t3 != t2[0:4]:
				k4 = k4 + 1
				if k4 >= len(gt):
					count = count + 1
					break
				gt1 = gt[k4].split()
				t3 = gt1[0].strip("\n")

			if count == 0:
				t4 = gt1[1].strip("\n")
				if t4 != "None":
					if float(t4) < res:
						pdb = t2
			k3 = k3 + 1

		if k1 == 0:
			list1 = list1 + "{}".format(pdb)
		else:
			list1 = list1 + ",{}".format(pdb)

		k1 = k1 + 1

	m.write("{}".format(list1))
	m.write("\n")

	k = k + 1

m.close()
	

























						
			
			
		
