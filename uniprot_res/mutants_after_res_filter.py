# TO FORM THE DISTINCT MUTANT FILE WITH ONLY HIGHEST RESOLUTION PDB'S

def res_filter():
	f = open("PDB_classifier.txt","r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
	
		if t1 not in d.keys():
			d["{}".format(t1)] = t2

		k = k + 1

	return(d)

def mutants():

	resolution = res_filter()
	f = open("number_of_distinct_mutants.csv","r")
	ft = f.readlines()
	f.close()

	g = open("PDB_classifier.txt","r")
	gt = g.readlines()
	g.close()

	h = open("seq_sim.txt","r")
	ht = h.readlines()
	h.close()

	m = open("number_of_distinct_mutants_after_res_filter1.csv","w")

	k = 6989
	while k < 9450: #len(ft):
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
				try:
					t4 = resolution["{}".format(t2[0:4])]
				except:
					t4 = "None"
				if t4 != "None":
					if float(t4) < res:
						res = float(t4)
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

def mutations():

	resolution = res_filter()

	f = open("number_of_distinct_mutations_pdb.csv","r")
	ft = f.readlines()
	f.close()

	g = open("PDB_classifier.txt","r")
	gt = g.readlines()
	g.close()

	h = open("seq_sim.txt","r")
	ht = h.readlines()
	h.close()

	m = open("number_of_distinct_mutations_after_res_filter_pdb.csv","w")

	k = 0
	while k < 9450: #len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split(",")
		k1 = 0
		list1 = ""
		while k1 < len(ft1):
			t1 = ft1[k1].strip("\n")
			if k1 == 0:
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
					try:
						t4 = resolution["{}".format(t2[0:4])]
					except:
						t4 = "None"
					if t4 != "None":
						if float(t4) < res:
							res = float(t4)
							pdb = t2
					k3 = k3 + 1
				list1 = list1 + "{}".format(pdb)
			else:
				list1 = list1 + ",{}".format(t1)

			k1 = k1 + 1

		m.write("{}".format(list1))
		m.write("\n")

		k = k + 1

	m.close()
	
mutations()
























						
			
			
		
