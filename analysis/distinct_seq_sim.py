# TO REMOVE REPEATS IN THE 100 PERCENT SIMILAR SEQUENCES

def all_to_all():
	f = open("seq_sim.txt","r")
	ft = f.readlines()
	f.close()

	g = open("distinct_seq_sim.csv","w")

	h = open("capri_fit_data.txt","w")
	h1 = open("capri_fit_data_sc.txt","w")

	k = 0
	list1 = []
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split(",")
		if len(ft1) > 1:
			t0 = ft1[0].strip("\n")
			k1 = 0
			count = 0
			while k1 < len(list1):
				if list1[k1] == t0:
					count = count + 1
					k1 = len(list1)
				k1 = k1 + 1

			if count == 0:
				g.write("{}".format(ft[k]))
				k1 = 0
				while k1 < len(ft1):
					t1 = ft1[k1].strip("\n")
					list1.append(t1)
					pdb1 = t1[0:4]
					chain1 = t1[5:len(t1)]
				
					k2 = k1 + 1
					while k2 < len(ft1):
						t2 = ft1[k2].strip("\n")
						list1.append(t2)
						pdb2 = t2[0:4]
						chain2 = t2[5:len(t2)]

						if k1 == 0:
							h1.write("{} {} {} {}".format(pdb1,chain1,pdb2,chain2))
							h1.write("\n")
					
						h.write("{} {} {} {}".format(pdb1,chain1,pdb2,chain2))
						h.write("\n")
						k2 = k2 + 1
					k1 = k1 + 1

		k = k + 1

	g.close()
	h.close()	

def dictionary_sim_seq():
	f = open("seq_sim.txt","r")
	ft = f.readlines()
	f.close()

	g = open("distinct_seq_sim.csv","w")

	h = open("capri_fit_data.txt","w")
	h1 = open("capri_fit_data_sc.txt","w")

	d = dict()
	d1 = dict()
	cid = 0
	k = 0
	while k < len(ft):
		ft1 = ft[k].split(",")
		if len(ft1) > 1:
			t1 = ft1[0].strip("\n")
			count = 0
			k1 = 0
			while k1 < len(ft1):
				t2 = ft1[k1].strip("\n")
				if t2 in d.keys():
					count = count + 1
					break
				k1 = k1 + 1
			
			if count == 0:
				d["{}".format(t1)] = []
				d1[cid] = [t1]
				k1 = 1
				while k1 < len(ft1):
					t3 = ft1[k1].strip("\n")
					d["{}".format(t1)].append(t3)
					d1[cid].append(t3)
					k1 = k1 + 1	
				cid = cid + 1
		k = k + 1

	k = 0
	while k < len(d1):
		if k > 0:
			g.write("\n")
		t1 = "{}".format(d1[k][0])
		k1 = 1
		while k1 < len(d1[k]):
			t1 = t1 + ",{}".format(d1[k][k1])
			k1 = k1 + 1	
		g.write("{}".format(t1))

		k1 = 0
		while k1 < len(d1[k]):
			t2 = d1[k][k1].strip("\n")
			pdb1 = t2[0:4]
			chain1 = t2[5:len(t2)]
				
			k2 = k1 + 1
			while k2 < len(d1[k]):
				t3 = d1[k][k2].strip("\n")
				pdb2 = t3[0:4]
				chain2 = t3[5:len(t3)]

				if k1 == 0:
					h1.write("{} {} {} {}".format(pdb1,chain1,pdb2,chain2))
					h1.write("\n")
					
				h.write("{} {} {} {}".format(pdb1,chain1,pdb2,chain2))
				h.write("\n")
				k2 = k2 + 1
			k1 = k1 + 1
		
		k = k + 1

	g.close()
	h.close()
	h1.close()


dictionary_sim_seq()	
				
		
		
























			 
				
		
