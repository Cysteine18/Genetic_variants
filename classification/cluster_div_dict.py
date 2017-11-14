def pdbchain_level(arg):
	f = open("{}".format(arg),"r")
	ft = f.readlines()
	f.close()

	k = 1
	list1 = []
	list2 = []
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = "{}_{}".format(ft1[0],ft1[1])
		list1.append("{}".format(t1))	# KEY
		list2.append(ft1[2])	# VALUE
		k = k + 1

	#listall = zip(list1,list2)
	d = dict()
	k = 0
	while k < len(list1):
		t1 = list1[k]
		t2 = list2[k]
		#print(t1)
		if t1 in d.keys():
			k1 = 0
			count = 0
			while k1 < len(d[list1[k]]):
				if "{}".format(d[list1[k]][k1]) == "{}".format(list2[k]):
					count = count + 1
				k1 = k1 + 1
			if count == 0:
				#print("{} {}".format(d[list1[k]][0],list2[k]))
				d[list1[k]].append(list2[k])
		else:
			d[list1[k]] = [list2[k]]
		k = k + 1

	return(d)

def mut_zone(arg1,arg2,arg3):

	# FORMING THE DICTIONARY ELEMENT OF MUTATION ZONE

	"""
	D[$pdbchain] = [($resname(wildtype),$resname(mutant),$resid(wildtype),$resid(mutant) .........]

	"""
	
	f = open("{}".format(arg1),"r")
	ft = f.readlines()
	f.close()

	m = open("{}".format(arg2),"r")
	mt = m.readlines()
	m.close()

	g = open("{}".format(arg3),"r")
	gt = g.readlines()
	g.close()
	
	h = open("zone_dict.txt","w")

	zone = dict()
	k = 0
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = ft1[0]

		mt1 = mt[k].split(",")
	
		#print("{} :: {} of {}".format(t1,k,len(ft)))

		k1 = 1
		while k1 < len(ft1):
			list1 = []
			t2 = ft1[k1].strip("\n")
			mut = t2[1:(len(t2)-1)]
	
			# WILDTYPE
			k2 = 0
			gt1 = gt[k2].split(",")
			t3 = gt1[0].strip("\n") # PDB CHAIN
			k2 = k2 + 1
			gt1 = gt[k2].split(",")
			t4 = gt1[0]	# MUTATION POSITION
			count1 = 0
			count2 = 0
			while t3 != t1 and t4 != mut:
				if k2 > len(gt1):
					count1 = count1 + 1
					break
				k2 = k2 + 1
				gt1 = gt[k2].split(",")
				z = gt1[0].strip("\n")
				k2 = k2 + 1
				gt1 = gt[k2].split(",")
				t3 = gt1[0].strip("\n")
				k2 = k2 + 1
				gt1 = gt[k2].split(",")
				t4 = gt1[0]

			# MUTANT

			M = mt1[k1].strip("\n")
			k3 = 0
			gt1 = gt[k3].split(",")
			t3 = gt1[0].strip("\n") # PDB CHAIN
			k3 = k3 + 1
			gt1 = gt[k3].split(",")
			t4 = gt1[0]	# MUTATION POSITION
			count2 = 0
			while t3 != M and t4 != mut:
				if k3 > len(gt1):
					count2 = count2 + 1
					break
				k3 = k3 + 1
				gt1 = gt[k3].split(",")
				z = gt1[0].strip("\n")
				k3 = k3 + 1
				gt1 = gt[k3].split(",")
				t3 = gt1[0].strip("\n")
				k2 = k3 + 1
				gt1 = gt[k3].split(",")
				t4 = gt1[0]


			if count1 == 0 and count2 == 0: 
				z1 = gt[k2 + 1].strip("\n") # WILDTYPE
				res1 = gt[k2].strip("\n")
				z2 = gt[k3 + 1].strip("\n")	# MUTANT
				res2 = gt[k3].strip("\n")
				if "{}".format(t1) in zone.keys():
					list1 = (z1,z2,res1,res2)
					zone["{}".format(t1)].append(list1)
				else:
					list1 = [(z1,z2,res1,res2)]
					zone["{}".format(t1)] = list1
			elif count1 == 0 and count2 != 0:
				z1 = gt[k2 + 1].strip("\n") # WILDTYPE
				res1 = gt[k2].strip("\n")
				z2 = "NOT_FOUND"	# MUTANT
				res2 = "NOT_FOUND"
				if "{}".format(t1) in zone.keys():
					list1 = (z1,z2,res1,res2)
					zone["{}".format(t1)].append(list1)
				else:
					list1 = [(z1,z2,res1,res2)]
					zone["{}".format(t1)] = list1
			elif count1 != 0 and count2 == 0:
				z1 = "NOT_FOUND"	# WILDTYPE
				res1 = "NOT_FOUND"
				z2 = gt[k3 + 1].strip("\n") # MUTANT
				res2 = gt[k3].strip("\n")
				if "{}".format(t1) in zone.keys():
					list1 = (z1,z2,res1,res2)
					zone["{}".format(t1)].append(list1)
				else:
					list1 = [(z1,z2,res1,res2)]
					zone["{}".format(t1)] = list1
			elif count1 != 0 and count2 != 0:
				z1 = "NOT_FOUND"	# WILDTYPE
				res1 = "NOT_FOUND"
				z2 = "NOT_FOUND"	# MUTANT
				res2 = "NOT_FOUND"
				if "{}".format(t1) in zone.keys():
					list1 = (z1,z2,res1,res2)
					zone["{}".format(t1)].append(list1)
				else:
					list1 = [(z1,z2,res1,res2)]
					zone["{}".format(t1)] = list1

			k1 = k1 + 1
		k = k + 1
	return(zone)

def cluster():

	# INPUT COLUMN DATA

	#############################################################################

	f = open("number_of_distinct_mutants.csv","r")
	ft = f.readlines()
	f.close()

	g = open("organism.txt","r") # PRE PROCESSED
	gt = g.readlines()
	g.close()

	gc = open("cluster_length.txt","r") # PRE PROCESSED
	gct = gc.readlines()
	gc.close()

	uniprot = pdbchain_level("pdb_chain_uniprot.csv") # UNIPROT DICTIONARY TO ACCOUNT FOR CHIMERA'S

	n = open("number_of_distinct_mutations.csv","r")
	nt = n.readlines()
	n.close()

	zone = mut_zone("number_of_distinct_mutations.csv","number_of_distinct_mutants.csv","zone.txt") # DICTIONARY FOR ZONE AROUND THE MUTATION WITH CLUSTER NAME AS THE DICTIONARY KEY

	z = open("zone_dict.txt","w")
	z.write("{}".format(zone))
	z.close()

	##############################################################################

	h = open("cluster.txt","w")

	cluster = dict()

	k = 0
	while k < len(ft):
		ft1=ft[k].split(",")
		t1 = ft1[0].strip("\n")
		print(t1)
		K = k + 1 #key 

		print ("CLUSTER {} OF {}".format(K,len(ft)))
		gct1 = gct[k].split() # CLUSTER HEADER
		gt1 = gt[k].split() # ORGANISM INFO

		KE = ["{} {} {} {}".format(gct1[0],gct1[1],gct1[2],gct1[3].strip("\n"))]
		cluster[K] = [KE]
	
		k1 = 0
		while k1 < len(ft1):
			t1 = "{}".format(ft1[k1].split("\n"))
			t2 = gt1[k1].strip("\n")
			t3 = "{}".format(uniprot["{}".format(ft1[k1].strip("\n"))])
			if k1 == 0:
				t4= zone["{}".format(ft1[0].strip("\n"))][0][0]
				#t4 = "{}".format(zone["{}[0][0]".format((ft1[k1].strip("\n")))])
				#t5 = "{}".format(zone["{}[0][2]".format((ft1[k1].strip("\n")))])
				t5= zone["{}".format(ft1[0].strip("\n"))][0][2]
				KE1 = "{} {} {} {} {}".format(t1,t2,t3,t4,t5)
			else: 
				t4= zone["{}".format(ft1[0].strip("\n"))][(k1-1)][1]
				print(t4,ft1[k1])
				#t4 = "{}".format(zone["{}[{}][1]".format((ft1[k1].strip("\n")),(k1-1))])
				#t5 = "{}".format(zone["{}[{}][3]".format((ft1[k1].strip("\n")),(k1-1))])
				t5= zone["{}".format(ft1[0].strip("\n"))][(k1-1)][2]

			cluster[K].append(KE1)
		
			k1 = k1 + 1
	
		k = k + 1

	h.write(cluster)

cluster()
	
	
	
	
	
	
	
	
