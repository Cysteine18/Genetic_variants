# UNIPROT DISTRIBUTION AMONG THE CLUSTERS OF SIZE GREATER THAN 20

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

def uniprot_dis():
	
	uniprot = pdbchain_level("pdb_chain_uniprot.csv") # UNIPROT DICTIONARY

	f = open("cluster_length.txt","r")
	ft = f.readlines()
	f.close()

	g = open("uniprot_dis_clusters.csv","w")

	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")

		if int(t3) > 20:
			try:
				if len(uniprot["{}".format(t1)]) > 1:
					t4 = "{}".format(uniprot["{}".format(t1)][0])
					k1 = 1
					while k1 < len(uniprot["{}".format(t1)]):
						t4 = t4 + ",{}".format(uniprot["{}".format(t1)][k1])
						k1 = k1 + 1
				else:
					t4 = "{}".format(uniprot["{}".format(t1)])
					t4 = t4.strip("'|[|]")
			except:
				t4 = "NA"
			g.write("{},{},{},{}".format(t1,t2,t3,t4))
			g.write("\n")

		k = k + 1

	g.close()

def uniprot_count():

	import numpy as np

	f = open("uniprot_dis_clusters.csv","r")
	ft = f.readlines()
	f.close()
	
	g = open("unique_uniprot.txt","w")

	k = 0
	N = 0
	unique_uniprot = []
	num = []
	while k < len(ft):
		ft1 = ft[k].split(",")
		t2 = int(ft1[1]) - 1
		k1 = 3
		while k1 < len(ft1):
			t1 = ft1[k1].strip("\n")
			k2 = 0
			count = 0
			while k2 < len(unique_uniprot):
				if unique_uniprot[k2] == t1:
					num[k2] = num[k2] + t2
					count = count + 1
				k2 = k2 + 1
			if count == 0:
				num.append(t2)
				unique_uniprot.append(t1)
				N = N + 1

			k1 = k1 + 1

		k = k + 1

	print(len(unique_uniprot))

	for x in range(0,N):
		g.write("{} {} {}".format(x+1,unique_uniprot[x],num[x]))
		g.write("\n")
	g.close()

#uniprot_dis()
uniprot_count()
