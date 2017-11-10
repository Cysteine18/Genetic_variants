# CHIMERIC PDBs IN RCSB DATABASE

def pdbchain_level():
	f = open("pdb_chain_uniprot.csv","r")
	ft = f.readlines()
	f.close()

	g = open("chimera_pdbchain_level.csv","w")
	k = 1
	list1 = []
	list2 = []
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = "{}{}".format(ft1[0],ft1[1])
		list1.append("{}".format(t1))	# KEY
		list2.append(ft1[2])	# VALUE
		k = k + 1

	#listall = zip(list1,list2)
	d = dict()
	k = 0
	while k < len(list1):
		t1 = list1[k]
		t2 = list2[k]
		print(t1)
		if t1 in d.keys():
			k1 = 0
			count = 0
			while k1 < len(d[list1[k]]):
				if "{}".format(d[list1[k]][k1]) == "{}".format(list2[k]):
					count = count + 1
				k1 = k1 + 1
			if count == 0:
				print("{} {}".format(d[list1[k]][0],list2[k]))
				d[list1[k]].append(list2[k])
		else:
			d[list1[k]] = [list2[k]]
		k = k + 1
	
	#d = dict(listall)
	#print(len(d['4qxeA']))

	for key in d:
		if len(d[key]) > 1:
			row1 = "{},{}".format(key[0:4],key[4:len(key)])
			k1 = 0
			while k1 < len(d[key]):
				row1 = row1 + ",{}".format(d[key][k1]) 
				k1 = k1 + 1
			g.write("{}".format(row1))
			g.write("\n")

	g.close()

def uniprot_level():
	f = open("chimera_pdbchain_level.csv","r")
	ft = f.readlines()
	f.close()

	g = open("chimera_uniprot_level.csv","w")

	# KEYS
	k = 0
	list1 = []
	list2 = []
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = "{}{}".format(ft1[0],ft1[1])
		k1 = 2
		while k1 < len(ft1):
			t2 = ft1[k1].strip("\n")
			list1.append(t2)	# KEY
			list2.append("{}".format(t1))	# VALUE
			k1 = k1 + 1
		k = k + 1

	d = dict()
	k = 0
	while k < len(list1):
		t1 = list1[k]
		t2 = list2[k]
		print(t2)
		if t1 in d.keys():
			if "{}".format(d[list1[k]][0]) != "{}".format(list2[k]):
				print("{} {}".format(d[list1[k]][0],list2[k]))
				d[list1[k]].append(list2[k])
		else:
			d[list1[k]] = [list2[k]]
		k = k + 1

	for key in d:
		row1 = "{}".format(key)
		k1 = 0
		while k1 < len(d[key]):
			row1 = row1 + ",{},{}".format(d[key][k1][0:4],d[key][k1][4:len(d[key][k1])]) 
			k1 = k1 + 1
		g.write("{}".format(row1))
		g.write("\n")

	g.close()


#pdbchain_level()
uniprot_level()
