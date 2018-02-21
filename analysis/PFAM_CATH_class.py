def PFAM(file1):

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close()

	# DICTIONARY OF ALL THE UNIQUE PDBS

	p = dict()
	k = 0
	while k < len(ft):
		#print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t1 = t1.lower()
		t2 = ft1[1].strip("\n")
		t3 = ft1[2].strip("\n")
		t4 = ft1[3].strip("\n")
		t5 = ft1[4].strip("\n")
		key1 = "{}_{}".format(t1,t2)
		list1 = (t3,t4,t5)
		if key1 not in p.keys():
			p["{}".format(key1)] = [list1]
		else:
			p["{}".format(key1)].append(list1)
		k = k + 1

def cath(file1):

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close()

	# DICTIONARY OF ALL THE UNIQUE PDBS

	p = dict()
	k = 0
	while k < len(ft):
		#print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		temp = ft1[0].strip("\n")
		t1 = temp[0:4]
		t2 = temp[4:(len(temp)-2)]
		t3 = temp[len(temp):(len(temp)-2)]
		key1 = "{}_{}".format(t1,t2)
		list1 = [t3]
		if key1 not in p.keys():
			p["{}".format(key1)] = list1
		else:
			p["{}".format(key1)].append(list1)
		k = k + 1

	p1 = dict()
	for x in p.keys():
		t1 = len(p["{}".format(x)])
		p1["{}".format(x)] = t1

PFAM("pdb_pfam_mapping.txt")
cath("cath_domain.txt")

			
			
