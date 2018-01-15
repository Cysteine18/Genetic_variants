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

def main_function():
	import sys

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("failed_uniprot_test_sim_struct_rows.txt","w")

	g1 = open("failed_uniprot_test_sim_struct_pairwise.txt","w")

	uniprot = pdbchain_level("pdb_chain_uniprot.csv") # UNIPROT DICTIONARY TO ACCOUNT FOR CHIMERA'S

	k = 0
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split(",")
		k1 = 0

		# ALL AGAINST ALL COMPARISON

		while k1 < len(ft1):
			t1 = ft1[k1].strip("\n")
			k2 = 0
			while k2 < len(ft1):
				if k1 != k2:
					t2 = ft1[k2].strip("\n")
					try:
						ut1 = uniprot["{}".format(t1)]
						ut2 = uniprot["{}".format(t2)]
					except:
						ut1 = "NF"
						ut2 = "NF"
	
					if ut1 != ut2:
						g1.write("{} {}".format(t1,t2))
						g1.write("\n")

				k2 = k2 + 1
			k1 = k1 + 1
		k = k + 1

	g.close()
				
main_function()











			
