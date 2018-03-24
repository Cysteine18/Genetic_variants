import sys
pathfiles = "/Users/tarunkhanna/Documents/Bioinformatics/file_links/top_8000"

def top_8000():
	num_files = 4
	diff_hom = ["50","70","90","95"]

	d = dict()
	for x in range(0,num_files):
		f = open("{}/Top8000-best_hom{}_pdb_chain_list.csv".format(pathfiles,diff_hom[x]),"r")
		ft = f.readlines()
		f.close()

		k = 1
		while k < len(ft):
			ft1 = ft[k].split(",")
			t1 = ft1[0]
			t2 = ft1[1].strip("\n")
			key1 = "{}_{}".format(t1,t2)

			if key1 not in d.keys():
				d["{}".format(key1)] = "IN"
			k = k + 1

	return(d)

def our_data():

	top = top_8000()

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("mod_mut_data.txt","w")

	k = 0
	count = 0
	d2 = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0]
		t2 = ft1[1]

		count1 = 0
		if t1 not in d2.keys():
			d2["{}".format(t1)] = "IN"
			if t1 in top.keys():
				g.write("{}".format(ft[k]))
				count1 = count1 + 1
				count = count + 1

		if t2 not in d2.keys():
			d2["{}".format(t2)] = "IN"
			if t2 in top.keys():
				count = count + 1
				if count1 == 0:
					g.write("{}".format(ft[k]))

		k = k + 1

	print(count)
	print(len(ft))

	g.close()

our_data()






