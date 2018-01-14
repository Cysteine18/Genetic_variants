def res_filter():
	f = open("PDB_classifier.txt","r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		print(t1)
		t2 = ft1[1].strip("\n")
	
		if t1 not in d.keys():
			d["{}".format(t1)] = t2

		k = k + 1

	print(d["101m"])

res_filter()

