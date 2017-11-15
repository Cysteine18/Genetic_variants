# DIVITION OF THE WHOLE DATA INTO LENGTH AND ORGANISM

f = open("number_of_distinct_mutants.csv","r")
ft = f.readlines()
f.close()

g = open("pdb_seqres.txt","r")
gt = g.readlines()
g.close()

h = open("cluster_length.txt","w")
h1 = open("organism.txt","w")

k = 0
while k < len(ft):
	s1 = "{}".format(ft[k])
	ft1 = s1.split(",")
	count = 0
	k2=0
	org_str = ""
	while k2 < len(ft1):
		t1 = ft1[k2].strip("\n")
		print("### {} :: {} of {} ###".format(t1,k,len(ft)))
		k1 = 0
		gt1 = gt[k1].split()
		t2 = gt1[0].strip(">")
		while t2 != t1:
			k1 = k1 + 2
			gt1 = gt[k1].split()
			t2 = gt1[0].strip(">")

		if len(gt1) > 3:
			org = gt1[3]
			k3 = 4
			while k3 < len(gt1):
				org = org + "_{}".format(gt1[k3])
				k3 = k3 + 1
		else:
			org = "NF"

		org_str = org_str + "{} ".format(org)
		if k2 == 0:
			old_org = org
			id_clus = k1
			rep_res = t1
		else:
			count1 = 0
			if old_org in org:
				count1 = count1 + 1
			if org in old_org:
				count1 = count1 + 1
			if org.lower() == old_org:  
				count1 = count1 + 1
			if old_org.lower() == org:  
				count1 = count1 + 1
			if old_org.upper() == org:  
				count1 = count1 + 1
			if org.upper() == old_org:  
				count1 = count1 + 1
			if org == old_org:
				count1 = count1 + 1
			if count1 == 0:
				count = count + 1
		k2 = k2+ 1

	t1len = len(gt[id_clus+1])
	t1len = t1len - 1
	if count == 0:
		h.write("{} {} {} {}".format(rep_res,len(ft1),t1len,old_org))
	else:
		h.write("{} {} {} NC".format(rep_res,len(ft1),t1len,org_str))
	h.write("\n")
	h1.write("{}".format(org_str))
	h1.write("\n")

	k = k + 1

h.close()
h1.close()
