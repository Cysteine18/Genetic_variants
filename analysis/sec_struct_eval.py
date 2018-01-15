# THIS CODE EVALUATES THE CHANGE IN THE SECONDARY STRUCTURE UPON MUTATION

def sec_struct():
	f = open("distinct_mutants_pdb.txt","r")
	ft = f.readlines()
	f.close()

	g = open("distinct_mutants_struct.txt","r")
	gt = g.readlines()
	g.close()

	d = dict()
	k = 0
	k1 = 0
	while k < len(ft):
		ft1 = ft[k].split(", ")
		t1 = ft1[0].strip("[|'|]|\n")
		d["{}".format(t1)] = []
		oldk1 = k1
		gt1 = gt[k1].split(", ")
		t2 = gt1[0].strip("[|'|]|\n")
		count = 0
		while t1 != t2:
			k1 = k1 + 1
			if k1 > len(gt):
				count = count + 1
				k1 = oldk1
				break
			gt1 = gt[k1].split(", ")
			t2 = gt1[0].strip("[|'|]|\n")

		check1 = gt1[1].strip("[|'|]|\n")
		if count == 0 and check1 != "DATA_NOT_FOUND":
			k2 = 1
			while k2 < len(ft1):
				t2 = ft1[k2].strip("[|'|]|\n")
				t3 = gt1[k2].strip("[|'|]|\n")
				
				list1 = (t2,t3)
				d["{}".format(t1)].append(list1)
				k2 = k2 + 1
		k = k + 1
	
	return(d)

def mutants():
	import sys

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("secondary_struct_var.txt","w")
	
	# 0 IF THE STRUCTURE REMAIN SAME, 1 FOR CHANGE IN SECONDARY STRUCTURE AND 2 LARGE CHANGE

	# LARGER CLASSES DICTIONARY
	lc = dict()
	lc["G"] = 1
	lc["H"] = 1 
	lc["I"] = 1
	lc["B"] = 2
	lc["E"] = 2
	lc["T"] = 3
	lc["S"] = 3
	lc[" "] = 3
	
	SS = sec_struct()

	d = dict()		# MUTANTS DICTIONARY
	k = 0
	cid = 1
	nmut = 0
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		try:
			t2 = ft1[1].strip("\n")
		except:
			t2 = "NF"
		if t1 == "#" and t2 == "{}".format(cid):
			#print("LOOKING FOR MUTANTS IN CLUSTER {} WITH CLUSTER ID {}".format(cid,t1))	
			cid = cid + 1
			k = k + 1
			ft1 = ft[k].split()
			wt = ft1[0].strip("('|',")
			k = k + 1
			ft1 = ft[k].split()
			t1 = ft1[0].strip("\n")
			while t1 != "#" and k < (len(ft)-1):
				mut = ft1[0].strip("('|',")
				check = ft1[(len(ft1)-1)].strip("'|',|,)|\n")
				check1 = ft1[(len(ft1)-2)].strip("'|',|")
				if check != "NOT_FOUND" and check != "NA" and check1 != "NOT_FOUND" and check1 != "NA":
					pos = ft1[(len(ft1)-1)].split(",")
					pos1 = pos[0].strip("'|')")
					d[nmut] = []
					list1 = (wt,mut,pos1)
					d[nmut].append(list1)
					nmut = nmut + 1
					k = k + 1
					ft1 = ft[k].split()
					t1 = ft1[0].strip("\n")
				else:
					k = k + 1
					ft1 = ft[k].split()
					t1 = ft1[0].strip("\n")
		else:
			k = k + 1

	k = 0
	c0 = 0
	c1 = 0
	c2 = 0
	while k < len(d):
		print("ANALYSING MUTANT {} OF {}".format(k,len(d)))
		k1 = 0
		while k1 < len(d[k]):
			wt = d[k][k1][0]
			mut = d[k][k1][1]
			pos = d[k][k1][2]
		
			ss_wt = SS["{}".format(wt)]
			ss_mut = SS["{}".format(mut)]

			if len(ss_wt) > 0 and len(ss_mut) > 0:
				# FOR WILD_TYPE
				k2 = 0
				pos1 = ss_wt[k2][0]
				while pos != pos1:
					k2 = k2 + 1
					pos1 = ss_wt[k2][0]
				wt_ss = ss_wt[k2][1]

				# FOR MUTANT
				k2 = 0
				pos1 = ss_mut[k2][0]
				while pos != pos1:
					k2 = k2 + 1
					pos1 = ss_mut[k2][0]
				mut_ss = ss_mut[k2][1]

				if mut_ss == wt_ss:
					g.write("{} {} {} {} {} {} 0".format((k+1),wt,mut,wt_ss,mut_ss,pos))
					c0 = c0 + 1
					g.write("\n")
				else:
					# CHECKING THE DIFFERENT LARGER CLASSES
					if lc["{}".format(mut_ss)] == lc["{}".format(wt_ss)]:			
						g.write("{} {} {} {} {} {} 1".format((k+1),wt,mut,wt_ss,mut_ss,pos))
						c1 = c1 + 1
					else:
						g.write("{} {} {} {} {} {} 2".format((k+1),wt,mut,wt_ss,mut_ss,pos))
						c2 = c2 + 1
					g.write("\n")
			k1 = k1 +1
		k = k + 1

	print("\n")
	print("#### SECONDARY STRUCTURE :: SAME = {} AND DIFFERENT = {},{} ####".format(c0,c1,c2))
	g.write("#### SECONDARY STRUCTURE :: SAME = {} AND DIFFERENT = {},{} ####".format(c0,c1,c2))
	g.close()
	
mutants()




































