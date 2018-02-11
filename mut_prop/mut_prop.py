# ASSIGNMENT OF PROPERTIES TO EACH MUTATIONS

def res_filter():
	f = open("PDB_classifier.txt","r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
	
		if t1 not in d.keys():
			d["{}".format(t1)] = t2

		k = k + 1

	return(d)

def R_free():
	# THIS FUNCTIO CALCULATES THE AVERAGE B FACTOR AND R FREE FOR THE DIFFERENT PDB'S

	f = open("r_factor.csv","r")
	ft = f.readlines()
	f.close()

	k = 0
	r_free = dict()
	abfactor = dict()

	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = ft1[0].strip("\n|\"")
		t1 = t1.lower()
		t2 = ft1[4].strip("\n|\"")
		t3 = ft1[5].strip("\n|\"")

		if t1 not in r_free.keys():
			if t2 == "":
				r_free["{}".format(t1)] = "NA"
			else:
				r_free["{}".format(t1)] = t2
			if t3 == "":
				abfactor["{}".format(t1)] = "NA"
			else:
				abfactor["{}".format(t1)] = t3

		k = k + 1

	return(r_free,abfactor)

def temp_factor(file1,var1,var2):
	# file1 = file and var1 = Number of rows and var2 = row of the property 

	f = open("{}".format(file1),"r")
	ft = f.readlines()
	f.close()

	avg = dict()
	max1 = dict()
	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		p = ft1[0].strip("\n")
		pos = ft[k+1].split(",")
		p1 = pos[0].strip("\n")		
		t1 = ft[(k+(int(var2)-1))].split(",")
		check = t1[0].strip("\n")
		if check != "NA":
			k1 = 0
			m = 0.0
			a = 0.0
			while k1 < len(t1):
				t2 = t1[k1].strip("\n")
				if float(t2) > m:
					m = float(t2)
				a = a + float(t2)
				k1 = k1 + 1

			a = a / len(t1)
			a = round(a,2)

			v = "({},{})".format(p,p1)
			if v not in max1.keys():
				max1["{}".format(v)] = m
				avg["{}".format(v)] = a
			
		k = k + int(var1)

	return(avg,max1)
		
def mut_prop():

	import sys

	# THIS CODE ROSS AND SANDER DEFINATION OF MAXIMUM ACCESSIBILITY TO CALCULATE THE BURIED AND EXPOSED RESIDUES

	list1 = ['A','B','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','X','Y','Z','a']
	list2 = [106,160,135,163,194,197,84,184,169,205,164,188,157,136,198,248,130,142,142,227,180,222,196,135]

	acc = dict()
	for x in range(0,len(list1)):
		acc["{}".format(list1[x])] = list2[x]

	# DIVISON OF AMINO ACIDS BASED ON THEIR NATURE

	tup1 = ('R','K','D','E')	# CHARGED SIDE CHAIN ; KEY 1
	tup2 = ['Q','N','H','S','T','Y','C','W','a']	# POLAR AMINO ACIDS ; KEY 2
	tup3 = ['A','I','L','M','F','V','P','G']	# HYDROPHOBIC AMINO ACIDS	; KEY 3

	nature = dict()
	nature[1] = tup1
	nature[2] = tup2
	nature[3] = tup3

	NC = 0
	NNC = 0


	# BINARY MODEL 'B' 'E'
	modB1 = 0.16
	BN = 0
	BY = 0
	BB = 0
	BE = 0
	
	# TERNARY MODEL 'B' 'I' 'E'
	modT1 = 0.09
	modT2 = 0.36
	TN = 0
	TY = 0
	TB = 0
	TI = 0
	TE = 0

	res = res_filter()		# FOR RESOLUTION OF EACH PDB FILE

	rf = R_free() # R FACTOR FOR EACH PDB FILE ALONG WITH AVERAGE B FACTOR

	f = open("solvent_ass.csv","r")
	ft = f.readlines()
	f.close()

	g = open("mut_prop.csv","w")
	g.write("#wt,chwt,mut,chmut,reswt,resmut,pos,wt_acc,mut_acc,B/E_WT,B/E_MUT,change,B/I/E_WT,B/I/E_MUT,change,nature_change,type_of_change,ontology,local_rmsd,mut_localisation,pos_seqres, seq_length,WT_sec_struct,mut_sec_struct,c_alpha_wt,c_alpha_mut,res_wt,res_mut,r_free_wt,r_free_mut,bfactor_avg_wt, bfactor_mut_avg,avg_bfactor_ms_wt,avg_bactor_mut,max_bfactor_ms_wt,max_bfactor_ms_mut,avg_AVG_factor_10Azone_WT,avg_AVG_bfactor_10AzoneMUT, max_AVG_bfactor_10Azone_WT,max_AVG_bfactor_10Azone_MUT,avg_MAX_bfactor_10AzoneMUT,avg_MAX_bfactor_10AzoneMUT, max_MAX_bfactor_10Azone_WT,max_MAX_bfactor_10Azone_MUT")

	g.write("\n")

	k = int(sys.argv[1])
	end = sys.argv[2]
	if end == "END" or end == "end":
		end = len(ft)
	end = int(end)
	while k < end: #len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split(",")
		check = ft1[4].strip("\n")
		if check != "ERROR":
			count = 0
			t1 = ft1[0].strip("\n")
			t2 = ft1[1].strip("\n")
			t3 = ft1[2].strip("\n")
			t4 = ft1[3].strip("\n")
			t5 = ft1[4].strip("\n")	# WILDTYPE RESIDUE
			t6 = ft1[5].strip("\n")	# MUTANT RESIDUE
			t7 = ft1[6].strip("\n")
			t8 = ft1[7].strip("\n")
			t9 = ft1[8].strip("\n")

			# DETERMING THE NATURE OF AMINO ACIDS

			# FOR WILDTYPE
			count1 = 0
			for x in range(1,4):
				k1 = 0	
				while k1 < len(nature[x]):
					if nature[x][k1] == t5:
						nwt = x
						k1 = len(nature[x])
						x = 3
						count1 = count1 + 1
					k1 = k1 + 1
			if count1 == 0:
				nwt = "ERROR"

			# FOR MUTANT
			count1 = 0
			for x in range(1,4):
				k1 = 0	
				while k1 < len(nature[x]):
					if nature[x][k1] == t6:
						nmut = x
						k1 = len(nature[x])
						x = 3
						count1 = count1 + 1
					k1 = k1 + 1
			if count1 == 0:
				nmut = "ERROR"

			# COMPARISON
				
			if nwt != "ERROR" and nmut != "ERROR":
				t17 = "{}->{}".format(nwt,nmut)
				if nwt == nmut:
					t16 = "NO"
					NNC = NNC + 1
				else:
					t16 = "YES"
					NC = NC + 1
			else:
				t17 = "{}->{}".format(nwt,nmut)
				t16 = "ERROR"
	
			# COLUMN FOR ONTOLOGY

			o = open("TM_cluster_name.txt","r")
			ot = o.readlines()
			o.close()

			ot1 = ot[0].split(", ")
			k1 = 0
			count1 = 0
			while k1 < len(ot1):
				oc = ot1[k1].strip("'|[|['|']")
				oc = oc.lower()
				if oc == t1:
					count1 = count1 + 1
					k1 = len(ot1)
				k1 = k1 + 1

			if count1 == 0:
				t18 = "G"
			else:
				t18 = "TM"

			# COLUMN FOR LOCAL RMSD

			l = open("local_rmsd_clustal.out","r")
			lt = l.readlines()
			l.close()

			l1 = "{}_{}".format(t1,t2)
			l2 = "{}_{}".format(t3,t4)

			k1 = 0
			lt1 = lt[k1].split()
			l1c = lt1[1]
			l2c = lt1[2]
			count1 = 0
			while l1c != l1 or l2c != l2:
				k1 = k1 + 1
				if k1 >= len(lt):
					count1 = count1 + 1
					break
				lt1 = lt[k1].split()
				l1c = lt1[1]
				l2c = lt1[2]

			if count1 == 0:
				t19 = lt1[3].strip("\n")
			else:
				t19 = "ERROR"

			# COLUMN FOR MUTATION LOCALISATION

			m = open("mut_coverage.txt","r")
			mt = m.readlines()
			m.close()

			k1 = 0
			mt1 = mt[k1].split()
			m1 = mt1[0]
			m2 = mt1[1]
			count1 = 0
			while m1 != l1 or m2 != l2:
				k1 = k1 + 1
				if k1 >= len(mt):
					count1 = count1 + 1
					break
				mt1 = mt[k1].split()
				m1 = mt1[0]
				m2 = mt1[1]

			if count1 == 0:
				t20 = mt1[2].strip("\n")
				t21 = mt1[3].strip("\n")
				t22 = mt1[4].strip("\n")
			else:
				t20 = "ERROR"
				t21 = "ERROR"
				t22 = "ERROR"

			# COLUMNS FOR SECONDARY STRUCTURE ASSIGNMENT

			s = open("secondary_struct_var.txt","r")
			st = s.readlines()
			s.close()

			k1 = 0
			st1 = st[k1].split()
			s1 = st1[1]
			s2 = st1[2]
			count1 = 0
			while s1 != l1 or s2 != l2:
				k1 = k1 + 1
				if k1 >= len(st):
					count1 = count1 + 1
					break
				st1 = st[k1].split()
				s1 = st1[1]
				s2 = st1[2]
			if count1 == 0:
				t23 = st1[3]
				t24 = st1[4]
			else:
				t23 = "ERROR"
				t24 = "ERROR"

			# COLUMN FOR C-ALPHA DENSITY BASED ON 5A SPHERE AROUND THE MUTATION SITE

			c = open("c_alpha_density.txt","r")
			ct = c.readlines()
			c.close()

			k1 = 0
			ct1 = ct[k1].split()
			c1 = ct1[1]
			c2 = ct1[2]
			count1 = 0
			while c1 != l1 or c2 != l2:
				k1 = k1 + 1
				if k1 >= len(ct):
					count1 = count1 + 1
					break
				ct1 = ct[k1].split()
				c1 = ct1[1]
				c2 = ct1[2]

			if count1 == 0:
				t25 = ct1[3].strip("\n")
				t26 = ct1[4].strip("\n")
			else:
				t25 = "ERROR"
				t26 = "ERROR"

			# COLUMN FOR RESOLUTION

			try:
				t27 = res["{}".format(t1)]
				t28 = res["{}".format(t3)]
			except:
				t27 = "ERROR"
				t28 = "ERROR"

			# COLUMN FOR R FREE AND AVERAGE B FACTOR

			try:
				t29 = rf[0]["{}".format(t1)]
				t31 = rf[1]["{}".format(t1)]
				t30 = rf[0]["{}".format(t3)]
				t32 = rf[1]["{}".format(t3)]
			except:
				t29 = "ERROR"
				t30 = "ERROR"
				t31 = "ERROR"
				t32 = "ERROR"

			# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE MUTATION SITE

			tf = temp_factor("TF_mut_site.txt",3,3)

			try:
				tf1 = "({},{})".format(l1,t7)
				tf2 = "({},{})".format(l2,t7)
				t33 = tf[0]["{}".format(tf1)]
				t35 = tf[1]["{}".format(tf1)]
				t34 = tf[0]["{}".format(tf2)]
				t36 = tf[1]["{}".format(tf2)]
			except:
				t33 = "NA"
				t35 = "NA"
				t34 = "NA"
				t36 = "NA"

			# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE AVERAGE ZONE MUTATION SITE

			tf = temp_factor("B_factor.txt",4,3)

			try:
				tf1 = "({},{})".format(l1,t7)
				tf2 = "({},{})".format(l2,t7)
				t37 = tf[0]["{}".format(tf1)]
				t39 = tf[1]["{}".format(tf1)]
				t38 = tf[0]["{}".format(tf2)]
				t40 = tf[1]["{}".format(tf2)]
			except:
				t37 = "NA"
				t38 = "NA"
				t39 = "NA"
				t40 = "NA"

			# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE MAXIMUM ZONE MUTATION SITE

			tf = temp_factor("B_factor.txt",4,4)

			try:
				tf1 = "({},{})".format(l1,t7)
				tf2 = "({},{})".format(l2,t7)
				t41 = tf[0]["{}".format(tf1)]
				t43 = tf[1]["{}".format(tf1)]
				t42 = tf[0]["{}".format(tf2)]
				t44 = tf[1]["{}".format(tf2)]
			except:
				t41 = "NA"
				t42 = "NA"
				t43 = "NA"
				t44 = "NA"
			
			try:
				wtra = int(t8) / acc["{}".format(t5)]	# WILDTYPE
				mutra = int(t9) / acc["{}".format(t6)]	# MUTANT
			except:
				count = count + 1

			if count == 0:
			
				# BINARY ASSIGNMENT

				if wtra < modB1:
					t10 = "B"
					BB = BB + 1
				else:
					t10 = "E"
					BE = BE + 1
			
				if mutra < modB1:
					t11 = "B"
					BB = BB + 1
				else:
					t11 = "E"
					BE = BE + 1

				if t10 == t11:
					t12 = "NO"
					BN = BN + 1
				else:
					t12 = "YES"
					BY = BY + 1

				# TERNARY ASSIIGNMENT

				if wtra < modT1:
					t13 = "B"
					TB = TB + 1
				elif wtra > modT1 and wtra < modT2:
					t13 = "I"
					TI = TI + 1
				elif wtra > modT2:
					t13 = "E"
					TE = TE + 1

				if mutra < modT1:
					t14 = "B"
					TB = TB + 1
				elif mutra > modT1 and mutra < modT2:
					t14 = "I"
					TI = TI + 1
				elif mutra > modT2:
					t14 = "E"
					TE = TE + 1

				if t13 == t14:
					t15 = "NO"
					TN = TN + 1
				else:
					t15 = "YES"
					TY = TY + 1
				
				g.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14,t15,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44))
				g.write("\n")
			else:
				g.write("{},{},{},{},{},{},{},{},{},ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(t1,t2,t3,t4,t5,t6,t7,t8,t9,t16,t17,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44))
				g.write("\n")
		else:
				t1 = ft1[0].strip("\n")
				t2 = ft1[1].strip("\n")
				t3 = ft1[2].strip("\n")
				t4 = ft1[3].strip("\n")

				# COLUMN FOR ONTOLOGY

				o = open("TM_cluster_name.txt","r")
				ot = o.readlines()
				o.close()

				ot1 = ot[0].split(", ")
				k1 = 0
				count1 = 0
				while k1 < len(ot1):
					oc = ot1[k1].strip("'|[|['|']")
					oc = oc.lower()
					if oc == t1:
						count1 = count1 + 1
						k1 = len(ot1)
					k1 = k1 + 1

				if count1 == 0:
					t18 = "G"
				else:
					t18 = "TM"

				# COLUMN FOR LOCAL RMSD

				l = open("local_rmsd_clustal.out","r")
				lt = l.readlines()
				l.close()

				l1 = "{}_{}".format(t1,t2)
				l2 = "{}_{}".format(t3,t4)

				k1 = 0
				lt1 = lt[k1].split()
				l1c = lt1[1]
				l2c = lt1[2]
				count1 = 0
				while l1c != l1 or l2c != l2:
					k1 = k1 + 1
					if k1 >= len(lt):
						count1 = count1 + 1
						break
					lt1 = lt[k1].split()
					l1c = lt1[1]
					l2c = lt1[2]

				if count1 == 0:
					t19 = lt1[3].strip("\n")
				else:
					t19 = "ERROR"

				# COLUMN FOR MUTATION LOCALISATION

				m = open("mut_coverage.txt","r")
				mt = m.readlines()
				m.close()

				k1 = 0
				mt1 = mt[k1].split()
				m1 = mt1[0]
				m2 = mt1[1]
				count1 = 0
				while m1 != l1 or m2 != l2:
					k1 = k1 + 1
					if k1 >= len(mt):
						count1 = count1 + 1
						break
					mt1 = mt[k1].split()
					m1 = mt1[0]
					m2 = mt1[1]

				if count1 == 0:
					t20 = mt1[2].strip("\n")
					t21 = mt1[3].strip("\n")
					t22 = mt1[4].strip("\n")
				else:
					t20 = "ERROR"
					t21 = "ERROR"
					t22 = "ERROR"

				# COLUMNS FOR SECONDARY STRUCTURE ASSIGNMENT

				s = open("secondary_struct_var.txt","r")
				st = s.readlines()
				s.close()

				k1 = 0
				st1 = st[k1].split()
				s1 = st1[1]
				s2 = st1[2]
				count1 = 0
				while s1 != l1 or s2 != l2:
					k1 = k1 + 1
					if k1 >= len(st):
						count1 = count1 + 1
						break
					st1 = st[k1].split()
					s1 = st1[1]
					s2 = st1[2]
				if count1 == 0:
					t23 = st1[3]
					t24 = st1[4]
				else:
					t23 = "ERROR"
					t24 = "ERROR"

				# COLUMN FOR C-ALPHA DENSITY BASED ON 5A SPHERE AROUND THE MUTATION SITE

				c = open("c_alpha_density.txt","r")
				ct = c.readlines()
				c.close()

				k1 = 0
				ct1 = ct[k1].split()
				c1 = ct1[1]
				c2 = ct1[2]
				count1 = 0
				while c1 != l1 or c2 != l2:
					k1 = k1 + 1
					if k1 >= len(ct):
						count1 = count1 + 1
						break
					ct1 = ct[k1].split()
					c1 = ct1[1]
					c2 = ct1[2]

				if count1 == 0:
					t25 = ct1[3].strip("\n")
					t26 = ct1[4].strip("\n")
				else:
					t25 = "ERROR"
					t26 = "ERROR"

				# COLUMN FOR RESOLUTION

				try:
					t27 = res["{}".format(t1)]
					t28 = res["{}".format(t3)]
				except:
					t27 = "ERROR"
					t28 = "ERROR"

				# COLUMN FOR R FREE AND AVERAGE B FACTOR

				try:
					t29 = rf[0]["{}".format(t1)]
					t31 = rf[1]["{}".format(t1)]
					t30 = rf[0]["{}".format(t3)]
					t32 = rf[1]["{}".format(t3)]
				except:
					t29 = "ERROR"
					t30 = "ERROR"
					t31 = "ERROR"
					t32 = "ERROR"


				# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE MUTATION SITE

				tf = temp_factor("TF_mut_site.txt",3,3)

				try:
					tf1 = "({},{})".format(l1,t7)
					tf2 = "({},{})".format(l2,t7)
					t33 = tf[0]["{}".format(tf1)]
					t35 = tf[1]["{}".format(tf1)]
					t34 = tf[0]["{}".format(tf2)]
					t36 = tf[1]["{}".format(tf2)]
				except:
					t33 = "NA"
					t35 = "NA"
					t34 = "NA"
					t36 = "NA"

				# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE AVERAGE ZONE MUTATION SITE

				tf = temp_factor("B_factor.txt",4,3)

				try:
					tf1 = "({},{})".format(l1,t7)
					tf2 = "({},{})".format(l2,t7)
					t37 = tf[0]["{}".format(tf1)]
					t38 = tf[1]["{}".format(tf1)]
					t39 = tf[0]["{}".format(tf2)]
					t40 = tf[1]["{}".format(tf2)]
				except:
					t37 = "NA"
					t38 = "NA"
					t39 = "NA"
					t40 = "NA"

				# COLUMN FOR AVERAGE AND MAX TEMP FACTOR FOR THE MAXIMUM ZONE MUTATION SITE

				tf = temp_factor("B_factor.txt",4,4)

				try:
					tf1 = "({},{})".format(l1,t7)
					tf2 = "({},{})".format(l2,t7)
					t41 = tf[0]["{}".format(tf1)]
					t43 = tf[1]["{}".format(tf1)]
					t42 = tf[0]["{}".format(tf2)]
					t44 = tf[1]["{}".format(tf2)]
				except:
					t41 = "NA"
					t42 = "NA"
					t43 = "NA"
					t44 = "NA"
			

				g.write("{},{},{},{},ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,ERROR,{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}".format(t1,t2,t3,t4,t18,t19,t20,t21,t22,t23,t24,t25,t26,t27,t28,t29,t30,t31,t32,t33,t34,t35,t36,t37,t38,t39,t40,t41,t42,t43,t44))
				g.write("\n")

		k = k + 1

	g.close()
	h = open("buried_exposed_stats","w")
	h.write("# BINARY ASSIGNMENT\n")
	h.write("BURIED :: {}\n".format(BB))
	h.write("EXPOSED :: {}\n".format(BE))
	h.write("CHANGE :: {}\n".format(BY))
	h.write("NO_CHANGE :: {}\n".format(BN))
	nf = len(ft) - int(BY) - int(BN)
	h.write("CAN'T ASSIGN :: {}\n".format(nf))
	h.write("\n")

	h.write("# TERNARY ASSIGNMENT\n")
	h.write("BURIED :: {}\n".format(TB))
	h.write("INTERMEDIATE :: {}\n".format(TI))
	h.write("EXPOSED :: {}\n".format(TE))
	h.write("CHANGE :: {}\n".format(TY))
	h.write("NO_CHANGE :: {}\n".format(TN))
	nf = len(ft) - int(TY) - int(TN)
	h.write("CAN'T ASSIGN :: {}\n".format(nf))
	h.write("\n")

	h.write("# CHANGE OF NATURE\n")
	h.write("CHANGE :: {}\n".format(NC))
	h.write("NO_CHANGE :: {}\n".format(NNC))
	nf = len(ft) - int(NC) - int(NNC)
	h.write("CAN'T ASSIGN :: {}\n".format(nf))

	h.close()
	

mut_prop()













				
			


