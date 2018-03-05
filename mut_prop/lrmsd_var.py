# THIS CODE READS THE COLUMNS OF THE MUTATION PROPERTIES CSV FILE TO REDUCE DOWN THE SAMPLE 

def dis(var1,var2,var3):

	import sys

	# THIS FUNCTION CALCULATES THE DISTRIBUTION OF ANY COLUMN OF THE DATA DIVIDING INTO 10 BINS
	binsize = 15.0

	# MUTATION PEOPERTY FILE

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	min1 = 1000.0
	max1 = 0.0
	k = 1
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = ft1[int(var1)-1].strip("\n")
		if t1 != "NA" and t1 != "ERROR":
			if float(t1) > max1:
				max1 = float(t1)
			if float(t1) < min1:
				min1 = float(t1)
		k = k + 1

	max1 = round(max1,2)
	min1 = round(min1,2)
	bin1 = (var3 - var2) / binsize
	bin1 = round(bin1,2) 
	k = var2
	list1 = []
	list2 = []
	print("MINIMUM = {} , MAXIMUM = {}".format(min1,max1))
	while k < var3:
		list1.append(k)
		k1 = 1
		count = 0
		while k1 < len(ft):
			ft1 = ft[k1].split(",")
			t1 = ft1[int(var1)-1].strip("\n")
			if t1 != "NA" and t1 != "ERROR":
				if float(t1) >= k and float(t1) < (k + bin1):
					count = count + 1
			k1 = k1 + 1
		list2.append(count)
		k = k + bin1
		k = round(k,2)

	k1 = 1
	count = 0
	list1.append(k)
	while k1 < len(ft):
		ft1 = ft[k1].split(",")
		t1 = ft1[int(var1)-1].strip("\n")
		if t1 != "NA" and t1 != "ERROR":
			if float(t1) >= k:
				count = count + 1
		k1 = k1 + 1
	list2.append(count)

	return(list1,list2)

def main_func():
	import sys

	# MUTATION PEOPERTY FILE
	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("{}.txt".format(sys.argv[2]),"w")

	g1 = open("mut_prop_{}.csv".format(sys.argv[2]),"w")

	g1.write("#wt,chwt,mut,chmut,reswt,resmut,pos,wt_acc,mut_acc,B/E_WT,B/E_MUT,change,B/I/E_WT,B/I/E_MUT,change,nature_change,type_of_change,ontology,local_rmsd,mut_localisation,pos_seqres, seq_length,WT_sec_struct,mut_sec_struct,sec_struct_change,c_alpha_wt,c_alpha_mut,res_wt,res_mut,r_free_wt,r_free_mut,bfactor_avg_wt, bfactor_mut_avg,avg_bfactor_ms_wt,avg_bactor_mut,max_bfactor_ms_wt,max_bfactor_ms_mut,avg_AVG_factor_10Azone_WT,avg_AVG_bfactor_10AzoneMUT, max_AVG_bfactor_10Azone_WT,max_AVG_bfactor_10Azone_MUT,avg_MAX_bfactor_10AzoneMUT,avg_MAX_bfactor_10AzoneMUT, max_MAX_bfactor_10Azone_WT,max_MAX_bfactor_10Azone_MUT,author_WT,author_MUT,IF_change_author,TM_GRMSD,TM_LRMSD,TM_LRMSD_SC,WT_aligned_ratio,MUT_aligned_ratio, deposotion_year_WT,deposition_year_MUT,exp_type_WT,exp_type_MUT,potential_bound_unbound,cath_domains_WT,cath_domains_MUT,PFAM_acc_WT, PFAM_acc_MUT,DISULPHIDE_BOND,SALT_BRIDGE,BIOLOGICAL_UNIT\n")

	# COL1 = WT, COL2 = WT_CHAIN, COL3 = MUT, COL4 = MUT_CHAIN, COL5 = RES_WT, COL6 = RES_MUT, COL7 = POS, COL8 = WT_SASA, COL9 = MUT_SASA, COL10 = B/E_WT, COL11 = B/E_MUT, COL12 = CHANGE, COL13 = B/I/E_WT, COL14 = B/I/E_MUT, COL15 = CHANGE, COL16 = IF CHANGE IN NATURE, COL17 = TYPE OF CHANGE, COL18 = ONTOLOGY, COL19 =  LRMSD_CLUSTAL, COL20 = MUTATION LOCALISATION, COL21 = POSITION SEQRES, COL22 = SEQUENCE LENGTH, COL23 = SEC_STRUCT_WT, COL24 = SEC_STRUCT_MUT, COL25 = SEC_STRCUT_CHANGE, COL26 = C_ALPHA_WT, COL27 = C_ALPHA_MUT, COL28 = RES_WT, COL29 = RES_MUT, COL30 = RFREE WT, COL31 = RFREE MUT, COL32 = BFACTOR_WT, COL33 = BFACTOR_MUT, COL34 = AVG_BFACTOR_MUT_SITE WT, COL35 = AVG_BFACTOR_MUT_SITE MUT, COL36 = MAX_BFACTOR_MUT_SITE WT, COL37 = MAX_BFACTOR_MUT_SITE MUT, COL38 = AVG AVG_BFACTOR_10A_ZONE WT, COL39 = AVG AVG_BFACTOR_10A_ZONE MUT, COL40 = MAX AVG_BFACTOR_10A_ZONE WT, COL41 = MAX AVG_BFACTOR_10A_ZONE MUT, COL42 = AVG MAX_BFACTOR_10A_ZONE WT, COL43 = AVG MAX_BFACTOR_10A_ZONE WT, COL44 = MAX MAX_BFACTOR_10A_ZONE WT, COL45 = MAX MAX_BFACTOR_10A_ZONE MUT, COL46 = GRMSD_TM, COL47 = LRMSD_TM, COL48 = LRMSDSC_TM

	# PROPERTIES LOCAL

	# COL12 = CHANGE, COL15 = CHANGE, COL16 = IF CHANGE IN NATURE, COL17 = TYPE OF CHANGE, COL19 =  LRMSD_CLUSTAL, COL21 = POSITION SEQRES, COL22 = SEQUENCE LENGTH, COL23 = SEC_STRUCT_WT, COL24 = SEC_STRUCT_MUT, COL25 = SEC_STRCUT_CHANGE, COL34 = AVG_BFACTOR_MUT_SITE WT, COL35 = AVG_BFACTOR_MUT_SITE MUT, COL36 = MAX_BFACTOR_MUT_SITE WT, COL37 = MAX_BFACTOR_MUT_SITE MUT, COL38 = AVG AVG_BFACTOR_10A_ZONE WT, COL39 = AVG AVG_BFACTOR_10A_ZONE MUT, COL40 = MAX AVG_BFACTOR_10A_ZONE WT, COL41 = MAX AVG_BFACTOR_10A_ZONE MUT, COL42 = AVG MAX_BFACTOR_10A_ZONE WT, COL43 = AVG MAX_BFACTOR_10A_ZONE WT, COL44 = MAX MAX_BFACTOR_10A_ZONE WT, COL45 = MAX MAX_BFACTOR_10A_ZONE MUT

	# PROPERTIES GLOBAL

	# COL28 = RES_WT, COL29 = RES_MUT, COL30 = RFREE WT, COL31 = RFREE MUT, COL32 = BFACTOR_WT, COL33 = BFACTOR_MUT

	col1 = [10,11,12,13,14,15,16,17,19,21,22,23,24,25,28,29,30,31,32,33,26,27,5,6,49,50,51,48,52,53,54,55,56,57,58,59,60,63,64,65,34,35,36,37,38,39,40,41,42]

	valuecol1 = [["=/!=","B/E/ERROR"],["=/!=","B/E/ERROR"],["=","YES/NO/ERROR"],["=/!=","B/I/E/ERROR"],["=/!=","B/I/E/ERROR"],["!","YES/NO/ERROR"],["=", "YES/NO/ERROR"],["=/!=", "1->2,1->3,2->3"],[">/</=/<>","float number1/float number1 float number2/ERROR"],[">/</=/<>","NUMBER1/NUMBER1 NUMBER2"],[">/</=/<>","NUMBER1/NUMBER1 NUMBER2"],["=/!=","H/G/I/B/E/T/S/O/ERROR"],["=/!=","H/G/I/B/E/T/S/O/ERROR"],["=/!=","0/1/2"],[">/</=/<>", "float num1/float num1 float num2/None"],[">/</=/<>","float num1/float num1 float num2/None"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","int num1/int num1 int num2/ERROR"],[">/</=/<>","int num1/int num1 int num2/ERROR"],["=/!=","AMINO_ACID"],["=/!=","AMINO_ACID"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],["=","SAME/DIFFERENT/NA"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],[">/</=/<>","year/year1 year2/ERROR/NA"],[">/</=/<>","year/year1 year2/ERROR/NA"],["=/!=","diffraction,nmr,crystallography"],["=/!=","diffraction,nmr,crystallography"],["=/!=","YES/NO"],["=/!=","1/2.../NA"],["=/!=","1/2.../NA"],["=/!=","YES/NO/NA"],["=/!=","YES/NO/NA"],["=/!=","YES/NO"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"]]

	# DICTIONARY WITH KEY AS THE COLUMN

	NC = 65
	d = dict()
	for x in range(0,NC):
		d[(x+1)] = []

	k = 1
	while k < len(ft):
		ft1 = ft[k].split(",")
		k1 = 0
		while k1 < len(ft1):
			t1 = ft1[k1].strip("\n")
			d[(k1+1)].append("{}".format(t1))
			k1 = k1 + 1
		k = k + 1

	print("DO YOU WANT TO CALCULATE THE DISTRIBUTION OF ANY OF THE COLUMNS?(y/n)")
	input1 = input()
	
	if input1 == "Y" or input1 == "y":
		print("ENTER THE COLUMN NUMBER")
		input2 = input()
		print("ENTER THE MINIMUM VALUE OF THE PROPERTY")
		min1 = input()
		print("ENTER THE MAXIMUM VALUE OF THE PROPERTY")
		max1 = input()
		d = dis(int(input2),float(min1),float(max1))
		k = 0
		while k < len(d[0]):
			g.write("{} {}\n".format(d[0][k],d[1][k]))
			k = k + 1
	else:
		g.write("#number row_number wt wt_chain mut mut_chain pos lrmsd grmsd_TM lrmsd_TM lrmsdSC_TM\n")
		print("ENTER THE NUMBER OF COLUMNS YOU WANT THE CRITERIA TO BE BASED ON?)")
		input3 = input()
		list1 = []
		list2 = []
		for x in range(0,int(input3)):
			print("ENTER THE COLUMN {}".format((x+1)))
			inp = input()
	
			k1 = 0
			count = 0
			while k1 < len(col1):
				if int(inp) == col1[k1]:
					count = count + 1
					break
				k1 = k1 + 1

			if count != 0:
				p = valuecol1[k1]
			else:
				quit()
			
			print("ENTER THE VALUE COLUMN {} CAN TAKE SEPARATED BY SPACE".format(inp))
			print("HINT COL{} CAN TAKE {} VALUES".format(inp,p))
			value = input()
			list1.append(inp)
			list2.append(value)

		print("ENTER THE NUMBER OF NUMBER OF CROSS COLUMN CRITERIAS?)")
		input3 = input()
		list3 = []
		for x in range(0,int(input3)):
			print("ENTER THE CRITERIA {}".format((x+1)))
			print("FORMAT :: $col1 <,>,=,!= $col2/$col1 +,-,*,/,diff $col2  <,> value")
			inp = input()
			list3.append(inp)

		k = 0
		nterms = 1
		while k < len(d[1]):
			k1 = 0
			count = 0
			while k1 < len(list1):
				t1 = int(list1[k1])
				t2 = list2[k1].split()
				prop = d[t1][k]
				if prop != "ERROR" and prop != "NA" and prop != "nan":
					t = t2[0]
					k2 = 1
					if len(t2) < 3 or t == "<>":
						if t == "=":
							if prop != t2[k2]:
								count = count + 1
								k1 = len(list1)
						if t == "<":
							if prop != "None":
								if float(prop) > float(t2[k2]):
									count = count + 1
									k1 = len(list1)
							else:
								count = count + 1
								k1 = len(list1)
						if t == ">":
							if prop != "None":
								if float(prop) < float(t2[k2]):
									count = count + 1
									k1 = len(list1)
							else:
								count = count + 1
								k1 = len(list1)
						if t == "!=":
							if prop == t2[k2]:
								count = count + 1
								k1 = len(list1)
						if t == "<>":
							if float(prop) < float(t2[k2]) and float(prop) > float(t2[k2+1]):
								count = count + 1
								k1 = len(list1)
					else:
						# FOR MULTIPLE ENTRIES
						count1 = 0
						while k2 < len(t2):
							if t == "=":
								if prop == t2[k2]:
									count1 = count1 + 1
							if t == "!=":
								if prop != t2[k2]:
									count1 = count1 + 1
							k2 = k2 + 1
						if count1 == 0:
							count = count + 1
							k1 = len(list1)
							
					k1 = k1 + 1
				else:
					count = count + 1
					k1 = len(list1)
			
			if count == 0:
				k1 = 0
				while k1 < len(list3):
					c = list3[k1].split()
					t1 = int(c[0])
					t = c[1]
					t2 = int(c[2])
					prop1 = d[t1][k]
					prop2 = d[t2][k]
					if prop1 != "ERROR" and prop1 != "NA" and prop1 != "nan" and prop2 != "ERROR" and prop2 != "NA" and prop2 != "nan":
						if t == "=":
							if prop1 != prop2:
								count = count + 1
								k1 = len(list3)
						if t == "!=":
							if prop1 == prop2:
								count = count + 1
								k1 = len(list3)
						if t == ">":
							if prop1 != "None" and prop2 != "None":
								if float(prop1) <= float(prop2):
									count = count + 1 
									k1 = len(list3)
							else:
								count = count + 1
								k1 = len(list3)
						if t == "<":
							if prop != "None" and prop2 != "None":
								if float(prop1) >= float(prop2):
									count = count + 1
									k1 = len(list3)
							else:
								count = count + 1
								k1 = len(list3) 
						if t == "+":
							tt = c[3]
							t3 = float(c[4])
							if tt == ">":
								if (float(prop1) + float(prop2)) <= t3 :
									count = count + 1
									k1 = len(list3)
							if tt == "<":
								if (float(prop1) + float(prop2)) >= t3 :
									count = count + 1
									k1 = len(list3)
						if t == "-":
							tt = c[3]
							t3 = float(c[4])
							if tt == ">":
								if (float(prop1) - float(prop2)) <= t3 :
									count = count + 1
									k1 = len(list3)
							if tt == "<":
								if (float(prop1) - float(prop2)) >= t3 :
									count = count + 1
									k1 = len(list3)
						if t == "diff":
							tt = c[3]
							t3 = float(c[4])
							if tt == ">":
								if abs((float(prop1) - float(prop2))) <= t3 :
									count = count + 1
									k1 = len(list3)
							if tt == "<":
								if abs((float(prop1) - float(prop2))) >= t3 :
									count = count + 1
									k1 = len(list3)
						if t == "*":
							tt = c[3]
							t3 = float(c[4])
							if tt == ">":
								if (float(prop1) * float(prop2)) <= t3 :
									count = count + 1
									k1 = len(list3)
							if tt == "<":
								if (float(prop1) * float(prop2)) >= t3 :
									count = count + 1
									k1 = len(list3)
						if t == "/":
							tt = c[3]
							t3 = float(c[4])
							if tt == ">":
								if (float(prop1) / float(prop2)) <= t3 :
									count = count + 1
									k1 = len(list3)
							if tt == "<":
								if (float(prop1) / float(prop2)) >= t3 :
									count = count + 1
									k1 = len(list3)
						k1 = k1 + 1
					else:
						count = count + 1
						k1 = len(list3)


			if count == 0:
				g.write("{} {} {} {} {} {} {} {} {} {} {}\n".format(nterms,(k+2),d[1][k],d[2][k],d[3][k],d[4][k],d[7][k],d[19][k],d[49][k],d[50][k],d[51][k]))
				nterms = nterms + 1
				kk = 1
				st1 = ""
				while kk <= NC:
					if kk == 1:
						st1 = st1 + "{}".format(d[kk][k])
					else:
						st1 = st1 + ",{}".format(d[kk][k])
					kk = kk + 1
				
				g1.write("{}\n".format(st1)) 	

			k = k + 1

	g.close()

main_func()
			
	
	 
	












