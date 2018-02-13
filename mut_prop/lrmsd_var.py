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

	g = open("{}".format(sys.argv[2]),"w")

	# COL1 = WT, COL2 = WT_CHAIN, COL3 = MUT, COL4 = MUT_CHAIN, COL5 = RES_WT, COL6 = RES_MUT, COL7 = POS, COL8 = WT_SASA, COL9 = MUT_SASA, COL10 = B/E_WT, COL11 = B/E_MUT, COL12 = CHANGE, COL13 = B/I/E_WT, COL14 = B/I/E_MUT, COL15 = CHANGE, COL16 = IF CHANGE IN NATURE, COL17 = TYPE OF CHANGE, COL18 = ONTOLOGY, COL19 =  LRMSD_CLUSTAL, COL20 = MUTATION LOCALISATION, COL21 = POSITION SEQRES, COL22 = SEQUENCE LENGTH, COL23 = SEC_STRUCT_WT, COL24 = SEC_STRUCT_MUT, COL25 = SEC_STRCUT_CHANGE, COL26 = C_ALPHA_WT, COL27 = C_ALPHA_MUT, COL28 = RES_WT, COL29 = RES_MUT, COL30 = RFREE WT, COL31 = RFREE MUT, COL32 = BFACTOR_WT, COL33 = BFACTOR_MUT, COL34 = AVG_BFACTOR_MUT_SITE WT, COL35 = AVG_BFACTOR_MUT_SITE MUT, COL36 = MAX_BFACTOR_MUT_SITE WT, COL37 = MAX_BFACTOR_MUT_SITE MUT, COL38 = AVG AVG_BFACTOR_10A_ZONE WT, COL39 = AVG AVG_BFACTOR_10A_ZONE MUT, COL40 = MAX AVG_BFACTOR_10A_ZONE WT, COL41 = MAX AVG_BFACTOR_10A_ZONE MUT, COL42 = AVG MAX_BFACTOR_10A_ZONE WT, COL43 = AVG MAX_BFACTOR_10A_ZONE WT, COL44 = MAX MAX_BFACTOR_10A_ZONE WT, COL45 = MAX MAX_BFACTOR_10A_ZONE MUT, COL46 = GRMSD_TM, COL47 = LRMSD_RM, COL48 = LRMSDSC_RM

	# PROPERTIES LOCAL

	# COL12 = CHANGE, COL15 = CHANGE, COL16 = IF CHANGE IN NATURE, COL17 = TYPE OF CHANGE, COL19 =  LRMSD_CLUSTAL, COL21 = POSITION SEQRES, COL22 = SEQUENCE LENGTH, COL23 = SEC_STRUCT_WT, COL24 = SEC_STRUCT_MUT, COL25 = SEC_STRCUT_CHANGE, COL34 = AVG_BFACTOR_MUT_SITE WT, COL35 = AVG_BFACTOR_MUT_SITE MUT, COL36 = MAX_BFACTOR_MUT_SITE WT, COL37 = MAX_BFACTOR_MUT_SITE MUT, COL38 = AVG AVG_BFACTOR_10A_ZONE WT, COL39 = AVG AVG_BFACTOR_10A_ZONE MUT, COL40 = MAX AVG_BFACTOR_10A_ZONE WT, COL41 = MAX AVG_BFACTOR_10A_ZONE MUT, COL42 = AVG MAX_BFACTOR_10A_ZONE WT, COL43 = AVG MAX_BFACTOR_10A_ZONE WT, COL44 = MAX MAX_BFACTOR_10A_ZONE WT, COL45 = MAX MAX_BFACTOR_10A_ZONE MUT

	# PROPERTIES GLOBAL

	# COL28 = RES_WT, COL29 = RES_MUT, COL30 = RFREE WT, COL31 = RFREE MUT, COL32 = BFACTOR_WT, COL33 = BFACTOR_MUT

	col1 = [10,11,12,15,16,17,19,21,22,23,24,25,28,29,30,31,32,33,26,27,5,6,46,47,48]

	valuecol1 = [["=/!=","B/E/ERROR"],["=/!=","B/E/ERROR"],["=","YES/NO/ERROR"],["!","YES/NO/ERROR"],["=", "YES/NO/ERROR"],["=", "1->2,1->3,2->3"],[">/</=/<>","float number1/float number1 float number2/ERROR"],[">/</=/<>","NUMBER1/NUMBER1 NUMBER2"],[">/</=/<>","NUMBER1/NUMBER1 NUMBER2"],["=/!=","H/G/I/B/E/T/S/O/ERROR"],["=/!=","H/G/I/B/E/T/S/O/ERROR"],["=/!=","0/1/2"],[">/</=/<>", "float num1/float num1 float num2/None"],[">/</=/<>","float num1/float num1 float num2/None"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","float num1/float num1 float num2/NA"],[">/</=/<>","int num1/int num1 int num2/ERROR"],[">/</=/<>","int num1/int num1 int num2/ERROR"],["=/!=","AMINO_ACID"],["=/!=","AMINO_ACID"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"],[">/</=/<>","float number1/float number1 float number2/ERROR/NA"]]

	# DICTIONARY WITH KEY AS THE COLUMN

	NC = 48
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
		g.write("#number row_number wt wt_chain mut mut_chain pos lrmsd lrmsd_TM lrmsdSC_TM\n")
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
							if prop > t2[k2]:
								count = count + 1
								k1 = len(list1)
						if t == ">":
							if prop < t2[k2]:
								count = count + 1
								k1 = len(list1)
						if t == "!=":
							if prop == t2[k2]:
								count = count + 1
								k1 = len(list1)
						if t == "<>":
							if prop < t2[k2] and prop > t2[k2+1]:
								count = count + 1
								k1 = len(list1)
					else:
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
				g.write("{} {} {} {} {} {} {} {} {} {}\n".format(nterms,(k+2),d[1][k],d[2][k],d[3][k],d[4][k],d[7][k],d[19][k],d[47][k],d[48][k]))
				nterms = nterms + 1 	

			k = k + 1

	g.close()

main_func()
			
	
	 
	












