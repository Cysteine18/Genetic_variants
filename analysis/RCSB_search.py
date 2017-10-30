# THIS SMALL PYTHON CODE HELPS LOOK AT THE FASTA FILE FORM RCSB AND TRY TO FIND THE RELATIONSHIP BETWEEN SEQUENCE SIMILARITY AND MUTANTS

def match(list1,list2,mv):
	# THIS FUNCTION MATCHS THE STRING 1,str1 TO STRING 2, str2 WITH THE RESOLUTION DEFINED BY mv 
	#print "*** {} {}".format(str1,str2)

	str1 = "{}".format(list1[0])
	str2 = "{}".format(list2[0])

	if str1 != str2:
		k2 = 0
		count = 0
		while k2 < len(str1):
			t1 = list1[0][k2]
			t2 = list2[0][k2]
			if t1 != t2:
				count = count + 1
				diffid = k2 + 1
				diffseq1 = t1
				diffseq2 = t2
			if count > mv:
				return[-1,"NA","NA","NA"]
				k2 = len(str1)

			k2 = k2 + 1
		
		if count <= mv:
			return[0,diffid,diffseq1,diffseq2]

	else:
		return[-1,"NA","NA","NA"]

def exact_match(list1,list2):
	str1 = "{}".format(list1[0])
	str2 = "{}".format(list2[0])

	if str1 == str2:
		return[0,"NA","NA","NA"]
	else :
		return[1,"NA","NA","NA"]


f = open("pdb_seqres.txt","r")
ft = f.readlines()
f.close()

res = 1		# SEQUENCE SIMILARITY VARIABLE

# DETERMINING THE PROTEIN DATA INDEX

k = 0
gt=ft[k].split()
R = gt[1][4:len(gt[1])]


while R == "protein":
	gt=ft[k].split()
	R = gt[1][4:len(gt[1])]
	k=k+2

proindex = k
searchindex = proindex		# VARIABLE TO CONTROL THE SEARCH - END
startindex = 0			# START

#g = open("mutants.txt","w")
g = open("seq_sim.txt","w")
#h = open("mutations.txt","w")

k = startindex
while k < searchindex:
	gt=ft[k].split()
	if gt[0][0] == ">":
		pdb = gt[0][1:len(gt[0])]
		length = gt[2][7:len(gt[2])]
		#name = gt[3]

		perdone = (100 * k) / searchindex
		print ("		CHECKING PDB {}".format(pdb))
		print ("	{} PERCENT DONE".format(perdone))

		k = k + 1
		gt=ft[k].split() 		# SEQUENCE TO BE COMPARED
		
		mutant = "{}".format(pdb)
		mutation = "{}".format(pdb)


		# MATCHING WITH THE REST OF THE FILE

		k1 = 0
		while k1 < proindex:
			ct=ft[k1].split()
			if ct[0][0] == ">" and k1 != (k - 1):
				pdbc = ct[0][1:len(ct[0])]
				lengthc = ct[2][7:len(ct[2])]
				#namec = ct[3]
			
				#print pdbc
				k1 = k1 + 1

				ct=ft[k1].split() 		# SEQUENCE TO BE COMPARED WITH
				
				# MATCHING
				if length == lengthc:
					#M = match(gt, ct, res) 
					M = exact_match(gt,ct)
					if M[0] == 0:
						#print "PDB {} AND {} ARE POTENTIAL MUTANTS".format(pdb,pdbc)
						mutant = mutant + ",{}".format(pdbc)
						#mutation = mutation + ",{}{}{}".format(M[2],M[1],M[3])
						#g.write("{} {}".format(pdb,pdbc))
					
			k1 = k1 + 1

	g.write(mutant)
	g.write("\n")
	#h.write(mutation)
	#h.write("\n")
	k = k + 1

g.close()
#h.close()

