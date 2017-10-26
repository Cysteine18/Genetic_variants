# TO REDUCE THE NUMBER OF ENTRIES FOR STRUCTURAL ANALYSIS AND TO MAKE THE BEST CHOICE FOR THE 100 SEQUENCIALLY IDENTICAL STRUCTURES

def res_classifier():

  # CLASSIFICATION BASED ON THE RESOLUTION OF THE PDB FILE

	f=open("seq_sim.txt","r")
	ft=f.readlines()
	f.close()

	g=open("PDB_classifier.txt","r")
	gt=g.readlines()
	g.close()

	h=open("unique_PDB_RES.txt","w")

	k=0
	resdone = []
	while k < len(ft):
		seq=ft[k].split(",")

		# QUICK CHECK IF THE RESIDUE IS ALREADY COVERED

		k1 = 0
		while k1 < len(seq):
			lpdb=seq[k1].strip("\n")
			k2 = 0
			count = 0
			while k2 < len(resdone):
				if resdone[k2] == lpdb:
					count = count + 1
					k2 = len(resdone)
				k2 = k2 + 1
			k1 = k1 + 1

		if count == 0:
			k1 = 0
			low=100.0	
			lpdb=seq[0].strip("\n")
			print("{} ({} OF {})".format(lpdb,k,len(ft)))
			while k1 < len(seq):
				pdb = seq[k1][0:4]
				k2 = 0
				count = 0
				while k2 < len(gt):
					cseq=gt[k2].split()
					cpdb = cseq[0].strip("\n")
					if cpdb == pdb:
						count = count + 1
						if cseq[1] != "None":
							if float(low) > float(cseq[1]):
								low = float(cseq[1])
								lpdb = seq[k1].strip("\n")
					k2 = k2 + 1
				k1 = k1 + 1
			if count != 0:
				resdone.append(lpdb)
				h.write("{}".format(lpdb))
				h.write("\n")
		k = k + 1
	h.close()

res_classifier()
	

		


