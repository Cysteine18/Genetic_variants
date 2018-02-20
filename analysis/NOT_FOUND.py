def main_func():
	import sys

	f = open("{}".format(sys.argv[1]),"r")		# CLUSTER
	ft = f.readlines()
	f.close()

	h = open("{}".format(sys.argv[2]),"r")  # MUTATION SEQRES
	ht = h.readlines()
	h.close()

	m = open("{}".format(sys.argv[3]),"r")		# MUTANTS
	mt = m.readlines()
	m.close()

	g1 = open("NOT_FOUND.txt","w")
	g = open("need_review.txt","w")

	k = 0
	cid = 1
	while k < len(ft):
		print("{} of {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		try:
			t2 = ft1[1].strip("\n")
		except:
			t2 = "NA"
		if t1 == "#" and t2 == "{}".format(cid):
			cid = cid + 1
			k = k + 1
			ft1 = ft[k].split()
			wt = ft1[0].strip("('|',|\n")
			k = k + 1
			ft1 = ft[k].split()
			t1 = ft1[0].strip("('|',|\n")
			while t1 != "#" and (k+1) < len(ft):
				mut = t1
				t2 = ft1[(len(ft1)-1)].strip("')|'|\n|',")
				t3 = ft1[(len(ft1)-2)].strip("')|'|\n|',")
				t4 = ft1[(len(ft1)-3)].strip("')|'|\n|',")
				t5 = ft1[(len(ft1)-4)].strip("')|'|\n|',")

				if t2 == "NOT_FOUND" or t3 == "NOT_FOUND" or t4 == "NOT_FOUND" or t5 == "NOT_FOUND" or t4 == "NA" or t5 == "NA":
					k1 = 0
					mt1 = mt[k1].split(",")
					t1c = mt1[0].strip("\n")
					while t1c != wt:
						k1 = k1 + 1
						if k1 > len(mt):
							print("ERROR")
							quit()
						mt1 = mt[k1].split(",")
						t1c = mt1[0].strip("\n")
					wtp = k1

					k1 = 1
					count = 0
					while k1 < len(mt1):
						t1c = mt1[k1].strip("\n")
						#print(t1c)
						if t1c == mut:
							count = count + 1
							break
						k1 = k1 + 1
					if count != 0:
						mutp = k1
					else:
						print("MUT ERROR {}".format(mut))
						quit()

					ht1 = ht[wtp].split(",")
					cmut = ht1[mutp].strip("\n")
					mpos = cmut[1:(len(cmut)-1)]
				
					g.write("{} {} {}\n".format(wt,mut,mpos))
					g1.write("{} {} {}\n".format(wt,mut,mpos))
				else:
					k1 = 0
					mt1 = mt[k1].split(",")
					t1c = mt1[0].strip("\n")
					while t1c != wt:
						k1 = k1 + 1
						if k1 > len(mt):
							print("ERROR")
							quit()
						mt1 = mt[k1].split(",")
						t1c = mt1[0].strip("\n")

					wtp = k1

					k1 = 1
					count = 0
					while k1 < len(mt1):
						t1c = mt1[k1].strip("\n")
						#print(t1c)
						if t1c == mut:
							count = count + 1
							break
						k1 = k1 + 1
					if count != 0:
						mutp = k1
					else:
						print("MUT ERROR {}".format(mut))
						quit()

					# CHECKING

					ht1 = ht[wtp].split(",")
					cmut = ht1[mutp].strip("\n")
					rwt = cmut[0:1]
					rmut = cmut[(len(cmut)):(len(cmut)-1)]
					mpos = cmut[1:(len(cmut)-1)]
				
					crmut = t4[0:1]
					crwt = t5[0:1]
					if crwt != rwt and crmut != rmut:
						g.write("{} {} {}\n".format(wt,mut,mpos))
				
				k = k + 1
				ft1 = ft[k].split()
				t1 = ft1[0].strip("('|',|\n")
		else:
			k = k + 1

	g.close()
	g1.close()

def mutpos_dict(file1,file2,file3):
	
	f = open("{}".format(file1),"r")	# SEQRES
	ft = f.readlines()
	f.close()

	g = open("{}".format(file2),"r")	# PDB
	gt = g.readlines()
	g.close()

	m = open("{}".format(file3),"r")	# need review FILE
	mt = m.readlines()
	m.close()

	h = open("ALL_REVIEW.txt","w")

	k = 0
	nf = 0
	while k < len(mt):
		mt1 = mt[k].split()
		wt = mt1[0].strip("\n")
		mut = mt1[1].strip("\n")
		posseq = mt1[2].strip("\n")
		print("FINDING POSITION FOR MUTATION {} OF {} WITH MUTANTS {} AND {} AT POSITION {}".format(k,len(mt),wt,mut,posseq))
	
		k1 = 0
		ft1 = ft[k1].split(", ")
		t1 = ft1[0].strip("['|']|\n")
		while wt != t1:
			k1 = k1 + 1
			if k1 >= len(ft):
				print("WT ERROR")
				quit()
			ft1 = ft[k1].split(", ")
			t1 = ft1[0].strip("['|']|\n")
		wseq = k1
	
		# SEQRES FILE

		k1 = 0
		ft1 = ft[k1].split(", ")
		t1 = ft1[0].strip("['|']|\n")
		while mut != t1:
			k1 = k1 + 1
			if k1 >= len(ft):
				print("MUT ERROR")
				quit()
			ft1 = ft[k1].split(", ")
			t1 = ft1[0].strip("['|']|\n")
		mseq = k1

		k1 = 1
		ft1 = ft[wseq].split(", ")
		wpos = ft1[k1].strip("['|']|\n")
		wseq = k1
		while wpos != posseq:
			k1 = k1 + 1
			wpos = ft1[k1].strip("['|']|\n")
			wseq = k1
	
		k1 = 1
		ft2 = ft[mseq].split(", ")
		mpos = ft2[k1].strip("['|']|\n")
		mseq = k1
		while mpos != posseq:
			k1 = k1 + 1
			mpos = ft2[k1].strip("['|']|\n")
			mseq = k1

		# PDB FILE		
		count = 0

		k1 = 0
		gt1 = gt[k1].split(", ")
		t1 = gt1[0].strip("['|']|\n")
		while wt != t1:
			k1 = k1 + 1
			if k1 >= len(gt):
				print("WT ERROR PDB")
				count = count + 1
				break
			gt1 = gt[k1].split(", ")
			t1 = gt1[0].strip("['|']|\n")
		wpdb = k1
	
		k1 = 0
		gt1 = gt[k1].split(", ")
		t1 = gt1[0].strip("['|']|\n")
		while mut != t1:
			k1 = k1 + 1
			if k1 >= len(gt):
				print("MUT ERROR PDB")
				count = count + 1
				break
			gt1 = gt[k1].split(", ")
			t1 = gt1[0].strip("['|']|\n")
		mpdb = k1

		if count == 0:
			gt1 = gt[wpdb].split(", ")
			poswt = gt1[wseq].strip("['|']|\n")

			gt1 = gt[mpdb].split(", ")
			posmut = gt1[mseq].strip("['|']|\n")

			h.write("{} {} {} {}\n".format(wt,mut,poswt,posmut))
		else:
			nf = nf + 1

		k = k + 1

	print("{} MUTATIONS NOT FOUND".format(nf))

	h.close()


#main_func()
mutpos_dict("distinct_mutants_only_cluster_seqres.txt","distinct_mutants_only_cluster_PDB.txt","need_review.txt")










		
