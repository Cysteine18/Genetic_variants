# THIS CODE TRIES TO FIND THE MUTATIONS WHICH CAUSES LARGE SCALE CHANGE IN THE 10A ZONE AROUND THE MUTATION SITE

# FILTERING VARIABLE
fil = 20

# VOLUME OF THE 10A SPHERE

#vol = (4 * 3.14 * 1.0 * 1.0 * 1.0)/3.0 # IN nm^3
vol = 1

f = open("cluster.txt","r")
ft = f.readlines()
f.close()

g = open("del_mut.txt","w")
g1 = open("ins_mut.txt","w")
g2 = open("del_ins_mut.txt","w")
g3 = open("c_alpha_density.txt","w")

# DETERMING THE NUMBER OF DISTINCT CLUSTERS

k = 0
cid = 1
nc = 0
i = 1
while k < len(ft):
	ft1 = ft[k].split()
	if ft1[0] == "#" and ft1[1] == "{}".format(cid):
		print("CLUSTER {}".format(cid))
		cid = cid + 1
		fcheck = ft1[4].strip("'|,")
		if int(fcheck) > fil:
			k = k + 1
			ft1 = ft[k].split()
			nwt = ft1[0].strip("(|'|,")
			k = k + 1
			ft1 = ft[k].split()
			while ft1[0] != "#" and k < (len(ft) - 1):
				nmut = ft1[0].strip("(|'|,")
				wt = ft1[(len(ft1)-2)].strip("'|',")
				mut = ft1[(len(ft1)-1)].strip("'|\n|)")
				# DIVIDING THE MUTATIONS INTO THREE SETS OF DELETION, INSERTION AND DELETION-INSERTION MUTATIONS
				wt1 = wt.split(",")
				mut1 = mut.split(",")
				if len(wt1) > 1 and len(mut1) > 1:
					wt1 = wt.split(",")
					mut1 = mut.split(",")
					k1 = 0
					list1 = []
					while k1 < len(mut1):
						t1 = mut1[k1]
						k2 = 0
						count = 0
						while k2 < len(wt1):
							t2 = wt1[k2]
							if t1 == t2:
								count = count + 1
								k2 = len(wt1)
							k2 = k2 + 1
						if count == 0:
							list1.append("{}".format(t1))
						k1 = k1 + 1
					if len(mut1) == len(wt1):
						if wt1 != mut1:
							g1.write("{} {} {} {} {}".format(nwt,nmut,len(list1),len(wt1),len(mut1)))
							g1.write("\n")
					else:
						if len(list1) > 0:
							g2.write("{} {} {} {} {}".format(nwt,nmut,len(list1),len(wt1),len(mut1)))
							g2.write("\n")
						else:
							g.write("{} {} {} {}".format(nwt,nmut,len(wt1),len(mut1)))
							g.write("\n")
				if len(wt1) > 1 and len(mut1) > 1:
					pdwt = round((len(wt1)/vol),3)
					pdmut = round((len(mut1)/vol),3)
					diff = abs(pdmut - pdwt)
					g3.write("{} {} {} {} {} {}".format(i,nwt,nmut,pdwt,pdmut,diff))
					g3.write("\n")
					i = i + 1
				k = k + 1
				ft1 = ft[k].split()
	else:
		k = k + 1

g.close()
g1.close()
g2.close()
g3.close()

