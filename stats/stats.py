# TO COUNT THE NUMBER OF DIFFERENT UNIQUE MUTANTS AND MUTATIONS PRODUCES A FILE NAMED STATS.txt WHICH CONTAINS THE STATS FOR DISTINCT MUTANTS AND MUTATIONS

f = open("all_mutants.txt","r")
ft = f.readlines()
f.close()

g = open("all_mutations.txt","r")
gt = g.readlines()
g.close()

dis = open("distinct_mutants","w")
k = 0
while k < len(ft):
	dis_mut = []
	mu=ft[k].split(',')
	mut=gt[k].split(',')

	if len(mu) > 1:			
		pdb = mu[0].strip("\n")
		dis_mut.append(pdb)
		k1 = 1
		while k1 < len(mut):
			dumstr2 = mut[k1].strip("\n")
			mut_res=dumstr2[1:len(dumstr2)-1]
			k2 = 1
			count = 0
			while k2 < len(dis_mut):
				if dis_mut[k2] == mut_res:
					count = count + 1
				k2 = k2 + 1
			if count == 0:
				dis_mut.append(mut_res)
			k1 = k1 + 1
		dis.write('{}'.format(dis_mut))
		dis.write('\n')
	k = k + 1
				
dis.close()

# CHECKING THE MULTIPLICITY

print ("CHECKING THE MULTIPLICITY")

dis = open("distinct_mutants","r")
ht = dis.readlines()
dis.close()

s = open("stats.txt","w")
k = 0
while k < len(ht):
	print (k)
	dis_list = []
	l = ht[k].split()
	res = l[0].strip('[|,|]|\'')
	dis_list.append(res)
	k1 = 1
	while k1 < len(l):
		resid = l[1].strip('[|,|]|\'')
		dis_list.append(resid)
		count = 0
		k2 = 0
		while k2 < len(gt):
			mut = gt[k2].split(",")
			if mut[0] == res:
				k3 = 1
				while k3 < len(mut):
					dumstr2 = mut[k3].strip("\n")
					mut_res=dumstr2[1:len(dumstr2)-1]
					if mut_res == resid:
						count = count + 1
					k3 = k3 + 1
				dis_list.append(count)
				k1 = len(gt)
			k2 = k2 + 1
		k1 = k1 + 1
	s.write("{}".format(dis_list))
	k = k + 1

s.close()
	






















