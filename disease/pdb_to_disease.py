# MAPPING THE DISEASE RELATED UNIPROT ID TO THE PDB FILE AND THE SPECIFIC MUTATION

from Bio.PDB.Polypeptide import three_to_one as tto

f = open("omim_sprot.csv","r")
ft = f.readlines()
f.close()

g = open("uniprot_pdb.csv","r")
gt = g.readlines()
g.close()

h = open("disease_causing_mutations.txt","w")

k = 0
while k < len(ft):
	ft1 = ft[k].split(",")
	uid = ft1[2]
	res1 = ft1[3]
	res1 = tto("{}".format(res1))
	res2 = ft1[5]
	res2 = tto("{}".format(res2))
	pos = ft1[7]
	dis = ft1[8].strip("\n")

	print("{} of {}".format(k,len(ft)))

	k1 = 2
	gt1 = gt[k1].split(",")
	count = 0
	while gt1[0] != uid:
		if k1 >= len(gt):
			pdb = "NA"
			count = count + 1
			break
		gt1 = gt[k1].split(",")
		k1 = k1 + 1
	if count == 0:
		pdb = gt1[1].strip("\n")
	
	h.write("{},{},{}{}{},{}".format(uid,pdb,res1,pos,res2,dis))
	h.write("\n")

	k = k + 1

h.close()
