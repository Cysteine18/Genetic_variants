# THIS CODE MAPS THE MUTATED RESIDUE ONTO THE CRYSTAL STRUCTURE

f = open("distinct_mutants_seqres.txt","r")
ft = f.readlines()
f.close()

h = open("ss.txt","r") # 76 column file
ht = h.readlines()
h.close

g = open("distinct_mutants_struct.txt","w")

k = 0
while k < len(ft):
	struct = ""
	ft1 = ft[k].split(", ")
	t1 = ft1[0].strip("'|\n|[|]")
	pdb = t1[0:4]
	pdb = pdb.upper()
	chain = t1[5:len(ft1[0])]

	struct = ["{}".format(t1)]

	print("{} :: {} of {}".format(pdb,k,len(ft)))

	# DETERMINING THE WILDTYPE

	k2 = 0
	ht1 = ht[k2].split(":")
	count = 0
	while ht1[0].strip(">") != pdb or ht1[1] != chain or ht1[2].strip("\n") != "secstr":
		k2 = k2 + 1
		if k2 >= len(ht):
			print("NOT FOUND")
			count = count + 1
			break
		ht1 = ht[k2].split(":")

	if count == 0:
		k1 = 1
		while k1 < len(ft1):
			pos = k2
			mutid = ft1[k1].strip("'|\n|[|]")
			mutid = int(mutid)

			i = mutid / 76
			pos = k2 + 1 + int(i)
			sc = ht[pos].strip("\n")
			start = mutid - 1 - (76*int(i))
			end = mutid - (76*int(i))
			s = sc[start:end]
			struct.append("{}".format(s))

			k1 = k1 + 1
	else:
		struct.append("DATA_NOT_FOUND")

	g.write("{}".format(struct))
	g.write("\n")
	k = k + 1

g.close()
		
		
