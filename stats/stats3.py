# NUMBER OF DIFFERENT MUTATIONS BASED ON SEQUENCE SIMILARITY

f = open("seq_sim.txt","r")
ft = f.readlines()
f.close()

k=0
blacklist = []
nent = len(ft)
while k < len(ft):
	gt=ft[k].split(",")
	t1 = gt[0].strip("\n")
	print ("READING ENTRY {} OF {}".format(k,int(nent)))
	count = 0
	k1=0
	while k1 < len(blacklist):
		if t1 == blacklist[k1]:
			count = count + 1
			k1 = len(blacklist)
		k1 = k1 + 1
	if count == 0:
		k2=1
		while k2 < len(gt):
			t1 = gt[k2].strip("\n")
			blacklist.append(t1)
			k2 = k2 + 1
	k = k + 1

g = open("all_mutants.txt","r")
gt = g.readlines()
g.close()

k=0
count = 0
while k < len(gt):
	ht=gt[k].split(",")
	t1 = ht[0].strip("\n")
	print ("CHECKING ENTRY {} OF {}".format(k,int(nent)))
	k1 = 0
	tempcount = 0
	if len(ht) > 1:
		while k1 < len(blacklist):
			if t1 == blacklist[k1]:
				tempcount = tempcount + 1
				k1 = len(blacklist)
			k1 = k1 + 1

	if tempcount == 0:
		count = count + len(ht) - 1

	k = k + 1

nmut = count / 2
print (int(nmut))
