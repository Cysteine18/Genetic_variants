# FILTER BASED ON THE LENGTH OF THE CLUSTER

f = open("cluster_length.txt","r")
ft = f.readlines()
f.close()

g = open("cluster_sizes.txt","w")

# determining the maximum and the minimum cluster size
k = 0
mincs = 1000
maxcs = 0
x = 0		# total number of pdb chains represnted
while k < len(ft):
	ft1 = ft[k].split()
	t1 = int(ft1[2])
	t2 = int(ft1[1])
	x = x + t2
	g.write("{}".format(t1))
	g.write("\n")
	
	if t1 < mincs:
		mincs = t1
		minc = ft1[0]
	if t1 > maxcs:
		maxcs = t1
		maxc = ft1[0]

	k = k + 1

sf = 20
k = 0
count = 0
y = 0		# total number of pdb chains below sf
while k < len(ft):
	ft1 = ft[k].split()
	t1 = int(ft1[2])
	t2 = int(ft1[1])
	if t1 < sf:
		y = y + t2
		count = count + 1
	k = k + 1

print("{} CLUSTERS OUT OF {} CLUSTERS ARE BELOW {} RESIDUES".format(count,len(ft),sf))
print("THIS REPRESENT {} PDBs OUT OF {}".format(y,x))
g.close()
