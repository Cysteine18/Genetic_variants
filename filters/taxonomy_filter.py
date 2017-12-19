# FILTER BASED ON THE TAXONOMY OF THE PDB's

import sys

inp = len(sys.argv)

# CLUSTER FILE

f = open("cluster.txt","r")
ft = f.readlines()
f.close()

# INPUT SHOULD BE A .csv FILE

var1 = dict()
var2 = dict()
for x in range(1,inp):
	file1 = "{}".format(sys.argv[x])
	lf = file1.split(".")
	var1["{}".format(x)] = []
	var2["{}".format(x)] = lf[0]

# OPENING THE CLUSTER FILE AND DETERMINING THE TAXONOMY DISTRIBUTION

k=0
cid = 1
pdb = []
while k < len(ft):
	ft1 = ft[k].split()
	if ft1[0] == "#" and ft1[1] == "{}".format(cid):
		print("APPLYING TAXONOMY FILTER TO CLUSTER {}".format(cid))
		k = k + 1
		ft1 = ft[k].split()
		while ft1[0] != "#" and k < (len(ft)-1):
			t1 = ft[k].split(", ")
			t2 = t1[0].strip("'|(")
			t3 = t2[0:4]
			t4 = t3.upper()
		
			count = 0
			k1 = 0
			while k1 < len(pdb):
				if pdb[k1] == t4:
					count = count + 1
					k1 = len(pdb)
				k1 = k1 + 1
			
			if count == 0:
				pdb.append(t4)

				# CHECKING THE INPUT TAXONOMY FILES

				k1 = 1
				while k1 < inp:
					g = open("{}".format(sys.argv[k1]),"r")
					gt = g.readlines()
					g.close()
					gt1 = gt[0].split(", ")
					k2 = 0
					while k2 < len(gt1):
						t4c = gt1[k2].strip("\n")
						if t4c == t4:
							var1["{}".format(k1)].append("{}".format(t4))
							k2 = len(gt1)
							k1 = inp
						k2 = k2 + 1
					k1 = k1 + 1

			k = k + 1

			ft1 = ft[k].split()
		cid = cid + 1
	else:
		k = k + 1

# STORING THE FILES

h = open("taxonomy.txt","w")
h1 = open("pdb_tax_div.txt","w")

for x in range(1,inp):
	v1 = var2["{}".format(x)]
	v2 = len(var1["{}".format(x)])

	h.write("{} {} {}".format(x,v1,v2))
	h1.write("{}".format(var1["{}".format(x)]))
	h.write("\n")
	h1.write("\n")

h.close()

















						
