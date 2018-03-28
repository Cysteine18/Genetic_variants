import sys
import math

f = open("{}".format(sys.argv[1]),"r")		# RMSD SCWRL
ft = f.readlines()
f.close()

h = open("SCWRL_compare.csv","w")
h.write("#CRY_MUT,CHAIN,PRED_MUT,CHAIN,LRMSD_SC(0),LRMSD_SC(4),LRMSD_SC(5),LRMSD_SC(6), LRMSD_SC(7),LRMSD_SC(8),LRMSD_SC(9),LRMSD_SC(10)\n")

h1 = open("mean_value.txt","w")
d = dict()
k = 0
while k < len(ft):
	print("{} of {}".format(k,len(ft)))
	ft1 = ft[k].split(",")
	wt = ft1[0]
	wtc = ft1[1]
	mut = ft1[2]
	mutc = ft1[3]
	lenwt = ft1[4]
	lenmut = ft1[5]
	allen = ft1[6]

	temp = mut.split("-")
	wt1 = temp[0][0:4]
	wtc1 = temp[0][4:(len(temp[0]))]
	mut1 = temp[1][0:4]
	mutc1 = temp[1][4:(len(temp[0]))]

	list1 = []
	for x in range(0,8):
		ft1 = ft[(k+x)].split(",")
		t1 = ft1[9]
		if t1 != "NA":
			if x in d.keys():
				d[x].append(t1)
			else:
				d[x] = [t1]
			list1.append(t1)
		if len(list1) == 8:
			h.write("{},{},{}-{},{},{},{},{},{},{},{},{},{}\n".format(wt,wtc,wt1,mut1,mutc,list1[0],list1[1],list1[2],list1[3],list1[4],list1[5], list1[6],list1[7]))
	
	k = k + 8
	

# MEAN FOR EACH ZONE

list1 = [1,4,5,6,7,8,9,10]
for x in range(0,8):
	k = 0
	mean = 0.0
	while k < len(d[x]):
		t1 = float(d[x][k])
		mean = mean + t1
		k = k + 1
	mean = mean / len(d[x])
	SD = 0.0
	k = 0
	while k < len(d[x]):
		t1 = float(d[x][k])
		temp = t1 - mean
		temp = temp * temp
		SD = SD + temp
		k = k + 1
	SD = SD / len(d[x])
	SD = math.sqrt(SD)
	mean = round(mean,2)
	SD = round(SD,2)
	h1.write("{}	{}	{}\n".format(list1[x],mean,SD))

h1.close()
h.close()

