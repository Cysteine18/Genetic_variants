import sys
import math

f = open("{}".format(sys.argv[1]),"r")
ft = f.readlines()
f.close()

g = open("{}.txt".format(sys.argv[2]),"w")

g1 = open("NUM_{}.txt".format(sys.argv[2]),"w")

d = dict()
d1 = dict()
list1 = ["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V"]
k = 0
while k < len(list1):
	t1 = list1[k]
	k1 = 0
	while k1 < len(list1):
		t2 = list1[k1]
		key1 = (t1,t2)
		d[key1] = []
		d1[key1] = []
		k1 = k1 + 1
	k = k + 1

k = 1
while k < len(ft):
	ft1 = ft[k].split(",")

	t1 = ft1[4]
	t2 = ft1[5]
	t3 = ft1[49]
	t4 = ft1[50]

	key1 = (t1,t2)

	if key1 in d.keys() and t1 != t2:
		d[key1].append(t3)
		d1[key1].append(t4)
	else:
		print("{} ROW = {} COMBINATION NOT FOUND".format(key1,k+1))
	
	k = k + 1

# MEAN VALUE FOR EACH PAIR

for x in d.keys():

	# LOCAL RMSD

	if len(d[x]) > 0:
		k = 0
		mean = 0.0
		while k < len(d[x]):
			mean = float(d[x][k]) + mean
			k = k + 1
		mean = mean / len(d[x])
		mean = round(mean,2)

		SD = 0.0
		k = 0
		while k < len(d[x]):
			temp = mean - float(d[x][k])
			temp = temp * temp
			SD = temp + SD
			k = k + 1

		SD = SD / len(d[x])
		SD = math.sqrt(SD)
		SD = round(SD,2)

		g.write("{}	{}	{}	{}".format(x[0],x[1],mean,SD))

		# LOCAL RMSD SIDE CHAIN
		k = 0
		mean = 0.0
		while k < len(d1[x]):
			mean = float(d1[x][k]) + mean
			k = k + 1
		mean = mean / len(d1[x])
		mean = round(mean,2)

		SD = 0.0
		k = 0
		while k < len(d1[x]):
			temp = mean - float(d1[x][k])
			temp = temp * temp
			SD = temp + SD
			k = k + 1

		SD = SD / len(d1[x])
		SD = math.sqrt(SD)
		SD = round(SD,2)

		g.write("	{}	{}\n".format(mean,SD))
		g1.write("{}	{}	{}\n".format(x[0],x[1],len(d[x])))

	else:
		g.write("{}	{}	0.0	0.0	0.0	0.0\n".format(x[0],x[1]))
		g1.write("{}	{}	{}\n".format(x[0],x[1],len(d[x])))

g.close()

	








