# THIS CODE CREATES THE INITIAL FILES FOR GENERATING THE CLUSTER FILE

f = open("ALL_REVIEW.txt","r")
ft = f.readlines()
f.close()

g = open("number_of_distinct_mutants.csv","w")

g1 = open("number_of_distinct_mutations.csv","w")

d = dict()
d1 = dict()

k = 0
while k < len(ft):
	ft1 = ft[k].split()
	t1 = ft1[0]
	t2 = ft1[1]
	t3 = ft1[2]
	t4 = ft1[3].strip("\n")
	t5 = "{}-{}".format(t3,t4)

	if t1 in d.keys():
		d["{}".format(t1)].append(t2)
		d1["{}".format(t1)].append(t5)
	else:
		d["{}".format(t1)] = [t2]
		d1["{}".format(t1)] = [t5]

	k = k + 1

for x in d.keys():
	k = 0
	str1 = "{}".format(x)
	str2 = "{}".format(x)
	while k < len(d["{}".format(x)]):
		t1 = ",{}".format(d["{}".format(x)][k])
		t2 = ",{}".format(d1["{}".format(x)][k])

		str1 = str1 + t1
		str2 = str2 + t2

		k = k + 1

	g.write("{}\n".format(str1))
	g1.write("{}\n".format(str2))

g.close()
g1.close()

		
