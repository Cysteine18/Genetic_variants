# THIS CODE SHOWS THE DISTRIBUTION OF VARIOUS AMINO ACIDS AMONG THE MUTANTS

f = open("number_of_distinct_mutations.csv","r")
ft = f.readlines()
f.close()

num = dict()
g = open("AA_dis.txt","w")

k = 0
while k < 9449: #len(ft):
	ft1 = ft[k].split(",")
	k1 = 1
	while k1 < len(ft1):
		t1 = ft1[k1].strip("\n")
		res1 = t1[0:1]
		res2 = t1[(len(t1)-1):len(ft)]

		if "{}".format(res1) in num.keys():
			num["{}".format(res1)].append(1)
		else:
			num["{}".format(res1)] = [1]		
		
		if "{}".format(res2) in num.keys():
			num["{}".format(res2)].append(1)
		else:
			num["{}".format(res2)] = [1]

		k1 = k1 + 1

	k = k + 1

i = 1
for x in num.keys():
	number = len(num["{}".format(x)])
	g.write("{}	{} {}".format(i,x,number))
	g.write("\n")
	i = i + 1

g.close()

