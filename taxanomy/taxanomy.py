# HUMAN AND NON_HUMAN MUTANTS

def distinct_mutants ():
	f = open("humans_all.txt","r")
	ft = f.readlines()
	f.close()

	f1 = open("number_of_distinct_mutants.txt","r")
	ft1 = f1.readlines()
	f1.close()

	g = open("humans_mutants.txt","w")

	h = open ("non_human_mutants.txt","w")

	k = 0
	count = 0
	while k < len(ft1):
		ft11 = ft1[k].split()
		t1 = ft11[0].strip("{")
		t1 = t1[0:4]
		t1 = t1.upper()
		print("{} :: {} of {}".format(t1,k,len(ft1)))
		k2 = 0
		t2 = ft[0].split(", ")
		while k2 < len(t2):
			count1 = 0
			if t1 == t2[k2]:
				count = count + 1
				count1 = count1 + 1
				k2 = len(t2)
			k2 = k2 + 1
		if count1 != 0:
			g.write("{},".format(t1))
		else:
			h.write("{},".format(t1))
		k = k + 1

	print ("HUMANS::{}".format(count))	
	g.close()
	h.close()

def unique_distinct_mutants (var1):
	f = open("{}".format(var1),"r")
	ft = f.readlines()
	f.close()

	g = open("unique_{}".format(var1),"w")

	t1 = ft[0].split(",")
	k = 0
	list1 = []
	while k < len(t1):
		k2 = 0
		count = 0
		while k2 < len(list1):
			if list1[k2] == t1[k]:
				count = count + 1
			k2 = k2 + 1
		if count == 0:
			list1.append(t1[k])
			if len(list1) > 1:
				g.write("\n")
			g.write("{}".format(t1[k].lower()))

		k = k + 1

	g.close()

unique_distinct_mutants ("humans_mutants.txt")
