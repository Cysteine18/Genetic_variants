# THIS CODE READS THE ENRIES.idx FILE FROM RCSB DATABASE AND LOOKS FOR KEYWORDS IN THE TITLE TO DETERMINE BOUND/UNBOUND PDB's

def HETATM_KW():

	keywords = ["BOUND","COMPLEXED","COMPLEX"]
	# IF THE KEYWORD IS "COMPLEX" THEN IT HAS TO BE FOLLOWED BY "IN"

	f = open("entries.idx","r")
	ft = f.readlines()
	f.close()

	g = open("bound.txt","w")

	d = dict()
	k = 0
	while k < len(ft):
		#print("CHECKING ROW {} OF {}".format(k,end))
		ft1 = ft[k].split()
		k1 = 0
		count = 0
		while k1 < len(ft1):
			t2 = ft1[k1].strip("\n")
			k2 = 0
			while k2 < len(keywords):
				if t2 == keywords[k2]:
					if t2 == "COMPLEX":
						if ft1[(k1-2)] == "IN":
							g.write("{} IN {}\n".format(ft1[0],t2))
							d["{}".format(ft1[0].lower())] = "BOUND"
							count = count + 1
							k1 = len(ft1)
					else:
						g.write("{} {}\n".format(ft1[0],t2))
						d["{}".format(ft1[0].lower())] = "BOUND"
						count = count + 1
						k1 = len(ft1)
				k2 = k2 + 1
			k1 = k1 + 1
		
		if count == 0:
			if len(ft1[0]) == 4:
				d["{}".format(ft1[0].lower())] = "UNBOUND"

		k = k + 1

	g.close()
	return(d)

def res_filter():
	f = open("PDB_classifier.txt","r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
	
		if t1 not in d.keys():
			d["{}".format(t1)] = t2

		k = k + 1

	return(d)

def main_func():

	import sys

	d = HETATM_KW()

	res = res_filter()

	f = open("cluster_input.txt","r")
	ft = f.readlines()
	f.close()

	h = open("seq_sim.txt","r")
	ht = h.readlines()
	h.close()

	g = open("potential_bound_unbound.txt","w")
	
	g1 = open("potential_substitutes.txt","w")

	start = "{}".format(sys.argv[1])
	end = "{}".format(sys.argv[2])
	if end == "END" or end == "end":
		end = len(ft)

	start = int(start)
	end = int(end)
	k = start

	while k < end:
		print("CHECKING {} OF {}".format(k,len(ft)))
		ft1 = ft[k].split()
		t1 = ft1[0]
		cw = ft1[1]
		t2 = ft1[2]
		cm = ft1[3].strip("\n")

		try:
			if d["{}".format(t1)] != d["{}".format(t2)]:
				g.write("{} {}\n".format(t1,t2))

				if d["{}".format(t1)] == "BOUND":
					t1c = "{}_{}".format(t1,cw)
				else:
					t1c = "{}_{}".format(t2,cm)
				# LOOKING FOR ANY SUBSITUTE FOR THE BOUND FORM
				k1 = 0
				ht1 = ht[k1].split(",")
				tc = ht1[0].strip("\n")
				count = 0
				while tc != t1c:
					if k1 >= (len(ht)-1):
						count = count + 1
						break
					k1 = k1 + 1
					ht1 = ht[k1].split(",")
					tc = ht1[0].strip("\n")

				list1 = []
				min1 = 100.0
				psub = "NF"
				r = "NF"
				if count == 0:
					k2 = 1
					while k2 < len(ht1):
						tc = ht1[k2][0:4]
						if d["{}".format(tc)] != "BOUND":
							k3 = 0
							count1 = 0
							while k3 < len(list1):
								if list1[k3] == tc:
									count1 = count1 + 1
									k3 = len(list1)
								k3 = k3 + 1
							if count1 == 0:
								r = res["{}".format(tc)]
								if r != "None":
									if float(r) < min1:
										min1 = float(r)
										psub = ht1[k2].strip("\n")
								list1.append("{}".format(tc))
						k2 = k2 + 1
					r1 = res["{}".format(t1c[0:4])]
					g1.write("{} {} {} {}\n".format(t1c,psub,r1,min1))
				#else:
					#g1.write("{} NOT_FOUND\n".format(t1c))

		except:
			print("ERROR WITH {} {}".format(t1,t2))
		k = k + 1
		
	g.close()


main_func()
