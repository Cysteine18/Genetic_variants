# THIS CODE SHOWS THE DISTRIBUTION OF THE MUTATION ALONG THE PROTEIN LENGTH
# THIS CODE TAKS "CLUSTER LENGTH" AND "DISTINCT MUTATIONS" AS INPUT

def frange(start,stop,step):
	i = start
	while i < stop:
		yield i
		i = i + step
		i = round(i,2)

def function1():
	import sys

	# CLUSTER LENGTH FILTER
	lenfilter = 20

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	g = open("{}".format(sys.argv[2]),"r")
	gt = g.readlines()
	g.close()

	h = open("mut_coverage.txt","w")

	k = 0
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = ft1[0].strip("\n")

		k1 = 0
		gt1 = gt[k1].split()
		t2 = gt1[0].strip("\n")
		print(k,t1)
		while t1 != t2:
			k1 = k1 + 1
			if k1 >= len(gt):
				print("ERROR")
				quit()
			gt1 = gt[k1].split()
			t2 = gt1[0].strip("\n")
	
		length = gt1[2].strip("\n")

		if int(length) > lenfilter:
			k1 = 1
			while k1 < len(ft1):
				t3 = ft1[k1].strip("\n")
				resid = t3[1:(len(t3)-1)]
		
				pos = int(resid)/int(length)
				pos = round(pos,3)

				h.write("{}	{}".format(t1,pos))
				h.write("\n")
				k1 = k1 + 1

		k = k + 1

	h.close()

def function2():

	# BIN SIZE
	binsize = 20.0
	

	h = open("mut_coverage.txt","r")
	ht = h.readlines()
	h.close()

	m = open("coverage_stats.txt","w")

	step = 1.0/binsize
	step = round(step,2)
	for x in frange(0,1,step):
		k = 0
		count = 0
		while k < len(ht):
			ht1 = ht[k].split()
			t1 = ht1[1].strip("\n")
			if float(t1) > x and float(t1) < (x+step):
				count = count + 1
			k = k + 1
		m.write("{} {}".format(x,count))
		m.write("\n")

	m.close()
			

function2()

	
	












