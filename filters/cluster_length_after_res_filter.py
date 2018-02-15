# THIS CODE CHANGES THE CLUSTER LENGTH FILE FORMED FOR THE CLUSTER DEFINATION TO INCLUDE THE HIGH RESULTION STRUCTURE FOR EACH CLUSTER

f = open("cluster_length.txt","r")
ft = f.readlines()
f.close()

g = open("PDB_classifier.txt","r")
gt = g.readlines()
g.close()

h = open("seq_sim.txt","r")
ht = h.readlines()
h.close()

m = open("cluster_length_after_resfilter.txt","w")

k = 0
k1 = 0
while k < len(ft):
	print("#### CHECKING CLUSTER {} OF {} ####".format(k,len(ft)))
	ft1 = ft[k].split()
	t1 = ft1[0]
	t2 = ft1[1]
	t3 = ft1[2]
	t4 = ft1[3].strip("\n")

	ht1 = ht[k1].split(",")
	t1c = ht1[0].strip("\n")
	 
	while t1c != t1:
		k1 = k1 + 1
		if k1 >= len(ht):
			print("ERROR")
			quit()
		ht1 = ht[k1].split(",")
		t1c = ht1[0].strip("\n")

	pdb = t1
	k2 = 0
	low = 100.0
	while k2 < len(ht1):
		t1c = ht1[k2].strip("\n")
		k3 = 0
		count = 0
		gt1 = gt[k3].split()
		t1cc = gt1[0].strip("\n")
		while t1cc != t1c[0:4]:
			k3 = k3 + 1
			if k3 >= len(gt):
				count = count + 1
				break
			gt1 = gt[k3].split()
			t1cc = gt1[0].strip("\n")
			
		if count == 0:
			res = gt1[1].strip("\n")
			#print(res)
			if res != "None":
				if float(res) < low:
					low = float(res)
					pdb = t1c
		k2 = k2 + 1

	print("CHANGED {} TO {} WITH RESOLUTION {}".format(t1,pdb,res))
	m.write("{} {} {} {}".format(pdb,t2,t3,t4))
	m.write("\n")

	k = k + 1

m.close()
					
			







		
				
