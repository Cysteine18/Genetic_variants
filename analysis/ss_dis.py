f = open("cluster_length.txt","r")
ft = f.readlines()
f.close()

g = open("distinct_mutants_struct.txt","r")
gt = g.readlines()
g.close()

h = open("ss_dis.txt","w")

# dssp nomenclature for various structure elements

H = 0
B = 0
E = 0
G = 0
I = 0
T = 0
S = 0

k = 0
k1 = 0
while k < len(ft):
	ft1 = ft[k].split()
	t1 = ft1[0]
	print("CLUSTER {} OF {}".format(k,len(ft)))

	gt1 = gt[k1].split(", ")
	t2 = gt1[0].strip("[|'|]")
	while t2 != t1:
		k1 = k1 + 1
		if k1 >= len(gt):
			print("ERROR")
			quit()
		gt1 = gt[k1].split(", ")
		t2 = gt1[0].strip("[|'|]")
	
	k2 = 1
	while k2 < len(gt1):
		t3 = gt1[k2].strip("[|'|]|\n")

		if t3 == "H":
			H = H + 1
		if t3 == "B":
			B = B + 1
		if t3 == "E":
			E = E + 1
		if t3 == "G":
			G = G + 1
		if t3 == "I":
			I = I + 1
		if t3 == "T":
			T = T + 1
		if t3 == "S":
			S = S + 1

		k2 = k2 + 1

	k = k + 1

h.write("H {}" .format(H))
h.write("\n")
h.write("B {}" .format(B))
h.write("\n")
h.write("E {}" .format(E))
h.write("\n")
h.write("G {}" .format(G))
h.write("\n")
h.write("I {}" .format(I))
h.write("\n")
h.write("T {}" .format(T))
h.write("\n")
h.write("S {}" .format(S))
h.write("\n")

h.close()

print("ALPHA HELIX          = {}".format(H))
print("ISOLATED BETA-BRIDE  = {}".format(B))
print("EXTEND STRAND        = {}".format(E))
print("3/10 HELIX           = {}".format(G))
print("PI HELIX             = {}".format(I))
print("HYDROGEN BONDED TURN = {}".format(T))
print("BEND                 = {}".format(S))
