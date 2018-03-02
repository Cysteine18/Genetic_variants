import sys

pathentries = "/Users/tarunkhanna/Documents/Bioinformatics/file_links"

f = open("{}/entries.idx".format(pathentries),"r")
ft = f.readlines()
f.close()

inp = sys.argv[1]
inp1 = inp.lower()
inp2 = inp.upper()

k = 0
while k < len(ft):
	ft1 = ft[k].split("\t")
	t1 = ft1[0]

	if t1 == inp1 or t1 == inp2:
		print(len(ft1))
		t2 = ft1[3]
		t3 = ft1[4]
		print("**\n")
		print(t2,t3)
		print("\n**\n")

	k = k + 1

