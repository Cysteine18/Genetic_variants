# STEP 2: GOING FROM SEQUENCE TO STRUCTURE

import tarfile
import gzip
import zipfile
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=0)

# PATH FOR THE PDB/mmCIF FILES

pathmmcif = "/bmm/data/rcsb/data/structures/all/mmCIF"
pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"

f = open("mutants.txt","r")
ft = f.readlines()
f.close()

g = open("mutations.txt","r")
gt = g.readlines()
g.close()

h = open("structure","w")

# DUMMY VARIABLE TO REMOVE REPETITION IN THE CALCULATIONS

mutant = ""

k=1
while k < len(f):
	mu=ft[k].split(',')
	mut=gt[k].split(',')
	dumstr = mu
	
	k1=1
	while k1 < len(dumstr):
		pdb=mu[k1][0:4]
		chain=mu[k1][4:5]
		dumstr2 = mut[k1]
		mut_res=dumstr2[1:len(dumstr2)-1]

		# PREVENTING FROM PROCESSING THE PDB AGAIN

		k2=0
		count=0
		while k2 < len(mutant):
			if mutant[k2] == pdb:
				count = count + 1
				pos = k2
			k2=k2+1

		if count == 0:

			# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

			pdbfile = "{}/pdb{}.ent.gz".format(pathPDB,pdb)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess.pdb","wb")
			out.write(tar.read())
			tar.close()
			out.close()

			mutant = mutant + pdb

				






















