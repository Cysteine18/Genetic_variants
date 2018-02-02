# THIS CODE CHECKS IF ANY OF THE MUTATIONS I PICKED UP ARE DUE TO DISORDERS

import sys
import re
import gzip
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser(QUIET=True)
#from Bio.PDB.PDBParser import PDBParser
#parser = PDBParser(PERMISSIVE=0,QUIET=True)

pathmmcif = "/bmm/data/rcsb/data/structures/all/mmCIF"

dis = open("distinct_mutants_only_cluster.txt","r")
ht = dis.readlines()
dis.close()

g = open("{}".format(sys.argv[3]),"w")

start = sys.argv[1]
end = int(sys.argv[2])
k=int(start)
while k < end: # end = len(ht)
	mutant = []
	mu=ht[k].split(',')

	pdbid=mu[0].strip('[|\,|\'|]')
	pdb=pdbid[0:4]		# PDB NAME
	C=pdbid[5:6]		# CHAIN
		
	print("*** {} :: {} of {} ***" .format(pdb,k,len(ht)))

	count = 0
	if count == 0:
	#try:
		pdbfile = "{}/{}.cif.gz".format(pathmmcif,pdb)
		tar = gzip.open("{}".format(pdbfile),"rb")
		out = open("pdbprocess{}.cif".format(start),"wb")
		out.write(tar.read())
		tar.close()
		out.close()

		structure_id = "{}".format(pdb)
		filename = "pdbprocess{}.cif".format(start)
		structure = parser.get_structure(structure_id,filename)

		model = structure[0]

		chain = model["{}".format(C)]
		c1 = chain.get_list()		# LIST ALL THE RESIDUES
		k2 = 1
		while k2 < len(mu):
			pos= mu[k2].strip(' |[|,|]|\'|\n')
			pos = int(pos)
			k1 = 0
			while k1 < len(c1):
				c2 = c1[k1].get_id()
				resid = c2[1]
				if c2[0] == " " and resid == pos:
					residue = chain[c2]
					check = residue.is_disordered()
					if check == 2:
						print("{} IS DISORDED AT POS {}".format(pdbid,pos))
						g.write("{} {}".format(pdbid,pos))
						g.write("\n")
					k1 = len(c1)
				k1 = k1 + 1
			k2 = k2 + 1
	#except:
		#print("FILE NOT FOUND")

	k = k + 1

g.close()
	
	
