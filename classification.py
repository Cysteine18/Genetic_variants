# VARIOUS CLASSIFIERS FOR THE PDB FILES

import gzip
import sys
from Bio.PDB.PDBParser import PDBParser
parser = PDBParser(PERMISSIVE=0,QUIET=True)

pathPDB = "/bmm/data/rcsb/data/structures/all/pdb"

# LOADING THE FILE WHICH WAS USED TO DEFINE THE MUTANTS

f = open("pdb_seqres.txt","r")
ft = f.readlines()
f.close()

k = 0
gt=ft[k].split()
R = gt[1][4:len(gt[1])]
while R == "protein":
	gt=ft[k].split()
	R = gt[1][4:len(gt[1])]
	k=k+2

proindex = k

# OUTPUT

g = open("PDB_classifier.txt","w")
h = open("struct_not_found.txt","w")

k = 0
oldname = "NA"
while k < proindex:
	gt = ft[k].split()
	pdbid = gt[0]
	pdbname = pdbid[1:5]
	if pdbname != oldname:
		print("{} ({} of {})".format(pdbname,k,proindex))
		oldname = pdbname
		try:
			pdbfile = "{}/pdb{}.ent.gz".format(pathPDB,pdbname)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess.pdb","wb")
			out.write(tar.read())
			tar.close()
			out.close()

			structure_id = "{}".format(pdbname)
			filename = "pdbprocess.pdb"
			structure = parser.get_structure(structure_id,filename)
			try:
				resolution = structure.header['resolution']
			except:
				resolution = "NA"
			try:
				date = structure.header['release_date']
			except:
				date = "NA"
			try:
				structure_method = structure.header['structure_method']
			except:
				structure_method = "NA"

			g.write("{} {} {} {}".format(pdbname,resolution,date,structure_method))
			g.write("\n")
		except:
			print ("FILE NOT FOUND")
			h.write("{} ".format(pdbname))
	
	k = k + 2

g.close()
h.close()
