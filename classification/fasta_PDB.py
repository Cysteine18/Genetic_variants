import sys
import re
import gzip
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser(QUIET=True)
from Bio.PDB.Polypeptide import three_to_one as tto
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

def seqres_atom_map(mmcif_dict,c):


	category = "_pdbx_poly_seq_scheme"
	seq_len = len(mmcif_dict[category + ".seq_id"])
	seqres = {}
	for i in range(seq_len):
		seqres_index = mmcif_dict["_pdbx_poly_seq_scheme.seq_id"][i]
		pdb_seq_id = int(mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"][i])
		chain = mmcif_dict["_pdbx_poly_seq_scheme.pdb_strand_id"][i]
		if chain == c:
			res = mmcif_dict["_pdbx_poly_seq_scheme.pdb_mon_id"][i]
			if res == "?":
				sres = "-"
			else:
				sres = tto("{}".format(res))
			key1 = (seqres_index,chain)

			seqres[key1] = sres

	return seqres

def main():

	pathmmcif = "/Volumes/BIOINFO/mmCIF"
	pathfasta = "/Users/tarunkhanna/Documents/Bioinformatics/file_links/"

	f = open("{}/pdb_seqres.txt".format(pathfasta),"r")
	ft = f.readlines()
	f.close()

	g = open("PDB_fasta.txt","w")

	start = sys.argv[1]
	end = sys.argv[2]
	if end == "END" or end == "end":
		end = len(ft)
	end = int(end)
	start = int(start)
	k = start
	while k < end:
		print(k,len(ft))
		ft1 = ft[k].split()
		t1 = ft1[0].strip(">")
		pdb = t1[0:4]
		chain = t1[5:len(t1)]
	
		k = k + 2

		#count = 0
		#if count == 0:
		try:
			fol = pdb[1:3]		
			pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess1{}.cif".format(start),"wb")
			out.write(tar.read())
			tar.close()
			out.close()

			mmcif = MMCIF2Dict("pdbprocess1{}.cif".format(start))
			idmap1 = seqres_atom_map(mmcif,chain)
			k1 = 1
			str1 = ""
			while k1 <= len(idmap1):
				t2 = "{}".format(k1)
				key1 = (t2,chain)
				res = idmap1[key1]
				if k%100 == 0:
					str1 = str1 + "{}\n".format(res)
				else:
					str1 = str1 + "{}".format(res)
				k1 = k1 + 1
			g.write(">{}\n".format(t1))
			g.write("{}\n".format(str1))
				
		except:
			print("FILE NOT_FOUND")

				
main()







