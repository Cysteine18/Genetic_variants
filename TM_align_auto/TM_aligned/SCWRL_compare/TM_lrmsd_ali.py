import sys
import os
import re
import gzip
from Bio.PDB.MMCIFParser import MMCIFParser
parser = MMCIFParser(QUIET=True)
from Bio.PDB.Polypeptide import three_to_one as tto
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

pathTMalign = "/Users/tarunkhanna/Documents/Bioinformatics/TM_align"
pathSCWRL = "/Volumes/BIOINFO/SCWRL_PREDICTION"
pathmmcif = "/Volumes/BIOINFO/mmCIF"

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
				try:
					sres = tto("{}".format(res))
				except:
					sres = "X"
			key1 = (seqres_index,chain)

			seqres[key1] = sres

	return seqres

def fasta(t1):

	pdb = t1[0:4]
	chain = t1[5:len(t1)]

	g = open("test.txt","w")

	#count = 0
	#if count == 0:
	try:
		fol = pdb[1:3]		
		pdbfile = "{}/{}/{}.cif.gz".format(pathmmcif,fol,pdb)
		tar = gzip.open("{}".format(pdbfile),"rb")
		out = open("pdbprocess1.cif","wb")
		out.write(tar.read())
		tar.close()
		out.close()

		mmcif = MMCIF2Dict("pdbprocess1.cif")
		g.write("{}".format(mmcif))
		idmap1 = seqres_atom_map(mmcif,chain)
		k1 = 1
		str1 = ""
		while k1 <= len(idmap1):
			t2 = "{}".format(k1)
			key1 = (t2,chain)
			res = idmap1[key1]
			if k1%100 == 0:
				str1 = str1 + "{}\n".format(res)
			else:
				str1 = str1 + "{}".format(res)
			k1 = k1 + 1

		return(str1)
				
	except:
		return("NA")

def fasta_seqres():
	
	pathfasta = "/Users/tarunkhanna/Documents/Bioinformatics/file_links/"
	
	f = open("{}/pdb_seqres.txt".format(pathfasta),"r")
	ft = f.readlines()
	f.close()

	d = dict()
	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip(">")
		k = k + 1
		ft1 = ft[k].split()
		t2 = ft1[0]
		str1 = ""
		if len(t2) > 100:
			num = len(t2) / 100
			num = int(num)
			rem = len(t2) - (num*100)
			for x in range(0,num):
				temp = t2[(x*100):((x+1)*100)]
				str1 = str1 + "{}\n".format(temp)
			temp = t2[((x+1)*100):(((x+1)*100)+rem)]
			str1 = str1 + "{}".format(temp)
		k = k + 1

		if t1 not in d.keys():
			d[t1] = str1

	return(d)
	
def main(): 

	import subprocess

	os.system("mkdir -p TM_align")

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	m = open("RMSD_TM_align1.csv","w")
	m.write("WT,WT_CHAIN,MUT,MUT_CHAIN,WT_LEN_PDB,MUT_LEN_PDB,ALIGNED_LENGTH,GRMSD,LRMSD_CA,LRMSD_SC,CHI1,CHI2,CHI3,CHI4,CHI5\n")

	start = int(sys.argv[2])
	end = sys.argv[3]
	if end == "END" or end == "end":
		end = len(ft)
	end = int(end)
	
	k = start

	while k < end:
		print(k,end)
		ft1 = ft[k].split()
		t1 = ft1[0]
		t2 = ft1[1]
		t3 = ft1[2]
		t4 = ft1[3]
		t5 = ft1[4].strip("\n")

		# JUST FOR SCWRL / FOLDX COMPARISON
		te = t3.split("-")
		te1 = te[0][0:4]
		te2 = te[0][4:len(te[0])]

		fasta1 = fasta("{}_{}".format(t1,t2))
		fasta2 = fasta("{}_{}".format(te1,te2))

		#fasta2 = fasta("{}_{}".format(t3,t4))

		h = open("aligned.txt","w")
		h.write("> chain1\n")
		h.write("{}\n".format(fasta1))
		h.write("> chain2\n")
		h.write("{}\n".format(fasta2))
		h.close()

		#count = 0
		#if count == 0:
		try:
			#subprocess.call(['python', "script2.py", "0", "{}".format(t1), "{}".format(t2), "{}".format(t3), "{}".format(t4)], stdout=FNULL, stderr=FERR, close_fds=True)

			#subprocess.call(['{}/./TMalign'.format(pathTMalign), "chain1.pdb", "chain2.pdb", "-o", "output", "-I", "aligned.txt", "|", "tee", "temp"], stdout=FNULL, stderr=FERR, close_fds=True)

			os.system("python script1.py 0 {} {}".format(t1,t2))
			os.system("cp {}/{}.pdb ./chain2.pdb".format(pathSCWRL,t3))
			os.system("{}/./TMalign chain1.pdb chain2.pdb -o output -I aligned.txt > temp".format(pathTMalign))
			os.system("python script1.py 2 {}".format(t5))

			#subprocess.call(['python', "script1.py", "2", "{}".format(t5)], stdout=FNULL, stderr=FERR, close_fds=True)

			g = open("results","r")
			gt = g.readlines()
			g.close()

			gt1 = gt[0].split()
			grmsd = gt1[0]
			lwt = gt1[1]
			lmut = gt1[2]
			coverage = gt1[3]
			gt1 = gt[1].split()
			lrmsd1 = gt1[0]
			lrmsd2 = gt1[1]
			chi1 = gt1[2]
			chi2 = gt1[3]
			chi3 = gt1[4]
			chi4 = gt1[5]
			chi5 = gt1[6]

			m.write("{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}\n".format(t1,t2,t3,t4,lwt,lmut,coverage,grmsd,lrmsd1,lrmsd2,chi1,chi2,chi3,chi4,chi5))

			os.system("mv file2.pdb ./TM_align/{}{}-{}{}.pdb".format(t1,t2,t3,t4))
			os.system("mv temp ./TM_align/{}{}-{}{}.log".format(t1,t2,t3,t4))
		except:
			print("ERROR")

		k = k + 1

	delete()

	m.close()

def delete():
	try:
		os.remove("output")
		os.remove("output_all")
		os.remove("output_atm")
		os.remove("output_all_atm")
		os.remove("output_all_atm_lig")
		os.remove("chain1.pdb")
		os.remove("chain2.pdb")
		os.remove("results")
		os.remove("file1.pdb")
		os.remove("file2.pdb")
		os.remove("temp")
		os.remove("pdbprocess1.cif")
		os.remove("pdbprocess.cif")
	except:
		 k = 0

if __name__ == "__main__":
	main()


