import sys
import os

pathTMalign = "/Users/tarunkhanna/Documents/Bioinformatics/TM_align"

def fasta():
	
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
		#print(t2)
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

	fa = fasta()

	os.system("mkdir -p TM_align")

	f = open("{}".format(sys.argv[1]),"r")
	ft = f.readlines()
	f.close()

	m = open("RMSD_TM_align.csv","w")
	m.write("WT,WT_CHAIN,MUT,MUT_CHAIN,WT_LEN_PDB,MUT_LEN_PDB,ALIGNED_LENGTH,GRMSD,LRMSD_CA,LRMSD_SC\n")

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

		fasta1 = fa["{}_{}".format(t1,t2)]
		try:
			fasta2 = fa["{}_{}".format(t3,t4)]
		except:
			fasta2 = fasta1

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

			os.system("python script1.py 0 {} {} {} {}".format(t1,t2,t3,t4))
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
	except:
		 k = 0

if __name__ == "__main__":
	main()


