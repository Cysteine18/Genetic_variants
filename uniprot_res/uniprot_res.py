#pathfiles = "/bmm/home/tkhanna1/Documents/uniprot_res/files"

#####  OUTPUTS THE HIGH RESOLUTION PDB FOR EACH UNIPROT ID IN INPUT FILE $input #####

# REQUIRES : uniprot_pdb.csv and PDB_classifier.txt FILE"]

# EXECUTION : python uniprot_res.py input.txt

# OUTPUT :: output.txt ; 2 COLUMNS - FIRST UNIPROT AND SECOND HIGH RES PDB

def uniprot(arg):

	# DICTIONARY ELEMENT FOR UNIPROT TO PDB MAPPING; KEY UNIPROT ID AND VALUE LISTS OF PDB's

	f = open("{}".format(arg),"r")
	ft = f.readlines()
	f.close()

	d = dict()
	k = 2
	while k < len(ft):
		ft1 = ft[k].split(",")
		t1 = ft1[0]

		tmp = ft1[1]
		t2 = tmp.split(";")
		if t1 not in d.keys():
			d[t1] = t2
		k = k + 1

	return(d)

def res_filter(arg):

	# DICTIONARY FOR RESOLUTION OF ALL THE PDB FILES
 
	f = open("{}".format(arg),"r")
	ft = f.readlines()
	f.close()

	k = 0
	d = dict()
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		t2 = ft1[1].strip("\n")
	
		if t1 not in d.keys():
			d[t1] = t2

		k = k + 1

	return(d)

def main():

	# MAIN FUNCTION

	import sys
	
	f = open("{}".format(sys.argv[1]),"r")		# A SIMPLE TEXT FILE WITH LIST OF UNIPROT IDS (1 COLUMN)
	ft = f.readlines()
	f.close()

	g = open("output.txt","w")	# A SIMPLE TAB-SEPARATED TEXT FILE AS OUTPUT WITH 2 COLUMNS OF UNIPROT AND LOW RES PDB

	up = uniprot("uniprot_pdb.csv")
	res = res_filter("PDB_classifier.txt")

	k = 0
	while k < len(ft):
		ft1 = ft[k].split()
		t1 = ft1[0].strip("\n")
		try:
			pdblist = up[t1]	# LIST OF PDB's FOR A UNIPROT
		except:
			pdblist = "NA"
		if pdblist != "NA":
			k1 = 0
			pdb = pdblist[0].strip("\n")
			r = 1000.0
			while k1 < len(pdblist):
				t2 = pdblist[k1].strip("\n")
				try:
					temp = res[t2]
				except:
					temp = 10000.0	

				if temp == "None":
					temp = 1000.0		# IGNORING THE NMR STRUCTURES
			
				if float(temp) < r:
					pdb = t2
					r = float(temp)				
				k1 = k1 + 1
			g.write("{}	{}\n".format(t1,pdb))
		else:
			g.write("{}	NOT_FOUND\n".format(t1))

		k = k + 1

	g.close()
			
if __name__ == "__main__":
	main()

			
