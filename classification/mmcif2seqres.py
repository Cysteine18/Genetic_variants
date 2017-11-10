
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import gzip

def seqres_atom_map(mmcif_dict):

    # FUNCTION FROM STEPHAN
    """
    The mapping between the canonical sequence of a protein structure
    (the ``SEQRES`` field) is given in the ``_pdbx_poly_seq_scheme``
    category of the mmCIF dictionary.
    
    This function parses the dictionary of mmCIF data returned by BioPython's
    MMCIF2Dict class, and returns a dictionary mapping the indices of the
    ``SEQRES`` sequence to the PDB ID of that residue. A mapping is returned
    rather than a list because the PDB structure may be missing some residues
    that are defined in the ``SEQRES``, so the list may not be continuous. The
    PDB residue IDs are given as a 2-tuple of residue ID and insertion code (or
    `None` if there is no insertion code).
    
    Example output:
    
        {
            1: (0, None),
            2: (1, None),
            # ...
        }
    """
    category = "_pdbx_poly_seq_scheme"
    seq_len = len(mmcif_dict[category + ".seq_id"])
    seqres = {}
    for i in range(seq_len):
        seqres_index = int(mmcif_dict["_pdbx_poly_seq_scheme.seq_id"][i])
        pdb_seq_id = int(mmcif_dict["_pdbx_poly_seq_scheme.pdb_seq_num"][i])
        icode = mmcif_dict["_pdbx_poly_seq_scheme.pdb_ins_code"][i]
        if icode == ".":
            icode = None
        seqres[seqres_index] = (pdb_seq_id, icode)
    return seqres



def all_seqres_pdb_map():

	# PATH FOR THE PDB/mmCIF FILES
	import gzip

	pathmmcif = "/bmm/data/pdbmmcif/data/structures/all/mmCIF"

	dis = open("distinct_mutants.txt","r")
	ht = dis.readlines()
	dis.close()

	h = open("distinct_mutants_pdb.txt","w")

	k=0
	while k < len(ht): # end = len(ht)
		mutant = []
		mu=ht[k].split(', ')

		pdbid=mu[0].strip('[|\,|\'|]')
		pdb=pdbid[0:4]		# PDB NAME
		C=pdbid[5:6]		# CHAIN
		
		print("*** {} :: {} of {} ***" .format(pdb,k,len(ht)))

		# EXCUTE THE CODE TO PICK UP THE DESIRED ZONE AROUD THE RESIDUE

		try:
			pdbfile = "{}/{}.cif.gz".format(pathmmcif,pdb)
			tar = gzip.open("{}".format(pdbfile),"rb")
			out = open("pdbprocess.cif","wb")
			out.write(tar.read())
			tar.close()
			out.close()
	
			mmcif = MMCIF2Dict("pdbprocess.cif")
			idmap = seqres_atom_map(mmcif)

			reslist = [pdbid]
			k1 = 1
			while k1 < len(mu):
				id1 = int(mu[k1].strip("[|'|]|\n"))
				id2 = idmap[id1][0]
				reslist.append("{}".format(id2))
				k1 = k1 + 1
			h.write("{}".format(reslist))
			h.write("\n")
		except:
			print("FILE NOT FOUND")
	
		k = k + 1

all_seqres_pdb_map()

