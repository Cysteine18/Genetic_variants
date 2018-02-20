set pathpdb "/Users/tarunkhanna/Documents/Bioinformatics/Office_comp/HIGH_RES_CLUSTER/TMa/TM_align"

set pdb1 [lindex $::argv 0]
set pdb2 [lindex $::argv 1]

exec /Applications/chimera.app//Contents/MacOS/chimera $pathpdb/$pdb1-$pdb2.pdb
