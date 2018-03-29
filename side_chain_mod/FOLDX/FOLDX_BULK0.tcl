# SCRIPT 1 FOR THE BULK EXECUTION OF FOLDX FOR THE CALCULATION OF VARIANT DATA

exec mkdir -p FOLDX_PREDICTION

set tlc { ALA ARG ASN ASP ASX CYS GLU GLN GLX GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL }

set slc { A R N D B C E Q Z G H I L K M F P S T W Y V }

# MUT_DATA FILE

set f [open "[lindex $::argv 0]" "r"]	
set data [read $f]
close $f

set g [open "FOLDX_dataset.txt" "w"]

set k [lindex $::argv 1]
set k [expr { $k * 5 }]
set end [lindex $::argv 2]
if { $end == "END" || $end == "end" } {
	set end [llength $data]
} else {
	set end [expr { $end * 5 }]
}
set actual_end [expr { $end / 5 }]	

while { $k < $end } {
	puts "[expr { $k / 5 }] of $actual_end"
	# WT
	set wt [lindex $data $k]
	set t1 [string range $wt 0 3]
	set t2 [string range $wt 5 [expr { [string length $wt] -1 }]]

	# MUT
	set mut [lindex $data [expr { $k + 1 }]]
	set t3 [string range $mut 0 3]
	set t4 [string range $mut 5 [expr { [string length $mut] -1 }]]

	set t5 [lindex $data [expr { $k + 2 }]]
	set t6 [lindex $data [expr { $k + 3 }]]
	set t7 [lindex $data [expr { $k + 4 }]]

	set k1 0
	while { [lindex $slc $k1] != $t7 } {
		incr k1
	}
	set t8 [lindex $tlc $k1]

	set k1 0
	while { [lindex $slc $k1] != $t6 } {
		incr k1
	}
	set t9 [lindex $tlc $k1]

	catch {
		exec python3 FOLDX_getpdb.py $t1 $t2
		set h [open "input.cfg" "w"]
		puts $h "command=PositionScan"
		puts $h "pdb=WT.pdb"
		puts $h "positions=$t6$t2$t5$t7"
		close $h
		exec foldx -f input.cfg 

		exec mv $t8$t5\_WT.pdb ./FOLDX_PREDICTION/$t1$t2-$t3$t4.pdb
		exec mv energies_$t5\_WT.txt ./FOLDX_PREDICTION/energies_$t1$t2-$t3$t4.txt

		puts $g "$t3 $t4 $t1$t2-$t3$t4 $t4 $t5"

		file delete $t8$t5\_WT.pdb
		file delete $t9$t5\_WT.pdb
		file delete PS_WT.fxout
		file delete PS_WT_scanning_output.txt
	}
	incr k 5
}
file delete pdbprocess.cif
file delete WT.pdb
file delete MUT.pdb
file delete 




	
	
