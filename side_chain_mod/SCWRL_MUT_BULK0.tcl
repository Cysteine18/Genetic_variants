# SCRIPT 1 FOR THE BULK EXECUTION OF SCWRL FOR THE CALCULATION OF VARIANT DATA

exec mkdir -p SCWRL_PREDICTION

set 	SCWRLpath "/Users/tarunkhanna/Documents/Bioinformatics/SCWRL4"

set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set g [open "SCWRL_dataset3.txt" "w"]

set k [lindex $::argv 1]
set k [expr { $k * 8 }]
set end [lindex $::argv 2]
if { $end == "END" || $end == "end" } {
	set end [llength $data]
} else {
	set end [expr { $end * 8 }]
}
set actual_end [expr { $end / 8 }]	
while { $k < $end } {
	puts "[expr { $k / 8 }] of $actual_end"
	set t1 [lindex $data $k]
	set t2 [lindex $data [expr { $k + 1 }]]
	set t3 [lindex $data [expr { $k + 2 }]]
	set t4 [lindex $data [expr { $k + 3 }]]
	set t5 [lindex $data [expr { $k + 4 }]]
	set t6 [lindex $data [expr { $k + 5 }]]
	set t7 [lindex $data [expr { $k + 6 }]]
	set t8 [lindex $data [expr { $k + 7 }]]

	catch {
		exec python3 SCWRL_MUT_BULK.py $t1 $t2 $t6 $t7 $t8
		exec $SCWRLpath/./Scwrl4 -i MUT.pdb -o $t1$t2-$t3$t4.pdb -h

		exec mv $t1$t2-$t3$t4.pdb ./SCWRL_PREDICTION

		puts $g "$t3 $t4 $t1$t2-$t3$t4 $t4 $t5"
	}
	incr k 8
}
file delete pdbprocess.cif
file delete WT.pdb
file delete MUT.pdb




	
	
