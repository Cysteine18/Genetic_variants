# SCRIPT 1 FOR THE BULK EXECUTION OF SCWRL FOR THE CALCULATION OF VARIANT DATA

exec mkdir -p SCWRL_PREDICTION

set 	SCWRLpath "/Users/tarunkhanna/Documents/Bioinformatics/SCWRL4"

set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set g [open "SCWRL_dataset.txt" "w"]

set k [lindex $::argv 1]
set k [expr { $k * 14 }]
set end [lindex $::argv 2]
if { $end == "END" || $end == "end" } {
	set end [llength $data]
} else {
	set end [expr { $end * 14 }]
}
set actual_end [expr { $end / 14 }]	
while { $k < $end } {
	puts "[expr { $k / 14 }] of $actual_end"

	# WT
	set t1 [lindex $data $k]	
	# CHAIN
	set t2 [lindex $data [expr { $k + 1 }]]	
	# MUT
	set t3 [lindex $data [expr { $k + 2 }]]	
	# CHAIN
	set t4 [lindex $data [expr { $k + 3 }]]	
	# MUT POSITION
	set t5 [lindex $data [expr { $k + 4 }]]	
	# RESWT
	set t6 [lindex $data [expr { $k + 12 }]]
	# RESMUT
	set t7 [lindex $data [expr { $k + 13 }]]

	for {set i 0} {$i < 8} {incr i} {
		# ZONE
		set t8 [lindex $data [expr { $k + 4 + $i}]]	
		catch {
			exec python3 SCWRL_MUT_BULK.py $t1 $t2 $t5 $t6 $t7 $t8
			exec $SCWRLpath/./Scwrl4 -i MUT.pdb -o $t1$t2-$t3$t4-mod$i.pdb -h

			exec mv $t1$t2-$t3$t4-mod$i.pdb ./SCWRL_PREDICTION

			puts $g "$t3 $t4 $t1$t2-$t3$t4-mod$i $t4 $t8"
		}
	}
	incr k 14
}
file delete pdbprocess.cif
file delete WT.pdb
file delete MUT.pdb




	
	
