set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set ncol [lindex $::argv 1]

set h [open "RMSD_TM_align.csv" "w"]
puts $h "WT,WT_CHAIN,MUT,MUT_CHAIN,GRMSD,LRMSD_CA,LRMSD_SC"

set k [lindex $::argv 2]
set end [lindex $::argv 3]
while { $k < $end } {
	puts "$k of [llength $data]"
	set t1 [lindex $data $k]
	set t2 [lindex $data [expr { $k + 1 }]]
	set t3 [lindex $data [expr { $k + 2 }]]
	set t4 [lindex $data [expr { $k + 3 }]]

	if  { $ncol == 5 } {
	
		# FOR LOCAL RMSD CALCULATIONS

		set t5 [lindex $data [expr { $k + 4 }]]
		catch {
			exec python3 script1.py 0 $t1 $t2 $t3 $t4
			exec ./TMalign chain1.pdb chain2.pdb -o output
			exec python3 script1.py 2 $t5

			set g [open "results" "r"]
			set data1 [read $g]
			close $g

			set grmsd [lindex $data1 0]
			set lrmsd1 [lindex $data1 1]
			set lrmsd2 [lindex $data1 2]
			puts $h "$t1,$t2,$t3,$t4,$grmsd,$lrmsd1,$lrmsd2"
		}
	} else {
		catch {

			exec python3 script1.py 0 $t1 $t2 $t3 $t4
			exec ./TMalign chain1.pdb chain2.pdb -o output
			exec python3 script1.py 1

			set g [open "results" "r"]
			set data1 [read $g]
			close $g

			set grmsd [lindex $data1 0]
			puts $h "$t1,$t2,$t3,$t4,$grmsd"
		}
	}

	catch {
		file delete output
		file delete output_all
		file delete output_atm
		file delete output_all_atm
		file delete output_all_atm_lig
		file delete chain1.pdb
		file delete chain2.pdb
		file delete results
		file delete file1.pdb
		file delete file2.pdb
	}
	incr k 4
}
close $h
#file delete output
#file delete output_all
#file delete output_atm
#file delete output_all_atm
#file delete output_all_atm_lig
