set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set h [open "[lindex $::argv 1]" "w"]

puts $h "#WT	MUT	no_sim_wt	no_sim_mut grmsdwt	SDgrmsdwt	grmsdmut	SDgrmsdmut	lrmsdwt	SDlrmsdwt	lrmsdmut	SDlrmsdmut	lrmsdscwt	SDlrmsdscwt	lrmsdscmut	SDlrmsdscmut"

#exec mkdir -p TM_align

set k 11

while { $k < [llength $data] } {

	set wt ""
	set mut ""

	set t1 [lindex $data [expr { $k + 2 }]]
	set t2 [lindex $data [expr { $k + 3 }]]
	set wt [append wt $t1 "_" $t2 ]
	set t3 [lindex $data [expr { $k + 4 }]]
	set t4 [lindex $data [expr { $k + 5 }]]
	set mut [append mut $t3 "_" $t4 ]

	set pos [lindex $data [expr { $k + 6 }]]
	
	# FOR LOCAL RMSD CALCULATIONS

	puts "[expr { $k / 11 }] of [expr { [llength $data] / 11 }] $wt $mut"
	catch {
		# WILD TYPE
		exec python3 sim_struct.py $wt $pos "zone.txt"

		set g [open "zone.txt" "r"]
		set data2 [read $g]
		close $g

		set nw 0
		set grmsdwt "NA"
		set lrmsdwt "NA"
		set lrmsdscwt "NA"
		set grmsdmut "NA"
		set lrmsdmut "NA"
		set lrmsdscmut "NA"
		set sdgrmsdwt "NA"
		set sdlrmsdwt "NA"
		set sdlrmsdscwt "NA"
		set sdgrmsdmut "NA"
		set sdlrmsdmut "NA"
		set sdlrmsdscmut "NA"
		if { [llength $data2] > 1 } {
			set k1 0
			set sl 0.0
			set sg 0.0
			set ssc 0.0
			set term1 0
			set term2 0
			set term3 0
			while { $k1 < [llength $data2] } {
				puts "		## WILDTYPE :: [expr { ($k1 / 5) + 1 }] of [expr { [llength $data2] / 5 }] ##"
				set t1 [lindex $data2 [expr { $k1 + 0 }]]
				set t2 [lindex $data2 [expr { $k1 + 1 }]]
				set t3 [lindex $data2 [expr { $k1 + 2 }]]
				set t4 [lindex $data2 [expr { $k1 + 3 }]]
				set t5 [lindex $data2 [expr { $k1 + 4 }]]

				exec python3 script1.py 0 $t1 $t2 $t3 $t4
				exec ./TMalign chain1.pdb chain2.pdb -o output | tee temp
				exec python3 script1.py 2 $t5

				set g [open "results" "r"]
				set data1 [read $g]
				close $g

				set grmsd [lindex $data1 0]
				if { $grmsd != "NA" } {
					set sgl($term1) $grmsd
					incr term1
					set sg [expr { $grmsd + $sg }]
				}

				set lrmsd1 [lindex $data1 4]
				if { $lrmsd1 != "NA" } {
					set sll($term2) $lrmsd1
					incr term2
					set sl [expr { $lrmsd1 + $sl }]
				}

				set lrmsd2 [lindex $data1 5]
				if { $lrmsd2 != "NA" } {
					set sscl($term3) $lrmsd2
					incr term3
					set ssc [expr { $lrmsd2 + $ssc }]
				}

				incr k1 5
				incr nw
			}
			
			if { $term1 > 0 } {
				set grmsdwt [format "%.2f" [expr { $sg / $term1 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term1 } {
					set diff [expr { $grmsdwt - $sgl($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdgrmsdwt [expr { $sd / $term1 }]
				set sdgrmsdwt [format "%.2f" [expr { sqrt($sdgrmsdwt) }]]	
			} else {
				set grmsdwt "NA"
				set sdgrmsdwt "NA"
			}
			if { $term2 > 0 } {
				set lrmsdwt [format "%.2f" [expr { $sl / $term2 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term2 } {
					set diff [expr { $lrmsdwt - $sll($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdlrmsdwt [expr { $sd / $term2 }]
				set sdlrmsdwt [format "%.2f" [expr { sqrt($sdlrmsdwt) }]]
			} else {
				set lrmsdwt "NA"
				set sdlrmsdwt "NA"
			}
			if { $term3 > 0 } {
				set lrmsdscwt [format "%.2f" [expr { $ssc / $term3 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term3 } {
					set diff [expr { $lrmsdscwt - $sscl($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdlrmsdscwt [expr { $sd / $term3 }]
				set sdlrmsdscwt [format "%.2f" [expr { sqrt($sdlrmsdscwt) }]]
			} else {
				set lrmsdscwt "NA"
				set sdlrmsdscwt "NA"
			}
		}

		# MUTANT

		exec python3 sim_struct.py $mut $pos "zone.txt"

		set g [open "zone.txt" "r"]
		set data2 [read $g]
		close $g

		set nm 0
		if { [llength $data2] > 1 } {
			set k1 0
			set sl 0.0
			set sg 0.0
			set ssc 0.0
			set term1 0
			set term2 0
			set term3 0
			while { $k1 < [llength $data2] } {
				puts "		## MUTANT :: [expr { ($k1 / 5) + 1 }] of [expr { [llength $data2] / 5 }] ##"
				set t1 [lindex $data2 [expr { $k1 + 0 }]]
				set t2 [lindex $data2 [expr { $k1 + 1 }]]
				set t3 [lindex $data2 [expr { $k1 + 2 }]]
				set t4 [lindex $data2 [expr { $k1 + 3 }]]
				set t5 [lindex $data2 [expr { $k1 + 4 }]]

				exec python3 script1.py 0 $t1 $t2 $t3 $t4
				exec ./TMalign chain1.pdb chain2.pdb -o output | tee temp
				exec python3 script1.py 2 $t5


				set g [open "results" "r"]
				set data1 [read $g]
				close $g

				set grmsd [lindex $data1 0]
				if { $grmsd != "NA" } {
					set sgl($term1) $grmsd
					incr term1
					set sg [expr { $grmsd + $sg }]
				}

				set lrmsd1 [lindex $data1 4]
				if { $lrmsd1 != "NA" } {
					set sll($term2) $lrmsd1
					incr term2
					set sl [expr { $lrmsd1 + $sl }]
				}

				set lrmsd2 [lindex $data1 5]
				if { $lrmsd2 != "NA" } {
					set sscl($term3) $lrmsd2
					incr term3
					set ssc [expr { $lrmsd2 + $ssc }]
				}
				incr nm
				incr k1 5
			}

			if  { $term1 > 0 } {
				set grmsdmut [format "%.2f" [expr { $sg / $term1 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term1 } {
					set diff [expr { $grmsdmut - $sgl($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdgrmsdmut [expr { $sd / $term1 }]
				set sdgrmsdmut [format "%.2f" [expr { sqrt($sdgrmsdmut) }]]	
			} else {
				set grmsdmut "NA"
				set sdgrmsdmut "NA"
			}
			if { $term2 > 0 } {
				set lrmsdmut [format "%.2f" [expr { $sl / $term2 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term2 } {
					set diff [expr { $lrmsdmut - $sll($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdlrmsdmut [expr { $sd / $term2 }]
				set sdlrmsdmut [format "%.2f" [expr { sqrt($sdlrmsdmut) }]]
			} else {
				set lrmsdmut "NA"
				set sdlrmsdmut "NA"
			}
			if { $term3 > 0 } {
				set lrmsdscmut [format "%.2f" [expr { $ssc / $term3 }]]
				set k1 0
				set sd 0.0
				while { $k1 < $term3 } {
					set diff [expr { $lrmsdscmut - $sscl($k1) }]
					set diff [expr { $diff * $diff }]
					set sd [expr { $sd + $diff }]
					incr k1
				}
				set sdlrmsdscmut [expr { $sd / $term3 }]
				set sdlrmsdscmut [format "%.2f" [expr { sqrt($sdlrmsdscmut) }]]
			} else {
				set lrmsdscmut "NA"
				set sdlrmsdscmut "NA"
			}
		}			
		puts $h "$wt	$mut	$nw	$nm	$grmsdwt	$sdgrmsdwt	$grmsdmut	$sdgrmsdmut	$lrmsdwt	$sdlrmsdwt	$lrmsdmut	$sdlrmsdmut	$lrmsdscwt	$sdlrmsdscwt	$lrmsdscmut	$sdlrmsdscmut"
	}
	#catch {
		#exec mv file2.pdb ./TM_align/$t1$t2-$t3$t4.pdb
	#}

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
		file delete temp
		file delete pdbprocess.cif
	}
	incr k 11
}
close $h
#file delete output
#file delete output_all
#file delete output_atm
#file delete output_all_atm
#file delete output_all_atm_lig
