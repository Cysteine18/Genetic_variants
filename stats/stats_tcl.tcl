proc distinct_mutants {} {

	# TO REMOVE SIMILAR SEQUENCES IN THE MUTANT LIST

	set f [open "all_mutants_list.txt" "r"]
	set data [read $f]
	close $f

	set g [open "seq_sim_list.txt" "r"]
	set data1 [read $g]
	close $g

	set h [open "unique_mutants.txt" "w"]
	set k 0
	#for {set i 0} {$i < 1000000} {incr i} {
		#set dummy($i) ""
	#}
	set i 0
	set l 0
	while { $k < [llength $data] } {
		puts "			#### $k of [llength $data] ####"
		set k1 0
		set count1 0
		set t1 [lindex $data $k $k1]
		set count1 0
		for {set j 0} {$j < $i} {incr j} {
			if { $dummy($j) == $t1 } {
				incr count1
			}
		}
		incr k1
		set list1 ""
		if { $count1 == 0 } {
			set k2 0
			while { $k2 < [llength $data1] } {
				if { [lindex $data1 $k2 0] == $t1 } {
					for {set j 1} {$j < [llength [lindex $data1 $k2]]} {incr j} {
						set dummy($i) [lindex $data1 $k2 $j]
						incr i
					}
					set k2 [llength $data]
				}
				incr k2
			}
			lappend list1 $t1
			while { $k1 < [llength [lindex $data $k]] } {
				set t1 [lindex $data $k $k1]
				set count 0
				for {set j 0} {$j < $i} {incr j} {
					if { $dummy($j) == $t1 } {
						incr count
					}
				}
				if { $count == 0 } {
					set k2 0
					while { $k2 < [llength $data1] } {
						if { [lindex $data1 $k2 0] == $t1 } {
							for {set j 1} {$j < [llength [lindex $data1 $k2]]} {incr j} {
								set dummy($i) [lindex $data1 $k2 $j]
								incr i
							}
							set k2 [llength $data]
						}
						incr k2
					}
					lappend list1 $t1
				}
				incr k1
			}
		}
		if { $list1 != "" } {
			puts $h "\{$list1\}"
		}
		incr k
	}
	close $h
}


proc mutant_count {} {

	
	# COUNTING THE NUMBER OF MUTANTS

	set f [open "unique_mutants.txt" "r"]
	set data [read $f]
	close $f

	set g [open "number_of_distinct_mutants.txt" "w"]

	set k 0
	set count 0.0
	set nterms 0.0
	set i 0
	while { $k < [llength $data] } {
		set list2 ""
		puts "		*** $k of [llength $data] ***"
		set t1 [lindex $data $k 0]
		set len [llength [lindex $data $k]]
		if { $len > 1 } {
			lappend list2 $t1
			set k1 1
			while { $k1 < $len } {
				set t2 [lindex $data $k $k1]
				set count 0
				for {set j 0} {$j < $i} {incr j} {
					if { $t2 == $M1($j) && $t1 == $M2($j)} {
						incr count
					}
				}
				if { $count == 0 } {
					set M1($i) $t1
					set M2($i) $t2
					lappend list2 $t2
					incr i
				}
				incr k1
			}
			if { [llength $list2] > 1 } {	
				puts $g "{$list2}"
			}
		}
		incr k
	}
	puts "		#### NUMBER OF DISTINCT MUTANTS = $i ####"
	puts $g ""
	puts $g "NUMBER OF DISTINCT MUTANTS = $i"
}

mutant_count















