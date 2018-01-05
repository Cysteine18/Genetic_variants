# THIS CODE MAKE USE OF SUHAIL'S CODE getpdb AND capri-fit TO CALCULATE THE RMSD BETWEEN THE WIILDTYPE AND THE MUTANTS

set f [open "cluster.txt" "r"]
set data1 [read $f]
close $f
set data [split $data1 "\n"]

set g [open "mutations_rmsd_clustal.out" "w"]
set g1 [open "mutations_rmsd_mammoth.out" "w"]

set g2 [open "mutations_rmsd_errors.out" "w"]
set g3 [open "pdb_data.txt" "w"]

exec mkdir -p overlay_images

set k 0
set cid 1
set nmut 1
while { $k < [llength $data] } {
	if { [lindex $data $k 0] == "#" && [lindex $data $k 1] == $cid } {
		incr k
		set wt [lindex $data $k 0]
		set wt [string trim "$wt" "(|'|',"]
		puts "		### CLUSTER $cid WITH $wt AS SELECTED WILDTYPE ###"
		incr cid
		set wtpdb [string range $wt 0 3]
		set wtchain [string range $wt 5 [string length $wt]]
		# GETTING THE PDB FILE
		set exv 0
		catch {
			exec getpdb $wtpdb
			incr exv
		}
		if { $exv != 0 } {
			incr k
			while { [lindex $data $k 0] != "#" } {
				set mut [lindex $data $k 0]
				set mut [string trim "$mut" "(|'|',"]
				set mutpdb [string range $mut 0 3]
				set mutchain [string range $mut 5 [string length $mut]]
				incr k
				set exv 0
				catch {
					exec getpdb $mutpdb
					incr exv
				}
				if { $exv != 0 } {
					set exv 0
					catch {
						exec capri-fit -iref $wtpdb.pdb $wtchain -i $mutpdb.pdb -ic $mutchain -o rmsdC -clustal
						exec capri-fit -iref $wtpdb.pdb $wtchain -i $mutpdb.pdb -ic $mutchain -o rmsdM -mammoth
						incr exv
					}
					if { $exv != 0 } {

						# CLUSTAL VALUES

						set t1 [open "rmsdC.log" "r"]
						set data2 [read $t1]
						close $t1
	
						set data3 [split $data2 "\n"]

						set k1 0
						set count 0
						while { [lindex $data3 $k1 0] != "\$DATASUM" } {
							incr k1
							if { $k1 > [llength $data3] } {
								break
								incr count
							}
						}
						set value [lindex $data3 $k1 2]
				
						puts $g "$nmut	$wt	$mut	$value"
						
						# MAMMMOTH VALUES

						set t1 [open "rmsdM.log" "r"]
						set data2 [read $t1]
						close $t1
	
						set data3 [split $data2 "\n"]

						set k1 0
						while { [lindex $data3 $k1 0] != "\$DATASUM" } {
							incr k1
							if { $k1 > [llength $data3] } {
								break
								incr count
							}
						}
						if { $count == 0 } {
							set value [lindex $data3 $k1 2]
							puts $g1 "$nmut	$wt	$mut	$value"
						} else {
							puts $g2 "$wt $mut"
						}
						
						incr nmut
					} else {
						puts $g2 "$wt $mut"
					}
				} else {
					puts $g2 "$wt $mut"
				}
				puts $g3 "$wtpdb $wtchain $mutpdb $mutchain"
			}
		} else {
			puts $g2 "$wt $mut"
		}

		# DELETING THE FILES CREATED IN THE ABOVE PROCEDURE
		
		catch {
			set nv [expr { $cid - 1 }]
			exec mv rmsdC.pdb ./overlay_images/C.$nv.pdb
			exec mv rmsdM.pdb ./overlay_images/M.$nv.pdb
			file delete {*}[glob *.pdb]
			file delete {*}[glob *.log]
		}

	} else {
		incr k
	}
}				
close $g
close $g1
close $g2
close $g3
