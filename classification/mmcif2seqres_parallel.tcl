set f [open "[lindex $::argv 0]" "r"]
set data1 [read $f]
close $f

set data [split $data1 "\n"]
set nprocess [lindex $::argv 1]
set N 0
set N [expr { [llength $data] / $nprocess }]

set rem [expr { [llength $data] - ($N*$nprocess) }]

for {set i 1} {$i <= $nprocess} {incr i} {
	puts "WAIT"
	exec python mmcif2seqres.py [expr { ($i - 1)*$N }] [expr { $i * $N }] distant_mutants_pdb_$i.txt &
	after 1000
}
exec python mmcif2seqres.py [expr { ($i - 1)*$N }] [expr { (($i-1) * $N) + $rem }] distant_mutants_pdb_$i.txt 
puts "SETUP DONE"
