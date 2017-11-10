set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set data1 [split $data "\n"]
 
set code [lindex $::argv 1]
set out [lindex $::argv 2]
set nprocess [lindex $::argv 3]

set N [expr { [llength $data] / $nprocess }]
set rem [expr { [llength $data] - ($N*$nprocess) }]

for {set i 1} {$i <= $nprocess} {incr i} {
	puts "wait ........."
	exec python $code [expr { $N * ($i - 1) }] [expr { $N * $i }] $out.$i &
	after 1000
}  
exec python $code [expr { $N * $i }] [expr { $N * $i + $rem }] $out.$i
puts "DONE!!!!!"
file delete *.pdb
