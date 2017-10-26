# THIS CODE TRIES TO GET THE MUTANTS BASED ON THE SEQUENCE SIMILARITY APPROACH

# FASTA FILE INPUT

set fasta_file "[lindex $::argv 0]"

# NUMBER OF ALLOWED MUTATION

set nb [lindex $::argv 1]

set f [open "$fasta_file" "r"]
set data [read $f]
close $f

set data1 [string range $data 0 [expr { [string length $data] / 100 }]]

puts "hi"
puts [llength $data1]

set k 0

while { $k < [llength $data] } {
	set term [lindex $data $k]

	if { [string range 0 1 $term] == ">" } {
		#puts "			•••• CHECKING PDB [string range 1 expr [{ [string length $term] - 1 }]]
	}
}
