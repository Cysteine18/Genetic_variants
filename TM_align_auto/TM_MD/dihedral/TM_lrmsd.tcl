set pathSCWRL "/Volumes/BIOINFO/SCWRL_PREDICTION"
set pathTMalign "/Users/tarunkhanna/Documents/Bioinformatics/TM_align"
 
set f [open "[lindex $::argv 0]" "r"]
set data [read $f]
close $f

set ncol [lindex $::argv 1]

set sframe 95
set eframe 100

exec mkdir -p TM_align

set k [lindex $::argv 2]
set k [expr  { $k * 7 }]
if { [lindex $::argv 3] == "end" || [lindex $::argv 3] == "END" } {
	set end [llength $data]
} else {
	set end [expr { [lindex $::argv 3] * 7 }]
}

while { $k < $end } {
	puts "$k of [llength $data]"
	set t1 [lindex $data $k]
	set t2 [lindex $data [expr { $k + 1 }]]
	set t3 [lindex $data [expr { $k + 5 }]]
	set t4 [lindex $data [expr { $k + 3 }]]

	if  { $ncol == 7 } {
	
		# FOR LOCAL RMSD CALCULATIONS

		set t5 [lindex $data [expr { $k + 6 }]]
		catch {
			set j [expr { $k / 7 }]
			set h [open "RMSD_TM_align$j.csv" "w"]
			puts $h "WT,WT_CHAIN,MUT,MUT_CHAIN,WT_LEN_PDB,MUT_LEN_PDB,ALIGNED_LENGTH,GRMSD,LRMSD_CA,LRMSD_SC,CHI1,CHI2,CHI3,CHI4,CHI5"
			for {set i $sframe} {$i < $eframe} {incr i} {

				puts " ** $i of $eframe **"

				exec python3 script1.py 0 $t1 $t2 $t3 $t4 $i
				#exec cp $pathSCWRL/$t3.pdb ./chain2.pdb
				exec $pathTMalign/./TMalign chain1.pdb chain2.pdb -o output | tee temp
				exec python3 script1.py 2 $t5

				set g [open "results" "r"]
				set data1 [read $g]
				close $g

				set grmsd [lindex $data1 0]
				set lwt [lindex $data1 1]
				set lmut [lindex $data1 2]
				set coverage [lindex $data1 3]
				set lrmsd1 [lindex $data1 4]
				set lrmsd2 [lindex $data1 5]
				set chi1 [lindex $data1 6]
				set chi2 [lindex $data1 7]
				set chi3 [lindex $data1 8]
				set chi4 [lindex $data1 9]
				set chi5 [lindex $data1 10]
				#set lrmsd3 [lindex $data1 6]

				puts $h "$t1,$t2,$t3,$t4,$lwt,$lmut,$coverage,$grmsd,$lrmsd1,$lrmsd2,$chi1,$chi2,$chi3,$chi4,$chi5"
			}
			close $h
		}
	} else {
		catch {
			set j [expr { $k / 5 }]
			set h [open "RMSD_TM_align$j.csv" "w"]
			puts $h "WT,WT_CHAIN,MUT,MUT_CHAIN,WT_LEN_PDB,MUT_LEN_PDB,ALIGNED_LENGTH,GRMSD"
			for {set i $sframe} {$i < $eframe} {incr i} {
		
				puts " ** $i of $eframe (GRMSD CALCULATION)**"
				exec python3 script1.py 0 $t1 $t2 $t3 $t4 $i
				exec $pathTMalign/./TMalign chain1.pdb chain2.pdb -o output | tee temp
				exec python3 script1.py 1

				set g [open "results" "r"]
				set data1 [read $g]
				close $g

				set grmsd [lindex $data1 0]
				set lwt [lindex $data1 1]
				set lmut [lindex $data1 2]
				set coverage [lindex $data1 3]

				puts $h "$t1,$t2,$t3,$t4,$lwt,$lmut,$coverage,$grmsd"
			}
			close $h
		}
	}
	catch {
		exec mv file2.pdb ./TM_align/$t1$t2-$t3$t4.pdb
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
		file delete temp
	}
	incr k 7
}
#close $h
#file delete output
#file delete output_all
#file delete output_atm
#file delete output_all_atm
#file delete output_all_atm_lig
