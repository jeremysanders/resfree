proc ratio_scatter { dg outfile } {
    # get points with scatter
    newpar 50 1
    newpar 51 1
    tclout plot m x $dg
    set xvars $xspec_tclout
    tclout plot m y $dg
    set yscat $xspec_tclout

    # turn off scatter
    newpar 50 0
    newpar 51 0
    tclout plot m y $dg
    set ynoscat $xspec_tclout
    
    # write to file
    set f [open $outfile "w"]

    set nopts [llength $xvars]
    for { set i 0 } { $i < $nopts } { incr i } {
	if { [lindex $ynoscat $i] > 0. } {
	    set ratio [expr [lindex $yscat $i]/[lindex $ynoscat $i]]
	    puts $f "[lindex $xvars $i] $ratio"
	}
    }
    puts $f "log x on"
    puts $f "font roman"
    puts $f "time off"
    puts $f "la x Energy (keV)"
    puts $f "la y Ratio"
    puts $f "lw 2"
    puts $f "csize 1.3"
    puts $f "la f"
    close $f

    # restore scatter
    newpar 50 1
}

tclout datasets
set nodatasets $xspec_tclout

for {set ds 1} {$ds <= $nodatasets} {incr ds} {
    ratio_scatter $ds ratio_scatter_${ds}.qdp
}
