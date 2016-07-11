
package provide wccnt 0.1

proc wccnt { args } {

    proc usage {} {
	vmdcon -info {usage: wccnt <command> [args...]
	    
	    Build PSF/PDB/TOP files for CNT:
	     CNTbuild [options...]
	    
	    Fake oscillations for CNT:
	     CNToscillate [options...]

            Mean square displacement of water - for diffusion:
             MSDwater [options...]	    

	    Convert LAMMPS trj into pos- and vel- DCD
	     trj2dcd [options...]

	    Convert PSF/PDB and PAR/TOP into DATA/IN
	     ch2lmp [options...]

	}
	return
    }
    
    if { [llength $args] < 1 } then { usage; return }
    set command [ lindex $args 0 ]
    set args [lrange $args 1 end]
    set fullcommand "::wccnt::wccnt$command"
        
    if { [ string length [namespace which -command $fullcommand]] } then {
	eval "$fullcommand $args"
    } else { usage; return }

}

