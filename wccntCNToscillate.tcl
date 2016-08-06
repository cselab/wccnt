##############################################
#
# This script oscillates a CNT wall
#
# usage :
# wcccnt CNToscillate $pdbFile $numOsc $radialShifts $outName
#

package provide wccnt 0.1

namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccntCNToscillate { args } {

	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt CNToscillate [args...]
		-pdbFile      : coordinate file
		-numOsc       : number of oscillations [ default : 1 ]
		-radialShifts : list of radial displacement (in A)
		-outName      : output name		
	    }
	    return
	}
	if { [llength $args] < 1 } then { usage; return }
	
	
	# Set the defaults
	set numOsc 1;
	
	
	# Parse options
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	    set arg [ lindex $args $argnum ]
	    set val [ lindex $args [expr $argnum + 1]]
	    switch -- $arg {
		"-pdbFile"      { set pdbFile      $val; incr argnum; }
		"-numOsc"       { set numOsc       $val; incr argnum; }
		"-radialShifts" { set radialShifts $val; incr argnum; }
		"-outName"      { set outName      $val; incr argnum; }
		default { error "error: CNToscillate: unknown option: $arg" }
	    }
	}
	
	
	# check non-default variables    
	set checkPDBFILE      [ info exists pdbFile ];
	set checkRADIALSHIFTS [ info exists radialShifts ];
	set checkOUTNAME      [ info exists outName ];
	
	if { $checkPDBFILE < 1 } {
	    error "error: CNToscillate: need to define variable -pdbFile"
	}    
	if { $checkRADIALSHIFTS < 1 } {
	    error "error: CNToscillate: need to define variable -radialShifts"
	}    
	if { $checkOUTNAME < 1 } {
	    error "error: CNToscillate: need to define variable -outName"
	}
	
	# -----------------
	# start procedures
	# -----------------		

	proc shiftRad { pdbFile radShift outName } {
	    
	    # load molecule
	    mol new $pdbFile type pdb waitfor all;
	    set molID [ molinfo top ];
	    set nAtom [ molinfo $molID get numatoms ];
	    
	    set i 0;    
	    while { $i < $nAtom } {
		
		# current X Y positions
		set selOne [ atomselect $molID "index $i" ];
		set currX [ $selOne get x ];
		set currY [ $selOne get y ];
		
		# magnitude of vector
		set radMag [ expr sqrt( ($currX * $currX) + ($currY * $currY) ) ];
		
		# unit vector
		set uVecX [ expr $currX / $radMag ];
		set uVecY [ expr $currY / $radMag ];
		
		# displacement 
		set dispMag [ expr $radMag + $radShift ];
		set dispX   [ expr $dispMag * $uVecX ];
		set dispY   [ expr $dispMag * $uVecY ];
		
		# do the move
		$selOne set x $dispX;
		$selOne set y $dispY;
		$selOne delete;
		
		incr i;
	    }
	    
	    # output shifted coordinates
	    set selAll [ atomselect $molID all ];
	    $selAll writepdb $outName.pdb;
	    
	    mol delete $molID;
	}
	
	
	# 2.- MAIN
	# ----------
	
	# create temporary directory
	file mkdir $outName;
		
	# do oscillations
	set iOsc 0;
	set iPDB 0;
	
	while { $iOsc < $numOsc } {
	    
	    puts "MESSAGE :: working on oscillation $iOsc of $numOsc"
	    
	    foreach item $radialShifts {
		
		# do the expansion/contraction
		shiftRad $pdbFile $item $outName.$iPDB;
		
		# move PDB
		file rename $outName.$iPDB.pdb $outName/$outName.$iPDB.pdb
		
		incr iPDB;
	    }
	    
	    incr iOsc;
	}
	
	
	# join PDBs
	mol new $pdbFile;
	
	set iDCD 0;    
	while { $iDCD < $iPDB } {
	    mol addfile $outName/$outName.$iDCD.pdb;
	    file delete $outName/$outName.$iDCD.pdb;
	    incr iDCD;
	}
	
	
	# write trajectory
	#>>animate write lammpstrj $outName.lammpstrj beg 0 end -1 waitfor all;
	animate write dcd $outName.dcd beg 0 end -1 waitfor all;
	
	
	#clean
	file delete $outName;
	mol delete all;		
    }    
}



