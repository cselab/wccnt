############################################################
#
# This script converts LAMMPS trj into DCD
# position and velocities are saved in separate DCDs
#
# usage:
# set trjFile [path_to_file]/[name_of_file].lammpstrj;
# set outName [name]
# trj2dcd $trjFile $outName
#

package provide wccnt 0.1


namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccnttrj2dcd { args } {
	
	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt trj2dcd [args...]
		-trjFile  : trajectory in LAMMPS format
		-outName : output name		
	    }
	    return
	}
	if { [llength $args] < 1 } then { usage; return }
	
	
	# Parse options
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	    set arg [ lindex $args $argnum ]
	    set val [ lindex $args [expr $argnum + 1]]
	    switch -- $arg {
		"-trjFile"  { set trjFile  $val; incr argnum; }
		"-outName"  { set outName  $val; incr argnum; }
		default { error "error: trj2dcd: unknown option: $arg" }
	    }
	}

	
	# check non-default variables    
	set checkTRJFILE   [ info exists trjFile ];
	set checkOUTNAME   [ info exists outName ];
	
	if { $checkTRJFILE < 1 } {
	    error "error: trj2dcd: need to define variable -trjFile"
	}
	if { $checkOUTNAME < 1 } {
	    error "error: trj2dcd: need to define variable -outName"
	}


	########### MAIN ############

	# load molecule
	mol new $trjFile type lammpstrj first 0 last -1 step 1 waitfor all;
	set molID1   [ molinfo top ];
	
	# export original dcd file (containing positions as x y z)
	set selAll [ atomselect $molID1 all ];
	animate write dcd pos_$outName.dcd beg 0 end -1 sel $selAll waitfor all $molID1;
	$selAll delete;
	puts "DONE with position dump"
	
	# loop over frames to replace positions with velocities
	set numFrames [ molinfo $molID1 get numframes ];
	
	set i 0;
	while { $i < $numFrames } {
	    
	    # get positions in i time	
	    set selI [ atomselect $molID1 all frame $i ];
	    
	    # get velocities in i time	
	    set vxyzI [ $selI get { vx vy vz } ];
	    
	    # replace { x y z } with { vx vy vz } 
	    $selI set { x y z } $vxyzI
	    
	    # clean
	    $selI delete;	
	    unset vxyzI;
	    
	    incr i;
	} 
	
	# export velocity dcd file (containing velocities as x y z)
	set selAll [ atomselect $molID1 all ];
	animate write dcd vel_$outName.dcd beg 0 end -1 sel $selAll  waitfor all $molID1;
	$selAll delete;
	puts "DONE with velocity dump"
	
	# clean
	mol delete $molID1;    	
    }
}

