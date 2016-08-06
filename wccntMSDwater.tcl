############################################################
#
# This script calculates the MSD/N vs time
# where N = number of water oxygens
# the output is in nm^2 vs ns
#
# usage :
# wccnt MSDoxygen $psfFile $dcdFile $freqDCD $outName
#


package provide wccnt 0.1

namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccntMSDwater { args } {		   
	
	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt MSDwater [args...]
		-psfFile : structure file
		-dcdFile : trajectory file
		-freqDCD : frame frequency (in ns)
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
		"-psfFile"  { set psfFile $val; incr argnum; }
		"-dcdFile"  { set dcdFile $val; incr argnum; }
		"-freqDCD"  { set freqDCD $val; incr argnum; }
		"-outName"  { set outName $val; incr argnum; }
		default { error "error: MSDwater: unknown option: $arg" }
	    }
	}
	
	
	# check non-default variables    
	set checkPSF     [ info exists psfFile ];
	set checkDCD     [ info exists dcdFile ];
	set checkFREQDCD [ info exists freqDCD ];
	set checkOUTNAME [ info exists outName ];
	
	if { $checkPSF < 1 } {
	    error "error: MSDwater: need to define variable -psfFile"
	}    
	if { $checkDCD < 1 } {
	    error "error: MSDwater: need to define variable -dcdFile"
	}
	if { $checkFREQDCD < 1 } {
	    error "error: MSDwater: need to define variable -freqDCD"
	}    
	if { $checkOUTNAME < 1 } {
	    error "error: MSDwater: need to define variable -outName"
	}

		
	# -----
	#  run
	# -----
	
	# load molecule
	mol new $psfFile type psf waitfor all;
	set molID1   [ molinfo top ];
	mol addfile  $dcdFile type dcd molid $molID1 waitfor all;  
	
	# oxygen indexes
	set selO [ atomselect $molID1 "water and oxygen" ];
	set indexOlist [ $selO get index ];
	set lenOlist   [ llength $indexOlist ];
	$selO delete;
	
	# output file
	set outMSD    [ open $outName.msd w ];
	
	# loop over frames
	set numFrames [ molinfo $molID1 get numframes ];
	
	set i 0;
	while { $i < $numFrames } {
	    # initialize sum
	    set sumSquare 0;
	
	    foreach indexO $indexOlist {	
		# 0 time
		set sel0 [ atomselect $molID1 "index $indexO" frame 0 ];
		set xyz0 [ $sel0 get { x y z } ];
		set xyz0 [ lindex $xyz0 0 ];
		
		# i time	
		set selI [ atomselect $molID1 "index $indexO" frame $i ];
		set xyzI [ $selI get { x y z } ];
		set xyzI [ lindex $xyzI 0 ];
		
		# vector operations
		set vecDiff   [ vecsub $xyzI $xyz0 ];
		set magDiff2  [ veclength2 $vecDiff ];
		set sumSquare [ expr $sumSquare + $magDiff2 ];		
		
		# clean
		$sel0 delete;
		$selI delete;	
		unset xyz0;
		unset xyzI;
		unset vecDiff;
		unset magDiff2;
	    }
	    
	    # nm/ns units
	    set currTime [ expr $i * $freqDCD];
	    set msdNM2   [ expr $sumSquare / 100.0 ];
	    set msdNM2perAtom  [ expr $msdNM2 / $lenOlist ];
	    
	    puts "frame $i of $numFrames";
	    #puts $outMSD "$currTime $msdNM2 $msdNM2perAtom";
	    puts $outMSD "$currTime $msdNM2perAtom";
	    
	    incr i;
	}
	
	# clean
	close $outMSD;
	mol delete $molID1;
		
    }    
}


