############################################################
#
# This script creates a CNT : PSF, PDB and TOP files
#
# usage:
# CNTbuild $nIndex $mIndex $lengthNM $perZ $segName $resName $resID $outName 
#


package provide wccnt 0.1

namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccntCNTbuild { args } {

	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt CNTbuild [args...]
		-nIndex   : n chiral index
		-mIndex   : m chiral index
		-lengthNM : CNT length (in nm)
		-perZ     : periodic bonds along Z-axis. 0=no, 1=yes [ default : 0 ]
		-segName  : segname   [ default : CNT ]
		-resName  : resname   [ default : CNT ]
		-atomType : atom type [ default : CA ]
		-resID    : resid     [ default : 1 ]
		-outName  : output name		
	    }
	    return
	}
	if { [llength $args] < 1 } then { usage; return }


	# Set the defaults
	set perZ     0;
	set segName  CNT;
	set resName  CNT;
	set atomType CA;
	set resID    1;
	
	
	# Parse options
	for { set argnum 0 } { $argnum < [llength $args] } { incr argnum } {
	    set arg [ lindex $args $argnum ]
	    set val [ lindex $args [expr $argnum + 1]]
	    switch -- $arg {
		"-nIndex"   { set nIndex   $val; incr argnum; }
		"-mIndex"   { set mIndex   $val; incr argnum; }
		"-lengthNM" { set lengthNM $val; incr argnum; }
		"-perZ"     { set perZ     $val; incr argnum; }
		"-segName"  { set segName  $val; incr argnum; }
		"-resName"  { set resName  $val; incr argnum; }
		"-atomType" { set atomType $val; incr argnum; }
		"-resID"    { set resID    $val; incr argnum; }
		"-outName"  { set outName  $val; incr argnum; }
		default { error "error: CNTbuild: unknown option: $arg" }
	    }
	}
	
	
	# check non-default variables    
	set checkNINDEX    [ info exists nIndex ];
	set checkMINDEX    [ info exists mIndex ];
	set checkLENGTHNM  [ info exists lengthNM ];
	set checkOUTNAME   [ info exists outName ];
	
	if { $checkNINDEX < 1 } {
	    error "error: CNTbuild: need to define variable -nIndex"
	}    
	if { $checkMINDEX < 1 } {
	    error "error: CNTbuild: need to define variable -mIndex"
	}
	if { $checkLENGTHNM < 1 } {
	    error "error: CNTbuild: need to define variable -lengthNM"
	}    
	if { $checkOUTNAME < 1 } {
	    error "error: CNTbuild: need to define variable -outName"
	}
	
	
	# -----------------
	# start procedures
	# -----------------		

	proc CNTper { nIndex mIndex lengthNM perZ outName} {
	    
	    #
	    # This script creates a periodic CNT
	    # the periodicity is along the Z-direction
	    #
	    # nIndex   : n index for CNT - integer
	    # mIndex   : m index for CNT - integer
	    # lengthNM : length in nm - real
	    # perZ     :  0 for non periodic system; 1 for periodic system - integer
	    # outName  : output name
	    #
	    # test:
	    # -----
	    # set nIndex 24;
	    # set mIndex 24;
	    # set lengthNM 20;
	    # set perZ 1; 
	    # set outName test04per
	    # CNTper  $nIndex $mIndex $lengthNM $perZ $outName
	    #
	    
	    
	    ################## MAIN #################
	    
	    # 1.- create nanotube
	    # --------------------
	    
	    # create molecule
	    package require nanotube
	    nanotube -l $lengthNM -n $nIndex -m $mIndex; # generate improper infomation
	    
	    # output PSF/PDB
	    set molID [ molinfo top ];
	    set selAll [ atomselect $molID all ]
	    animate write psf $outName.NonPer.psf sel $selAll waitfor all $molID;
	    animate write pdb $outName.NonPer.pdb sel $selAll waitfor all $molID;
	    
	    # pbc information
	    set aPBC [ molinfo $molID get a ];
	    set bPBC [ molinfo $molID get b ];
	    set cPBC [ molinfo $molID get c ];
	    set alphaPBC [ molinfo $molID get alpha ];
	    set betaPBC  [ molinfo $molID get beta ];
	    set gammaPBC [ molinfo $molID get gamma ];
	    
	    # clean
	    $selAll delete;
	    mol delete $molID;
	    
	    	    
	    # 2.- remove improper terms
	    # --------------------------
	    
	    # load molecule
	    mol new $outName.NonPer.psf type psf waitfor all;
	    set molID2   [ molinfo top ];
	    mol addfile  $outName.NonPer.pdb type pdb molid $molID2 waitfor all;
	    set selAll   [ atomselect $molID2 all ];
	    
	    # remove improper terms
	    package require topotools
	    topo -molid $molID2 -sel $selAll clearimpropers
    
	    # output PSF/PDB
	    animate write psf $outName.NonImpr.psf sel $selAll waitfor all $molID2;
	    animate write pdb $outName.NonImpr.pdb sel $selAll waitfor all $molID2;	    
	    
	    # clean
	    $selAll delete;
	    mol delete $molID2;
	    file delete $outName.NonPer.psf;
	    file delete $outName.NonPer.pdb;
	    
	    
	    # 3.- periodic bonds
	    # -------------------
	    
	    if { $perZ > 0 } {
		
		# load molecule
		mol new $outName.NonImpr.psf type psf waitfor all;
		set molID3 [ molinfo top ];
		mol addfile  $outName.NonImpr.pdb type pdb molid $molID3 waitfor all;
		
		# up and down rings
		set selAll [ atomselect $molID3 all ];
		foreach { cenX cenY cenZ } [ measure center $selAll ] { break };
		$selAll delete;		
		set selRingUp   [ atomselect $molID3 "(numbonds == 2) and (z > $cenZ)" ];
		set selRingDown [ atomselect $molID3 "(numbonds == 2) and (z < $cenZ)" ];

		# move down ring close to up ring
		set zPer     [ molinfo $molID3 get c ];
		set moveUp   "0 0 $zPer";
		set moveDown [ vecscale -1 $moveUp ];
		$selRingDown moveby $moveUp;
		
		# list of bonds between up and down rings
		set indexRingUp   [ $selRingUp get index ];
		set indexRingDown [ $selRingDown get index ];
		
		set cutoffCNT 1.43;
		set perBonds "";
		
		foreach index $indexRingDown {
		    set nearSel [ atomselect $molID3 "(all within $cutoffCNT of index $index ) and (index $indexRingUp)" ];
		    set nearIndex [ $nearSel get index ];
		    
		    foreach index2 $nearIndex {
			lappend perBonds "$index $index2";
		    }
		    
		    $nearSel delete;
		}
		
		# move down ring back
		$selRingDown moveby $moveDown;
		
		# add bonds
		foreach eachPair $perBonds {
		    foreach { indexLeft indexRight } $eachPair { break };
		    topo -molid $molID3 addbond $indexLeft $indexRight
		}
		
		# output PSF/PDB
		set selAll [ atomselect $molID3 all ];
		animate write psf $outName.Per.psf sel $selAll waitfor all $molID3;
		animate write pdb $outName.Per.pdb sel $selAll waitfor all $molID3;
		
		# clean
		$selAll delete;
		$selRingUp delete;
		$selRingDown delete;
		unset perBonds;
		mol delete $molID3;
		file delete $outName.NonImpr.psf;
		file delete $outName.NonImpr.pdb;
		
	    } else {	
		file rename $outName.NonImpr.psf $outName.Per.psf;
		file rename $outName.NonImpr.pdb $outName.Per.pdb;
	    }
    

	    # 4.- regenerate angles/dihedrals
	    # --------------------------------
	    package require psfgen
	    
	    resetpsf    
	    readpsf  $outName.Per.psf;
	    coordpdb $outName.Per.pdb;    
	    regenerate angles dihedrals    
	    writepsf $outName.psf
	    writepdb $outName.pdb
	    
	    # clean
	    resetpsf
	    file delete $outName.Per.psf;
	    file delete $outName.Per.pdb;
	    
	    
    
	    # 5.- add PBC info back
	    # -----------------------
	    
	    # load molecule
	    mol new $outName.psf type psf waitfor all;
	    set molID4 [ molinfo top ];
	    mol addfile  $outName.pdb type pdb molid $molID4 waitfor all;
	    
	    # set values
	    molinfo $molID4 set a $aPBC;
	    molinfo $molID4 set b $bPBC;
	    molinfo $molID4 set c $cPBC;
	    molinfo $molID4 set alpha $alphaPBC;
	    molinfo $molID4 set beta  $betaPBC;
	    molinfo $molID4 set gamma $gammaPBC;
	    
	    # output molecule
	    set selAll [ atomselect $molID4 all ];
	    animate write pdb $outName.pdb sel $selAll waitfor all $molID4;

	    # clean
	    $selAll delete;	    
	    mol delete $molID4
	}
	
	
	
	proc CNTtop { psfFile pdbFile segName resName resID atomType outName } {
	    
	    #
	    # This script generates a CNT topology
	    # from a PSF/PDB structure to be used in CHARMM2LAMMPS
	    #
	    # psfFile : PSF file
	    # pdbFile : PDB file
	    # segName : segname - string up to 4 characters long
	    # outName : output name
	    #
	    # test:
	    # -----
	    # set psfFile test04per.psf
	    # set pdbFile test04per.pdb
	    # set segName A; # no more than 4 letters
	    # set resName CNT; # no more than 4 letters
	    # set resID   1;
	    # set outName test10;
	    # CNTtop $psfFile $pdbFile $segName $resName $resID $outName;
	    #
	    #
	    # NOTE: this script renames atoms up to 475253 
	    #       if you have more atoms, re-write procedure chainName3
	    #
	    
	    
	    # 0.- procedures and previous calculations
	    # ----------------------------------------
	    
	    # defaul values for resName and resID here!!!!!!
	    
	    
	    # procedure to produce string from A to ZZZ
	    proc chainName3 { numChain2 } {
		
		# 475253 = ZZZZ
		# 18277 = ZZZ
		# 701 = ZZ 
		# 26 = Z
		# more numbers require new design
		
		proc chainName { numChain } {
		    
		    # template names
		    set listTemplate "A B C D E F G H I J K L M N O P Q R S T U V W X Y Z";
		    set lengthTemplate 26;
		    
		    set chainName [ lindex $listTemplate  [ expr $numChain%26 ] ];
		    
		    set prefixIndex [ expr ($numChain/26) - 1 ];
		    set prefix [ lindex $listTemplate $prefixIndex ];
		    
		    return "$prefix$chainName"	
		}
		
		if { $numChain2 <= 701 } {
		    # two letters
		    chainName $numChain2;	
		} else {
		    # more than two letters
		    set num1 [ expr $numChain2 - 702 ];
		    set perNum [ expr 702 - 26 ];
		    set num2 [ expr $num1/$perNum ];
		    set num3  [ expr $num1 - ( $num2 * $perNum ) + 26 ];
	    
		    set preLeft [ chainName $num2 ];
		    set preRight [ chainName $num3 ];
		    
		    return "$preLeft$preRight";	
		}    
	    }
	    
	    
	    ############## MAIN ################
	    
	    # 1.- change molecule info
	    # -------------------------
	    
	    # load structure
	    mol new $psfFile type psf waitfor all;
	    set molID1  [ molinfo top ];
	    mol addfile $pdbFile type pdb molid $molID1 waitfor all;	    

	    # PBC information
	    set aPBC     [ molinfo $molID1 get a ];
	    set bPBC     [ molinfo $molID1 get b ];
	    set cPBC     [ molinfo $molID1 get c ];
	    set alphaPBC [ molinfo $molID1 get alpha ];
	    set betaPBC  [ molinfo $molID1 get beta ];
	    set gammaPBC [ molinfo $molID1 get gamma ];
	    
	    # rename 
	    set selAll [ atomselect $molID1 all ];
	    $selAll set segname $segName;
	    $selAll set segid   $segName;
	    $selAll set resname $resName;
	    $selAll set resid   $resID;
	    $selAll set type    $atomType;
	    $selAll set beta 0;
	    $selAll set occupancy 0;
	    
	    # rename atom names
	    set listIndex [ $selAll get index ];
	    
	    set i 0;
	    foreach index $listIndex {
		set selOne [ atomselect $molID1 "index $i" ];
		
		# new atom name
		set newName [ chainName3 $i ];
		$selOne set name $newName;
		unset newName;
		
		$selOne delete;
		incr i;
	    }

	    # output molecule with new info
	    animate write psf $outName.TMP.psf sel $selAll waitfor all $molID1;
	    animate write pdb $outName.TMP.pdb sel $selAll waitfor all $molID1;
	    $selAll delete;
	    mol delete $molID1;
	    unset listIndex;
	    	    
	    
	    # 2.- topology info
	    # ------------------
	    
	    # topology file
	    set outTOP [ open $outName.top w ];
	    
	    # output mass/resname
	    puts $outTOP "31  1";
	    puts $outTOP " ";
	    puts $outTOP "MASS     1 $atomType    12.01100 C";
	    puts $outTOP " ";
	    puts $outTOP "AUTO ANGLES DIHE";
	    puts $outTOP " ";
	    puts $outTOP "RESI $resName          0.00";
	    
	    # load molecule
	    mol new $outName.TMP.psf type psf waitfor all;
	    set molID2  [ molinfo top ];
	    mol addfile $outName.TMP.pdb type pdb molid $molID2 waitfor all;
	    
	    # output atom names
	    set selAll    [ atomselect $molID2 all ];
	    set listIndex [ $selAll get index ];
	    $selAll delete;
	    
	    foreach index $listIndex {
		set  selOne  [ atomselect $molID2 "index $index" ];
		set  nameOne [ $selOne get name ];
		puts $outTOP "ATOM $nameOne $atomType      0.00";
		$selOne delete;
	    }
	    
	    
	    # create bond list
	    set bondList "";
	    
	    foreach index $listIndex {
		# find neigh around each index
		set selOne [ atomselect $molID2 "index $index" ];
		set neighIndex [ $selOne getbonds ];
		set neighIndex [ lindex $neighIndex 0 ];
		$selOne delete;
		
		# create a list of all indexes
		foreach indexB $neighIndex { 
		    set selA [ atomselect $molID2 "index $index" ];
		    set selB [ atomselect $molID2 "index $indexB" ];
		    set nameA [ $selA get name ];
		    set nameB [ $selB get name ];
		    
		    if { $index < $indexB } {
			lappend bondList "$nameA $nameB";
		    } elseif { $index > $indexB } {
			lappend bondList "$nameB $nameA";
		    } else {
		    }
		    
		    $selA delete;
		    $selB delete;
		}
	    }
	    
	    # sort bond list
	    set bondList [ lsort -unique $bondList ];
	    
	    # output bonds
	    foreach valPAIRunique $bondList {
		set valA [  lindex $valPAIRunique 0 ];
		set valB [  lindex $valPAIRunique 1 ];
		puts $outTOP "BOND $valA $valB";
	    }
	    
	    # close topology file;
	    puts  $outTOP " ";
	    close $outTOP;
	    
	    # clear
	    mol delete $molID2;
	    unset bondList;
	    unset listIndex;
    
	    
	    
	    # 3.- create molecule from topology file
	    # ---------------------------------------
    
	    package require psfgen
	    resetpsf
	    psfcontext reset;
	    topology $outName.top;
	    segment  $segName { pdb $outName.TMP.pdb } 
	    coordpdb $outName.TMP.pdb $segName;
	    regenerate angles dihedrals;
	    writepsf $outName.psf;
	    writepdb $outName.pdb;
	    resetpsf;
	    psfcontext reset;
	    
	    file delete $outName.TMP.psf;
	    file delete $outName.TMP.pdb;
	    
    
	    
	    # 4.- add PBC info back
	    # -----------------------
	    
	    # load molecule
	    mol new $outName.psf type psf waitfor all;
	    set molID3 [ molinfo top ];
	    mol addfile  $outName.pdb type pdb molid $molID3 waitfor all;
	    
	    # set PBC values
	    molinfo $molID3 set a $aPBC;
	    molinfo $molID3 set b $bPBC;
	    molinfo $molID3 set c $cPBC;
	    molinfo $molID3 set alpha $alphaPBC;
	    molinfo $molID3 set beta  $betaPBC;
	    molinfo $molID3 set gamma $gammaPBC;
    
	    # output molecule
	    set selAll [ atomselect $molID3 all ];
	    animate write pdb $outName.pdb sel $selAll waitfor all $molID3;
	    
	    # clean
	    mol delete $molID3
	    $selAll delete;
	    
	}
	
	
	# ---------------
	# run procedures
	# ---------------
	
	# create CNT
	CNTper  $nIndex $mIndex $lengthNM $perZ $outName.NoTop
	
	# create topology
	CNTtop $outName.NoTop.psf $outName.NoTop.pdb $segName $resName $resID $atomType $outName;

	# clean
	file delete  $outName.NoTop.psf;
	file delete  $outName.NoTop.pdb;

	unset nIndex;
	unset mIndex;
	unset lengthNM;
	unset perZ;
	unset outName;
	unset segName;
	unset resName;
	unset resID;
	unset atomType;

	        
    }   
}

