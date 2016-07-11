###########################################################
#
# This script uses charmm2lammps.pl from LAMMPS
# to convert PSF/PDB and TOP/PAR files from charmm format into 
# DATA/IN files from LAMMPS format
#
# charmm2lammps.pl was written by Pieter J. in 't Veld (pjintve@sandia.gov)
# and Paul Crozier (pscrozi@sandia.gov), Sandia National Laboratories, 2005.
# The original charmm2lammps code is available in the LAMMPS repository 
# at https://github.com/lammps/lammps and used here under GLPv3 license
#

package provide wccnt 0.1

variable datadir $env(WCCNTDIR);


namespace eval ::wccnt:: {
    namespace export wccnt*

    proc wccntch2lmp { args } {

	global datadir; #directory containing WCCNT

	
	# Info
	proc usage {} {
	    vmdcon -info {usage: wccnt ch2lmp [args...]
		-psf     : structure file
		-pdb     : coordinate file	
		-par     : parameter file
		-top     : topology file
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
		"-psf"     { set psfFile $val; incr argnum; }
		"-pdb"     { set pdbFile $val; incr argnum; }
		"-par"     { set parFile $val; incr argnum; }
		"-top"     { set topFile $val; incr argnum; }
		"-outName" { set outName $val; incr argnum; }
		default { error "error: ch2lmp: unknown option: $arg" }
	    }
	}
	
	
	# check non-default variables    
	set checkPDBFILE      [ info exists pdbFile ];
	set checkPSFFILE      [ info exists psfFile ];
	set checkPARFILE      [ info exists parFile ];
	set checkTOPFILE      [ info exists topFile ];
	set checkOUTNAME      [ info exists outName ];
	
	if { $checkPDBFILE < 1 } {
	    error "error: ch2lmp: need to define variable -pdb"
	}    
	if { $checkPSFFILE < 1 } {
	    error "error: ch2lmp: need to define variable -psf"
	}    
	if { $checkPARFILE < 1 } {
	    error "error: ch2lmp: need to define variable -par"
	} 
	if { $checkTOPFILE < 1 } {
	    error "error: ch2lmp: need to define variable -top"
	}
	if { $checkOUTNAME < 1 } {
	    error "error: ch2lmp: need to define variable -outName"
	}
	

	# get PBC
	mol new  $psfFile type psf waitfor all;
	set molID  [ molinfo top ];
	mol addfile $pdbFile type pdb molid $molID waitfor all;
	
	set aPBC     [ molinfo $molID get a ];
	set bPBC     [ molinfo $molID get b ];
	set cPBC     [ molinfo $molID get c ];

	mol delete $molID;


	# create tmp files
	file copy $psfFile $outName.PDBPSF.psf
	file copy $pdbFile $outName.PDBPSF.pdb
	file copy $parFile par_$outName.TOPPAR.prm
	file copy $topFile top_$outName.TOPPAR.rtf
	
	
	# charmm2lammps
	set dirCH2LMP "$datadir/ch2lmp/"; # location of ch2lmp
	exec perl $dirCH2LMP/charmm2lammps.pl $outName.TOPPAR $outName.PDBPSF -lx=$aPBC -ly=$bPBC -lz=$cPBC; # execute charmm2lammps


	# clean 
	file rename -force $outName.PDBPSF.in   $outName.in;
	file rename -force $outName.PDBPSF.data $outName.data;
	
	file delete -force $outName.PDBPSF.psf;
	file delete -force $outName.PDBPSF.pdb;
	file delete -force par_$outName.TOPPAR.prm
	file delete -force top_$outName.TOPPAR.rtf
	file delete -force $outName.PDBPSF_ctrl.psf;
	file delete -force $outName.PDBPSF_ctrl.pdb;

    }
}

