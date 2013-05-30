# ### BiaxialConcrete_matTest.tcl
# ### single truss for testing of material stress-strain relation for the biaxial concrete and associated truss element
# ### Yuan Lu 05/22/2013

wipe;
model BasicBuilder -ndm 2 -ndf 3;

set tcl_precision 12; # set precision
 
set in 1.;
set kip 1.;
set sec 1.;
set kN [expr 0.2248089431*$kip];
set ksi [expr $kip/pow($in,2)];
set psi [expr $ksi/1000.];
set ft [expr 12.*$in];
set lb [expr $kip/1000.];
set pcf [expr $lb/pow($ft,3)];
set cm [expr $in/2.54];
set mm [expr $in/25.4];
set meter [expr 100.*$cm];
set MPa [expr 145*$psi];
set GPa [expr 1000*$MPa];
set PI [expr 2*asin(1.0)];
set g [expr 32.2*$ft/pow($sec,2)];
set U 1.e10;
set u [expr 1/$U];

# --- MATERIAL MODELS
set concHoriz	10;
set elasticMat	100;

# General Concrete stuff
set fc 		[expr -5.8*$ksi];		# CONCRETE Compressive Strength, ksi   (+Tension, -Compression)
set ec 		-0.002;
set Ec 		[expr {57000.*sqrt(-$fc/$psi)*$psi}];	# Concrete Elastic Modulus
set lambda 0.5;			# ratio between unloading slope at $eps2 and initial slope $Ec

set ft [expr {4*sqrt(-$fc/$psi)*$psi}];		# tensile strength +tension
set epscr [expr $ft/$Ec];
set epstu [expr 1.0]


set fc1  [expr $fc];
set eps1 [expr $ec];
set fc2 [expr 0.01*$fc];
set temp [expr 0.5*($Ec + $fc1/$eps1)];
set eps2 [expr ($eps1 + -0.002)];
set fc3 [expr 0.0*$fc];		# ultimate stress
set eps3 [expr -1.0];			# strain at ultimate stress

set ft1  [expr $ft];
set ft2  [expr 0.01*$ft];
set epst2  [expr 4.*$epscr];
set ft3  [expr 0.0*$ft];
set epst3  [expr 1.0];

uniaxialMaterial ConcretewBeta $concHoriz $fc1 $eps1 $fc2 $eps2 $fc3 $eps3 $ft1 $ft2 $epst2 $ft3 $epst3 -beta 0.3 [expr 0.01] 0.1 [expr 0.025] -alpha 0.5 -lambda $lambda -E $Ec;

uniaxialMaterial Elastic $elasticMat 400;


# # Trusses...
node 1   	0.0  0.0
node 2 	 	1.0  0.0
node 3   	0.5  0.0
node 4 	 	0.5  1.0

# set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
fix 1 1 1 1
fix 2 0 1 1
fix 3 1 1 1
fix 4 1 0 1

element	Truss2 1 1 2 3 4 1 $concHoriz;
element	truss 2 1 2 1 $elasticMat;
element	truss 3 3 4 1 $elasticMat;

# ### Analysis 1 - with beta = 1

set load 1;
set controlNode 2;
pattern Plain 1 Linear -factor 1 {
   load $controlNode $load 0 0;			
}

pattern Plain 2 Constant -factor 1 {
   load 4 0 6 0;			
}

# recorder Node -file "disp2_dof1.out" -time -node 2 -dof 1 disp
# recorder Node -file "disp4_dof2.out" -time -node 4 -dof 2 disp


# --- LOADING LOOP
set file [open "cyclicDisp.dat" r];
set forceFile [open "cyclicForce.dat" r];
set factor 1;
set Tol 1.0e-6;			# convergence tolerance for test

constraints Transformation
numberer RCM
system BandGeneral

set d1 0.;
set Nstep 5;

while {[gets $file temp] > 0} {
	gets $forceFile forceTemp;
	foreach d2 $temp forceVal $forceTemp {
		set d2 [expr {($factor*$d2)}];
		set DLoad [expr {($d2 - $d1)/$Nstep}];
		# puts [expr {$d1 + $Nstep*$DLoad}];
		
		for {set i 0} {$i < $Nstep} {incr i 1} {
			test NormDispIncr $Tol 10 0
			algorithm Newton
			integrator DisplacementControl $controlNode 1 $DLoad
			analysis Static
		    set ok [analyze 1]
			if {$ok != 0} {
				puts "Trying Newton with Initial Tangent.."
				test NormDispIncr $Tol 1000 0
				algorithm Newton -initial
				integrator DisplacementControl $controlNode 1 $DLoad
				analysis Static
				set ok [analyze 1]
			}
			if {$ok != 0} {
				remove recorders;
				return -1
		    }; # end if
		}; #end for
		
		set compForce [getTime];
		if { [expr abs($compForce - $forceVal)] <= 1.e-5} {
		   puts "PASSED"
		} else {
		  puts "FAILED ----"
		  remove recorders;
		  return;
		}
		
		set d1 $d2;
	}
}
puts "PASSED ALL TESTS ---"
remove recorders;
