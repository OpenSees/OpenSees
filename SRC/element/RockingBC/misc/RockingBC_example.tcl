# RockingBC Examples - Similar to the ones in Avgenakis E. and Psycharis I.N. (2020) “An integrated macroelement formulation for the dynamic response of inelastic deformable rocking bodies.” Earthquake Engineering and Structural Dynamics.

# ##### INITIALIZATION

# Model properties 

set B 1.0
set w 1.0
set rho 2.5
set mu 0
set g 9.81

set nu 0.0
set Nw 15

set a 0.2	
set zeta 0.05
set e0 5.0e-6 ;# These parameters are defined in: Avgenakis E. and Psycharis I.N. (2020) “An integrated macroelement formulation for the dynamic response of inelastic deformable rocking bodies.” Earthquake Engineering and Structural Dynamics.
set sr 5.0e-3

set H [expr $B/tan($a)]
set E [expr $rho*$g*$H/$e0] 
set sy [expr -$rho*$g*$H/$sr]

set R [expr sqrt($H*$H+$B*$B)/2.]
set p [expr sqrt(3.0*$g/4.0/$R)]
set total_mass [expr $H*$B*$w*$rho]
set rot_inertia [expr 1./12.*$total_mass*($B*$B+$H*$H)]

# RockingBC element Parameters

set convlim 1.0e-15
set useshear 0
set blevery 1
set useUelNM 1
set usecomstiff 0

set af 1.0
set aflim 0.4
set convlimmult 1.0
set maxtries 50

set NlimN 0.001
set NlimT 1.0
set DtMax 1.0e-5
set DtMin 1.0e-8
set errorifNexceeds 1 ;# A non-default value is used here, together with a timestep correction strategy at the end of the script for better accuracy during impacts

# Analysis Parameters

set Tol 1.0e-10
set NumIter 100
set nEigen 1

set DtfacL 0.5
set DtfacU 1.5
set incrafter 10
set Nlimratio 0.9

# Analysis cases

set THcase "Free" ; # CHANGE ANALYSIS TYPE AS FOLLOWS

set pi [expr acos(-1.0)]
if {$THcase=="Free"} {
	set ur 0.5
	set THousnerratio 3.0
	set th0 [expr $ur*$a]	
	set ang [expr 1.0/(1.0-$th0/$a)]
	set Thousner [expr 4.0/$p*log($ang+sqrt($ang*$ang-1))]
	puts "Thousner = $Thousner"
	set TmaxAnalysis [expr $THousnerratio*$Thousner]
} elseif {$THcase=="Sinepulse"} {
	set wpratio 6.0
	set apratio 2.0
	set Tpratio 18.0
	set wp [expr $wpratio*$p]
	set ap [expr $apratio*$a*$g]
	set Tp [expr 2.0*$pi/$wp]
	set TmaxAnalysis [expr $Tpratio*$Tp]	
} elseif {$THcase=="EQ"} {
	set EQfile "Northridge.dat"
	set skiprows 5
	set accmult 9.81
} else {
	error "Error: No such timehistory case"
}

set name "$THcase"
file mkdir $name

# ##### BUILD MODEL

wipe

model basic -ndm 2 -ndf 3

node 1   0.0  0.0
node 2	 0.0  [expr 0.5*$H]
node 3	 0.0 $H

element RockingBC 1 2 1 $Nw $E $nu $sy $B $w $mu -convlim $convlim -useshear $useshear -blevery $blevery -useUelNM $useUelNM -usecomstiff $usecomstiff -af $af -aflim $aflim -convlimmult $convlimmult -maxtries $maxtries -NlimN $NlimN -NlimT $NlimT -Dtlim $DtMin -errorifNexceeds $errorifNexceeds

geomTransf Corotational 1
element elasticBeamColumn 2 2 3 [expr $B*$w] $E [expr 1./12.*$w*$B*$B*$B] 1

fix 1 1 1 1

mass 2 $total_mass $total_mass $rot_inertia

pattern Plain 1 Linear {
    load 2 0.0 [expr -$total_mass*$g] 0.0
}

recorder Node -file $name\\topdisps.txt -time -nodes 3 -dof 1 2 3 disp
recorder Node -file $name\\centdisps.txt -time -nodes 2 -dof 1 2 3 disp
recorder Node -file $name\\centvel.txt -time -nodes 2 -dof 1 2 3 vel
recorder Node -file $name\\reaction.txt -time -nodes 1 -dof 1 2 3 reaction
recorder Element -file $name\\sL.txt -time -ele 1 sL ; # Sliding at the rocking end
recorder Element -file $name\\Nforces.txt -time -ele 1 basicForce ; # Corotational system forces
recorder Element -file $name\\forces.txt -time -ele 1 force ; # Global system forces
# recorder Element -file $name\\dummy -time -ele 1 $name\\Dist ; # Uncomment to record Stress and Plastic displacement distributions each step inside the folder given at the end of the command (comment the same command at the end of the script)

constraints Transformation
numberer RCM
system UmfPack
test EnergyIncr $Tol $NumIter
algorithm Newton
set NstepGravity 10
set DGravity [expr 1./$NstepGravity]
integrator LoadControl $DGravity
analysis Static
analyze $NstepGravity

puts "Dead load analysis complete"
loadConst -time 0.0

set lambdaN [eigen 1]
set lambdaI [lindex $lambdaN [expr $nEigen-1]]
set omegaI [expr pow($lambdaI,0.5)]
set zetab $zeta
set betaKcurr [expr 2.*$zeta/$omegaI]

set T1 [expr 2.*$pi/$omegaI]
puts "T1 = $T1"   

set infofile [open $name\\info.txt w+]   
puts $infofile "H $H"
puts $infofile "B $B"
puts $infofile "w $w"
puts $infofile "a $a"
puts $infofile "R $R"
puts $infofile "p $p"
puts $infofile ""
puts $infofile "rho $rho"
puts $infofile "total_mass $total_mass"
puts $infofile "rot_inertia $rot_inertia"
puts $infofile ""
puts $infofile "e0 $e0"
puts $infofile "sr $sr"
puts $infofile "E $E"
puts $infofile "nu $nu"
puts $infofile "sy $sy"
puts $infofile ""
puts $infofile "Nw $Nw"
puts $infofile "beta $betaKcurr"
puts $infofile "zeta $zetab"
puts $infofile "nEigen $nEigen"
puts $infofile "mu $mu"
puts $infofile ""
puts $infofile "T1 $T1"
puts $infofile "DtMax $DtMax"
puts $infofile "DtMin $DtMin"
puts $infofile "convlim $convlim"
puts $infofile "useshear $useshear"
puts $infofile "blevery $blevery"
puts $infofile "NlimN $NlimN"
puts $infofile "NlimT $NlimT"
puts $infofile ""
if {$THcase=="Free"} {
	puts $infofile "th0 $th0"
	puts $infofile "ur $ur"
	puts $infofile "Thousner $Thousner"
	puts $infofile "TmaxAnalysis $TmaxAnalysis"
} elseif {$THcase=="Sinepulse"} {
	puts $infofile "wp $wp"
	puts $infofile "ap $ap"
} elseif {$THcase=="EQ"} {
	puts $infofile "EQfile $EQfile"
	puts $infofile "accmult $accmult"
} else {
	error "Error: No such timehistory case"
}

close $infofile

# ##### PUSHOVER IF NEEDED

if {$THcase=="Free"} {

	pattern Plain 2 Linear {
		load 2 1.0 0.0 0.0
		}

	set ctrlNode 2
	set ctrlDOF 3
	set Dmax [expr -$th0]
	set Dincr [expr 0.005*$Dmax]

	constraints Transformation
	numberer RCM
	system UmfPack
	test EnergyIncr $Tol $NumIter
	algorithm Newton
	analysis Static

	set Dstep 0.0
	set ok 0
	set Nk 100
	
	while {$Dstep < 1.0 && $ok == 0} {

		set controlDisp [lindex [nodeDisp $ctrlNode $ctrlDOF] 0]
		set Dstep [expr $controlDisp/$Dmax]
		if {$Dmax>0} {
			set Dincrcur [expr min(abs($Dincr),abs($Dmax-$controlDisp))]
		} else {
			set Dincrcur [expr -1.0*min(abs($Dincr),abs($Dmax-$controlDisp))]
		}
		
		integrator DisplacementControl $ctrlNode $ctrlDOF $Dincrcur		
		set ok [analyze 1]

		if {$ok != 0} {
			incr ierrors
			set DincrReduced [expr $Dincrcur/$Nk];
			integrator DisplacementControl $ctrlNode $ctrlDOF $DincrReduced
			for {set ik 1} {$ik <=$Nk} {incr ik 1} {
				set ok [analyze 1]
			}
		}
	}

	set infofile [open $name\\success.txt w+]   
	puts $infofile "pushover success $ok"
	close $infofile
	
	puts "Pushover complete"
} else {
	set infofile [open $name\\success.txt w+]   
	puts $infofile "nopushover success 0"
	close $infofile
}

# ##### DYNAMIC ANALYSIS

setTime 0.0
record
rayleigh 0 $betaKcurr 0 0; # Use only stiffness-dependent damping for the model

set TS [list]
set SS [list]
set pi [expr acos(-1.0)]

if {$THcase=="Free"} {
	remove loadPattern 2 

	set Nsteps [expr int($TmaxAnalysis/$DtMax)]
	for {set i 0} {$i < $Nsteps} {incr i} {
		lappend TS [expr $i*$DtMax]
		lappend SS [expr 0.0]
	}
	
} elseif {$THcase=="Sinepulse"} {

	set Nstepspulse [expr int($Tp/$DtMax)]
	set Dtpulse [expr $Tp/$Nstepspulse]
	
	for {set i 0} {$i < $Nstepspulse} {incr i} {
		lappend TS [expr $i*$Dtpulse]
		lappend SS [expr $ap*sin($wp*$i*$Dtpulse)]
	}	
	while {[lindex TS end]<$TmaxAnalysis} {
		lappend TS [expr [lindex TS end]+$DtMax]
		lappend SS [expr 0.0]
	}

} else {

	set fp [open $EQfile r]
	set file_data [read $fp]
	close $fp
	set data [split $file_data "\n"]

	set ifd 0
	foreach line $data {
		if {$ifd<$skiprows} {
			set ifd [expr $ifd+1]
		} else {
			set lsp [split $line]
			lappend TS [lindex $lsp 0]
			lappend SS [expr $accmult*[lindex $lsp 1]]
		}
	}
	set TmaxAnalysis [lindex $TS end]

}

timeSeries Path 200 -time $TS -values $SS -factor 1.0
pattern UniformExcitation 300 1 -accel 200
constraints Transformation
numberer RCM
system UmfPack
algorithm Newton
integrator Newmark 0.5 0.25
analysis Transient
test RelativeNormUnbalance $Tol $NumIter
algorithm Newton

set ok 0

set curTime [getTime]
set Dtcur $DtMax
set igood 0

while {$curTime < $TmaxAnalysis && $ok == 0} {
	
	set ok [analyze 1 $Dtcur]
	
	if {$ok != 0 && $Dtcur==$DtMin} { # Since the element does not throw an error due to large axial force changes when $Dtcur==$DtMin, the problem lies elsewhere and the analysis is stopped prematurely
		break
	} elseif {$ok!=0} {
		set Dtcur [expr max($DtMin,$DtfacL*$Dtcur)]; # If an error is detected, possibly due to a large axial force change inside the RockingBC element, $Dt is lowered until the minimum value
		set ok 0
		set igood 0
	} else {
		incr igood
		if {$igood>=$incrafter} { # Check that some successful iterations have been performed
			if {[lindex [eleResponse 1 forceratioNmax] 0]<$Nlimratio*$NlimN && [lindex [eleResponse 1 forceratioTmax] 0]<$Nlimratio*$NlimT} { # Check that the maximum axial force ratios are somewhat lower than the respective limits
				set Dtcur [expr min($DtMax,$DtfacU*$Dtcur)] ;# The timestep may be gradually increased again until the maximum value
			}
		}
	}
	
	set curTime [getTime]
	
}

set infofile [open $name\\success.txt a+]  
puts $infofile "timehistory success $ok"
close $infofile

puts "Timehistory complete"

recorder Element -file $name\\dummy -time -ele 1 $name\\Dist ; #Record Stress and Plastic displacement distributions at the last step inside the folder given at the end of the command
record

wipe

