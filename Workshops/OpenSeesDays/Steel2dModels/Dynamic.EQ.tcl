# ------------------------------------------------------------------------------
# Dynamic Earthquake Analysis
# Created by:  Vesna Terzic, UC Berkeley, 2011                               
# execute this file after you have built the model, and after you apply gravity
# ------------------------------------------------------------------------------


# source in procedures
source ReadSMDfile.tcl;		# procedure for reading GM file and converting it to proper format

# Uniform Earthquake ground motion (uniform acceleration input at all support nodes)
set GMdirection 1;													# ground-motion direction
set GMfileH [format "Oak_%i_50_%i_%s" $Hazard $i $eqDir] ;			# ground-motion filenames: horizntal component
set GMfileV [format "Oak_%i_50_%i_UP" $Hazard $i] ;					# ground-motion filenames: vertical component
set GMfact 1.0;														# ground-motion scaling factor


#############################################################################
#                     Define & Apply Damping
#############################################################################
# RAYLEIGH damping parameters 
# D=$alphaM*M + $betaKcurr*Kcurrent + $betaKcomm*KlastCommit + $beatKinit*$Kinitial

if { $Hazard == 2 } {
    set xDamp 0.03; 			# damping ratio
    set Tcoeff 1.5;             # coefficient that modifies the period for which the Railigh coefficients are calculated
} elseif { $Hazard == 10 } {	 
    set xDamp 0.03
    set Tcoeff 1.5
} elseif { $Hazard == 50 } {	 
    set xDamp 0.03
    set Tcoeff 1.0
}	

# define RAYLEIGH coefficients for Ti=1.5*T1 (or Ti=T1) and Tj=T3			
set MpropSwitch 1.0;
set KcurrSwitch 0.0;
set KcommSwitch 1.0;
set KinitSwitch 1.0;

set lambda [eigen 3]


#puts "EIGENvalues: $lambda"

#columns
set omegaI [expr pow([lindex $lambda 0],0.5)/$Tcoeff];
set omegaJ [expr pow([lindex $lambda 2],0.5)];
set alphaM [expr $MpropSwitch*$xDamp*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	# M-prop. damping; D = alphaM*M
set betaKcurr [expr $KcurrSwitch*2.*$xDamp/($omegaI+$omegaJ)];         			# current-K;      +beatKcurr*KCurrent
set betaKcomm [expr $KcommSwitch*2.*$xDamp/($omegaI+$omegaJ)];   				# last-committed K;   +betaKcomm*KlastCommitt
set betaKinit [expr $KinitSwitch*2.*$xDamp/($omegaI+$omegaJ)];         			# initial-K;     +beatKinit*Kini


rayleigh $alphaM 0. 0. $betaKcomm

# Uniform EXCITATION: acceleration input
set IDloadTag 400;	# for uniformSupport excitation	
set inFileH $GMdir$GMfileH.acc
set inFileV $GMdir$GMfileV.acc
set outFileH $GMdir$GMfileH.g3
set outFileV $GMdir$GMfileV.g3
ReadSMDFile $inFileH $outFileH dt nPt;		# call procedure to convert the horizontal ground-motion file
ReadSMDFile $inFileV $outFileV dt nPt;		# call procedure to convert the vertical ground-motion file
set GMfatt [expr $g*$GMfact];				# data in input file is in g Unifts -- ACCELERATION TH

timeSeries Path 10 -dt $dt -filePath $outFileH -factor  $GMfatt;   # horizonatal time series information
timeSeries Path 11 -dt $dt -filePath $outFileV -factor  $GMfatt;   # veritcal time series information

pattern UniformExcitation $IDloadTag $GMdirection -accel 10;   # create Unifform excitation for horizontal GM
pattern UniformExcitation [expr $IDloadTag+1] 2 -accel 11;     # create Unifform excitation for horizontal GM 

# record max absolute floor accelerations
recorder EnvelopeNode -file [format "$dataDir/$subDir1/FloorAccEnv_%i_%s.out" $i $eqDir] -timeSeries 10 -node 11 12 13 14 -dof 1 accel
recorder Node -file [format "$dataDir/$subDir1/FloorAcc_%i_%s.out" $i $eqDir] -timeSeries 10 -node 11 12 13 14 -dof 1 accel	

recorder EnvelopeNode -file [format "$dataDir/$subDir1/FloorDispEnv_%i_%s.out" $i $eqDir] -node 11 12 13 14 64 74 -dof 1 disp
recorder Node -file [format "$dataDir/$subDir1/FloorDisp_%i_%s.out" $i $eqDir] -node 11 12 13 14 -dof 1 disp	

set ViewScale 15;
#DisplayModel2D DeformedShape $ViewScale ;	# display deformed shape, the scaling factor needs to be adjusted for each model

# Define analysis objects and performe analysis
set tFinal	[expr $dt*$nPt];	# maximum duration of ground-motion analysis
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1.0e-6 10 
algorithm Newton
integrator Newmark 0.5 0.25
analysis VariableTransient

set ok 0
set currentTime 0.0
while {$ok == 0 && $currentTime < $tFinal} {
    set ok [analyze 1 $dt]
    if {$ok != 0} {
	test NormDispIncr 0.5e-3 2000 1
	algorithm ModifiedNewton -initial
	set ok [analyze 1 0.005 0.0005 0.01 6]
	test NormDispIncr 1.0e-6 10 
	algorithm Newton
    } 
    set currentTime [getTime]
}

#free vibrations with 20 times bigger damping ratio
set MFD 20.; #multiplication factor for damping
set alphaM [expr $MpropSwitch*($xDamp*$MFD)*(2*$omegaI*$omegaJ)/($omegaI+$omegaJ)];	# M-prop. damping; D = alphaM*M
set betaKcurr [expr $KcurrSwitch*2.*($xDamp*$MFD)/($omegaI+$omegaJ)];         			# current-K;      +beatKcurr*KCurrent
set betaKcomm [expr $KcommSwitch*2.*($xDamp*$MFD)/($omegaI+$omegaJ)];   				# last-committed K;   +betaKcomm*KlastCommitt
set betaKinit [expr $KinitSwitch*2.*($xDamp*$MFD)/($omegaI+$omegaJ)];         			# initial-K;     +beatKinit*Kini


rayleigh $alphaM 0. 0. $betaKcomm	

set tFinal2 [expr $dt*500.+$currentTime]
while {$ok == 0 && $currentTime < $tFinal2} {
    set ok [analyze 1 $dt]
    if {$ok != 0} {
	test NormDispIncr 0.5e-3 2000 1
	algorithm ModifiedNewton –initial
	set ok [analyze 1 0.005 0.0005 0.01 6]
	test NormDispIncr 1.0e-6 10 
	algorithm Newton
    } 
    set currentTime [getTime]
    #puts "time :$currentTime"
}

puts "Ground Motion Done. End Time: [getTime]. tFinal: $tFinal. Run #$i"
