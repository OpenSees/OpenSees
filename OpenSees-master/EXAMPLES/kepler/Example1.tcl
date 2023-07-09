#input eigenvalue & peermotion

wipe

set PI 3.14159265
set g 386.4 
set M 1.0;
set A 1.0;
set L 1.0;
set loc 0;

set gMotion ABRUZZO/ATI-WE.AT2
set gMotionList [split $gMotion "/"]
set gMotionDir  [lindex $gMotionList end-1]
set gMotionNameInclAT2 [lindex $gMotionList end]
set gMotionName [string range $gMotionNameInclAT2 0 end-4 ]
     
set Gaccel "PeerDatabase $gMotionDir $gMotionName -accel $g -dT dT -nPts nPts"


model basic -ndm 1 -ndf 1
set K [expr $eigenvalue * $M];
set E [expr $K * $L / $A];
node 1 0.0
node 2 $L -mass $M
fix 1 1
uniaxialMaterial Elastic 1 $E
element truss 1 1 2 $A 1
pattern UniformExcitation 2 1 -accel $Gaccel
recorder EnvelopeNode -file envelope.out -node 2 -dof 1 disp
system ProfileSPD
constraints Plain
numberer Plain
algorithm Linear
integrator Newmark 0.5 0.25
analysis Transient
analyze 2000 $dT;
remove recorders
if [catch {open envelope.out r} InFileID] {
    puts "0.0"
}
set min [gets $InFileID]
set max [gets $InFileID]
set absMax [gets $InFileID]
close $InFileID
set last "$absMax"
