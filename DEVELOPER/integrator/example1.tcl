# Example: 1 element truss model subjcted to sine wave input motion at base
#
# model  Uf,Af --> |______/\/\/\/\/\_____**  M  -->U,V,A
#                  |          K          **
#
# subject to motion:  Uf = USinwt
#
#     | M  0 | A   + | K  -K|U  = 0
#     | 0  M2| Af    | -K  K|Uf   0
# 
# rewriting eqn 1 of equtaions:
#
#       MA + KU = KU Sinwt  subject to: U(0) = 0, V(0) = 0;
#
# Exact Solution:
#
#     U(t) = ____U_____ Sin(wt) -  _U_ (w/wn)_ Sin(wnt)
#            (1-(w/wn)^2)         (1-(w/wn)^2)
# 
# 
# Uniform Excitation Case: Af=-Uw^2 Sinwt   Vf(0) = -Uw
#
#     U(t) = Uf(t) + Ur(t)
#
#       MAr + KUr = -MAf = MUw^2 Sinwt, Ur(0) = 0; Vr(0) = -Uw 
#
# Exact Solution for Uniform Case:
#
#     Ur(t) = ___MUw^2____ Sin(wt) -  Uw Sin(wnt) -  _MUw^2 (w/wn)_ Sin(wnt)
#             K(1-(w/wn)^2)           wn              K(1-(w/wn)^2)
#
#           = _____Uw^2______ Sin(wt) -  Uw Sin(wnt) - __Uw^2_ (w/wn)__   Sin(wnt)
#             wn^2(1-(w/wn)^2)           wn             wn^2 (1-(w/wn)^2)
#
#           = _____Uw^2______ Sin(wt) -  ____ U w ____  Sin(wnt)
#             wn^2(1-(w/wn)^2)           wn(1-(w/wn)^2)
#
#     U(t) = Ur(t) + Uf(t)
#
#     U(t) = ____U_____ Sin(wt) -  _U_ (w/wn)_ Sin(wnt)
#            (1-(w/wn)^2)         (1-(w/wn)^2)            
#

proc doModel {K periodN U period tFinish inputT} {

    # some constants
    set g 386.4
    set PI 3.14159265

    # some clculated parameters for basic model
    set omega [expr 2*$PI/$period]
    set omegaN [expr 2*$PI/$periodN]
    set M [expr $K/($omegaN*$omegaN)]
    set factor $U

    puts "$period $periodN $omega $omegaN"

    #
    # Create Model
    #

    model BasicBuilder -ndm 1 -ndf 1;		# Define the model builder, ndm=#dimension, ndf=#dofs
    wipe
    node 1 0 			# node#, X, Y
    node 2 1.0 -mass $M
    uniaxialMaterial Elastic 1 $K
    element truss 1 1 2 1 1

    #
    # load pattern
    #

    if {$inputT == "Uniform" } {
	set vel0 [expr -1.0*$U*$omega]
#	set vel0 0.0
	puts "vel0 $vel0"
	set factor [expr -1.0*$U*$omega*$omega]
	set AccelSeries "Sine 0.0 $tFinish $period -factor $factor";
	pattern UniformExcitation 1 1 -accel  $AccelSeries -vel0 $vel0  ;		# create Unifform excitation
	fix 1 1
	recorder Node    -file DUniform.out -time -node 2 -dof 1 disp;		# displacements of free nodes
#	recorder plot DUniform.out  "DISP"  10 10 625 450 -dT [expr 10*0.01] -columns 1 2 -columns 1 3

    } else {

	set DispSeries "Sine 0.0 $tFinish $period -factor $factor";
	
	pattern MultipleSupport 1  {
	    groundMotion 1 Plain -disp  $DispSeries  
	    imposedMotion 1  1 1
	};    
	recorder Node -file DMultiSupport.out -time -node 1 2 -dof 1 disp;		# displacements of free nodes
#	recorder plot DMultiSupport.out  "DISP"  10 10 625 450 -dT [expr 10*0.01] -columns 1 2 -columns 1 3
    }

    #
    # Create analysis
    #

#    integrator HHT 1.0
#    integrator Newmark 0.25 0.5
    integrator Transient Trapezoidal
    system ProfileSPD
    algorithm Linear
    constraints Transformation
    analysis Transient

    # 
    # Perform Analysis
    #
    
    set dt [expr $periodN/100]
    set nSteps [expr int($tFinish/$dt)]

    analyze $nSteps $dt
    print node 1 2
}

#
# set the parameters & use the procedure
# 

set K 100.0
set periodN 2.0

set period 4.0
set U 10.0
set tFinish 6.0

set inputT MultiSupport

doModel $K $periodN $U $period $tFinish $inputT

set inputT Uniform

doModel $K $periodN $U $period $tFinish $inputT

#
# create Matlab script:
#

set matlabF [open matlab.m w]
puts $matlabF "load DUniform.out"
puts $matlabF "load DMultiSupport.out"
puts $matlabF "a=DMultiSupport;"
puts $matlabF "b=DUniform;"
puts $matlabF "w=2*pi/$period"
puts $matlabF "wn=2*pi/$periodN"
puts $matlabF "U=$U"
puts $matlabF "exact=(U/(1-(w/wn)*(w/wn)))*sin(w*a(:,1)) - ((U*w/wn)/(1-(w/wn)*(w/wn)))*sin(wn*a(:,1));";
puts $matlabF "plot(a(:,1),a(:,3),a(:,1),b(:,2),'*',a(:,1),a(:,2)+b(:,2),a(:,1),exact(:,1));"
puts $matlabF "grid on;"
puts $matlabF "legend(\"MultiSupport\",\"Uniform\",\"Ground+Uniform\",\"Exact\")";
close $matlabF;


exit

