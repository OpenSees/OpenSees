pattern Plain 3 Linear {
load 2 0. -100000. 0.

}

constraints Plain
numberer Plain
system BandGeneral
test NormDispIncr 1.e-8 6
algorithm Newton
integrator DisplacementControl 2 2 -0.01
analysis Static
analyze 50

print node 2