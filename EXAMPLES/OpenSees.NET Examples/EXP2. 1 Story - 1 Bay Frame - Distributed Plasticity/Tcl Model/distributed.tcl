pattern Plain 2 Linear {
eleLoad -ele 21 22 -type -beamUniform -5e4
}
constraints Plain
numberer Plain
system BandGeneral
test NormDispIncr 1.e-8 6
algorithm ModifiedNewton
integrator LoadControl 1
analysis Static
analyze 1
loadConst -time 0.0