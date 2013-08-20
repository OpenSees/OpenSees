# create the model builder
model Basic -ndm 1 -ndf 1

# create 2 nodes
node 1 0
node 2 1.0

# fix node 1
fix 1 1

# create material
set Fy 60.0
set E 30000.0
set b -0.1
uniaxialMaterial Steel01 1 $Fy $E $b

# create element
set A 1.0
element Truss 1 1 2 $A 1

# create time series and load pattern
timeSeries Linear 1
set P 100.0
pattern Plain 1 1 {
   load 2 $P
}

# create an analysis
constraints Plain
numberer RCM
test NormDispIncr 1.0e-12 6 0
algorithm Newton
system BandGen
integrator DisplacementControl 2 1 0.001
analysis Static

# perform the analysis & print the results
recorder Node -file node.out -time -node 2 -dof 1 disp
record 

foreach {numIter dU} {10 0.001 20 -0.001 10 0.001} {
    integrator DisplacementControl 2 1 $dU
    analyze $numIter
    set factor [getTime]
    puts "[expr $factor*$P]  [lindex [nodeDisp 2] 0]"
}

print node 2
