# create the model builder
model Basic -ndm 2 -ndf 3

# create 2 nodes
node 1   0.0   0.0
node 2 360.0   0.0
node 3   0.0 144.0
node 4 360.0 144.0

# fix node 1
fix 1 1 1 1
fix 2 1 1 1

geomTransf Linear 1
geomTransf Linear 2

source SteelWSections.tcl
set in 1.0
set E 30000

ElasticBeamWSection2d 1 1 3 W14X257 $E 1
ElasticBeamWSection2d 2 3 4 W24X68  $E 2
ElasticBeamWSection2d 3 2 4 W14X257 $E 1

# create time series and load pattern
set P 10000.0
timeSeries Constant 1
pattern Plain 1 1 {
   load 3 0.0 -$P 0.0
   load 4 0.0 -$P 0.0
}

# create an analysis
constraints Plain
numberer RCM
test NormDispIncr 1.0e-12 6 0
algorithm Newton
system ProfileSPD
integrator LoadControl 1.0
analysis Static

# perform the analysis
analyze 1

timeSeries Linear 2
pattern Plain 2 2 {
   load 3 1.0 0.0 0.0
   load 4 1.0 0.0 0.0
}

recorder Node -file node.out -time -node 3 -dof 1 disp

integrator DisplacementControl 3 1 0.1
analyze 100

print node 3