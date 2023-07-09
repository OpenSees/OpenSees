model Basic -ndm 2 -ndf 2

set thk 0.01
set L 5.0
set H 1.0

set P 100.
set E 200.e6
set nu 0.3

node 1 0. 0.
node 2 5. 0.
node 3 5. 1.
node 4 0. 1.
node 5 2.5 0.
node 6 5. .5
node 7 2.5 1.
node 8 0. .5
# node 9 2.5 .5; # comment for quad8n element

nDMaterial ElasticIsotropic 1 $E $nu

# element quad9n 1 1 2 3 4 5 6 7 8 9 $thk "PlaneStress" 1
element quad8n 1 1 2 3 4 5 6 7 8 $thk "PlaneStress" 1

fix 1 1 1
fix 4 1 0
fix 8 1 0

timeSeries Linear 1

pattern Plain 1 1 {
	load 2 $P 0.
	load 3 -$P 0.
}

analysis Static

analyze 1

print

# verification:
# tip vertical displacement (node 2 and 3) = 0.0075
# bottom Gauss Point stress_xx = 46475.8
# bottom extreme stress_xx (extrapolated) = 60000.0

exit
