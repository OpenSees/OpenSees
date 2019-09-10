wipe
model basic -ndm 3 -ndf 6

puts "node"
node 1          0.   1.0   0.
node 2          0.   0.5   0.
node 3          0.    0.   0.
node 4         1.0    1.   0.
node 5         1.0   0.5   0.
node 6         1.0    0.   0.
node 7         2.0    1.   0.
node 8         2.0   0.5   0.
node 9         2.0    0.   0.
node 10        3.0    1.   0.
node 11        3.0   0.5   0.
node 12        3.0    0.   0.
node 13        4.0    1.   0.
node 14        4.0   0.5   0.
node 15        4.0    0.   0.
node 16         5.    1.   0.
node 17         5.   0.5   0.
node 18         5.    0.   0.
node 19        6.0    1.   0.
node 20        6.0   0.5   0.
node 21        6.0    0.   0.
node 22        7.0    1.   0.
node 23        7.0   0.5   0.
node 24        7.0    0.   0.
node 25        8.0    1.   0.
node 26        8.0   0.5   0.
node 27        8.0    0.   0.
node 28        9.0    1.   0.
node 29        9.0   0.5   0.
node 30        9.0    0.   0.
node 31         10.    1.   0.
node 32         10.   0.5   0.
node 33         10.    0.   0.

fix 1 1 1 1 1 1 1;
fix 2 1 1 1 1 1 1;
fix 3 1 1 1 1 1 1;

puts "Material"
nDMaterial ElasticIsotropic 2 1.200E+006 0.0
nDMaterial PlateFiber 601 2
section PlateFiber 701 601 0.1
#section   LayeredShell       701       10     2 0.005 2 0.005 2 0.005 2 0.005 2 0.005 2 0.005 2 0.005 2 0.005 2 0.005 2 0.005      
puts "shell element"
element ShellNLDKGQ  1  2  5  4  1 701
element ShellNLDKGQ  2  3  6  5  2 701
element ShellNLDKGQ  3   5  8  7  4 701
element ShellNLDKGQ  4   6  9  8  5 701
element ShellNLDKGQ  5   8 11 10  7 701
element ShellNLDKGQ  6   9 12 11  8 701
element ShellNLDKGQ  7  11 14 13 10 701
element ShellNLDKGQ  8  12 15 14 11 701
element ShellNLDKGQ  9  14 17 16 13 701
element ShellNLDKGQ  10  15 18 17 14 701
element ShellNLDKGQ  11  17 20 19 16 701
element ShellNLDKGQ  12  18 21 20 17 701
element ShellNLDKGQ  13  20 23 22 19 701
element ShellNLDKGQ  14  21 24 23 20 701
element ShellNLDKGQ  15  23 26 25 22 701
element ShellNLDKGQ  16  24 27 26 23 701
element ShellNLDKGQ  17  26 29 28 25 701
element ShellNLDKGQ  18  27 30 29 26 701
element ShellNLDKGQ  19  29 32 31 28 701
element ShellNLDKGQ  20  30 33 32 29 701

puts "recorder"
recorder Node -file disp1.txt -time -nodeRange 1 33 -dof 1 disp;
recorder Node -file disp3.txt -time -nodeRange 1 33 -dof 3 disp;
recorder Node -file reaction1.txt -time -nodeRange 1 3 -dof 1 2 3 reaction;
puts "loading"
pattern Plain 1 Linear {
load 32 0.0  0.0 0.0 0 15.702963 0;
}
puts "analysis"
constraints Plain
numberer RCM
system BandGeneral
test NormDispIncr 1e-3 1000 2
algorithm KrylovNewton
integrator LoadControl 0.001
analysis Static
analyze 4000