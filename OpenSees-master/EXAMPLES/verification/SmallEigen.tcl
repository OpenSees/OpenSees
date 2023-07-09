
# Small 3dof eigen example
#
# 3 springs in series, K = 30,20,10, M = 11.1,22.2,33.3
#
# results as computed by matlab: 0.139300647653736 1.05865310768512 4.9582024008173

# using following 3 commands:
# format longG
# K=[[50 -20 0];[-20 30 -10];[0 -10 10]]
# M=[[11.1 0 0];[0 22.2 0];[0 0 33.3]]
# eig(K,M)

puts "SmallEigen.tcl: Small 3 dof 2 spring model: checking eigen results of smallest and largest values"
puts "  using masses defined through node and mass command, results compared for exact against Matlab\n"

set exactResults {0.139300647653736 1.05865310768512  4.9582024008173}

wipe
model Basic -ndm 1 -ndf 1
node 1 0
node 2 0 -mass 11.1
node 3 0 -mass 22.2
node 4 0 -mass 33.3
uniaxialMaterial Elastic 1 30
uniaxialMaterial Elastic 2 20
uniaxialMaterial Elastic 3 10
fix 1 1
element zeroLength 1 1 2 -mat 1 -dir 1
element zeroLength 2 2 3 -mat 2 -dir 1
element zeroLength 3 3 4 -mat 3 -dir 1

puts "Checking Computation of Smallest Eigenvalues:"
set testOK 0

set tol 1.0e-6
set numEigen 2
set eigenValues [eigen $numEigen]

set formatString {%15s%15s}
puts [format $formatString OpenSees Exact]
set formatString {%15.6f%15.6f}
for {set i 0} {$i<$numEigen} {incr i 1} {
    set OpenSeesR [lindex $eigenValues  $i]
    set exactR    [lindex $exactResults $i]
    puts [format $formatString $OpenSeesR  $exactR]

    if {[expr abs($OpenSeesR-$exactR)] > $tol} {
	set testOK -1;
	puts "failed-> [expr abs($OpenSeesR-$exactR)] $tol"
    }
}

puts "\nChecking Computation of Largest Eigenvalues:"
puts [eigen -findLargest $numEigen]
set eigenValues [eigen -findLargest $numEigen]

set formatString {%15s%15s}
puts [format $formatString OpenSees Exact]
set formatString {%15.6f%15.6f}
for {set i 0} {$i<$numEigen} {incr i 1} {
    set OpenSeesR [lindex $eigenValues  $i]
    set exactR    [lindex $exactResults [expr $i+1]]
    puts [format $formatString $OpenSeesR  $exactR]

    if {[expr abs($OpenSeesR-$exactR)] > $tol} {
	set testOK -1;
	puts "failed-> [expr abs($OpenSeesR-$exactR)] $tol"
    }
}


wipe
wipe
model Basic -ndm 1 -ndf 1
node 1 0
node 2 0
node 3 0
node 4 0
uniaxialMaterial Elastic 1 30
uniaxialMaterial Elastic 2 20
uniaxialMaterial Elastic 3 10
fix 1 1
element zeroLength 1 1 2 -mat 1 -dir 1
element zeroLength 2 2 3 -mat 2 -dir 1
element zeroLength 3 3 4 -mat 3 -dir 1

mass 2 11.1
mass 3 22.2
mass 4 33.3

puts "Checking Computation of Smallest Eigenvalues:"
set testOK 0

set tol 1.0e-6
set numEigen 2
set eigenValues [eigen $numEigen]

set formatString {%15s%15s}
puts [format $formatString OpenSees Exact]
set formatString {%15.6f%15.6f}
for {set i 0} {$i<$numEigen} {incr i 1} {
    set OpenSeesR [lindex $eigenValues  $i]
    set exactR    [lindex $exactResults $i]
    puts [format $formatString $OpenSeesR  $exactR]

    if {[expr abs($OpenSeesR-$exactR)] > $tol} {
	set testOK -1;
	puts "failed-> [expr abs($OpenSeesR-$exactR)] $tol"
    }
}

puts "\nChecking Computation of Largest Eigenvalues:"
set eigenValues [eigen -findLargest $numEigen]

set formatString {%15s%15s}
puts [format $formatString OpenSees Exact]
set formatString {%15.6f%15.6f}
for {set i 0} {$i<$numEigen} {incr i 1} {
    set OpenSeesR [lindex $eigenValues  $i]
    set exactR    [lindex $exactResults [expr $i+1]]
    puts [format $formatString $OpenSeesR  $exactR]

    if {[expr abs($OpenSeesR-$exactR)] > $tol} {
	set testOK -1;
	puts "failed-> [expr abs($OpenSeesR-$exactR)] $tol"
    }
}

set results [open results.out a+]
if {$testOK == 0} {
    puts "PASSED Verification Test SmallEigen.tcl \n\n"
    puts $results "PASSED : SmallEigen.tcl"
} else {
    puts "FAILED Verification Test SmallEigen.tcl \n\n"
    puts $results "FAILED : SmallEigen.tcl"
}
close $results
