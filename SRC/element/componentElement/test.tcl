set outFile [open res.out a]
setPrecision 10

set doSimple FALSE;

if {$doSimple == "TRUE"} {
    model Basic -ndm 2 
    node 1 0 0
    node 2 0 10
    node 3 20 0
    node 4 20 10.
    
    set A 10
    set E 100
    set I 1000.
    geomTransf Linear 1
    element elasticBeamColumn 1 1 2 $A $E $I 1 
    element elasticBeamColumn 2 3 4 $A $E $I 1 
    element elasticBeamColumn 3 2 4 $A $E $I 1 
    
    fix 1 1 1 1
    fix 3 1 1 1
    
    timeSeries Linear 1
    pattern Plain 1 1 {
	load 2 10 0 0
    }
    test NormDispIncr 1e-10 2 4
    algorithm Newton
    system FullGeneral
    integrator LoadControl 1.
    analysis Static
    analyze 1
    print node 2
    
    wipe
    model Basic -ndm 2 
    node 1 0 0
    node 2 0 10
    node 3 20 0
    node 4 20 10.
    
    uniaxialMaterial Elastic 1 1e10
    uniaxialMaterial Elastic 2 1e10
    
    set A 10
    set E 100
    set I 1000.
    geomTransf Linear 1
    element componentElement2d 1 1 2 $A $E $I 1 1 2
    element componentElement2d 2 3 4 $A $E $I 1 1 2
    element elasticBeamColumn 3 2 4 $A $E $I 1 
    
    fix 1 1 1 1
    fix 3 1 1 1
    
    timeSeries Linear 1
    pattern Plain 1 1 {
	load 2 10 0 0
    }
    test NormDispIncr 1e-6 3 4
    algorithm Newton
    system FullGeneral
    integrator LoadControl 1
    analysis Static
    analyze 100
    print node 
    
    wipe
exit
    model Basic -ndm 2 
    node 1 0 0
    node 2 0 10
    node 3 20 0
    node 4 20 10.
    
    uniaxialMaterial Elastic 1 1e10
    uniaxialMaterial Elastic 2 0

    set A 10
    set E 100
    set I 1000.
    geomTransf Linear 1
    element componentElement2d 1 2 1 $A $E $I 1 1 2
    element componentElement2d 2 4 3 $A $E $I 1 1 2
    element elasticBeamColumn 3 2 4 $A $E $I 1 
    
    fix 1 1 1 1
    fix 3 1 1 1
    
    timeSeries Linear 1
    pattern Plain 1 1 {
	load 2 10 0 0
    }
    test NormDispIncr 1e-6 3 4
    algorithm Newton
    system FullGeneral
    integrator LoadControl 1.
    analysis Static
    analyze 1
    print node 2
}

wipe
model Basic -ndm 2 
node 1 0 0
node 2 0 10
node 3 20 0
node 4 20 10.

node 21 0 10
node 41 20 10.
    
set A 10
set E 100
set I 1000.
geomTransf Linear 1
uniaxialMaterial Steel01 1 50 1e10 0.00000000000001
element elasticBeamColumn 1 1 21 $A $E $I 1 
element elasticBeamColumn 2 3 41 $A $E $I 1 
element elasticBeamColumn 3 2 4 $A $E $I 1 
element zeroLength 4 21 2 -mat 1 -dir 6
element zeroLength 5 41 4 -mat 1 -dir 6

fix 1 1 1 1
fix 3 1 1 1
equalDOF 2 21 1 2
equalDOF 4 41 1 2

timeSeries Linear 1
pattern Plain 1 1 {
    load 2 10 0 0
}
test NormDispIncr 1e-10 100 0
algorithm Newton -initial
system FullGeneral
integrator LoadControl .1
analysis Static
analyze 30
print node 2
print ele 1

puts $outFile [nodeDisp 2]

wipe
model Basic -ndm 2 
node 1 0 0
node 2 0 10
node 3 20 0
node 4 20 10.

set A 10
set E 100
set I 1000.
geomTransf Linear 1
uniaxialMaterial Elastic 1 1.0e10
uniaxialMaterial Steel01 2 50 1e10 0.00000000000001
#uniaxialMaterial Elastic 2 1.0e10
element componentElement2d 1 1 2 $A $E $I 1 1 2
element componentElement2d 2 3 4 $A $E $I 1 1 2
element elasticBeamColumn 3 2 4 $A $E $I 1 

fix 1 1 1 1
fix 3 1 1 1

timeSeries Linear 1
pattern Plain 1 1 {
    load 2 10 0 0
}
test NormDispIncr 1e-12 5 0
algorithm Newton 
system FullGeneral
integrator LoadControl 0.1
analysis Static
analyze 30
print node 2
print ele 1

puts $outFile [nodeDisp 2]
close $outFile
