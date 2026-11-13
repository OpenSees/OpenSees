wipe
if {[llength [info commands testUniaxialMaterial]]} {
    model BasicBuilder -ndm 1 -ndf 1
} else {
    model testUniaxial
}
uniaxialMaterial CFSTsteel 1 \
    206000.0 305.0 0.02 -342.07 26.58 -74.87 3.16 \
    4.07 6.38 3.12 2.36 5.82 1.0 1.0 1.0 1.0 1.0 \
    116.39 11.81 1 4467.68 38.39

if {[llength [info commands testUniaxialMaterial]]} {
    testUniaxialMaterial 1
    proc applyStrain {value} { setStrain $value }
    proc readStress {} { return [getStress] }
    proc readTangent {} { return [getTangent] }
} else {
    uniaxialTest 1
    proc applyStrain {value} { strainUniaxialTest $value }
    proc readStress {} { return [stressUniaxialTest] }
    proc readTangent {} { return [tangUniaxialTest] }
}
set output [open CFSTsteel_cyclic_tcl.csv w]
puts $output "strain,stress_MPa,tangent_MPa"
set last 0.0
foreach peak {0.0002 -0.0002 0.003 -0.004 0.01 -0.02 0.015 -0.03 0.02 -0.04 0.01 0.0} {
    for {set j 1} {$j <= 100} {incr j} {
        set strain [expr {$last + ($peak-$last)*$j/100.0}]
        applyStrain $strain
        puts $output [format "%.17g,%.17g,%.17g" $strain [readStress] [readTangent]]
    }
    set last $peak
}
close $output
puts "Wrote CFSTsteel_cyclic_tcl.csv (1200 strain steps)"
wipe
