# Reproducible METIS partition benchmark for OpenSeesMP.
# Run with: mpiexec -n 4 OpenSeesMP metisPartitionBenchmark.tcl

puts "METIS_BENCHMARK_START rank=[getPID] ranks=[getNP]"
flush stdout

model BasicBuilder -ndm 2 -ndf 2
uniaxialMaterial Elastic 1 2.0e11

set nx 40
set ny 20
set objective volume
set seed 42
if {$argc > 0} {
    set objective [lindex $argv 0]
}
if {$argc > 1} {
    set seed [lindex $argv 1]
}
set nodeTag 1
for {set j 0} {$j <= $ny} {incr j} {
    for {set i 0} {$i <= $nx} {incr i} {
        node $nodeTag [expr {double($i)}] [expr {double($j)}]
        incr nodeTag
    }
}

set eleTag 1
for {set j 0} {$j <= $ny} {incr j} {
    for {set i 0} {$i < $nx} {incr i} {
        set n1 [expr {$j * ($nx + 1) + $i + 1}]
        element truss $eleTag $n1 [expr {$n1 + 1}] 0.01 1
        incr eleTag
    }
}
for {set j 0} {$j < $ny} {incr j} {
    for {set i 0} {$i <= $nx} {incr i} {
        set n1 [expr {$j * ($nx + 1) + $i + 1}]
        element truss $eleTag $n1 [expr {$n1 + $nx + 1}] 0.01 1
        incr eleTag
    }
}

set totalBefore [getNumElements]
partition -objective $objective -ncuts 4 -niter 20 -ufactor 30 -seed $seed -info
puts "METIS_BENCHMARK rank=[getPID] ranks=[getNP] before=$totalBefore after=[getNumElements]"
