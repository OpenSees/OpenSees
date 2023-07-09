# TRUSS EXAMPLE
puts "Truss Example to Compare Truss, CorotTruss and PlanarTruss"


puts "\n      Element Type      Disp X         Disp Y"
puts "    -----------------------------------------"
set formatString {%15s%15.5f%15.5f}
foreach eleType {truss corotTruss planarTruss} {

    wipe
    model BasicBuilder -ndm 2 -ndf 2

    # Create nodes
    # ------------
    node 1   0.0  0.0
    node 2 144.0  0.0
    node 3 168.0  0.0
    node 4  72.0 96.0

    # Set the boundary conditions - command: fix nodeID xResrnt? yRestrnt?
    fix 1 1 1 
    fix 2 1 1
    fix 3 1 1

    # Define materials for truss elements
    # -----------------------------------
    
    # Create Elastic material prototype - command: uniaxialMaterial Elastic matID E
    uniaxialMaterial Elastic 1 3000
    section Elastic 1 3000 5.0 10000
    
    # Define elements
    # ---------------
    
    # Create truss elements - command: element truss trussID node1 node2 A matID
    if {$eleType == "planarTruss"} {
	element planarTruss 1 1 4 1 10.0
	element planarTruss 2 2 4 1 5.0 
	element planarTruss 3 3 4 1 5.0 
    } else {
	element $eleType 1 1 4 10.0 1
	element $eleType 2 2 4 5.0 1
	element $eleType 3 3 4 1
    }
    
    # Create a Plain load pattern with a linear TimeSeries
    pattern Plain 1 "Linear" {
	load 4 100 -50
    }
    
    system BandSPD
    numberer RCM
    constraints Plain
    integrator LoadControl 1.0
    test NormDispIncr 1.0e-10 10
    algorithm Newton
    analysis Static 
    
    recorder Node -file example.out -time -node 4 -dof 1 2 disp

    analyze 1

    set nodeDisp4 [nodeDisp 4]
    puts [format $formatString $eleType [lindex $nodeDisp4 0] [lindex $nodeDisp4 1]]
}
puts ""
