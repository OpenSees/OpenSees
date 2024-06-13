# ------------------------------
# Start of model generation
# ------------------------------

# Remove existing model
wipe

# Create ModelBuilder 
model BasicBuilder -ndm 1 -ndf 1

# Create nodes
# ------------
    
# Create nodes & add to Domain 
node 1 0.0 
node 2 0.0 

# Set the boundary conditions - command: fix nodeID xResrnt? 
fix 1 1  

# Define materials for truss elements
# -----------------------------------

# Create Elastic material prototype 

set Fy 200.0
set E 1000000.0
set fTravel 0.004
set fTravelInitial 0
set RatType 2

#set b -0.1
#uniaxialMaterial Steel01 1 $Fy $E $b
#uniaxialMaterial ElasticPPGap 1 $E $Fy 0.002 0 damage
#uniaxialMaterial GNG 1 $E $Fy $toothSize 0
uniaxialMaterial Ratchet 1 $E $fTravel $fTravelInitial $RatType
uniaxialMaterial Elastic 2 1

# Define elements
# ---------------

# Create element
element zeroLength 1 1 2 -mat 1 -dir 1
element zeroLength 2 1 2 -mat 2 -dir 1
    
# Define loads
# ------------

set P 300.0
#create a Linear TimeSeries (load factor varies linearly with time): command timeSeries Linear $tag
timeSeries Linear 1

# Create a Plain load pattern with a linear TimeSeries: command pattern Plain $tag $timeSeriesTag { $loads }
pattern Plain 1 1 {
    
    # Create the nodal load - command: load nodeID xForce
    load 2 $P
}
    
# ------------------------------
# Start of analysis generation
# ------------------------------


constraints Transformation
numberer RCM
test NormDispIncr 1.0e-6 6 0
algorithm ModifiedNewton
system BandGeneral
integrator DisplacementControl 2 1 0.001
analysis Static

# create a Recorder object for the nodal displacements at node 2
recorder Node -file RatchetTestOutput.out -time -node 2 -dof 1 disp

#foreach numIter {10 20 10 20} dU {0.001 -0.001 0.001 -0.001} {
#integrator DisplacementControl 2 1 $dU
#analyze $numIter
#set factor [getTime]
#puts "[expr $factor*$P] [lindex [nodeDisp 2] 0]"
#}

#foreach numIter {10 20 10 10 5 10} dU {0.001 -0.001 0.001 0.001 -0.001 0.001} {
#integrator DisplacementControl 2 1 $dU
#analyze $numIter
#set factor [getTime]
#puts "[expr $factor*$P] [lindex [nodeDisp 2] 0]"
#}

foreach numIter {10 20 10 20 10 20 10} dU {0.001 -0.001 0.001 -0.001 0.001 -0.001 0.001} {
integrator DisplacementControl 2 1 $dU
analyze $numIter
set factor [getTime]
puts "[expr $factor*$P] [lindex [nodeDisp 2] 0]"
}

print node 2
print ele
