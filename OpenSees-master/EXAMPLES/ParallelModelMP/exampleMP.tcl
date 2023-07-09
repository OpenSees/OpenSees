wipe

# some parameters
set numP [getNP]
set pid  [getPID]

source common.tcl

set dx [expr $d/($nx*1.0)]
set dy [expr $b/($ny*1.0)]
set dz [expr $L/($numP * $nz)]

# create the modelbuilder
model Basic -ndm 3 -ndf 3

# add a material
set E 30000.
set v 0.2
nDMaterial ElasticIsotropic   1   $E $v

# add the nodes
set counter 1
for {set i 0} {$i<=[expr $numP*$nz]} {incr i 1} {
    for {set j 0} {$j<=$ny} {incr j 1} {
	for {set k 0} {$k<=$nx} {incr k 1} {
	    if {($i >= [expr $pid*$nz]) && ($i <= [expr ($pid+1)*$nz])} {
		node $counter [expr $k*$dx] [expr $j*$dy] [expr $i*$dz]
	    }
	    incr counter
	}
    }
}

set numNode [expr ($nx+1)*($ny+1)*($numP*$nz+1)]

# add the elements
set counter 1
for {set i 0} {$i<$numP*$nz} {incr i 1} {
         for {set j 0} {$j<$ny} {incr j 1} {
          set iNode [expr $i*($nx+1)*($ny+1) + 1 + $j*($nx+1)]
          set jNode [expr $iNode+1]
          set lNode [expr $iNode+($nx+1)]
          set kNode [expr $lNode+1]
           
          set mNode [expr ($i+1)*($nx+1)*($ny+1) + 1 + $j*($nx+1)]
          set nNode [expr $mNode+1]
          set pNode [expr $mNode+($nx+1)]
          set oNode [expr $pNode+1]
  
	for {set k 0} {$k<$nx} {incr k 1} {
	           if {$i >= [expr $pid*$nz] && $i < [expr ($pid+1)*$nz]} {
		      element stdBrick $counter $iNode $jNode $kNode $lNode $mNode $nNode $oNode $pNode 1
	      }
	    incr counter
	    incr iNode 1
	    incr jNode 1
	    incr kNode 1
	    incr lNode 1
	    incr mNode 1
	    incr nNode 1
	    incr oNode 1
	    incr pNode 1
	}
    }
}

fixZ 0.0 1 1 1 

if {$pid == [expr $numP-1]} {
    timeSeries Linear 1
    pattern Plain 1 1 {
	load $numNode 10. 10. 0.
    }
}

integrator LoadControl 1.0
algorithm Linear
numberer Plain
constraints Plain
system Mumps
analysis Static

analyze 1
puts "HI"

if {$pid == [expr $numP-1]} {	
    set disp [nodeDisp $numNode 1]
    puts "$disp"
}
wipe
 
