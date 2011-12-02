# Define a procedure which will generate nodes and elements for
# a plane frame having absolute column line locations in the list
#'columnLine',  absolute girder line locations in the list 'girderLine',
# section IDs for the columns and girders, 'columnID' and 'girderID', and
# 'nIntPt' integration points along every member.
#
# Notes: automatically fixes the base nodes
#        only geneartes nonlinearBeamColumn elements
#        allows columns and girders to be spaced arbitrarily
#        does not add nodal masses or loads, but can be extended to do so
#        starts node numbering at 1
#
# Formal arguments
#    columnLine - a list of column line locations
#       The actual argument would be defined as so,
#           set columns {0 120 240 360}
#    girderLine - a list of grider line locations
#       The actual argument would be defined as so,
#           set girders {180 300 420}
#    columnID - an integer representing the section ID (tag) for
#               the columns in the frame
#    girderID - an integer representing the section ID (tag) for
#               the girders in the frame
#    nIntPt - an integer representing the number of integration points
#             along each member in the frame

proc genPlaneFrame {columnLine girderLine columnID girderID nIntPt} {

   # Node number counter 
   set n 1

   # Geometric transformation for all elements
   geomTransf Linear 1

   # For each column line
   foreach xLoc $columnLine {
      node $n $xLoc 0
      # Fix the base node
      fix $n 1 1 1
      incr n 1
      # For each girder line
      foreach yLoc $girderLine {
         node $n $xLoc $yLoc
         incr n 1
      }
   }

   # Useful variables
   set numCol [llength $columnLine]
   set numGir [llength $girderLine]

   # Element number counter
   set e 1

   # For each column line
   for {set i 1} {$i <= $numCol} {incr i 1} {
      # Node number at the base of this column line
      set bottom [expr ($i-1)*$numGir + $i]
      # Node number at the top of this column line
      set top [expr $i*$numGir + $i]
      # Travel up this column line creating elements
      for {set j $bottom} {$j <= [expr $top-1]} {incr j 1} {
         element nonlinearBeamColumn $e $j [expr $j+1] $nIntPt $columnID 1
         incr e 1
      }
   }

   # Difference in node numbers I and J for any girder in the frame
   set delta [expr $numGir+1]

   # For each girder line
   for {set j 1} {$j <= $numGir} {incr j 1} {
      # Node number at the left end of this girder line
      set left [expr $j+1]
      # Node number at the right end of this girder line
      set right [expr ($numCol-1)*$numGir + $numCol + $j]
      # Travel across this girder line creating elements
      for {set k $left} {$k < $right} {incr k $delta} {
         element nonlinearBeamColumn $e $k [expr $k+$delta] $nIntPt $columnID 1
         incr e 1
      }
   }
}
