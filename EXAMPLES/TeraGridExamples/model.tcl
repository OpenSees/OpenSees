
model basic -ndm 2 -ndf 3

set ok 0
if {$ok == 0} {

set width    360
set height   144
set P 180.0
set G 386.4
set m [expr $P/$G];       # expr command to evaluate an expression

node  1       0.0     0.0 
node  2    $width     0.0 
node  3       0.0 $height -mass $m $m 0.0
node  4    $width $height -mass $m $m 0.0

fix   1     1    1    1
fix   2     1    1    1

uniaxialMaterial Concrete01  1  -6.0  -0.004   -5.0     -0.014
uniaxialMaterial Concrete01  2  -5.0   -0.002   0.0     -0.006
set fy 60.0;      # Yield stress
set E 30000.0;    # Young's modulus
uniaxialMaterial Steel01  3  $fy $E 0.01

set colWidth 15
set colDepth 24 

set cover  1.5
set As    0.60;     # area of no. 7 bars

# some variables derived from the parameters
set y1 [expr $colDepth/2.0]
set z1 [expr $colWidth/2.0]

section Fiber 1 {

    # Create the concrete core fibers
    patch rect 1 10 1 [expr $cover-$y1] [expr $cover-$z1] [expr $y1-$cover] [expr $z1-$cover]

    # Create the concrete cover fibers (top, bottom, left, right)
    patch rect 2 10 1  [expr -$y1] [expr $z1-$cover] $y1 $z1
    patch rect 2 10 1  [expr -$y1] [expr -$z1] $y1 [expr $cover-$z1]
    patch rect 2  2 1  [expr -$y1] [expr $cover-$z1] [expr $cover-$y1] [expr $z1-$cover]
    patch rect 2  2 1  [expr $y1-$cover] [expr $cover-$z1] $y1 [expr $z1-$cover]

    # Create the reinforcing fibers (left, middle, right)
    layer straight 3 3 $As [expr $y1-$cover] [expr $z1-$cover] [expr $y1-$cover] [expr $cover-$z1]
    layer straight 3 2 $As 0.0 [expr $z1-$cover] 0.0 [expr $cover-$z1]
    layer straight 3 3 $As [expr $cover-$y1] [expr $z1-$cover] [expr $cover-$y1] [expr $cover-$z1]

}    


geomTransf Linear 1  

set np 5
element nonlinearBeamColumn  1   1   3   $np    1       1 
element nonlinearBeamColumn  2   2   4   $np    1       1 

geomTransf Linear 2  
element elasticBeamColumn   3   3   4    360    4030  8640    2

pattern Plain 1 "Linear" {
	load  3   0.0  [expr -$P] 0.0
	load  4   0.0  [expr -$P] 0.0
}


}
