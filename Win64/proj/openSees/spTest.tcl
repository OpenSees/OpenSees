wipeAnalysis

model BasicBuilder -ndm 2 -ndf 3

set L 100
set E 29000
set A 20
set I 1400

node 1 0 0
node 2 $L 0

fix 1 1 1 1
fix 2 1 1 1

geomTransf Linear 1

element elasticBeamColumn 1 1 2 $A $E $I 1

timeSeries Constant 1

# pattern Plain 1 1  {
#     sp 2 2 1.0
# }

pattern Plain 1 1  {
    sp 2 0 1 0 -constant 0
}
constraints Transformation
system UmfPack
analysis Static
analyze 1
reactions
nodeReaction 1
