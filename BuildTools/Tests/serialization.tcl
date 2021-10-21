source TestUtilities/Variables.tcl
source TestUtilities/MelolandModel.tcl

file mkdir $OUTDIR/serialization/
print -JSON -file $OUTDIR/serialization/Meloland.json

