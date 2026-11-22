# VTK recorder heap-overflow probe (see tests/vtk_recorder_valgrind.sh).
# A minimal stable truss model; every committed step makes the VTK
# recorder rebuild the .vtu filename with sprintf, which is where the
# two buffer overruns happened.
#
# The base name is long on purpose: it exercises both the recorder's
# stored-name copy and the per-step filename construction.

wipe
model basic -ndm 2 -ndf 2
node 1 0 0
node 2 1 0
fix 1 1 1
fix 2 0 1
uniaxialMaterial Elastic 1 1000.0
element truss 1 1 2 1.0 1
timeSeries Linear 1
pattern Plain 1 Linear {
    load 2 1.0 0.0
}
system BandGeneral
numberer Plain
constraints Plain
integrator LoadControl 0.2
algorithm Linear
analysis Static

recorder vtk vtk_overflow_probe disp

analyze 5
exit 0
