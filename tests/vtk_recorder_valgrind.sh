#!/bin/sh
# Heap-overflow probe for the VTK recorder (SRC/recorder/VTK_Recorder.cpp),
# run under valgrind. Linux only.
#
# Two bugs wrote past the end of heap blocks on every recorded step:
#
#   stored name:   new char[strlen(inputName+1)] -- strlen(inputName+1)
#                 measures the name MINUS its first character, and the +1
#                 for the NUL terminator is missing, so the strcpy that
#                 fills the buffer runs 2 bytes past the block
#   step names:   new char[2*strlen(name)+26] -- the '%s/%s%d%020d.vtu'
#                 format needs 2*strlen(name)+37 bytes ('/' + up to 11
#                 digits for the part number + 20 zero-padded timestep
#                 digits + '.vtu' + NUL), so every sprintf runs at least
#                 1 byte past, more with MPI partitions
#
# Most runs survive because the allocator pads small blocks. When they
# do not, the run aborts inside malloc/free long after the recorder
# line, or the .vtu filenames come out corrupted and ParaView cannot
# load the series. valgrind makes it deterministic. Measured with the
# 5-step probe in vtk_recorder_probe.tcl:
#
#   unfixed VTK_Recorder.cpp:  108 errors from 13 contexts,
#                              including 4 invalid writes
#   fixed:                      0 errors
#
# Usage: sh tests/vtk_recorder_valgrind.sh /path/to/OpenSees
# A relative path is accepted; it is resolved against the directory you
# run the script from (the script itself cd's into a scratch directory).
OPENSEES="${1:?usage: vtk_recorder_valgrind.sh /path/to/OpenSees}"
case "$OPENSEES" in
    /*) ;;
    *) OPENSEES="$PWD/$OPENSEES" ;;
esac
command -v valgrind >/dev/null 2>&1 || {
    echo "ERROR: valgrind not found in PATH (Linux tool, install your distro's valgrind)"
    exit 2
}
[ -x "$OPENSEES" ] || {
    echo "ERROR: not an executable: $1"
    exit 2
}
DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
SCRATCH=$(mktemp -d)
trap 'rm -rf "$SCRATCH"' EXIT
cd "$SCRATCH"

valgrind --error-exitcode=99 "$OPENSEES" "$DIR/vtk_recorder_probe.tcl" > valgrind.log 2>&1
rc=$?

grep -E "Invalid write|Invalid read|ERROR SUMMARY" valgrind.log
if [ $rc -eq 0 ]; then
    echo "PASS: VTK recorder runs clean under valgrind"
elif [ $rc -eq 99 ]; then
    echo "FAIL: valgrind found memory errors -- see the lines above"
else
    echo "ERROR: the run itself did not complete (exit $rc) -- this is not a memory verdict"
    tail -5 valgrind.log
fi
exit $rc
