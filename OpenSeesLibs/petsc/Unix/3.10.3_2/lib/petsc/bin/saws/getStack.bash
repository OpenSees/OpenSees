#!/bin/bash
#
#  Returns the SAWS published PETSc stack one function per line
#
which jshon > /dev/null
if [ $? != 0 ] ; then
  echo "You must install jshon before using this app"
  exit 1
fi

${PETSC_DIR}/bin/saws/getSAWs.bash PETSc/Stack | jshon -e directories -e Stack -e variables -e functions -e data 
#


