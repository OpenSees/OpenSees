#!/bin/bash
#
#  getSAWs.bash [dir1[/dir2[/variablename]]]
#
#
if [ "${SAWS_HOST}foo" == "foo" ]; then export SAWS_HOST=localhost; fi
if [ "${SAWS_PORT}foo" == "foo" ]; then export SAWS_PORT=8080; fi
if [ $# == 1 ]; then
  export mem=$1;
else 
  export mem="*";
fi

curl --silent --show-error "${SAWS_HOST}:${SAWS_PORT}/SAWs/${mem}"
