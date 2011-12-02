#!/bin/csh
# used to create version.txt file for OpenSees on UNIX
# Boris Jeremic@ucdavis.edu Oct2002
set VERSION = (`whoami` `uname -n` `uname -r` `uname -s` `date` `setenv  | grep CC+`) 
echo $VERSION
echo 'char version[132] = "compiled and linked by: '$VERSION'";' >! version.txt
