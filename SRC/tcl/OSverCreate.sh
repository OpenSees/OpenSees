#!/bin/bash
# used to create version.txt file for OpenSees on UNIX
# Boris Jeremic@ucdavis.edu Oct2002

VERSION="`whoami` `uname -n` `uname -r` `uname -s` `date`"
echo 'char version[132] = "compiled and linked by: '$VERSION'";' >> version.txt
