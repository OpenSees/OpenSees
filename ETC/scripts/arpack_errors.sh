#!/bin/bash

# Claudio Perez
#
# $1: ARPACK source directory

for fname in $1/*.f; do
  justfile=${fname##*/};
  printf "const char *ARPACK_${justfile%.f}(int code)\n{\n  switch (code) {\n";
  sed -nE 's/^c( *)= *(-?[0-9]*): *(.*)/  case \2: return "\3"/gp' $fname;
  #sed -nE 's/^c( *)= *(-?[0-9]*): *(.*)/  case \2: return "\3"/gp; tx; :x; s/^c *([A-z].*)/\n    "\1"/' $fname;
  printf "  default: return \"Unknown error code in ARPACK::${justfile%.f}\"\n  }\n}\n\n"
done

