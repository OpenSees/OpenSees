file(REMOVE_RECURSE
  "../../lib/libMUMPS.a"
  "../../lib/libMUMPS.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/MUMPS.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
