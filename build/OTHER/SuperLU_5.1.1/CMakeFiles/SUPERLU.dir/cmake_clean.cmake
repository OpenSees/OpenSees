file(REMOVE_RECURSE
  "../../lib/libSUPERLU.a"
  "../../lib/libSUPERLU.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/SUPERLU.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
