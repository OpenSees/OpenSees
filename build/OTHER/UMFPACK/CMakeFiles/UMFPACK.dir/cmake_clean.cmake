file(REMOVE_RECURSE
  "../../lib/libUMFPACK.a"
  "../../lib/libUMFPACK.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C)
  include(CMakeFiles/UMFPACK.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
