file(REMOVE_RECURSE
  "lib/libG3.a"
  "lib/libG3.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang C CXX)
  include(CMakeFiles/G3.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
