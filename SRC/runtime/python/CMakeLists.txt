#=====================================================================
# 
#   OpenSees -- Open System For Earthquake Engineering Simulation
#           Pacific Earthquake Engineering Research Center
#
#=====================================================================
#find_package(PythonExtensions REQUIRED)

# find_package(Python COMPONENTS Interpreter Development)
# find_package(pybind11 CONFIG)

# add_library(OpenSeesPyRT SHARED)
pybind11_add_module(OpenSeesPyRT)

#message(" :: ${pybind11_INCLUDE_DIRS}\n     ${pybind11_LIBRARIES}")

target_compile_definitions(OpenSeesPyRT PRIVATE USE_TCL_STUBS)

target_include_directories(OpenSeesPyRT 
  PRIVATE #${Python_INCLUDE_DIRS}
          ${pybind11_INCLUDE_DIRS}
          ${OPS_SRC_DIR}/domain/pattern/series
)

# set_target_properties(OpenSeesPyRT PROPERTIES CXX_VISIBILITY_PRESET hidden)
target_link_libraries(OpenSeesPyRT PRIVATE ${pybind11_LIBRARIES} ${TCL_STUB_LIBRARY})
target_link_libraries(OpenSeesPyRT PRIVATE OPS_Runtime OpenSeesRT)

target_sources(OpenSeesPyRT PRIVATE
  "OpenSeesPyRT.cpp"
)


