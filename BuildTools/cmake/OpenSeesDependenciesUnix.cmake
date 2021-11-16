#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#==============================================================================
#                            External Libraries
#
#==============================================================================
# Synopsis
# - opensees_load(<PACKAGE> [FLAGS] [<PATHS>])
#
# Flags:
# - FIND:   Use CMake to find library on system, fail if not found
#
# Keyword arguments:
#   Provide specific paths for library.
#
# - BUNDLED <path/to/OTHER/LIB/>
#   Provide path to OpenSees bundled library directory containing a 
#   CMakeLists.txt file.
#
# - LIBRARY <path/to/lib.a> INCLUDE <path/to/include/>
#
# - CONAN <conan-package/version>
#   Point to a conan package.
#
#----------------------------------------------------------------
#opensees_load(TCL CONAN tcl/8.6.10)
opensees_load(TCL                                            FIND)
set(TCL_LIBRARIES ${TCL_LIBRARY})

opensees_load(BLAS                                           FIND)

opensees_load(LAPACK                                         FIND)

opensees_load(SUPERLU
  BUNDLED "${OPS_BUNDLED_DIR}/SuperLU_5.1.1/"
  #LIBRARY "${OPS_BUNDLED_DIR}/bin/SuperLU_5.1.1/libSUPERLU.a"
)

opensees_load(ARPACK
  BUNDLED "${OPS_BUNDLED_DIR}/ARPACK/"
  #LIBRARY "${OPS_BUNDLED_DIR}/bin/ARPACK/libARPACK.a"
)

opensees_load(UMFPACK
  BUNDLED "${OPS_BUNDLED_DIR}/UMFPACK/" 
  #LIBRARY "${OPS_BUNDLED_DIR}/bin/UMFPACK/libUMFPACK.a" 
)

opensees_load(CSPARSE
  BUNDLED "${OPS_BUNDLED_DIR}/CSPARSE/"
  #LIBRARY "${OPS_BUNDLED_DIR}/bin/CSPARSE/libCSPARSE.a"
)

opensees_load(AMD
	BUNDLED "${OPS_BUNDLED_DIR}/AMD/"
  #LIBRARY "${OPS_BUNDLED_DIR}/bin/AMD/libAMD.a"
)

opensees_load(tet BUNDLED "${OPS_BUNDLED_DIR}/Tetgen/")
opensees_load(triangle BUNDLED "${OPS_BUNDLED_DIR}/Triangle/")

opensees_load(METIS                                          FIND)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)

find_package(Python COMPONENTS Development)

# Integrated exteral libraries
opensees_load(FEDEAS_Uniaxial
  LIBRARY FALSE
)


