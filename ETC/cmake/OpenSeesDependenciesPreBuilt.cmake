#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#==============================================================================
#                      External Libraries
#
#==============================================================================
# Synopsis
# - opensees_load(<PACKAGE> OPTIONS [<PATHS>])
#
# Options:
# - FIND:   Use CMake to find library, fail if not found
# - SEARCH: Try finding library with CMake, build OpenSees
#           Version if not found.
#
# Keyword arguments
#   Provide specific paths for library.
#
# - BUNDLED <path/to/OTHER/LIB/>
#   Provide path to OpenSees bundled library
#
# - LIBRARY <path/to/lib.a> INCLUDE <path/to/include/>
#
#----------------------------------------------------------------
opensees_load(TCL                                            FIND
)

set(TCL_LIBRARIES ${TCL_LIBRARY})

opensees_load(BLAS                                           FIND 
	#LIBRARY /home/user/lib/libBlas.a
)

opensees_load(LAPACK                                         FIND
  #LIBRARY /home/user/lib/libLapack.a
)

opensees_load(SUPERLU                                       #FIND 
  BUNDLED "${OPS_BUNDLED_DIR}/SuperLU_5.1.1/"
)

opensees_load(ARPACK                                        #FIND
  BUNDLED "${OPS_BUNDLED_DIR}/ARPACK/"
)

opensees_load(UMFPACK                                        #FIND
  BUNDLED "${OPS_BUNDLED_DIR}/UMFPACK/" 
)
opensees_load(CSPARSE                                        #FIND
  BUNDLED "${OPS_BUNDLED_DIR}/CSPARSE/"
)

opensees_load(AMD                                           #FIND
	BUNDLED "${OPS_BUNDLED_DIR}/AMD/" 
)

opensees_load(METIS                                          FIND)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)

# Integrated exteral libraries
opensees_load(FEDEAS_Uniaxial
  LIBRARY FALSE
)

