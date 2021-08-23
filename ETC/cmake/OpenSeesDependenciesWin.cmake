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
# Synopsis
# - opensees_load(<PACKAGE> [BUILD|FIND|SEARCH|PATHS] [<PATHS>])
#
# Options:
# - BUILD:  Build OpenSees provided library
# - FIND:   Use CMake to find library, fail if not found
# - SEARCH: Try finding library with CMake, build OpenSees
#           Version if not found.
# - BUNDLED:  Provide specific paths for library.
#
#==============================================================================
set(CONDA_DIR "C:/Users/claud/miniconda3")
set(CONDA_ENV "C:/Users/claud/miniconda3/envs/sim")
set(BUNDLE_LIBS "${PROJECT_SOURCE_DIR}/Win64/lib/debug/")

opensees_load(TCL 
    LIBRARY ${CONDA_DIR}/Library/lib/tcl86t.lib
    #LIBRARY "${BUNDLE_LIBS}/tcl.lib"
    INCLUDE ${CONDA_DIR}/Library/include 
)

set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIRS})
set(TCL_LIBRARY ${TCL_LIBRARIES})

message("TCL: ${TCL_INCLUDE_PATH}")

opensees_load(BLAS
    LIBRARY "${BUNDLE_LIBS}/blas.lib"
)

opensees_load(CBLAS
    LIBRARY "${BUNDLE_LIBS}/cblas.lib"
)

opensees_load(LAPACK
    LIBRARY "${BUNDLE_LIBS}/lapack.lib"
)

#set(ENV{SUPERLU_DIR})
opensees_load(SUPERLU
    BUNDLED ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
)

opensees_load(ARPACK                                       SEARCH
    BUNDLED ${OPS_BUNDLED_DIR}/ARPACK/
)

opensees_load(METIS                                        SEARCH)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          #FIND
    LIBRARY ${CONDA_ENV}/Library/lib/libmysql.lib
    INCLUDE ${CONDA_ENV}/Library/include/mysql
)

set(MYSQL_INCLUDE_DIR "${CONDA_ENV}/Library/include/mysql/")


