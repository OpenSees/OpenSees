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
#                             External Libraries
#
# - BLAS_LIBRARIES
# - BLAS_INCLUDE_DIRS
#
# - LAPACK_LIBRARIES
# - LAPACK_INCLUDE_DIRS
#
# - ARPACK_LIBRARIES
#
# - SUPERLU_LIBRARIES
# - SUPERLU_INCLUDE_DIRS
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
# - PATHS:  Provide specific paths for library.
#
#==============================================================================
set(CONDA_DIR "C:/Users/claud/miniconda3")
set(CONDA_ENV "C:/Users/claud/miniconda3/envs/sim")

opensees_load(TCL                                          #FIND
	LIBRARY ${CONDA_DIR}/Library/lib/tcl86t.lib
	INCLUDE ${CONDA_DIR}/Library/include 
)

set(TCL_INCLUDE_PATH ${TCL_INCLUDE_DIRS})
set(TCL_LIBRARY ${TCL_LIBRARIES})

message("TCL: ${TCL_INCLUDE_PATH}")

opensees_load(BLAS                                         #FIND
	LIBRARY ${CONDA_ENV}/Library/lib/blas.lib
	INCLUDE ${CONDA_ENV}/Library/include/
)

opensees_load(CBLAS                                         #FIND
	LIBRARY ${CONDA_ENV}/Library/lib/cblas.lib
	INCLUDE ${CONDA_ENV}/Library/include/
)

opensees_load(LAPACK                                       #FIND
	LIBRARY ${CONDA_ENV}/Library/lib/lapack.lib
	INCLUDE ${CONDA_ENV}/Library/include/
)

set(ENV{SUPERLU_DIR})
opensees_load(SUPERLU                                       #SEARCH
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



#==============================================================================
#                           Select Element Libraries
#
# Each element in this list ows and associated macro definition
#==============================================================================
set(OPS_Element_List
    #OPS_Element_PFEMElement
    #OPS_Element_beamWithHinges
    #OPS_Element_feap
    OPS_Element_LHMYS
    OPS_Element_PML
    OPS_Element_RockingBC
    OPS_Element_UP_ucsd
    OPS_Element_absorbentBoundaries
    OPS_Element_adapter
    OPS_Element_beam3d
    #OPS_Element_beam2d
    OPS_Element_catenaryCable
    OPS_Element_componentElement
    OPS_Element_dispBeamColumnInt
    OPS_Element_forceBeamColumn
    OPS_Element_elastomericBearing
    OPS_Element_frictionBearing
    OPS_Element_generic
    OPS_Element_gradientInelasticBeamColumn
    OPS_Element_joint
    OPS_Element_mixedBeamColumn
    OPS_Element_mvlem
    OPS_Element_pyMacro
    OPS_Element_shell
    OPS_Element_surfaceLoad
    OPS_Element_truss
    OPS_Element_updatedLagrangianBeamColumn
)


