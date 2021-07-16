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
# - opensees_load(<PACKAGE> [BUILD|FIND|SEARCH|PATHS] [<PATHS>])
#
# Options:
# - BUILD:  Build OpenSees provided library
# - FIND:   Use CMake to find library, fail if not found
# - SEARCH: Try finding library with CMake, build OpenSees
#           Version if not found.
#
# Keyword arguments
# - PATHS:  <PATH>..
#   Provide specific paths for library.
#
# - BUNDLED <PATH>
#
#----------------------------------------------------------------
opensees_load(TCL                                          FIND
	#LIBRARY /home/claudio/miniconda3/lib/libtcl8.6.so
	#INCLUDE /home/claudio/miniconda3/include 
	#AS TCL_LIBRARY TCL_INCLUDE_PATH
)
set(TCL_LIBRARIES ${TCL_LIBRARY})

opensees_load(BLAS                                           FIND 
	#LIBRARY /home/claudio/lib/libBlas.a
)

opensees_load(LAPACK                                         FIND
        #LIBRARY /home/claudio/lib/libLapack.a
)

opensees_load(SUPERLU                                       #FIND 
    BUNDLED ${OPS_BUNDLED_DIR}/SuperLU_5.1.1/
)

opensees_load(ARPACK                                        #FIND
    BUNDLED ${OPS_BUNDLED_DIR}/ARPACK/ 
)

opensees_load(AMD                                            FIND
	BUNDLED ${OPS_BUNDLED_DIR}/AMD/ 
)

opensees_load(METIS                                          FIND)

opensees_load(HDF5                                           FIND)

opensees_load(MySQL                                          FIND)


#==============================================================================
#                           Select Element Libraries
#
# Each element in this list ows and associated macro definition
#==============================================================================
set(OPS_Element_List

    OPS_Element_truss
    OPS_Element_beam3d
    #OPS_Element_beam2d
    OPS_Element_dispBeamColumnInt
    OPS_Element_forceBeamColumn
    OPS_Element_mixedBeamColumn

    #OPS_Element_beamWithHinges
    OPS_Element_LHMYS
    OPS_Element_PML
    OPS_Element_RockingBC
    OPS_Element_UP_ucsd
    OPS_Element_absorbentBoundaries
    OPS_Element_adapter
    OPS_Element_catenaryCable
    OPS_Element_componentElement

    OPS_Element_elastomericBearing
    OPS_Element_frictionBearing

    OPS_Element_generic
    OPS_Element_gradientInelasticBeamColumn
    OPS_Element_joint
    OPS_Element_mvlem
    OPS_Element_pyMacro
    OPS_Element_shell
    OPS_Element_surfaceLoad
    OPS_Element_updatedLagrangianBeamColumn
    #OPS_Element_feap
    #OPS_Element_PFEMElement
)


