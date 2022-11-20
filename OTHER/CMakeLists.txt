#==============================================================================
# 
#        OpenSees -- Open System For Earthquake Engineering Simulation
#                Pacific Earthquake Engineering Research Center
#
#     (c) Copyright 1999-2021 The Regents of the University of California
#                             All Rights Reserved
# (Copyright and Disclaimer @ http://www.berkeley.edu/OpenSees/copyright.html)
#
#------------------------------------------------------------------------------
# Pre-compile OpenSees numeric libraries
#==============================================================================

cmake_minimum_required(VERSION 3.16)
project(OpenSeesNumPack)

set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

add_subdirectory(AMD/)
add_subdirectory(UMFPACK/)
add_subdirectory(ITPACK/)
add_subdirectory(ARPACK/)
add_subdirectory(SuperLU_5.1.1/)
add_subdirectory(CSPARSE/)

if(MPI_FOUND)
message("adding MUMPS and SUPERLU_DIST")
add_subdirectory(MUMPS_5.4.1/)
add_subdirectory(SUPERLU_DIST_4.3/SRC)
else()
message("NOT adding MUMPS and SUPERLU_DIST")
endif()

