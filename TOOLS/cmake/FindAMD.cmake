# - Find the AMD includes and library
#
# This module defines
#  AMD_INCLUDE_DIR, where to find umfpack.h, etc.
#  AMD_LIBRARIES, the libraries to link against to use AMD.
#  AMD_FOUND, If false, do not try to use AMD.
# also defined, but not for general use are
#  AMD_LIBRARY, where to find the AMD library.
# None of the above will be defined unless UFconfig can be found.
# AMD depends on  UFConfig

#=============================================================================
# Copyright 2010, Martin Koehler
# http://www-user.tu-chemnitz.de/~komart/
#
# Distributed under the OSI-approved BSD License (the "License");
# see accompanying file Copyright.txt for details.
#
# This software is distributed WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the License for more information.
#=============================================================================
# Changelog:
#    - Apr. 1, 2011 	Martin Koehler
#    - June 25, 2013 	Martin Koehler, add _libdir and _incdir 

  if (WIN32)
    set(_libdir ENV LIB)
    set(_liblist $ENV{LIB})
    set(_incdir) 
    foreach ( dir ${_liblist}) 
    	set(_incdir ${_incdir} "${dir}/../include")
    endforeach() 
    set(_incdir "${_incdir}" ENV INC ENV INCLUDE ENV CPATH)
  elseif (APPLE)
    set(_libdir ENV DYLD_LIBRARY_PATH)
    string(REPLACE ":" ";" _liblist $ENV{DYLD_LIBRARY_PATH} "")
    set(_incdir) 
    foreach ( dir ${_liblist}) 
    	set(_incdir ${_incdir} "${dir}/../include")
    endforeach() 
    set(_incdir "${_incdir}" ENV INC ENV INCLUDE ENV CPATH)
  else ()
    set(_libdir ENV LD_LIBRARY_PATH)
    string(REPLACE ":" ";" _liblist $ENV{LD_LIBRARY_PATH} "")
    set(_incdir) 
    foreach ( dir ${_liblist}) 
    	set(_incdir ${_incdir} "${dir}/../include")
    endforeach() 
    set(_incdir "${_incdir}" ENV INC ENV INCLUDE ENV CPATH)
  endif ()


find_path(AMD_AMD_INCLUDE_DIR amd.h
	HINTS ${SUITESPARSE}/AMD/Include  #Local Setup
	${SUITESPARSE}/include
	${_incdir}
	/usr/local/include/suitesparse 	#FreeBSD
	/usr/include/suitesparse	#Debian
	/opt/local/include/ufsparse	#Macports
	NO_DEFAULT_PATH
  )
find_path(AMD_AMD_INCLUDE_DIR amd.h) 


set(AMD_NAMES ${AMD_NAMES} libamd amd)
set(AMD_PATH
	${SUITESPARSE}/AMD/Lib 
	${SUITESPARSE}/lib
 	/opt/local/lib	# Macports
	${_libdir} 
)
find_library(AMD_LIBRARY NAMES  ${AMD_NAMES} PATHS ${AMD_PATH} NO_DEFAULT_PATH )
#MESSAGE(STATUS  "AMD_LIB: ${AMD_LIBRARY}") 
find_library(AMD_LIBRARY NAMES  ${AMD_NAMES} )

#MESSAGE(STATUS  "AMD_INC: ${AMD_AMD_INCLUDE_DIR}" )
#MESSAGE(STATUS  "AMD_LIB: ${AMD_LIBRARY}") 

if (AMD_LIBRARY AND AMD_AMD_INCLUDE_DIR)
      SET(AMD_INCLUDE_DIR ${AMD_AMD_INCLUDE_DIR} )
      SET(AMD_LIBRARIES ${AMD_LIBRARY} )
endif (AMD_LIBRARY AND AMD_AMD_INCLUDE_DIR)


# handle the QUIETLY and REQUIRED arguments and set AMD_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(AMD  DEFAULT_MSG AMD_LIBRARY AMD_AMD_INCLUDE_DIR)

mark_as_advanced(AMD_AMD_INCLUDE_DIR AMD_LIBRARY )


