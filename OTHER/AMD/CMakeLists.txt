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
#
#  NOTE: THIS BUILD SOLUTION HAS BEEN DEVELOPED FOR USE WITH THE OPENSEES 
#     FRAMEWORK AND MAY DIFFER CONSIDERABLY FROM THE UPSTREAM LIBRARY.
# 
#           AMD Version 2.4, Copyright (c) 1996-2013 by Timothy A. Davis
#                  people.engr.tamu.edu/davis/suitesparse.html
#==============================================================================
# Claudio Perez

project(AMD)

add_library(AMD)

target_compile_options(AMD PRIVATE -fexceptions -fPIC -I.)

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------
target_sources(AMD PRIVATE
  amd_aat.c
  amd_1.c
  amd_2.c
  amd_dump.c
  amd_postorder.c
  amd_post_tree.c
  amd_defaults.c
  amd_order.c
  amd_control.c
  amd_info.c
  amd_valid.c
  amd_preprocess.c
  SuiteSparse_config.c
)


#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

#amd_i_%.o: amd_%.c $(INC)
#	$(C) -DDINT -c $< -o $@
#
## amd_l_%.o: amd_%.c $(INC)
## 	$(C) -DDLONG -c $< -o $@
#
#$(AMD_LIBRARY): $(AMDI) SuiteSparse_config.o
#	@$(RM) $(RMFLAGS) $(AMD_LIBRARY)
#	@$(AR) $(ARFLAGS) $(AMD_LIBRARY) $(AMDI) SuiteSparse_config.o
#	@$(RANLIB) $(AMD_LIBRARY)


