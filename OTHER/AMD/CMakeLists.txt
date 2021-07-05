project(AMD)

add_library(AMD)


add_compiler_flags(-fexceptions -fPIC -I.)

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------
target_sources(AMD
    amd_aat.c amd_1.c amd_2.c amd_dump.c amd_postorder.c amd_post_tree.c amd_defaults.c \
    amd_order.c amd_control.c amd_info.c amd_valid.c amd_preprocess
)

amd.h amd_internal.h SuiteSparse_config.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

#AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
#AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

#-------------------------------------------------------------------------------
# compile each int and long routine (with no real/complex version)
#-------------------------------------------------------------------------------

amd_i_%.o: amd_%.c $(INC)
	$(C) -DDINT -c $< -o $@

# amd_l_%.o: amd_%.c $(INC)
# 	$(C) -DDLONG -c $< -o $@

$(AMD_LIBRARY): $(AMDI) SuiteSparse_config.o
	@$(RM) $(RMFLAGS) $(AMD_LIBRARY)
	@$(AR) $(ARFLAGS) $(AMD_LIBRARY) $(AMDI) SuiteSparse_config.o
	@$(RANLIB) $(AMD_LIBRARY)


