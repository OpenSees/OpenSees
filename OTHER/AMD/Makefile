include ../../Makefile.def

all: $(AMD_LIBRARY)

C = $(CC) $(OPT_FLAG) -fexceptions -fPIC -I.

#-------------------------------------------------------------------------------
# source files
#-------------------------------------------------------------------------------

AMD = amd_aat amd_1 amd_2 amd_dump amd_postorder amd_post_tree amd_defaults \
	amd_order amd_control amd_info amd_valid amd_preprocess

INC = amd.h amd_internal.h SuiteSparse_config.h

#-------------------------------------------------------------------------------
# object files for each version
#-------------------------------------------------------------------------------

AMDI = $(addsuffix .o, $(subst amd_,amd_i_,$(AMD)))
AMDL = $(addsuffix .o, $(subst amd_,amd_l_,$(AMD)))

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

clean:
	@$(RM) $(RMFLAGS) *.o

wipe: clean
