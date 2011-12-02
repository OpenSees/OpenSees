include ../../Makefile.def

AMD	= amd.o \
	amd_1.o \
	amd_2.o \
	amd_defaults.o \
	amd_preprocess.o \
	amd_dump.o \
	amd_order.o \
	amd_valid.o \
	amd_aat.o \
	amd_global.o \
	amd_post_tree.o \
	amdbar.o \
	amd_control.o \
	amd_info.o \
	amd_postorder.o

all: $(AMD) 
	@$(RM) $(RMFLAGS) $(AMD_LIBRARY)
	@$(AR) $(ARFLAGS) $(AMD_LIBRARY) $(AMD)
	@$(RANLIB) $(AMD_LIBRARY)

clean:
	@$(RM) $(RMFLAGS) *.o

wipe: clean
