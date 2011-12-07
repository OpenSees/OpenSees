include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:  $(LIBRARY_F)

LIBRARY_F    	= elasticPPf.so

OBJS_F 		= elasticPPF.o 


$(LIBRARY_F): $(OBJS_F)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)  elasticPPF.o -lc

else

LIBRARY_F    	= elasticPPf.dylib

OTHER_OBJS = dynamicLoadFunctions.o
OBJS =  materialtypes.mod materialapi.mod 


all: $(OBJS)
	$(FC) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace elasticPPf.f $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY_F)

materialapi.mod: materialAPI.f materialtypes.mod
	$(FC)  -c  materialAPI.f

materialtypes.mod: materialTypes.f
	$(FC)  -c  materialTypes.f

endif

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so *.dylib

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
