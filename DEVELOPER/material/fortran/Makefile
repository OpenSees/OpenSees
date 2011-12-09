include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:  $(LIBRARY_F)

LIBRARY    	= elasticPPf.so
OBJS 		= materialtypes.mod materialapi.mod elasticPPf.o 


all: $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) elasticPPf.o -lc

elasticPPf.o:
	$(FC) -fPIC -c elasticPPf.f 

materialapi.mod: materialAPI.f materialtypes.mod
	$(FC)  -c  materialAPI.f

materialtypes.mod: materialTypes.f
	$(FC)  -c  materialTypes.f

else

LIBRARY    	= elasticPPf.dylib

OTHER_OBJS = dynamicLoadFunctions.o
OBJS =  materialtypes.mod materialapi.mod 


all: $(OBJS)
	$(FC) $(OSLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace elasticPPf.f $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)

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
