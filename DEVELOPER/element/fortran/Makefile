include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY    	= trussf.so
OBJS 		= elementtypes.mod elementapi.mod trussF.o 

all: $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)  trussF.o -lc



elementapi.mod: elementAPI.f elementtypes.mod
	$(FC)  -c  elementAPI.f

elementtypes.mod: elementTypes.f
	$(FC)  -c  elementTypes.f

trussF.o:
	$(FC) -fPIC -c trussF.f 


else

LIBRARY    	= trussf.dylib

OTHER_OBJS = dynamicLoadFunctions.o
OBJS = elementtypes.mod elementapi.mod 


all: $(OBJS)
	$(FC) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace trussf.f $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY_F)

elementapi.mod: elementAPI.f elementtypes.mod
	$(FC)  -c  elementAPI.f

elementtypes.mod: elementTypes.f
	$(FC)  -c  elementTypes.f

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
