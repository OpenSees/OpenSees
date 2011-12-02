include ../../../Makefile.def


OS_FLAG ?= -D_LINUX

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY         = trussf.so

all: $(LIBRARY)

OBJS = elementtypes.mod elementapi.mod trussF.o 


$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)  trussF.o -lc

else

LIBRARY         = trussf.dylib

OUTSIDE_OBJS = $(FE)/../PACKAGES/dynamicLoadFunctions.o

OBJS = elementtypes.mod elementapi.mod 

all: $(OBJS)
	$(FC) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace trussf.f $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)

endif

elementapi.mod: elementAPI.f elementtypes.mod
	$(FC)  -c  elementAPI.f

elementtypes.mod: elementTypes.f
	$(FC)  -c  elementTypes.f

trussF.o: trussF.f elementapi.mod
	$(FC) -shared -fPIC $(INCLUDES)  -c -Wall trussF.f

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core *.out *.so

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) $(PROGRAM) fake core $(LIBRARY)

# DO NOT DELETE THIS LINE -- make depend depends on it.
