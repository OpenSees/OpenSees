include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:   $(LIBRARY_F)

LIBRARY_F    	= trussf.so

OBJS_F 		= elementtypes.mod elementapi.mod trussF.o 

$(LIBRARY_F): $(OBJS_F)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)  trussF.o -lc

else

LIBRARY_F    	= trussf.dylib

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
