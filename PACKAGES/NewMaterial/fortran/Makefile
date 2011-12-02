include ../../../Makefile.def

OS_FLAG ?= -D_LINUX

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY         = elasticppf.so

all: $(LIBRARY)

OBJS = materialtypes.mod materialapi.mod elasticPPf.o 


$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY)   elasticPPf.o -lc

else

LIBRARY         = elasticPPf.dylib

OUTSIDE_OBJS = $(FE)/../PACKAGES/dynamicLoadFunctions.o

OBJS = materialtypes.mod materialapi.mod 

all: $(OBJS)
	$(FC) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace elasticPPf.f $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)


endif

materialapi.mod: materialAPI.f materialtypee.mod
	$(FC)  -c  materialAPI.f

materialtypes.mod: materialTypes.f
	$(FC)  -c  materialTypes.f

elasticPPf.o: elasticPPC.f materialapi.mod
	$(FC) -shared -fPIC $(INCLUDES)  -c -Wall elasticPPf.f

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
