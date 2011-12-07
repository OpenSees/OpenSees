include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:   $(LIBRARY_CPP) $(LIBRARY_C) $(LIBRARY_F)

LIBRARY_C    	= elasticPPc.so

OBJS_C 		= elasticPPc.o elasticPlaneStress.o

$(LIBRARY_C): $(OBJS_C)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY_CPP) $(OBJS_CPP) -lc

elasticPPC.o: elasticPPC.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall elasticPPC.c

elasticPlaneStressC.o: elasticPlaneStressC.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall elasticPlaneStress.c

else

LIBRARY_C    	= elasticPPC.dylib


all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace elasticPPC.c $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o elasticPPC.dylib
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace elasticPlaneStressC.c $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o elasticPlaneStressC.dylib

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
