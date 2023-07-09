include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)


OBJS 		= elasticPPc.o elasticPlaneStressC.o

all: $(OBJS)
	$(CC++) -shared -Wl,-soname,elasticPPc.so -o elasticPPc.so elasticPPc.o -lc
	$(CC++) -shared -Wl,-soname,elasticPPc.so -o elasticPlaneStressC.so elasticPlaneStressC.o -lc

elasticPPc.o: elasticPPc.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall elasticPPc.c

elasticPlaneStressC.o: elasticPlaneStressC.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall elasticPlaneStressC.c

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
