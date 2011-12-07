include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:   $(LIBRARY_C) 

LIBRARY_C    	= trussC.so

OBJS_C 		= TrussC.o 

$(LIBRARY_C): $(OBJS_C)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY_CPP) $(OBJS_CPP) -lc

trussC.o: trussC.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall trussC.c

else

LIBRARY_C    	= trussC.dylib

all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace trussC.c $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY_C)

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
