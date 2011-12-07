include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:   $(LIBRARY_CPP) $(LIBRARY_C) $(LIBRARY_F)

LIBRARY_CPP 	= Truss2D.so

OBJS_CPP	= Truss2D.o 

$(LIBRARY_CPP): $(OBJS_CPP)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY_CPP) $(OBJS_CPP) -lc

Truss2D.o: Truss2D.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall Truss2D.cpp

else

LIBRARY_CPP 	= Truss2D.dylib

all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace Truss2D.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY_CPP)

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


