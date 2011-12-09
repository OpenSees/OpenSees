include ../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

LIBRARY         = Trapezoidal.so

OBJS = Trapezoidal.o 

all: $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

Trapezoidal.o: Trapezoidal.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall Trapezoidal.cpp

else

LIBRARY         = Trapezoidal.dylib


all: 
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace Trapezoidal.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY)

endif

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
