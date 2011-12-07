include ../../Makefile.def

ifeq ($(OS_FLAG), -D_LINUX)

all:   $(LIBRARY_CPP) $(LIBRARY_C) $(LIBRARY_F)

LIBRARY_CPP 	= ElasticPPcpp.so
OBJS_CPP	= ElasticPPcpp.o 

$(LIBRARY_CPP): $(OBJS_CPP)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY_CPP) $(OBJS_CPP) -lc

ElasticPPcpp.o: ElasticPPcpp.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall ElasticPPcpp.cpp

elasticPPC.o: elasticPPC.c
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall elasticPPC.c

else

LIBRARY_CPP 	= ElasticPPCPP.dylib


all: $(OBJS)
	$(CC++) $(OS_FLAG) -dynamiclib $(INCLUDES) -Wl,-undefined,suppress,-flat_namespace ElasticPPcpp.cpp $(OUTSIDE_OBJS)  -current_version 1.0 -compatibility_version 1.0 -fvisibility=hidden -o $(LIBRARY_CPP)


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
