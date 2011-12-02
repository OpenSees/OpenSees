include ../../Makefile.def

LIBRARY         = MyElasticPP.so

all: $(LIBRARY)

OBJS = MyElasticPP.o TclMyElasticPP.o 

$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

MyElasticPP.o: MyElasticPP.cpp
	$(CC++) -fPIC $(INCLUDES)  -c -Wall MyElasticPP.cpp

TclMyElasticPP.o: TclMyElasticPP.cpp
	$(CC) -fPIC $(INCLUDES)  -c -Wall TclMyElasticPP.cpp


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
