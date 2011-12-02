include ../../Makefile.def

LIBRARY         = myTruss.so

all:   $(LIBRARY)

OBJS = MyTruss.o TclMyTruss.o 

$(LIBRARY): $(OBJS)
	$(CC++) -shared -Wl,-soname,$(LIBRARY) -o $(LIBRARY) $(OBJS) -lc

MyTruss.o: MyTruss.cpp
	$(CC++) -fPIC $(INCLUDES) -g -c -Wall MyTruss.cpp

TclMyTruss.o: TclMyTruss.cpp
	$(CC) -fPIC $(INCLUDES) -g -c -Wall TclMyTruss.cpp

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
