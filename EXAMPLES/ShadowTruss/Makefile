include ../../Makefile.def

SERVER_PROGRAM         = serverTruss

CLIENT_PROGRAM	       = example1

OBJS = ActorTruss.o ShadowTruss.o ClientMain.o ServerMain.o

all: $(OBJS)
	@$(LINKER) $(LINKFLAGS) ClientMain.o  ShadowTruss.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) -o $(CLIENT_PROGRAM)
	@$(LINKER) $(LINKFLAGS) ServerMain.o  ActorTruss.o\
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) -o $(SERVER_PROGRAM)

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core

spotless: clean
	@$(RM) $(RMFLAGS) fake core

wipe: spotless
	@$(RM) $(RMFLAGS) fake core $(CLIENT_PROGRAM) $(SERVER_PROGRAM)

# DO NOT DELETE THIS LINE -- make depend depends on it.
