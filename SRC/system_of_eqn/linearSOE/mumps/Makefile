# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM         = test

OBJS       = MumpsSOE.o \
	MumpsSolver.o 

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = MumpsSOE.o \
	MumpsSolver.o \
	MumpsParallelSOE.o \
	MumpsParallelSolver.o 

endif

all:         $(OBJS)

test: example.cpp $(OBJS)
	@$(CC) $(CFLAGS) $(INCLUDES) -c example.cpp
	@$(LINKER) $(LINKFLAGS) example.o $(OBJS) \
	-L$(HOME)/MUMPS_4.6/lib \
	-ldmumps -lpord  -lessl -lpesslsmp -lblacssmp \
	$(HPM_LIB) \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) \
	$(MACHINE_SPECIFIC_LIBS) -o $(PROGRAM)

test1: c_example.c
	@$(CC) $(CFLAGS) $(INCLUDES) -c c_example.c   
	@$(LINKER) $(LINKFLAGS) c_example.o  \
	-L$(HOME)/MUMPS_4.6/lib \
	-ldmumps -lpord  -lessl -lpesslsmp -lblacssmp \
	$(HPM_LIB) \
	-lm \
	$(MACHINE_SPECIFIC_LIBS) -o $(PROGRAM)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
