# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM         = test

OBJS = 

ifdef MUMPS

OBJS       =


ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS       = MumpsSOE.o \
	MumpsSolver.o \
	MumpsParallelSOE.o \
	MumpsParallelSolver.o 

endif

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

OBJS       = MumpsSOE.o \
	MumpsSolver.o \
	MumpsParallelSOE.o \
	MumpsParallelSolver.o 

endif
endif



GRAPHIC_LIBS_TCL = $(FE)/tcl/TclFeViewer.o \
	$(FE)/renderer/PlainMap.o 

all:         $(OBJS)

test: $(OBJS) example.o
	@$(CC) $(CFLAGS) $(INCLUDES) -c example.cpp
	@$(LINKER) $(LINKFLAGS) example.o  \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) $(GRAPHIC_LIBS_TCL) \
	$(MACHINE_SPECIFIC_LIBS) -o test

test1: 
	@$(CC) $(CFLAGS) $(INCLUDES) -I$(HOME)/MUMPS_4.6.4/include -c c_example.c   
	@$(LINKER) $(LINKFLAGS) c_example.o  mumps.o \
	-L$(HOME)/MUMPS_4.6.4/lib \
	-ldmumps -lpord  -lessl -lpesslsmp -lblacssmp \
	$(HPM_LIB) \
	-lm \
	$(MACHINE_SPECIFIC_LIBS) -o test1

test2: c_example1.o
	@$(LINKER) $(LINKFLAGS) c_example1.o  mumps.o \
	-L$(HOME)/MUMPS_4.6.4/lib \
	-ldmumps -lpord  -lessl -lpesslsmp -lblacssmp \
	$(HPM_LIB) \
	$(FE)/actor/channel/MPI_Channel.o \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) $(PARALLEL_LIB) \
	$(MACHINE_NUMERICAL_LIBS) $(TCL_LIBRARY) $(GRAPHIC_LIBS_TCL) \
	$(MACHINE_SPECIFIC_LIBS) -o test2

#	-lm \
#	$(MACHINE_SPECIFIC_LIBS) -o test2

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
