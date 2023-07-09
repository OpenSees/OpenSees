include ../../../Makefile.def

OBJS = MachineBroker.o 

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS = MachineBroker.o MPI_MachineBroker.o

endif

ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

OBJS = MachineBroker.o MPI_MachineBroker.o

endif




all: $(OBJS)

mpi: MPI_MachineBroker.o

gexec: GEXEC_MachineBroker.o

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) fake

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.



