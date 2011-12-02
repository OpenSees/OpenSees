include ../../../Makefile.def

OBJS = MachineBroker.o AlphaMachineBroker.o DecMachineBroker.o MillMachineBroker.o

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



