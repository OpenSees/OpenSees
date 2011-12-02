include ../../../Makefile.def

OBJS	=   ChannelAddress.o 



all:            $(OBJS)

mpi:	 MPI_ChannelAddress.o

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.



