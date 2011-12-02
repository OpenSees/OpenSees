include ../../../Makefile.def

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS	=	Channel.o TCP_Socket.o UDP_Socket.o MPI_Channel.o

else

OBJS	=	Channel.o TCP_Socket.o UDP_Socket.o 

endif


# Compilation control

all:       $(OBJS)

mpi: MPI_Channel.o

tcp: TCP_Socket.o UDP_Socket.o

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.



