include ../../../Makefile.def

OBJS	=	Channel.o TCP_Socket.o UDP_Socket.o Socket.o HTTP.o 

ifeq ($(PROGRAMMING_MODE), PARALLEL)

OBJS	=	Channel.o TCP_Socket.o UDP_Socket.o MPI_Channel.o HTTP.o Socket.o

endif


ifeq ($(PROGRAMMING_MODE), PARALLEL_INTERPRETERS)

OBJS	=	Channel.o TCP_Socket.o UDP_Socket.o MPI_Channel.o HTTP.o Socket.o

endif


# Compilation control

all:       $(OBJS)

mpi: MPI_Channel.o

tcp: TCP_Socket.o UDP_Socket.o

test: Test.o HTTP.o Socket.o	
	$(LINKER) Test.o Socket.o HTTP.o $(FE)/utility/NeesCentral.o -l ssl -o a.out

# Miscellaneous
tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean:  tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.



