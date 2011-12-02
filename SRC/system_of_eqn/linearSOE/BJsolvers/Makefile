# Makefile for fe objects

include ../../../../Makefile.def

PROGRAM = test

OBJS       = skymatr.o  profmatr.o


################### TARGETS ########################
all: $(OBJS) 


libBJsolvers.a: $(OBJS)
	$(AR) -rv $(HOME)/lib/libBJsolvers.a $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean
	@$(RM) $(RMFLAGS) $(PROGRAM) 

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
