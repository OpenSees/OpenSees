include ../../../Makefile.def


OBJS       =  

#OBJS       = beam2d02.o beam2d03.o beam2d04.o ElasticBeam2d.o
#OBJS       =  beam2d03.o beam2d04.o beam2d02.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
