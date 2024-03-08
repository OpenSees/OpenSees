include ../../../../Makefile.def

OBJS       = ReinforcedConcreteLayeredMembraneSection.o\
             LayeredMembraneSection.o\
			 ElasticMembraneSection.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.








