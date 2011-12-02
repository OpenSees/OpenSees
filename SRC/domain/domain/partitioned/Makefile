include ../../../../Makefile.def

OBJS       = PartitionedDomainEleIter.o PartitionedDomainSubIter.o \
	PartitionedDomain.o	

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
