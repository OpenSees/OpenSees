include ../../Makefile.def

OBJS       = etime.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy

	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



