include ../../Makefile.def

OBJS       = G3string.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core go

clean: tidy

	@$(RM) $(RMFLAGS) $(OBJS) *.o test go

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



