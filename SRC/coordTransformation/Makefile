include ../../Makefile.def

OBJS       = CrdTransf.o CrdTransf2d.o LinearCrdTransf2d.o \
                         CrdTransf3d.o LinearCrdTransf3d.o 

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.






