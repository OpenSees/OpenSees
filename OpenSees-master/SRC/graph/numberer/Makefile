include ../../../Makefile.def

OBJS       = RCM.o \
	AMDNumberer.o \
	SimpleNumberer.o \
	GraphNumberer.o \
	MyRCM.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless


# DO NOT DELETE THIS LINE -- make depend depends on it.


