include ../../../Makefile.def

OBJS       = AC3D8HexWithSensitivity.o \
	ASI3D8QuadWithSensitivity.o \
	AV3D4QuadWithSensitivity.o \
	VS3D4QuadWithSensitivity.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
