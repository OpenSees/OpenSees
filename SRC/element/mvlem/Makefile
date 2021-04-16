include ../../../Makefile.def

OBJS       = MVLEM.o \
	SFI_MVLEM.o \
	MVLEM_3D.o \
	SFI_MVLEM_3D.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
