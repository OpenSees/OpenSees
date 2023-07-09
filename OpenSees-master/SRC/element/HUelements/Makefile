include ../../../Makefile.def

OBJS       = MultipleShearSpring.o \
	MultipleNormalSpring.o \
	DBESI0.o \
	DBESI1.o \
	YamamotoBiaxialHDR.o \
	KikuchiBearing.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
