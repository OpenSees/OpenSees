include ../../../../Makefile.def

OBJS       = PySimple1.o \
	QzSimple1.o \
	TzSimple1.o \
	PySimple2.o \
	PySimple3.o \
	QzSimple2.o \
	TzSimple2.o \
	TclPyTzQzMaterialCommand.o \
	ShallowFoundationGen.o \
	PyLiq1.o \
	TzLiq1.o \
	PySimple1Gen.o \
	TzSimple1Gen.o 

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) 
	@$(RM) $(RMFLAGS) *.o
	@$(RM) $(RMFLAGS) .inslog*
	@$(RM) $(RMFLAGS) tca.map

spotless: clean					   

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
