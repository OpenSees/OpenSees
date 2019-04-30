############################################################################
#
#  Program:  OpenSees
#
#  Purpose:  A Top-level Makefile to create the libraries needed
#	     to use the OpenSees framework.
#
#  Written: fmk 
#  Created: 10/99
#
#  Send bug reports, comments or suggestions to fmckenna@ce.berkeley.edu
#
############################################################################


include Makefile.def


############################################################################
#
#  First, modify the definitions in Makefile.def to match 
#  your library archiver, compiler, and the options to be used.
#
#  Sample Makefile.def's can be found in the directory MAKES
#  
#  Then to create or add to the libraries needed, enter make after having
#  making the modifications to this file.
#
#  The name of the libraries created and their placement are defined 
#  in the file called Makefile.def.
#
#  To remove the object files after the libraries and testing executables
#  are created, enter
#       make clean
#  To remove the object files and the libraries specified in WIPE_LIBS, enter
#       make wipe
#  To just make the libs, enter 
#	make lib
#  To just build the interpreter type
#	make OpenSees
############################################################################

all: 
ifdef MKDIR
	$(MKDIR) $(HOME)/bin
	$(MKDIR) $(HOME)/lib
endif
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );
	@$(ECHO) LIBRARIES BUILT ... NOW LINKING OpenSees PROGRAM;
	@$(CD) $(FE)/tcl; $(MAKE) tcl;
	@$(CD) $(FE)/modelbuilder/tcl;  $(MAKE) tcl;

OpenSees: tcl

tcl:
ifdef MKDIR
	$(MKDIR) $(HOME)/bin
	$(MKDIR) $(HOME)/lib
endif
	@$(ECHO) Building OpenSees Program ..;
	@$(CD) $(FE)/tcl;  $(MAKE) tcl;
	@$(CD) $(FE)/modelbuilder/tcl;  $(MAKE) tcl;

OpenSeesTk: tk

tk:
	@$(ECHO) Building OpenSees Program ..;
	@$(CD) $(FE)/tcl;  $(MAKE) tk;
	@$(CD) $(FE)/modelbuilder/tcl;  $(MAKE) tk;

OpenSeesPy: python

python: 
ifdef MKDIR
	$(MKDIR) $(HOME)/bin
	$(MKDIR) $(HOME)/lib
endif
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );
	@$(ECHO) LIBRARIES BUILT ... NOW LINKING OpenSeesPy Module;
	@$(CD) $(FE)/interpreter; $(MAKE) pythonmodule;

libs:
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(MAKE); \
		$(CD) ..; \
	done );

clean:
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(ECHO) Making lib in $$f; \
		$(MAKE) clean; \
		$(CD) ..; \
	done );
	@$(RM) $(RMFLAGS) *.o *~ core
	@$(CD) $(FE)/../EXAMPLES;  $(MAKE) wipe;

wipe: 
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(ECHO) Making lib in $$f; \
		$(MAKE) wipe; \
		$(MAKE) clean; \
		$(CD) ..; \
	done );
	@$(RM) $(RMFLAGS) $(WIPE_LIBS) *.o *~ core 
	@$(CD) $(FE)/../EXAMPLES;  $(MAKE) wipe;

wipeall: 
	@( \
	for f in $(DIRS); \
	do \
		$(CD) $$f; \
		$(ECHO) Making lib in $$f; \
		$(MAKE) wipe; \
		$(CD) ..; \
	done );
	@$(RM) $(RMFLAGS) $(WIPE_LIBS) *.o *~ core
	@$(CD) $(FE)/../EXAMPLES;  $(MAKE) wipe
	@$(RM) $(RMFLAGS) $(OpenSees_PROGRAM);

help:
    @$(ECHO) "usage: make ?"






