#
#  This file is part of MUMPS 5.4.1, released
#  on Tue Aug  3 09:49:43 UTC 2021
#
topdir = .
libdir = $(topdir)/lib

default: d

.PHONY: default all s d c z prerequisites libseqneeded clean

all: prerequisites
	cd src; $(MAKE) all
	cd examples; $(MAKE) all

s: prerequisites
	cd src; $(MAKE) s
	cd examples; $(MAKE) s

d: prerequisites
	cd src; $(MAKE) d
	cd examples; $(MAKE) d

c: prerequisites
	cd src; $(MAKE) c
	cd examples; $(MAKE) c

z: prerequisites
	cd src; $(MAKE) z
	cd examples; $(MAKE) z


# Is Makefile.inc available ?
Makefile.inc:
	@echo "######################################################################"
	@echo "# BEFORE COMPILING MUMPS, YOU MUST HAVE AN APPROPRIATE Makefile.inc"
	@echo "# FILE AVAILABLE. PLEASE CHECK THE DIRECTORY ./Make.inc FOR EXAMPLES"
	@echo "# OF Makefile.inc FILES, AND USE Make.inc/Makefile.inc.generic IF YOU"
	@echo "# NEED TO BUILD A NEW ONE. SEE ALSO THE README AND INSTALL FILES"
	@echo "######################################################################"
	@exit 1

include Makefile.inc

prerequisites: Makefile.inc $(LIBSEQNEEDED) $(libdir)/libpord$(PLAT)$(LIBEXT)

# dummy MPI library (sequential version)

libseqneeded:
	(cd libseq; $(MAKE))

# Build the libpord.a library and copy it into $(topdir)/lib
$(libdir)/libpord$(PLAT)$(LIBEXT):
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); \
	  $(MAKE) CC="$(CC)" CFLAGS="$(OPTC)" AR="$(AR)" RANLIB="$(RANLIB)" OUTC="$(OUTC)" LIBEXT=$(LIBEXT); \
	fi;
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cp $(LPORDDIR)/libpord$(LIBEXT) $@; \
	fi;

clean:
	(cd src; $(MAKE) clean)
	(cd examples; $(MAKE) clean)
	(cd $(libdir); $(RM) *$(PLAT)$(LIBEXT))
	(cd libseq; $(MAKE) clean)
	if [ "$(LPORDDIR)" != "" ] ; then \
	  cd $(LPORDDIR); $(MAKE) realclean; \
        fi;

