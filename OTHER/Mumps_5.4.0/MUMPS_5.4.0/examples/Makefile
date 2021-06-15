#
#  This file is part of MUMPS 5.4.0, released
#  on Tue Apr 13 15:26:30 UTC 2021
#
topdir = ..
libdir = $(topdir)/lib

default: d

.PHONY: default all s d c z multi clean
.SECONDEXPANSION:

all:	c z s d multi

c:	csimpletest csimpletest_save_restore
z:	zsimpletest zsimpletest_save_restore
s:	ssimpletest ssimpletest_save_restore
d:	dsimpletest dsimpletest_save_restore c_example_save_restore c_example
multi:	multiple_arithmetics_example


include $(topdir)/Makefile.inc

LIBMUMPS_COMMON = $(libdir)/libmumps_common$(PLAT)$(LIBEXT)


LIBSMUMPS = $(libdir)/libsmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

ssimpletest:  $(LIBSMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)


LIBDMUMPS = $(libdir)/libdmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

dsimpletest: $(LIBDMUMPS)  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest.o  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)


LIBCMUMPS = $(libdir)/libcmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

csimpletest: $(LIBCMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)


LIBZMUMPS = $(libdir)/libzmumps$(PLAT)$(LIBEXT) $(LIBMUMPS_COMMON)

zsimpletest: $(LIBZMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest.o  $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

c_example:	$(LIBDMUMPS) $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)


multiple_arithmetics_example:	$(LIBSMUMPS) $(LIBDMUMPS) $(LIBCMUMPS) $(LIBZMUMPS) $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBSMUMPS) $(LIBDMUMPS) $(LIBCMUMPS) $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

ssimpletest_save_restore:  $(LIBSMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) ssimpletest_save_restore.o  $(LIBSMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

dsimpletest_save_restore: $(LIBDMUMPS)  $$@.o 
	$(FL) -o $@ $(OPTL) dsimpletest_save_restore.o  $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

csimpletest_save_restore: $(LIBCMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) csimpletest_save_restore.o  $(LIBCMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

zsimpletest_save_restore: $(LIBZMUMPS)  $$@.o
	$(FL) -o $@ $(OPTL) zsimpletest_save_restore.o  $(LIBZMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

c_example_save_restore:	$(LIBDMUMPS) $$@.o
	$(FL) -o $@ $(OPTL) $@.o $(LIBDMUMPS) $(LORDERINGS) $(LIBS) $(LIBBLAS) $(LIBOTHERS)

.SUFFIXES: .c .F .o
.F.o:
	$(FC) $(OPTF) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.F $(OUTF)$*.o
.c.o:
	$(CC) $(OPTC) $(CDEFS) -I. -I$(topdir)/include -I$(topdir)/src $(INCS) -c $*.c $(OUTC)$*.o


$(libdir)/libsmumps$(PLAT)$(LIBEXT):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

$(libdir)/libdmumps$(PLAT)$(LIBEXT):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

$(libdir)/libcmumps$(PLAT)$(LIBEXT):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

$(libdir)/libzmumps$(PLAT)$(LIBEXT):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

$(LIBMUMPS_COMMON):
	@echo 'Error: you should build the library' $@ 'first'
	exit 1

clean:
	$(RM) *.o [sdcz]simpletest c_example multiple_arithmetics_example ssimpletest_save_restore dsimpletest_save_restore csimpletest_save_restore zsimpletest_save_restore c_example_save_restore
