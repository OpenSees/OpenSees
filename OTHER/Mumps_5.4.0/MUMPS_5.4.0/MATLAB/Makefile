# Please only change make.inc, not this Makefile
include make.inc

# MUMPS include files
INCMUMPS = -I$(MUMPS_DIR)/include 

# MUMPS libraries
LIBMUMPS = -L$(MUMPS_DIR)/lib -l$(ARITH)mumps -lmumps_common

# Stub MPI/BLACS/ScaLAPACK
INCSEQ     = -I$(MUMPS_DIR)/libseq
LIBSEQ     = -L$(MUMPS_DIR)/libseq -lmpiseq

# MUMPS includes
INC = $(INCMUMPS) $(IORDERINGS) $(INCSEQ)
LIB = $(LIBMUMPS) $(LORDERINGS) $(LIBSEQ) $(LIBBLAS) $(LIBFORT)

all:	d z

d:
	$(MAKE) ARITH=d dmumpsmex.stamp
z:
	$(MAKE) ARITH=z zmumpsmex.stamp

clean:
	rm -f dmumpsmex.* zmumpsmex*

$(ARITH)mumpsmex.stamp: mumpsmex.c
	cp -f mumpsmex.c $(ARITH)mumpsmex.c
	$(MEX) $(OPTC) $(ARITH)mumpsmex.c -DMUMPS_ARITH=MUMPS_ARITH_$(ARITH) $(INC) $(LIB)
	rm -f $(ARITH)mumpsmex.c
	touch $@

