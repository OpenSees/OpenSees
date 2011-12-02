#  Written: Boris Jeremic@ucdavis.edu
#  Created: 10/97

include ../../Makefile.def

NDARRAY  = BJtensor.o \
           BJmatrix.o \
	   BJvector.o \
	   nDarray.o  \
	   basics.o   \
	   stresst.o  \
	   straint.o  \
	   Cosseratstresst.o  \
	   Cosseratstraint.o  

OBJS       = $(NDARRAY)

################### TARGETS ########################
all: $(OBJS)
	@$(RM) $(RMFLAGS) $(NDARRAY_LIBRARY);
	$(AR) $(ARFLAGS) $(NDARRAY_LIBRARY) $(NDARRAY)

libndarray.a: $(NDARRAY)
	$(AR) $(ARFLAGS) $(NDARRAY_LIBRARY) $(NDARRAY)


ndarray: libndarray.a

math_tst: $(NDARRAY) math_tst.o
	$(LINKER) $(LINKFLAGS) math_tst.o $(NDARRAY) -o math_tst

elast_tst: $(NDARRAY) elast_tst.o
	$(LINKER) $(LINKFLAGS) elast_tst.o $(NDARRAY) -o elast_tst

el-pl-test: $(NDARRAY) el-pl-test.o
	$(LINKER) $(LINKFLAGS) el-pl-test.o $(NDARRAY) -o el-pl-test

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core test
	@$(RM) $(RMFLAGS) math_tst elast_tst 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o core
	@$(RM) $(RMFLAGS) math_tst  
	@$(RM) $(RMFLAGS) elast_tst 
	@$(RM) $(RMFLAGS) el-pl-test 
	@$(RM) $(RMFLAGS) test.\$\$\$ 
	@$(RM) $(RMFLAGS) *.ti 
	@$(RM) $(RMFLAGS) tca.map 
	@$(RM) $(RMFLAGS) *.da 
	@$(RM) $(RMFLAGS) *.gcov
	@$(RM) $(RMFLAGS) *.bb
	@$(RM) $(RMFLAGS) *.bbg
	@$(RM) $(RMFLAGS) gmon.out
	@$(RM) $(RMFLAGS) .inslog*

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
