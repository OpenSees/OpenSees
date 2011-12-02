include ../../Makefile.def

OBJS       = OPS_Stream.o \
	StandardStream.o \
	FileStream.o \
	XmlFileStream.o \
	DataFileStream.o \
	DatabaseStream.o \
	DummyStream.o

TEST_OBJS = $(OBJS) TestDataOutputStreamHandler.o  TestDataOutputFileHandler.o TestDataOutputDatabaseHandler.o
# Compilation control

all:         $(OBJS)

test:  $(TEST_OBJS)
	$(LINKER) $(LINKFLAGS) TestDataOutputStreamHandler.o $(OBJS) $(FE_LIBRARY) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o testDataStreamHandler
	$(LINKER) $(LINKFLAGS) TestDataOutputFileHandler.o $(OBJS) $(FE_LIBRARY) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o testDataFileHandler

#	$(LINKER) $(LINKFLAGS) TestDataOutputDatabaseHandler.o $(OBJS) $(FE_LIBRARY) \
#	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
#	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
#	 -o testDataDatabaseHandler

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o test* unitTest* UnitTest*

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
