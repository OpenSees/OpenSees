#  Written: fmk 
#  Created: 10/97

include ../../Makefile.def

OBJS       = FE_Datastore.o \
	FileDatastore.o \
	TclDatabaseCommands.o \
	NEESData.o

# Compilation control

all:         $(OBJS)


MySQL_LIBRARY = MySQL.so 

mysql: MySqlDatastore.o TclMySQL.o
	$(CC++) -shared -Wl,-soname,$(MySQL_LIBRARY) -o $(MySQL_LIBRARY) \
	TclMySQL.o MySqlDatastore.o /usr/local/mysql/lib/mysql/libmysqlclient.a -lc

MySqlDatastore.o: MySqlDatastore.cpp MySqlDatastore.h
	$(CC++) -fPIC $(INCLUDES)  -c -Wall MySqlDatastore.cpp

TclMySQL.o: TclMySQL.cpp
	$(CC) -fPIC $(INCLUDES)  -c -Wall TclMySQL.cpp


BerkeleyDB_LIBRARY = BerkeleyDB.so 

berkeleydb: BerkeleyDbDatastore.o TclBerkeleyDB.o
	$(CC++) -shared -Wl,-soname,$(BerkeleyDB_LIBRARY) -o $(BerkeleyDB_LIBRARY) \
	TclBerkeleyDB.o BerkeleyDbDatastore.o /usr/lib/libdb-4.2.a -lc

BerkeleyDbDatastore.o: BerkeleyDbDatastore.cpp BerkeleyDbDatastore.h
	$(CC++) -fPIC -D_BerkeleyDB41orGreater $(INCLUDES)  -c -Wall BerkeleyDbDatastore.cpp

TclBerkeleyDB.o: TclBerkeleyDB.cpp
	$(CC) -fPIC $(INCLUDES)  -c -Wall TclBerkeleyDB.cpp



test:  NeesDataTest.o 
	$(LINKER) $(LINKFLAGS) NeesDataTest.o NEESData.o $(FE_LIBRARY) \
	$(FE_LIBRARY) $(MACHINE_LINKLIBS) \
	$(MACHINE_NUMERICAL_LIBS) $(MACHINE_SPECIFIC_LIBS) $(METIS_LIBRARY) \
	 -o test

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o *.so

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



