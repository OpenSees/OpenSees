/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2002, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.4 $
// $Date: 2003-10-16 00:06:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/BerkeleyDbDatastore.cpp,v $

#include <BerkeleyDbDatastore.h>
#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <errno.h>
#include <unistd.h>
#include <sys/stat.h>

BerkeleyDbDatastore::BerkeleyDbDatastore(const char *projectName, 
					 Domain &theDomain,
					 FEM_ObjectBroker &theObjectBroker,
					 char *dbType)
  :FE_Datastore(theDomain, theObjectBroker), dbTag(0),
   connection(true)
{
  // copy the projectName
  project = new char[strlen(projectName) + 1];
  strcpy(project, projectName);

  // create the directory if the database directory does not exist
  struct stat sb;
  if (stat(project, &sb) != 0) {

    // Create the directory, read/write/access owner only.
    if (mkdir(project, S_IRWXU) != 0) {
      opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed mkdir: ";
      opserr << project << " " << strerror(errno);
      connection = false;
      return;
    }
  }

  // create the environment handle
  int ret = 0;
  if ((ret = db_env_create(&dbenv, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed db_env_create: ";
    opserr << db_strerror(ret) << endln;
    connection = false;
  }

  // Set up error handling.
  dbenv->set_errpfx(dbenv, project);

  // now create the enviroment
  if ((ret = dbenv->open(dbenv, project,
			 DB_CREATE | DB_INIT_LOG | DB_PRIVATE | DB_INIT_MPOOL, 
			 S_IRUSR | S_IWUSR)) != 0) {
			 
    opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed db_env_create: ";
    opserr << db_strerror(ret) << endln;
    connection = false;
    return;
  }


  // initialise the db structuress
  int result;
  if ((result = db_create(&dbMatrix, dbenv, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to init dbMatrix";
    opserr << db_strerror(result) << endln;
    connection = false;
    return;
  }
  if ((result = db_create(&dbVector, dbenv, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to init dbVector";
    opserr << db_strerror(result) << endln;
    connection = false;
    return;
  }
  if ((result = db_create(&dbID, dbenv, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to init dbID";
    opserr << db_strerror(result) << endln;
    connection = false;
    return;
  }


  // set the database storage type
  DBTYPE type = DB_HASH;
  if (dbType != NULL) {
    if ((strcmp(dbType, "hash") == 0) || (strcmp(dbType, "HASH") == 0))
      type = DB_HASH;
    else if ((strcmp(dbType, "queue") == 0) || (strcmp(dbType, "QUEUE") == 0))
      type = DB_QUEUE;	 
    else if ((strcmp(dbType, "btree") == 0) || (strcmp(dbType, "BTREE") == 0))
      type = DB_BTREE;
    else if ((strcmp(dbType, "recno") == 0) || (strcmp(dbType, "RECNO") == 0))
      type = DB_RECNO;
  }

  // open the database files .. use DB_UNKNOWN first to see if database already exists
#ifdef _BerkeleyDB41orGreater

  if ((result = dbMatrix->open(dbMatrix, NULL, "Matrices.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    //type = DB_UNKNOWN;
    if ((result = dbMatrix->open(dbMatrix, NULL, "Matrices.db", NULL, type, 
				 DB_CREATE, 0664)) != 0) {
#else
  if ((result = dbMatrix->open(dbMatrix, "Matrices.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    //type = DB_UNKNOWN;
    if ((result = dbMatrix->open(dbMatrix, "Matrices.db", NULL, type, 
				 DB_CREATE, 0664)) != 0) {
#endif

      opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to open dbMatrix\n";
      opserr << db_strerror(result) << endln;
      connection = false;     
      return;
    } 
  }


#ifdef _BerkeleyDB41orGreater

  if ((result = dbVector->open(dbVector, NULL, "Vectors.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    //type = DB_UNKNOWN;
    if ((result = dbVector->open(dbVector, NULL, "Vectors.db", NULL, type, 
				 DB_CREATE, 0664)) != 0) {

#else

  if ((result = dbVector->open(dbVector, "Vectors.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    //type = DB_UNKNOWN;
    if ((result = dbVector->open(dbVector, "Vectors.db", NULL, type, 
				 DB_CREATE, 0664)) != 0) {

#endif
      opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to open dbVector\n";
      opserr << db_strerror(result) << endln;
      connection = false;     
      return;
    } 
  }

#ifdef _BerkeleyDB41orGreater
  if ((result = dbID->open(dbID, NULL, "IDs.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    type = DB_UNKNOWN;
    if ((result = dbID->open(dbID, NULL, "IDs.db", NULL, type, 
			     0, 0664)) != 0) {
#else
  if ((result = dbID->open(dbID, "IDs.db", NULL, type, 
			       DB_CREATE | DB_EXCL, 0664)) != 0) {
    type = DB_UNKNOWN;
    if ((result = dbID->open(dbID, "IDs.db", NULL, type, 
			     0, 0664)) != 0) {
#endif

      opserr << "BerkeleyDbDatastore::BerkeleyDbDatastore - failed to open dbID\n";
      opserr << db_strerror(result) << endln;
      connection = false;     
      return;
    }  else
      result = dbID->get_type(dbID, &type);
  }
}


BerkeleyDbDatastore::~BerkeleyDbDatastore()
{
  opserr << "CLOSING DATABASE\n";

  if (connection == true) 
    dbMatrix->close(dbMatrix, 0);
  if (connection == true) 
    dbVector->close(dbVector, 0);
  if (connection == true) 
    dbID->close(dbID, 0);
  if (connection == true) 
    dbenv->close(dbenv, 0);

  if (project != 0)
    delete [] project;
}


int
BerkeleyDbDatastore::getDbTag(void)
{
  dbTag++;
  return dbTag;
}


int 
BerkeleyDbDatastore::sendMsg(int dataTag, int commitTag, 
		       const Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "BerkeleyDbDatastore::sendMsg() - not yet implemented\n";
  return -1;
}		       

int 
BerkeleyDbDatastore::recvMsg(int dataTag, int commitTag, 
		       Message &, 
		       ChannelAddress *theAddress)
{
  opserr << "BerkeleyDbDatastore::recvMsg() - not yet implemented\n";
  return -1;
}		       

int 
BerkeleyDbDatastore::sendMatrix(int dbTag, int commitTag, 
				const Matrix &theMatrix, 
				ChannelAddress *theAddress)
{
  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theMatrix.dataSize);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  
  
  // set up the data structure
  data.data = (void *)theMatrix.data;
  data.size = theMatrix.dataSize * sizeof(double);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  int ret;
  if ((ret = dbMatrix->put(dbMatrix, NULL, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::sendMatrix() - failed to send the Matrix to database\n";
    dbMatrix->err(dbMatrix, ret, "DB->put");
    return -2;      
  }

  return 0;
}
 

int 
BerkeleyDbDatastore::recvMatrix(int dbTag, int commitTag, 
				Matrix &theMatrix,
				ChannelAddress *theAddress)
{
  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theMatrix.dataSize);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  
  
  // set up the data structure
  data.data = (void *)theMatrix.data;
  data.size = theMatrix.dataSize * sizeof(double);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  int ret;
  if ((ret = dbMatrix->get(dbMatrix, NULL, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::recvMatrix() - failed to get the Matrix from database\n";
    dbMatrix->err(dbMatrix, ret, "DB->get");
    return -2;      
  }

  return 0;
}


int 
BerkeleyDbDatastore::sendVector(int dbTag, int commitTag, 
				const Vector &theVector, 
				ChannelAddress *theAddress)
{
  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theVector.sz);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  
  
  // set up the data structure
  data.data = (void *)theVector.theData;
  data.size = theVector.sz * sizeof(double);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  int ret;
  if ((ret = dbVector->put(dbVector, NULL, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::sendVector() - failed to send the Vector to database\n";
    dbVector->err(dbVector, ret, "DB->put");
    return -2;      
  }

  return 0;
}
 

int 
BerkeleyDbDatastore::recvVector(int dbTag, int commitTag, 
				Vector &theVector,
				ChannelAddress *theAddress)
{
  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theVector.sz);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  
  
  // set up the data structure
  data.data = (void *)theVector.theData;
  data.size = theVector.sz * sizeof(double);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  int ret;
  if ((ret = dbVector->get(dbVector, NULL, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::recvVector() - failed to get the Vector from database\n";
    dbVector->err(dbVector, ret, "DB->get");
    return -2;      
  }

  return 0;
}


int 
BerkeleyDbDatastore::sendID(int dbTag, int commitTag, 
			    const ID &theID, 
			    ChannelAddress *theAddress)
{
  int ret = 0;

  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theID.sz);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  
  
  // set up the data structure
  data.data = (void *)theID.data;
  data.size = theID.sz * sizeof(int);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  if ((ret = dbID->put(dbID, 0, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::sendID() - failed to send the ID to database\n";
    dbID->err(dbID, ret, "DB->put");
    return -2;      
  
  }
  opserr << dbTag << " " << commitTag << " " << theID.Size() << " " << query << endln;

  return 0;
}
 

int 
BerkeleyDbDatastore::recvID(int dbTag, int commitTag, 
			    ID &theID,
			    ChannelAddress *theAddress)
{
  // check that we have a connection
  if (connection == false)
    return -1;

  // set up the key structure
  sprintf(query,"%d-%d-%d", dbTag, commitTag, theID.sz);
  key.data = query;
  key.size = strlen(query);
  key.ulen = key.size;
  key.doff = 0;
  key.flags = DB_DBT_USERMEM;  

  opserr << dbTag << " " << commitTag << " " << theID.Size() << " " << query << endln;


  
  // set up the data structure
  data.data = (void *)theID.data;
  data.size = theID.sz * sizeof(int);
  data.ulen = data.size;
  data.doff = 0;
  data.flags = DB_DBT_USERMEM;

  int ret;
  if ((ret = dbID->get(dbID, NULL, &key, &data, 0)) != 0) {
    opserr << "BerkeleyDbDatastore::recvID() - failed to get the ID from database\n";
    dbID->err(dbID, ret, "DB->get");
    return -2;      
  }

  return 0;
}
