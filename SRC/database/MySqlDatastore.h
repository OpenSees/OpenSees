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
                                                                        
// $Revision: 1.5 $
// $Date: 2005-11-07 21:34:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/MySqlDatastore.h,v $

#ifndef MySqlDatastore_h
#define MySqlDatastore_h

// Written:  fmk
// Created:  02/02
//
// Description: This file contains the class definition for MySqlDatastore.
// MySqlDatastore is a concrete subclas of FE_Datastore. A MySqlDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored in a MySQL Database which can be local or remote.
//
// What: "@(#) MySqlDatastore.h, revA"

#include <FE_Datastore.h>
#include <mysql.h>

class MySqlDatastore: public FE_Datastore
{
 public:
  MySqlDatastore(const char *databaseName,
		 Domain &theDomain, 
		 FEM_ObjectBroker &theBroker,
		 int dbRun = 0);    

  MySqlDatastore(const char *databaseName,
		 const char *host,
		 const char *user,
		 const char *passwd,
		 unsigned int port,
		 const char *socket,
		 unsigned int clientFlag,
		 Domain &theDomain, 
		 FEM_ObjectBroker &theBroker,
		 int dbRun = 0);    
  
  ~MySqlDatastore();

  // methods for sending and recieving matrices, vectors and id's
  int sendMsg(int dbTag, int commitTag, 
	      const Message &, 
	      ChannelAddress *theAddress =0);    
  int recvMsg(int dbTag, int commitTag, 
	      Message &, 
	      ChannelAddress *theAddress =0);        

  int sendMatrix(int dbTag, int commitTag, 
		 const Matrix &theMatrix, 
		 ChannelAddress *theAddress =0);
  int recvMatrix(int dbTag, int commitTag, 
		 Matrix &theMatrix, 
		 ChannelAddress *theAddress =0);

  int sendVector(int dbTag, int commitTag, 
		 const Vector &theVector, 
		 ChannelAddress *theAddress =0);
  int recvVector(int dbTag, int commitTag, 
		 Vector &theVector, 
		 ChannelAddress *theAddress =0);
    
  int sendID(int dbTag, int commitTag, 
	     const ID &theID, 
	     ChannelAddress *theAddress =0);
  int recvID(int dbTag, int commitTag, 
	     ID &theID, 
	     ChannelAddress *theAddress =0);    

  int createTable(const char *tableName, int numColumns, char *columns[]);
  int insertData(const char *tableName, char *columns[], int commitTag, const Vector &data);
  int getData(const char *tableName, char *columns[], int commitTag, Vector &data);

  int setDbRun(int run);
  int getDbRun(void);

 protected:
  int createOpenSeesDatabase(const char *projectName);
  int execute(const char *query);

 private:
  int dbTag;
  int dbRun;
  
  bool connection;
  MYSQL mysql;

  char *query;
  int sizeQuery;
  int sizeColumnString;
};

#endif
