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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-15 16:49:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/BerkeleyDbDatastore.h,v $

#ifndef BerkeleyDbDatastore_h
#define BerkeleyDbDatastore_h

// Written:  fmk
// Created:  02/02
//
// Description: This file contains the class definition for BerkeleyDbDatastore.
// BerkeleyDbDatastore is a concrete subclas of FE_Datastore. A BerkeleyDbDatastore 
// object is used in the program to store/restore the geometry and state 
// information in a domain at a particular instance in the analysis. The
// information is stored using the BerkeleyDB Database.
//
// What: "@(#) BerkeleyDbDatastore.h, revA"

#include <FE_Datastore.h>
#include <db.h>

class BerkeleyDbDatastore: public FE_Datastore
{
 public:
  BerkeleyDbDatastore(const char *projectName, 
		      Domain &theDomain, 
		      FEM_ObjectBroker &theBroker,
		      char *type = NULL);

  ~BerkeleyDbDatastore();

  // method to get a database tag
  int getDbTag(void);
  
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

 protected:

 private:
  int dbTag;

  DB_ENV *dbenv;      // database enviroment - directory containing files
  DB *dbMatrix;       // file holding the data for the Matrices
  DB *dbVector;       // file holding the data for the Vectors
  DB *dbID;           // file holding the data for the ID's

  bool connection;

  DBT key, data;  // structures for the key data pair
  char query[128];

  char *project;
};

#endif
