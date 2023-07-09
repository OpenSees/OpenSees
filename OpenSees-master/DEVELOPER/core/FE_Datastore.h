/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
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
                                                                        
// $Revision: 1.8 $
// $Date: 2005-11-08 00:04:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/FE_Datastore.h,v $
                                                                        
                                                                        
#ifndef FE_Datastore_h
#define FE_Datastore_h

// File: ~/database/FE_Datastore.h
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class definition for FE_Datastore.
// FE_Datastore is an abstract base class. An FE_datastore object is used
// in the program to store/restore the geometry and state information in
// a domain at a particular instance in the analysis.
//
// What: "@(#) FE_Datastore.h, revA"

#include <Channel.h>

class Domain;
class FEM_ObjectBroker;

class FE_Datastore: public Channel
{
  public:
    FE_Datastore(Domain &theDomain, FEM_ObjectBroker &theBroker);    
    virtual ~FE_Datastore();

    // methods defined in the Channel class interface which mean nothing for a Datastore
    char *addToProgram(void);
    int setUpConnection(void);
    int setNextAddress(const ChannelAddress &otherChannelAddress);
    ChannelAddress *getLastSendersAddress(void);

    // methods defined in the Channel class interface which can be implemented here    
    int sendObj(int commitTag,
		MovableObject &theObject, 
		ChannelAddress *theAddress =0);

    int recvObj(int commitTag,
		MovableObject &theObject, 
		FEM_ObjectBroker &theBroker,
		ChannelAddress *theAddress =0);
		
    // pure virtual functions in addition to those defined
    // in the ModelBuilder and Channel classes for database applications
	int getDbTag(void);
    virtual int isDatastore(void);

    virtual int commitState(int commitTag);    
    virtual int restoreState(int commitTag);        

    virtual int createTable(const char *tableName, int numColumns, char *columns[]);
    virtual int insertData(const char *tableName, char *columns[], 
			   int commitTag, const Vector &data);
    virtual int getData(const char *tableName, char *columns[], 
			int commitTag, Vector &data);
			
  protected:
    FEM_ObjectBroker *getObjectBroker(void);
    
  private:
    FEM_ObjectBroker *theObjectBroker;
    Domain *theDomain;
    static int lastDbTag;

};


#endif

