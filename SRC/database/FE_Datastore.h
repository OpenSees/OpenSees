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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:17 $
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

    // pure virtual functions in addition to those defined
    // in the ModelBuilder and Channel classes

    virtual int getDbTag(void) =0;
    virtual int commitState(int commitTag);    
    virtual int restoreState(int commitTag);        
    
  protected:
    FEM_ObjectBroker *getObjectBroker(void);
    
  private:
    FEM_ObjectBroker *theObjectBroker;
    Domain *theDomain;
};


#endif

