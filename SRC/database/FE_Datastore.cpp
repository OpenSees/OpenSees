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
// $Source: /usr/local/cvs/OpenSees/SRC/database/FE_Datastore.cpp,v $
                                                                        
                                                                        
// File: ~/database/FE_Datastore.C
//
// Written: fmk 
// Created: 10/98
// Revision: A
//
// Description: This file contains the class implementation for FE_Datastore.
// FE_Datastore is an abstract base class. An FE_Datastore object is used
// in the program to store/restore the geometry and state information in
// a domain at a particular instance in the analysis.
//
// What: "@(#) FE_Datastore.C, revA"

#include "FE_Datastore.h"
#include <FEM_ObjectBroker.h>
#include <Domain.h>
#include <G3Globals.h>



// FE_Datastore(int tag, int noExtNodes);
// 	constructor that takes the FE_Datastore's unique tag and the number
//	of external nodes for the FE_Datastore.

FE_Datastore::FE_Datastore(Domain &thDomain, FEM_ObjectBroker &theBroker) 
  :theObjectBroker(&theBroker), theDomain(&thDomain)
{

}


FE_Datastore::~FE_Datastore() 
{
    // does nothing
}




int
FE_Datastore::commitState(int commitTag)
{
  // invoke sendSelf on the domain object with this as an arg
  int res = 0;
  if (theDomain != 0) {
    res = theDomain->sendSelf(commitTag, *this);
    if (res < 0) {
      g3ErrorHandler->warning("FE_Datastore::commitState - domain failed to sendSelf\n");
    }
  }
    
  return res;
}



int
FE_Datastore::restoreState(int commitTag)
{
  // invoke sendSelf on the domain object with this as an arg
  int res = 0;
  if (theDomain != 0) {
    res = theDomain->recvSelf(commitTag, *this, *theObjectBroker);
    if (res < 0) {
      g3ErrorHandler->warning("FE_Datastore::restoreState - domain failed to recvSelf\n");
    }
  }
    
  return res;
}




