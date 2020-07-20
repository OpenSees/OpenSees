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
                                                                        
// $Revision: 1.3 $
// $Date: 2005-11-30 23:47:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/subdomain/ActorSubdomain.h,v $
                                                                        
                                                                        
#ifndef ActorSubdomain_h
#define ActorSubdomain_h

// Written: fmk 
// Revision: A
//
// Description: This file contains the class definition for ActorSubdomain.
// ActorSubdomain is a container class. The class is responsible for holding
// and providing access to the Elements, Nodes, LoadCases, SP_Constraints 
// and MP_Constraints that have been added to the ActorSubdomain.
//
// What: "@(#) ActorSubdomain.h, revA"

#include "Subdomain.h"
#include <Actor.h>

class ActorSubdomain: public Subdomain, public Actor
{
  public:
    ActorSubdomain(Channel &theChannel,
		   FEM_ObjectBroker &theBroker);
    
    virtual ~ActorSubdomain();
    
    virtual int run(void);
    virtual const Vector &getLastExternalSysResponse(void);

    virtual int  update(void);    
    virtual int  updateTimeDt(void);    
    virtual int  barrierCheck(int res);    

  protected:

    
  private:
    ID msgData;
    Vector *lastResponse;
};
	
		   
#endif
