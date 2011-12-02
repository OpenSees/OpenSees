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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-15 00:38:07 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ActorTruss.cpp,v $
                                                                        
// Written: fmk 
// Created: 08/03
//
// Description: This file contains the implementation for the ActorTruss class.
//
// What: "@(#) ActorTruss.C, revA"


// we specify what header files we need
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include <OPS_Globals.h>

#include "ActorTruss.h"
#include "ShadowActorTruss.h"

#include <Information.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <UniaxialMaterial.h>


// initialise the class wide variables
Matrix ActorTruss::trussK(4,4);
Vector ActorTruss::trussR(4);

// typical constructor
ActorTruss::ActorTruss(Channel &theChannel, 
		       FEM_ObjectBroker &theObjectBroker)
:Actor(theChannel, theObjectBroker, 0),     
 msgData(2),
 trans(1, 4), L(0.0), A(0.0), theMaterial(0)
{	

}



//  destructor - provided to clean up any memory
ActorTruss::~ActorTruss()
{
  if (theMaterial != 0)
    delete theMaterial;
}

int
ActorTruss::run(void)
{
  bool exitYet = false;
  while (exitYet == false) {

    this->recvID(msgData);
    int action = msgData(0);

    
    switch (action) {
      
    case ShadowActorTruss_setMaterial:
      this->setMaterial();
      break;

    case ShadowActorTruss_setDomain:
      this->setDomain();
      break;

    case ShadowActorTruss_commitState:
      this->commitState();
      break;

    case ShadowActorTruss_revertToLastCommit:
      this->revertToLastCommit();
      break;

    case ShadowActorTruss_revertToStart:
      this->revertToStart();
      break;

    case ShadowActorTruss_update:
      this->update();
      break;

    case ShadowActorTruss_getTangentStiff:
      this->getTangentStiff();
      break;

    case ShadowActorTruss_getInitialStiff:
      this->getInitialStiff();
      break;

    case ShadowActorTruss_getResistingForce:
      this->getResistingForce();
      break;


    case ShadowActorTruss_DIE:
      exitYet = true;
      break;

    default:
      opserr << "ActorTruss::invalid action " << action << "received\n";
      msgData(0) = -1;
    }
  }
  return 0;
}

int
ActorTruss::setMaterial(void) {
  // recv area & material from shadow
  static Vector data(1);
  this->recvVector(data);
  A = data(0);

  theMaterial = theBroker->getNewUniaxialMaterial(msgData(1));
  this->recvObject(*theMaterial);
  return 0;
}


int
ActorTruss::setDomain(void)
{
    // recv coord data from shadow
    static Vector coords(4);
    this->recvVector(coords);

    // determin transformation & length
    double dx = coords(2)-coords(0);
    double dy = coords(3)-coords(1);
    
    L = sqrt(dx*dx + dy*dy);
    
    if (L == 0.0) {
      opserr << "WARNING ActorTruss::setDomain() - ActorTruss has zero length\n";
      return -1;  // don't go any further - otherwise divide by 0 error
    }
	
    double cs = dx/L;
    double sn = dy/L;

    trans(0,0) = -cs;
    trans(0,1) = -sn;    
    trans(0,2) = cs;
    trans(0,3) = sn;
    return 0;
}   	 


int
ActorTruss::commitState()
{
    return theMaterial->commitState();
}

int
ActorTruss::revertToLastCommit()
{
  return theMaterial->revertToLastCommit();
}

int
ActorTruss::revertToStart()
{
  return theMaterial->revertToStart();
}

int
ActorTruss::update()
{
  // determine the current strain given trial displacements at nodes
  static Vector data(1);
  this->recvVector(data);
  double strain = data(0);

  // set the strain in the materials
  theMaterial->setTrialStrain(strain);
  return 0;
}

int
ActorTruss::getTangentStiff(void)
{
  if (L == 0.0) { // length = zero - problem in setDomain() warning message already printed
    trussK.Zero();

    return this->sendMatrix(trussK);
  }
  
  // get the current E from the material for the last updated strain
  double E = theMaterial->getTangent();
  
  // form the tangent stiffness matrix
  trussK = trans^trans;
  trussK *= A*E/L;  
  
  // return the matrix
  return this->sendMatrix(trussK);
}

int
ActorTruss::getInitialStiff(void)

{  
  if (L == 0.0) { // length = zero - problem in setDomain() warning message already printed
    trussK.Zero();
    return this->sendMatrix(trussK);
  }
  
  // get the current E from the material for the last updated strain
  double E = theMaterial->getInitialTangent();
  
  // form the tangent stiffness matrix
  trussK = trans^trans;
  trussK *= A*E/L;  
  
  // return the matrix
  return this->sendMatrix(trussK);
}



int
ActorTruss::getResistingForce()
{	
  if (L == 0.0) { // if length == 0, problem in setDomain()
    trussR.Zero();
    return this->sendVector(trussR);
  }

  // want: R = Ku - Pext
  
  // force = F * transformation 
  double force = A*theMaterial->getStress();
  for (int i=0; i<4; i++)
    trussR(i) = trans(0,i)*force;
  
  return this->sendVector(trussR);
}



