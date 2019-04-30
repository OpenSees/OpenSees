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
                                                                        
// $Revision: 1.31 $
// $Date: 2010-02-04 00:39:05 $
// $URL: /usr/local/cvs/OpenSees/SRC/domain/node/Node.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// This file contains the implementation of the Node class
//
// What: "@(#) Node.h, revA"
   
#include <Node.h>
#include <stdlib.h>

#include <Element.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <DOF_Group.h>
#include <Renderer.h>
#include <string.h>
#include <Information.h>
#include <Parameter.h>

// AddingSensitivity:BEGIN //////////////////////////
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
// AddingSensitivity:END ////////////////////////////

#include <NodalLoad.h> 
//Added by Liming Jiang for link NodalLoadPtr, [SIF]

#include <OPS_Globals.h>
#include <elementAPI.h>

Matrix **Node::theMatrices = 0;
int Node::numMatrices = 0;

int OPS_Node()
{
    Domain* theDomain = OPS_GetDomain();
    int ndm = OPS_GetNDM();
    int ndf = OPS_GetNDF();
    
    if(theDomain == 0) {
	opserr<<"WARNING: domain is not defined\n";
	return -1;
    }
    if(ndm<=0 || ndf<=0) {
	opserr<<"WARNING: system ndm and ndf are zero\n";
	return -1;
    }

    if(OPS_GetNumRemainingInputArgs() < 1+ndm) {
	opserr<<"insufficient number of arguments\n";
	return -1;
    }

    // get tag
    int tag = 0;
    int numData = 1;
    if(OPS_GetIntInput(&numData, &tag) < 0) {
	opserr<<"WARNING tag is not integer\n";
	return -1;
    }

    // get crds
    Vector crds(ndm);
    if(OPS_GetDoubleInput(&ndm, &crds(0)) < 0) {
	opserr<<"WARNING noda coord are not double\n";
	return -1;
    }

    // check options
    Vector disp,vel,mass,dispLoc;
    Matrix ndmass;
    while(OPS_GetNumRemainingInputArgs() > 0) {
	const char* type = OPS_GetString();
	
	if(strcmp(type,"-disp")==0 || strcmp(type,"-Disp")==0) {
	    if(OPS_GetNumRemainingInputArgs() < ndf) {
		opserr<<"incorrect number of nodal disp terms\n";
		return -1;
	    }
	    disp.resize(ndf);
	    if(OPS_GetDoubleInput(&ndf, &disp(0)) < 0) {
		opserr << "WARNING: failed to read disp\n";
		return -1;
	    }

	} else if(strcmp(type,"-vel")==0 || strcmp(type,"-Vel")==0) {
	    if(OPS_GetNumRemainingInputArgs() < ndf) {
		opserr<<"incorrect number of nodal vel terms\n";
		return -1;
	    }
	    vel.resize(ndf);
	    if(OPS_GetDoubleInput(&ndf, &vel(0)) < 0) {
		opserr << "WARNING: failed to read vel\n";
		return -1;
	    }

	} else if(strcmp(type,"-mass")==0 || strcmp(type,"-Mass")==0) {
	    if(OPS_GetNumRemainingInputArgs() < ndf) {
		opserr<<"incorrect number of nodal mass terms\n";
		return -1;
	    }
	    Vector data(ndf);
	    if(OPS_GetDoubleInput(&ndf, &data(0)) < 0) {
		opserr << "WARNING: failed to read mass\n";
		return -1;
	    }
	    ndmass.resize(ndf,ndf);
	    ndmass.Zero();
	    for(int i=0; i<ndf; i++) {
		ndmass(i,i) = data(i);
	    }

	} else if(strcmp(type,"-dispLoc")==0 || strcmp(type,"-dispLoc")==0) {
	    if(OPS_GetNumRemainingInputArgs() < ndm) {
		opserr<<"incorrect number of nodal mass terms\n";
		return -1;
	    }
	    dispLoc.resize(ndm);
	    if(OPS_GetDoubleInput(&ndm, &dispLoc(0)) < 0) {
		opserr << "WARNING: failed to read dispLoc\n";
		return -1;
	    }

	} else if(strcmp(type,"-ndf")==0 || strcmp(type,"-NDF")==0) {
	    if(OPS_GetNumRemainingInputArgs() < 1) {
		opserr<<"incorrect number for ndf\n";
		return -1;
	    }
	    int numdata = 1;
	    if(OPS_GetIntInput(&numdata, &ndf) < 0) {
		opserr << "WARNING: failed to read ndf\n";
		return -1;
	    }

	}
    }
    
    // create node
    Node* theNode = 0;
    if(ndm == 1) {
	theNode = new Node(tag,ndf,crds(0));
    } else if(ndm == 2) {
	theNode = new Node(tag,ndf,crds(0),crds(1));
    } else {
	theNode = new Node(tag,ndf,crds(0),crds(1),crds(2));
    }

    if(theNode == 0) {
	opserr<<"run out of memory for node "<<tag<<"\n";
	return -1;
    }

    // set node
    if(disp.Size() == ndf) {
	theNode->setTrialDisp(disp);
    }
    if(vel.Size() == ndf) {
	theNode->setTrialVel(vel);
    }
    if(ndmass.noRows() == ndf) {
	theNode->setMass(ndmass);
    }
    if(dispLoc.Size() == ndm) {
	theNode->setDisplayCrds(dispLoc);
    }
    theNode->commitState();

    // add node to domain
    if(theDomain->addNode(theNode) == false) {
	opserr<<"WARNING: failed to add node to domain\n";
	delete theNode;
	return -1;
    }

    return 0;
}

// for FEM_Object Broker to use
Node::Node(int theClassTag)
:DomainComponent(0,theClassTag), 
 numberDOF(0), theDOF_GroupPtr(0), 
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0), 
 incrDeltaDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0), 
 index(-1), reaction(0), displayLocation(0)
{
  // for FEM_ObjectBroker, recvSelf() must be invoked on object

  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]
}    


Node::Node(int tag, int theClassTag)
:DomainComponent(tag,theClassTag), 
 numberDOF(0), theDOF_GroupPtr(0), 
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0), 
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0), 
 index(-1), reaction(0), displayLocation(0)
{
  // for subclasses - they must implement all the methods with
  // their own data structures.
  
  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]
}

Node::Node(int tag, int ndof, double Crd1, Vector *dLoc)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0), 
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0), 
 index(-1), reaction(0), displayLocation(0)
{
  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]
  
  Crd = new Vector(1);
  (*Crd)(0) = Crd1;

  if (dLoc != 0) {
    displayLocation = new Vector(*dLoc);
  }
  
  index = -1;
}


//  Node(int tag, int ndof, double Crd1, double yCrd);
//	constructor for 2d nodes
Node::Node(int tag, int ndof, double Crd1, double Crd2, Vector *dLoc)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0), 
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 reaction(0), displayLocation(0)
{
  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(2);
  (*Crd)(0) = Crd1;
  (*Crd)(1) = Crd2;

  if (dLoc != 0) {
    displayLocation = new Vector(*dLoc);
  }
  
  index = -1;
}


//  Node(int tag, int ndof, double Crd1, double Crd2, double zCrd);
//	constructor for 3d nodes

 Node::Node(int tag, int ndof, double Crd1, double Crd2, double Crd3, Vector *dLoc)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0), 
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
 reaction(0), displayLocation(0)
{
  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]
  
  Crd = new Vector(3);
  (*Crd)(0) = Crd1;
  (*Crd)(1) = Crd2;
  (*Crd)(2) = Crd3;    

  if (dLoc != 0) {
    displayLocation = new Vector(*dLoc);
  }
  
  index = -1;
}


// used for domain decomposition & external nodes
//  copy everything but the mass 
//  we should really set the mass to 0.0
Node::Node(const Node &otherNode, bool copyMass)
  :DomainComponent(otherNode.getTag(),otherNode.getClassTag()), 
 numberDOF(otherNode.numberDOF), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 incrDeltaDisp(0), 
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
 R(0), mass(0), unbalLoadWithInertia(0), alphaM(0.0), theEigenvectors(0),
   reaction(0), displayLocation(0)
{
  // AddingSensitivity:BEGIN /////////////////////////////////////////
  dispSensitivity = 0;
  velSensitivity = 0;
  accSensitivity = 0;
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////////

  theNodalThermalActionPtr = 0;//Added by Liming for initializing NodalLoadPointer, [SIF]

  Crd = new Vector(otherNode.getCrds());
  if (Crd == 0) {
    opserr << " FATAL Node::Node(node *) - ran out of memory for Crd\n";
    exit(-1);
  }

  if (otherNode.displayLocation != 0) {
    displayLocation = new Vector(*(otherNode.displayLocation));
  }

  if (otherNode.commitDisp != 0) {
    if (this->createDisp() < 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for displacement\n";
      exit(-1);
    }
    for (int i=0; i<4*numberDOF; i++)
      disp[i] = otherNode.disp[i];
  }    
  
  if (otherNode.commitVel != 0) {
    if (this->createVel() < 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for velocity\n";
      exit(-1);
    }
    for (int i=0; i<2*numberDOF; i++)
      vel[i] = otherNode.vel[i];
  }    
  
  if (otherNode.commitAccel != 0) {
    if (this->createAccel() < 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for acceleration\n";
      exit(-1);
    }
    for (int i=0; i<2*numberDOF; i++)
      accel[i] = otherNode.accel[i];
  }    
  
  
  if (otherNode.unbalLoad != 0){
    unbalLoad = new Vector(*(otherNode.unbalLoad));    
    if (unbalLoad == 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for Load\n";
      exit(-1);
    }
    unbalLoad->Zero();
  }    

  if (otherNode.mass != 0 && copyMass == true) {
    mass = new Matrix(*(otherNode.mass)) ;
    if (mass == 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for mass\n";
      exit(-1);
    }
  }

  if (otherNode.R != 0) {
    R = new Matrix(*(otherNode.R));
    if (R == 0) {
      opserr << " FATAL Node::Node(node *) - ran out of memory for R\n";
      exit(-1);
    }
  }

  index = -1;
}


// ~Node():
// 	destructor

Node::~Node()
{
    // delete anything that we created with new
    if (Crd != 0)
	delete Crd;

    if (commitDisp != 0)
	delete commitDisp;

    if (commitVel != 0)
	delete commitVel;

    if (commitAccel != 0)
	delete commitAccel;

    if (trialDisp != 0)
	delete trialDisp;

    if (trialVel != 0)
	delete trialVel;

    if (trialAccel != 0)
	delete trialAccel;

    if (incrDisp != 0)
	delete incrDisp;
    
    if (incrDeltaDisp != 0)
	delete incrDeltaDisp;    
    
    if (unbalLoad != 0)
	delete unbalLoad;
    
    if (disp != 0)
	delete [] disp;

    if (vel != 0)
	delete [] vel;

    if (accel != 0)
	delete [] accel;

    if (mass != 0)
	delete mass;
    
    if (R != 0)
	delete R;

    if (unbalLoadWithInertia != 0)
      delete unbalLoadWithInertia;

    if (theEigenvectors != 0)
      delete theEigenvectors;

    // AddingSensitivity:BEGIN ///////////////////////////////////////
    if (dispSensitivity != 0)
      delete dispSensitivity;
    if (velSensitivity != 0)
	delete velSensitivity;
    if (accSensitivity != 0)
      delete accSensitivity;
    // AddingSensitivity:END /////////////////////////////////////////

    if (reaction != 0)
      delete reaction;


    if (displayLocation != 0)
      delete displayLocation;

    if (theDOF_GroupPtr != 0)
      theDOF_GroupPtr->resetNodePtr();
}


int
Node::getNumberDOF(void) const
{
    // return the number of dof
    return  numberDOF;
}


void
Node::setDOF_GroupPtr(DOF_Group *theDOF_Grp)
{
    // set the DOF_Group pointer
    theDOF_GroupPtr = theDOF_Grp;
}


DOF_Group *
Node::getDOF_GroupPtr(void)
{
    // return the DOF_Group pointer
    return theDOF_GroupPtr;
}


const Vector &
Node::getCrds() const
{
    // return the vector of nodal coordinates
    return *Crd;
}




const Vector &
Node::getDisp(void) 
{
    // construct memory and Vectors for trial and committed
    // displacement on first call to this method, getTrialDisp()
    // setTrialDisp() or incrTrialDisp()
    if (commitDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::getDisp() -- ran out of memory\n";
	    exit(-1);
	}
    }
    
    // return the committed disp
    return *commitDisp;
}

const Vector &
Node::getVel(void) 
{
    // construct memory and Vectors for trial and committed
    // velocity on first call to this method, getTrialVel()
    // setTrialVel() or incrTrialVel()    
    if (commitVel == 0) {
	if (this->createVel() < 0) {
	    opserr << "FATAL Node::getVel() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    // return the velocity
    return *commitVel;    
}


const Vector &
Node::getAccel(void) 
{
    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialAccel()
    // setTrialAccel() or incrTrialAccel()        
    if (commitAccel == 0) {
	if (this->createAccel() < 0) {
	    opserr << "FATAL Node::getAccel() -- ran out of memory\n";
	    exit(-1);
	}
    }

    return *commitAccel;    
}


/* *********************************************************************
**
**   Methods to return the trial response quantities similar to committed
**
** *********************************************************************/

const Vector &
Node::getTrialDisp(void) 
{
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *trialDisp;
}



const Vector &
Node::getTrialVel(void) 
{
    if (trialVel == 0) {
	if (this->createVel() < 0) {
	    opserr << "FATAL Node::getTrialVel() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *trialVel;
}



const Vector &
Node::getTrialAccel(void) 
{
    if (trialAccel == 0) {
	if (this->createAccel() < 0) {
	    opserr << "FATAL Node::getTrialAccel() - ran out of memory\n";
	    exit(0);
	}
    }    
    return *trialAccel;
}

const Vector &
Node::getIncrDisp(void) 
{
    if (incrDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *incrDisp;
}

const Vector &
Node::getIncrDeltaDisp(void) 
{
    if (incrDeltaDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *incrDeltaDisp;
}


int
Node::setTrialDisp(double value, int dof)
{
    // check vector arg is of correct size
    if (dof < 0 || dof >=  numberDOF) {
      opserr << "WARNING Node::setTrialDisp() - incompatable sizes\n";
      opserr << "node: " << this->getTag() << endln;
      return -2;
    }    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialDisp(),
    // getDisp(), or incrTrialDisp()        
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::setTrialDisp() - ran out of memory\n";
	    exit(-1);
	}    
    }

    // perform the assignment .. we dont't go through Vector interface
    // as we are sure of size and this way is quicker
    double tDisp = value;
    disp[dof+2*numberDOF] = tDisp - disp[dof+numberDOF];
    disp[dof+3*numberDOF] = tDisp - disp[dof];	
    disp[dof] = tDisp;

    return 0;
}

int
Node::setTrialDisp(const Vector &newTrialDisp)
{
    // check vector arg is of correct size
    if (newTrialDisp.Size() != numberDOF) {
      opserr << "WARNING Node::setTrialDisp() - incompatable sizes\n";
      opserr << "node: " << this->getTag() << endln;
      return -2;
    }    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialDisp(),
    // getDisp(), or incrTrialDisp()        
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::setTrialDisp() - ran out of memory\n";
	    exit(-1);
	}    
    }

    // perform the assignment .. we dont't go through Vector interface
    // as we are sure of size and this way is quicker
    for (int i=0; i<numberDOF; i++) {
        double tDisp = newTrialDisp(i);
	disp[i+2*numberDOF] = tDisp - disp[i+numberDOF];
	disp[i+3*numberDOF] = tDisp - disp[i];	
	disp[i] = tDisp;
    }

    return 0;
}

int
Node::setTrialVel(const Vector &newTrialVel)
{
    // check vector arg is of correct size
    if (newTrialVel.Size() != numberDOF) {
	    opserr << "WARNING Node::setTrialVel() - incompatable sizes\n";
	    return -2;
    }    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialVEl(),
    // getVEl(), or incrTrialVel()        
    if (trialVel == 0) {
	if (this->createVel() < 0) {
	    opserr << "FATAL Node::setTrialVel() - ran out of memory\n";
	    exit(-1);
	}
    }      
    
    // set the trial quantities
    for (int i=0; i<numberDOF; i++)
	vel[i] = newTrialVel(i);
    return 0;
}


int
Node::setTrialAccel(const Vector &newTrialAccel)
{
    // check vector arg is of correct size
    if (newTrialAccel.Size() != numberDOF) {
	    opserr << "WARNING Node::setTrialAccel() - incompatable sizes\n";
	    return -2;
    }    

    // create a copy if no trial exists    
    if (trialAccel == 0) {
	if (this->createAccel() < 0) {
	    opserr << "FATAL Node::setTrialAccel() - ran out of memory\n";
	    exit(-1);
	}
    }        
    
    // use vector assignment otherwise        
    for (int i=0; i<numberDOF; i++)
	accel[i] = newTrialAccel(i);

    return 0;
}

int
Node::incrTrialDisp(const Vector &incrDispl)
{
    // check vector arg is of correct size
    if (incrDispl.Size() != numberDOF) {
	opserr << "WARNING Node::incrTrialDisp() - incompatable sizes\n";
	return -2;
    }    

    // create a copy if no trial exists andd add committed
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    opserr << "FATAL Node::incrTrialDisp() - ran out of memory\n";
	    exit(-1);
	}    
	for (int i = 0; i<numberDOF; i++) {
	  double incrDispI = incrDispl(i);
	  disp[i] = incrDispI;
	  disp[i+2*numberDOF] = incrDispI;
	  disp[i+3*numberDOF] = incrDispI;
	}
	return 0;
    }

    // otherwise set trial = incr + trial
    for (int i = 0; i<numberDOF; i++) {
	  double incrDispI = incrDispl(i);
	  disp[i] += incrDispI;
	  disp[i+2*numberDOF] += incrDispI;
	  disp[i+3*numberDOF] = incrDispI;
    }

    return 0;
}


int
Node::incrTrialVel(const Vector &incrVel)
{
    // check vector arg is of correct size
    if (incrVel.Size() != numberDOF) {
	opserr << "WARNING Node::incrTrialVel() - incompatable sizes\n";
	return -2;
    }    

    // create Vectors and array if none exist and set trial
    if (trialVel == 0) {
	if (this->createVel() < 0) {
	    opserr << "FATAL Node::incrTrialVel - ran out of memory\n";
	    exit(-1);
	}    
	for (int i = 0; i<numberDOF; i++)
	    vel[i] = incrVel(i);

	return 0;
    }

    // otherwise set trial = incr + trial
    for (int i = 0; i<numberDOF; i++)
	vel[i] += incrVel(i);    

    return 0;
}


int
Node::incrTrialAccel(const Vector &incrAccel)
{
    // check vector arg is of correct size
    if (incrAccel.Size() != numberDOF) {
	opserr << "WARNING Node::incrTrialAccel() - incompatable sizes\n";
	return -2;
    }    

    // create a copy if no trial exists andd add committed    
    if (trialAccel == 0) {
	if (this->createAccel() < 0) {
	    opserr << "FATAL Node::incrTrialAccel() - ran out of memory\n";
	    exit(-1);
	}    
	for (int i = 0; i<numberDOF; i++)
	    accel[i] = incrAccel(i);

	return 0;
    }

    // otherwise set trial = incr + trial
    for (int i = 0; i<numberDOF; i++)
	accel[i] += incrAccel(i);    

    return 0;
}


void 
Node::zeroUnbalancedLoad(void)
{
    if (unbalLoad != 0)
	unbalLoad->Zero();
}

int
Node::addUnbalancedLoad(const Vector &add, double fact)
{
    // check vector arg is of correct size
    if (add.Size() != numberDOF) {
	opserr << "Node::addunbalLoad - load to add of incorrect size ";
	opserr << add.Size() << " should be " <<  numberDOF << endln;	
	return -1;
    }

    // if no load yet create it and assign
    if (unbalLoad == 0) {
	unbalLoad = new Vector(add); 
	if (unbalLoad == 0) {
	    opserr << "FATAL Node::addunbalLoad - ran out of memory\n";
	    exit(-1);
	}
	if (fact != 1.0)
	    (*unbalLoad) *= fact;
	return 0;
    }

    // add fact*add to the unbalanced load
    unbalLoad->addVector(1.0, add,fact);

    return 0;
}



int
Node::addInertiaLoadToUnbalance(const Vector &accelG, double fact)
{
  // simply return if node has no mass or R matrix
  if (mass == 0 || R == 0)
    return 0;

  // otherwise we must determine MR accelG
  if (accelG.Size() != R->noCols()) {
    opserr << "Node::addInertiaLoadToUnbalance - accelG not of correct dimension";
    return -1;
  }

  // if no load yet create it and assign
  if (unbalLoad == 0) {
      unbalLoad = new Vector(numberDOF); 
      if (unbalLoad == 0 || unbalLoad->Size() != numberDOF) {
	  opserr << "FATAL Node::addunbalLoad - ran out of memory\n";
	  exit(-1);
      }  
  }

  // form - fact * M*R*accelG and add it to the unbalanced load
  //(*unbalLoad) -= ((*mass) * (*R) * accelG)*fact;

  Matrix MR(mass->noRows(), R->noCols());
  MR.addMatrixProduct(0.0, *mass, *R, 1.0);
  unbalLoad->addMatrixVector(1.0, MR, accelG, -fact);

  return 0;
}



int
Node::addInertiaLoadSensitivityToUnbalance(const Vector &accelG, double fact, bool somethingRandomInMotions)
{
  // simply return if node has no mass or R matrix
  if (mass == 0 || R == 0)
    return 0;

  // otherwise we must determine MR accelG
  if (accelG.Size() != R->noCols()) {
    opserr << "Node::addInertiaLoadToUnbalance - accelG not of correct dimension";
    return -1;
  }

  // if no load yet create it and assign
  if (unbalLoad == 0) {
      unbalLoad = new Vector(numberDOF); 
      if (unbalLoad == 0 || unbalLoad->Size() != numberDOF) {
	  opserr << "FATAL Node::addunbalLoad - ran out of memory\n";
	  exit(-1);
      }  
  }

  // form - fact * M*R*accelG and add it to the unbalanced load
  //(*unbalLoad) -= ((*mass) * (*R) * accelG)*fact;


	Matrix massSens(mass->noRows(),mass->noCols());
	massSens = this->getMassSensitivity();

	Matrix MR(mass->noRows(), R->noCols());

	if (somethingRandomInMotions) {
	  MR.addMatrixProduct(0.0, *mass, *R, 1.0);
	}
	else {
	  MR.addMatrixProduct(0.0, massSens, *R, 1.0);
	}
	unbalLoad->addMatrixVector(1.0, MR, accelG, -fact);

  return 0;
}



const Vector &
Node::getUnbalancedLoad(void) 
{
    // make sure it was created before we return it
    if (unbalLoad == 0) {
	unbalLoad = new Vector(numberDOF);
	if (unbalLoad == 0 || unbalLoad->Size() != numberDOF) {
	    opserr << "FATAL Node::getunbalLoad() -- ran out of memory\n";
	    exit(-1);
	}
    }

    // return the unbalanced load

    return *unbalLoad;
}


    
const Vector &
Node::getUnbalancedLoadIncInertia(void) 
{
    // make sure it was created before we return it
    if (unbalLoadWithInertia == 0) {
	unbalLoadWithInertia = new Vector(this->getUnbalancedLoad());
	if (unbalLoadWithInertia == 0) {
	    opserr << "FATAL Node::getunbalLoadWithInertia -- ran out of memory\n";
	    exit(-1);
	}
    } else
      (*unbalLoadWithInertia) = this->getUnbalancedLoad();

    if (mass != 0) {

      const Vector &theAccel = this->getTrialAccel(); // in case accel not created
      unbalLoadWithInertia->addMatrixVector(1.0, *mass, theAccel, -1.0);

      if (alphaM != 0.0) {
	const Vector &theVel = this->getTrialVel(); // in case vel not created
	unbalLoadWithInertia->addMatrixVector(1.0, *mass, theVel, -alphaM);
      }
    } 

    return *unbalLoadWithInertia;
}



int
Node::commitState()
{
    // check disp exists, if does set commit = trial, incr = 0.0
    if (trialDisp != 0) {
      for (int i=0; i<numberDOF; i++) {
	disp[i+numberDOF] = disp[i];  
        disp[i+2*numberDOF] = 0.0;
        disp[i+3*numberDOF] = 0.0;
      }
    }		    
    
    // check vel exists, if does set commit = trial    
    if (trialVel != 0) {
      for (int i=0; i<numberDOF; i++)
	vel[i+numberDOF] = vel[i];
    }
    
    // check accel exists, if does set commit = trial        
    if (trialAccel != 0) {
      for (int i=0; i<numberDOF; i++)
	accel[i+numberDOF] = accel[i];
    }

    // if we get here we are done
    return 0;
}



int
Node::revertToLastCommit()
{
    // check disp exists, if does set trial = last commit, incr = 0
    if (disp != 0) {
      for (int i=0 ; i<numberDOF; i++) {
	disp[i] = disp[i+numberDOF];
	disp[i+2*numberDOF] = 0.0;
	disp[i+3*numberDOF] = 0.0;
      }
    }
    
    // check vel exists, if does set trial = last commit
    if (vel != 0) {
      for (int i=0 ; i<numberDOF; i++)
	vel[i] = vel[numberDOF+i];
    }

    // check accel exists, if does set trial = last commit
    if (accel != 0) {    
      for (int i=0 ; i<numberDOF; i++)
	accel[i] = accel[numberDOF+i];
    }

    // if we get here we are done
    return 0;
}


int
Node::revertToStart()
{
    // check disp exists, if does set all to zero
    if (disp != 0) {
      for (int i=0 ; i<4*numberDOF; i++)
	disp[i] = 0.0;
    }

    // check vel exists, if does set all to zero
    if (vel != 0) {
      for (int i=0 ; i<2*numberDOF; i++)
	vel[i] = 0.0;
    }

    // check accel exists, if does set all to zero
    if (accel != 0) {    
      for (int i=0 ; i<2*numberDOF; i++)
	accel[i] = 0.0;
    }
    
    if (unbalLoad != 0) 
	(*unbalLoad) *= 0;




// AddingSensitivity: BEGIN /////////////////////////////////
	if (dispSensitivity != 0) 
		dispSensitivity->Zero();
	
	if (velSensitivity != 0) 
		velSensitivity->Zero();
	
	if (accSensitivity != 0) 
		accSensitivity->Zero();
// AddingSensitivity: END ///////////////////////////////////



    // if we get here we are done
    return 0;
}


const Matrix &
Node::getMass(void) 
{
    if (index == -1) {
	setGlobalMatrices();
    }
    
    // make sure it was created before we return it
    if (mass == 0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else 
      return *mass;
}


int 
Node::setRayleighDampingFactor(double alpham) {
  alphaM = alpham;
  return 0;
}


const Matrix &
Node::getDamp(void) 
{
    if (index == -1) {
	setGlobalMatrices();
    }
    
    // make sure it was created before we return it
    if (mass == 0 || alphaM == 0.0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else {
      Matrix &result = *theMatrices[index];
      result = *mass;
      result *= alphaM;
      return result;
    } 
}


const Matrix &
Node::getDampSensitivity(void) 
{
    if (index == -1) {
	setGlobalMatrices();
    }
    
    // make sure it was created before we return it
    if (mass == 0 || alphaM == 0.0) {
      theMatrices[index]->Zero();
      return *theMatrices[index];
    } else {
      Matrix &result = *theMatrices[index];
	  result.Zero();
      //result = *mass;
      //result *= alphaM;
      return result;
    } 
}


int
Node::setMass(const Matrix &newMass)
{
    // check right size
    if (newMass.noRows() != numberDOF || newMass.noCols() != numberDOF) {
	opserr << "Node::setMass - incompatable matrices\n";
	return -1;
    }	

    // create a matrix if no mass yet set
    if (mass == 0) {
	mass = new Matrix(newMass);
	if (mass == 0 || mass->noRows() != numberDOF) {
	    opserr << "FATAL Node::setMass - ran out of memory\n";
	    return -1;
	}
	return 0;
    }

    // assign if mass has already been created
    (*mass) = newMass;
    
    return 0;
}



int 
Node::setNumColR(int numCol)
{
  if (R != 0) {
    if (R->noCols() != numCol) {
      delete R;
      R = new Matrix(numberDOF, numCol);
    }
  } else 
    R = new Matrix(numberDOF, numCol);

  if (R == 0 || R->noRows() != numberDOF) {
      opserr << "FATAL Node::setNumColR() - out of memory\n";
      exit(-1);
  }

  R->Zero();
  return 0;
}  

int 
Node::setR(int row, int col, double Value)
{
  // ensure R had been set
  if (R == 0) {
    opserr << "Node:setR() - R has not been initialised\n";
    return -1;
  }
  
  // ensure row, col in range (matrix assignment will catch this - extra work)
  if (row < 0 || row > numberDOF || col < 0 || col > R->noCols()) {
    opserr << "Node:setR() - row, col index out of range\n";
    return -1;
  }
  
  // do the assignment
  (*R)(row,col) = Value;
  
  /*
  // to test uniform excitation pattern with consistent mass matrices:
  // found that the static application of a unit ground displacement
  // needs to also be applied to the constrained DOFs
  Domain *theDomain = this->getDomain();
  SP_ConstraintIter &theSPs = theDomain->getSPs();
  SP_Constraint *theSP;
  // assign zero if there is a homogeneous SP
  while ((theSP = theSPs()) != 0) {
      if (theSP->getNodeTag() == this->getTag() &&
          theSP->getDOF_Number() == row &&
          theSP->isHomogeneous()) {
              (*R)(row,col) = 0.0;
      }
  }
  */
  
  return 0;
}



const Vector &
Node::getRV(const Vector &V)
{
    // we store the product of RV in unbalLoadWithInertia

    // make sure unbalLoadWithInertia was created, if not create it
    if (unbalLoadWithInertia == 0) {
	unbalLoadWithInertia = new Vector(numberDOF);
	if (unbalLoadWithInertia == 0) {
	    opserr << "Node::getunbalLoadWithInertia -- ran out of memory\n";
	    exit(-1);
	}
    } 
    
    // see if quick return , i.e. R == 0
    if (R == 0) {
	unbalLoadWithInertia->Zero();
	return *unbalLoadWithInertia;
    }
    
    // check dimesions of R and V
    if (R->noCols() != V.Size()) {
	opserr << "WARNING Node::getRV() - R and V of incompatable dimesions\n";
	opserr << "R: " << *R << "V: " << V;
	unbalLoadWithInertia->Zero();
	return *unbalLoadWithInertia;
    }    

    // determine the product
    unbalLoadWithInertia->addMatrixVector(0.0, *R, V, 1.0);
    return *unbalLoadWithInertia;
}


int 
Node::setNumEigenvectors(int numVectorsToStore)
{
  // ensure a positive number of vectors
  if (numVectorsToStore <= 0) {
    opserr << "Node::setNumEigenvectors() - " << numVectorsToStore << " < 0\n";
    return -1;
  }    

  // if matrix not yet assigned or not of correct size delete old and create new
  if (theEigenvectors == 0 || theEigenvectors->noCols() != numVectorsToStore) {
    if (theEigenvectors != 0)
      delete theEigenvectors;


    theEigenvectors = new Matrix(numberDOF, numVectorsToStore);
    if (theEigenvectors == 0 || theEigenvectors->noCols() != numVectorsToStore) {
      opserr << "Node::setNumEigenvectors() - out of memory\n";
      return -2;
    }
  } else
    // zero the eigenvector matrix
    theEigenvectors->Zero();
  
    return 0;
}
int 
Node::setEigenvector(int mode, const Vector &eigenVector)
{
  if (theEigenvectors == 0 || theEigenvectors->noCols() < mode) {
    opserr << "Node::setEigenvectors() - mode " << mode << " invalid\n";
    return -1;
  }

  if (eigenVector.Size() != numberDOF) {
      opserr << "Node::setEigenvectors() - eigenvector of incorrect size\n";
      return -2;
  }
  // set the values
  for (int i=0; i<numberDOF; i++)
    (*theEigenvectors)(i, mode-1) = eigenVector(i);

  return 0;
}
const Matrix &
Node::getEigenvectors(void)
{
  // check the eigen vectors have been set
	if (theEigenvectors == 0) {
    opserr << "Node::getEigenvectors() - eigenvectors have not been set\n";
	exit(0);
	}
  
  return *theEigenvectors;
}


int 
Node::sendSelf(int cTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();

    ID data(14);
    data(0) = this->getTag(); 
    data(1) = numberDOF; 
    
    // indicate whether vector quantaties have been formed
    if (disp == 0)       data(2) = 1; else data(2) = 0;
    if (vel == 0)        data(3) = 1; else data(3) = 0;
    if (accel == 0)      data(4) = 1; else data(4) = 0;
    if (mass == 0)       data(5) = 1; else data(5) = 0;
    if (unbalLoad  == 0) data(6) = 1; else data(6) = 0;    
    if (R == 0) 	 
	data(12) = 1; 
    else {
	data(12) = 0;        
	data(13) = R->noCols();
    }

    data(7) = Crd->Size();

    if (dbTag1 == 0)
      dbTag1 = theChannel.getDbTag();
    if (dbTag2 == 0)
      dbTag2 = theChannel.getDbTag();
    if (dbTag3 == 0)
      dbTag3 = theChannel.getDbTag();
    if (dbTag4 == 0)
      dbTag4 = theChannel.getDbTag();

    data(8) = dbTag1;
    data(9) = dbTag2;
    data(10) = dbTag3;
    data(11) = dbTag4;

    int res = 0;

    res = theChannel.sendID(dataTag, cTag, data);
    if (res < 0) {
      opserr << " Node::sendSelf() - failed to send ID data\n";
      return res;
    }

    res = theChannel.sendVector(dataTag, cTag, *Crd);
    if (res < 0) {
      opserr << " Node::sendSelf() - failed to send Vecor data\n";
      return res;
    }

    if (commitDisp != 0) {
	res = theChannel.sendVector(dbTag1, cTag, *commitDisp);	
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send Disp data\n";
	  return res;
	}
    }

    if (commitVel != 0) {
	res = theChannel.sendVector(dbTag2, cTag, *commitVel);		
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send Vel data\n";
	  return res;
	}
    }

    if (commitAccel != 0) {
	res = theChannel.sendVector(dbTag3, cTag, *commitAccel); 
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send Accel data\n";
	  return res;
	}
    }

    if (mass != 0) {
	res = theChannel.sendMatrix(dataTag, cTag, *mass);
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send Mass data\n";
	  return res;
	}
    }
    
    if (R != 0) {
	res = theChannel.sendMatrix(dataTag, cTag, *R);
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send R data\n";
	  return res;
	}
    }    

    if (unbalLoad  != 0) {
	res = theChannel.sendVector(dbTag4, cTag, *unbalLoad);	
	if (res < 0) {
	  opserr << " Node::sendSelf() - failed to send Load data\n";
	  return res;
	}
    }

    // if get here succesfull
    return 0;
}

int 
Node::recvSelf(int cTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker)
{
    int res = 0;
    int dataTag = this->getDbTag();

    
    ID data(14);
    res = theChannel.recvID(dataTag, cTag, data);
    if (res < 0) {
      opserr << "Node::recvSelf() - failed to receive ID data\n";
      return res;
    }

    this->setTag(data(0));
    numberDOF = data(1);
    int numberCrd = data(7);

    dbTag1 = data(8);
    dbTag2 = data(9);
    dbTag3 = data(10);
    dbTag4 = data(11);
    

    // create a Vector to hold coordinates IF one needed
    if (Crd == 0) {
      Crd = new Vector(numberCrd);
    }

    // check we did not run out of memory
    if (Crd == 0) {
      opserr << "Node::recvSelf() - out of memory creating Coordinate vector\n";
      return -1;
    }

    if (theChannel.recvVector(dataTag, cTag, *Crd) < 0) {
      opserr << "Node::recvSelf() - failed to receive the Coordinate vector\n";
      return -2;
    }

    if (data(2) == 0) {
      // create the disp vectors if node is a total blank
      if (commitDisp == 0)
	this->createDisp();

      // recv the committed disp
      if (theChannel.recvVector(dbTag1, cTag, *commitDisp) < 0) {
	opserr << "Node::recvSelf - failed to receive Disp data\n";
	return res;
      }

      // set the trial quantities equal to committed
      for (int i=0; i<numberDOF; i++)
	disp[i] = disp[i+numberDOF];  // set trial equal commited

    } else if (commitDisp != 0) {
      // if going back to initial we will just zero the vectors
      commitDisp->Zero();
      trialDisp->Zero();
    }
      
    
    if (data(3) == 0) {
      // create the vel vectors if node is a total blank
      if (commitVel == 0) 
	this->createVel();

      // recv the committed vel
      if (theChannel.recvVector(dbTag2, cTag, *commitVel) < 0) {
	opserr << "Node::recvSelf - failed to receive Velocity data\n";
	return -3;
      }

      // set the trial quantity
      for (int i=0; i<numberDOF; i++)
	vel[i] = vel[i+numberDOF];  // set trial equal commited
    }

    if (data(4) == 0) {
      // create the vel vectors if node is a total blank
      if (commitAccel == 0) 
	this->createAccel();

      // recv the committed accel
      if (theChannel.recvVector(dbTag3, cTag, *commitAccel) < 0) {
	opserr << "Node::recvSelf - failed to receive Acceleration data\n";
	return -4;
      }
      
      // set the trial values
      for (int i=0; i<numberDOF; i++)
	accel[i] = accel[i+numberDOF];  // set trial equal commited
    }

    if (data(5) == 0) {
      // make some room and read in the vector
      if (mass == 0) {
	mass = new Matrix(numberDOF,numberDOF);
	if (mass == 0) {
	  opserr << "Node::recvData -- ran out of memory\n";
	  return -5;
	}
      }
      if (theChannel.recvMatrix(dataTag, cTag, *mass) < 0) {
	opserr << "Node::recvSelf() - failed to receive Mass data\n";
	return -6;
      }
    }            
    
    if (data(12) == 0) {
      // create a matrix for R
      int noCols = data(13);
      if (R == 0) {
	R = new Matrix(numberDOF, noCols);
	if (R == 0) {
	  opserr << "Node::recvData -- ran out of memory\n";
	  return -1;
	}
      }
      // now recv the R matrix
      if (theChannel.recvMatrix(dataTag, cTag, *R) < 0) {
	opserr << "Node::recvSelf() - failed to receive R data\n";
	return res;
      }
    }


    if (data(6) == 0) {
      // create a vector for the load
      if (unbalLoad == 0) {
	unbalLoad = new Vector(numberDOF);
	if (unbalLoad == 0) {
	  opserr << "Node::recvData -- ran out of memory\n";
	  return -10;
	}
      }
      if (theChannel.recvVector(dbTag4, cTag, *unbalLoad) < 0) {
	opserr << "Node::recvSelf() - failed to receive Load data\n";
	return res;
      }
    }        


  index = -1;
  if (numMatrices != 0) {
    for (int i=0; i<numMatrices; i++)
      if (theMatrices[i]->noRows() == numberDOF) {
	index = i;
	i = numMatrices;
      }
  }
  if (index == -1) {
    Matrix **nextMatrices = new Matrix *[numMatrices+1];
    if (nextMatrices == 0) {
      opserr << "Element::getTheMatrix - out of memory\n";
      exit(-1);
    }
    for (int j=0; j<numMatrices; j++)
      nextMatrices[j] = theMatrices[j];
    Matrix *theMatrix = new Matrix(numberDOF, numberDOF);
    if (theMatrix == 0) {
      opserr << "Element::getTheMatrix - out of memory\n";
      exit(-1);
    }
    nextMatrices[numMatrices] = theMatrix;
    if (numMatrices != 0) 
      delete [] theMatrices;
    index = numMatrices;
    numMatrices++;
    theMatrices = nextMatrices;
  }

  return 0;
}



void
Node::Print(OPS_Stream &s, int flag)
{
    if (flag == OPS_PRINT_CURRENTSTATE) { // print out everything
        s << "\n Node: " << this->getTag() << endln;
        s << "\tCoordinates  : " << *Crd;
        if (commitDisp != 0)
            s << "\tDisps: " << *trialDisp;
        if (commitVel != 0)
            s << "\tVelocities   : " << *trialVel;
        if (commitAccel != 0)
            s << "\tcommitAccels: " << *trialAccel;
        if (unbalLoad != 0)
            s << "\t unbalanced Load: " << *unbalLoad;
        if (reaction != 0)
            s << "\t reaction: " << *reaction;
        if (mass != 0) {
            s << "\tMass : " << *mass;
            s << "\t Rayleigh Factor: alphaM: " << alphaM << endln;
            s << "\t Rayleigh Forces: " << *this->getResponse(RayleighForces);
        }
        if (theEigenvectors != 0)
            s << "\t Eigenvectors: " << *theEigenvectors;
        if (theDOF_GroupPtr != 0)
            s << "\tID : " << theDOF_GroupPtr->getID();
        s << "\n";
    }
    
    else if (flag == 1) { // print out: nodeId displacements
        s << this->getTag() << "  " << *commitDisp;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"ndf\": " << numberDOF << ", ";
        s << "\"crd\": [";
        int numCrd = Crd->Size();
        for (int i = 0; i < numCrd - 1; i++)
            s << (*Crd)(i) << ", ";
        s << (*Crd)(numCrd - 1) << "]";
        if (mass != 0) {
            s << ", \"mass\": [";
            for (int i = 0; i < numberDOF - 1; i++)
                s << (*mass)(i, i) << ", ";
            s << (*mass)(numberDOF - 1, numberDOF - 1) << "]";
        }
        s << "}";
    }
}
  
int
Node::displaySelf(Renderer &theRenderer, int displayMode, float fact)
{

  if (displayMode == 0)
    return 0;

//  const Vector &theDisp = this->getDisp();
  static Vector position(3);

  this->getDisplayCrds(position, fact, displayMode);
  
  if (displayMode == -1) { 
    // draw a text string containing tag
    static char theText[20];
    sprintf(theText,"%d",this->getTag());
    return theRenderer.drawText(position, theText, (int) strlen(theText));

  } else if (displayMode > 0) {
    // draw a point - pixel size equals displayMode tag
    return theRenderer.drawPoint(position, 0.0, this->getTag(), 0, displayMode);
  }


  return 0;
}


// createDisp(), createVel() and createAccel():
// private methods to create the arrays to hold the disp, vel and acceleration
// values and the Vector objects for the committed and trial quantaties.

int
Node::createDisp(void)
{
  // trial , committed, incr = (committed-trial)
  disp = new double[4*numberDOF];
    
  if (disp == 0) {
    opserr << "WARNING - Node::createDisp() ran out of memory for array of size " << 2*numberDOF << endln;
			    
    return -1;
  }
  for (int i=0; i<4*numberDOF; i++)
    disp[i] = 0.0;
    
  commitDisp = new Vector(&disp[numberDOF], numberDOF); 
  trialDisp = new Vector(disp, numberDOF);
  incrDisp = new Vector(&disp[2*numberDOF], numberDOF);
  incrDeltaDisp = new Vector(&disp[3*numberDOF], numberDOF);
  
  if (commitDisp == 0 || trialDisp == 0 || incrDisp == 0 || incrDeltaDisp == 0) {
    opserr << "WARNING - Node::createDisp() " <<
      "ran out of memory creating Vectors(double *,int)";
    return -2;
  }
    
  return 0;
}


int
Node::createVel(void)
{
    vel = new double[2*numberDOF];
    
    if (vel == 0) {
      opserr << "WARNING - Node::createVel() ran out of memory for array of size " << 2*numberDOF << endln;
      return -1;
    }
    for (int i=0; i<2*numberDOF; i++)
      vel[i] = 0.0;
    
    commitVel = new Vector(&vel[numberDOF], numberDOF); 
    trialVel = new Vector(vel, numberDOF);
    
    if (commitVel == 0 || trialVel == 0) {
      opserr << "WARNING - Node::createVel() %s" <<
	"ran out of memory creating Vectors(double *,int) \n";
      return -2;
    }
    
    return 0;
}

int
Node::createAccel(void)
{
    accel = new double[2*numberDOF];
    
    if (accel == 0) {
      opserr << "WARNING - Node::createAccel() ran out of memory for array of size " << 2*numberDOF << endln;
      return -1;
    }
    for (int i=0; i<2*numberDOF; i++)
	accel[i] = 0.0;
    
    commitAccel = new Vector(&accel[numberDOF], numberDOF);
    trialAccel = new Vector(accel, numberDOF);
    
    if (commitAccel == 0 || trialAccel == 0) {
      opserr << "WARNING - Node::createAccel() ran out of memory creating Vectors(double *,int)\n";
      return -2;
    }

    return 0;
}


// AddingSensitivity:BEGIN ///////////////////////////////////////

Matrix
Node::getMassSensitivity(void)
{
    if (index == -1) {
	setGlobalMatrices();
    }
    
	if (mass == 0) {
		theMatrices[index]->Zero();
		return *theMatrices[index];
	} 
	else {
		Matrix massSens(mass->noRows(),mass->noCols());
		if ( (parameterID == 1) || (parameterID == 2) || (parameterID == 3) ) {
			massSens(parameterID-1,parameterID-1) = 1.0;
		}
		if (parameterID == 7) {
		  massSens(0,0) = 1.0;
		  massSens(1,1) = 1.0;
		}
		if (parameterID == 8) {
		  massSens(0,0) = 1.0;
		  massSens(1,1) = 1.0;
		  massSens(2,2) = 1.0;
		}
		return massSens;
	}
}


int
Node::getCrdsSensitivity(void)
{
	if ( (parameterID == 4) || (parameterID == 5) || (parameterID == 6) ) {
		return (parameterID-3);
	}
	else {
		return 0;
	}
}


int
Node::setParameter(const char **argv, int argc, Parameter &param)
{
  // The following parameterID map is being used:
  // 1: nodal mass in direction 1	
  // 2: nodal mass in direction 2
  // 3: nodal mass in direction 3
  // 4: coordinate in direction 1
  // 5: coordinate in direction 2
  // 6: coordinate in direction 3

  if (argc < 2)
    return -1;

  if ((strstr(argv[0],"mass") != 0) || (strstr(argv[0],"-mass") != 0)) { 
    int direction = 0; // atoi(argv[1]);
    if ((strcmp(argv[1],"x") == 0)||(strcmp(argv[1],"X") == 0)||(strcmp(argv[1],"1") == 0)) {
      direction = 1;
      if (mass != 0)
	param.setValue((*mass)(0,0));
    }
    else if ((strcmp(argv[1],"y") == 0)||(strcmp(argv[1],"Y") == 0)||(strcmp(argv[1],"2") == 0)) {
      direction = 2;
      if (mass != 0)
	param.setValue((*mass)(1,1));
    }
    else if ((strcmp(argv[1],"z") == 0)||(strcmp(argv[1],"Z") == 0)||(strcmp(argv[1],"3") == 0)) {
      direction = 3;
      if (mass != 0)
	param.setValue((*mass)(2,2));
    }
    else if ((strcmp(argv[1],"xy") == 0)||(strcmp(argv[1],"XY") == 0)) {
      direction = 7;
      if (mass != 0)
	param.setValue((*mass)(0,0));
    }
    else if ((strcmp(argv[1],"xyz") == 0)||(strcmp(argv[1],"XYZ") == 0)) {
      direction = 8;
      if (mass != 0)
	param.setValue((*mass)(0,0));
    }
    
    if ((direction >= 1 && direction <= 3) || direction == 7 || direction == 8)
      return param.addObject(direction, this);
  }
  else if (strstr(argv[0],"coord") != 0) {
    int direction = atoi(argv[1]);
    if (direction >= 1 && direction <= 3) {
      if (Crd != 0)
	param.setValue((*Crd)(direction-1));
      return param.addObject(direction+3, this);
    }
  }
  else
    opserr << "WARNING: Could not set parameter in Node. " << endln;
  
  return -1;
}



int
Node::updateParameter(int pparameterID, Information &info)
{
  if (pparameterID >= 1 && pparameterID <= 3)
    (*mass)(pparameterID-1,pparameterID-1) = info.theDouble;

  else if (pparameterID == 7) {
    (*mass)(0,0) = info.theDouble;
    (*mass)(1,1) = info.theDouble;
  } else if (pparameterID == 8) {
    (*mass)(0,0) = info.theDouble;
    (*mass)(1,1) = info.theDouble;
    (*mass)(2,2) = info.theDouble;
  }

  else if (pparameterID >= 4 && pparameterID <= 6) {

    if ( (*Crd)(pparameterID-4) != info.theDouble) {

      // Set the new coordinate value
      (*Crd)(pparameterID-4) = info.theDouble;
      
      // Need to "setDomain" to make the change take effect. 
      Domain *theDomain = this->getDomain();
      ElementIter &theElements = theDomain->getElements();
      Element *theElement;
      while ((theElement = theElements()) != 0) {
	theElement->setDomain(theDomain);
      }
    }
    else {
      // No change in nodal coordinate
    }
  }
  
  return -1;
}




int
Node::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

int 
Node::saveDispSensitivity(const Vector &v, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (dispSensitivity == 0) {
    dispSensitivity = new Matrix( numberDOF, numGrads );
  } 

  if (dispSensitivity->noRows() != numberDOF ||
      dispSensitivity->noCols() != numGrads) {
    delete dispSensitivity;
    dispSensitivity = new Matrix( numberDOF, numGrads );
  }

  //opserr << "Node::saveDispSens " << dispSensitivity->noRows() << ' ' << dispSensitivity->noCols() << endln;
  for (int i=0; i<numberDOF; i++ )
    (*dispSensitivity)(i,gradIndex) = v(i);

  return 0;
}

int 
Node::saveVelSensitivity(const Vector &vdot, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (velSensitivity == 0) {
    velSensitivity = new Matrix( numberDOF, numGrads );
  } 

  for (int i=0; i<numberDOF; i++ )
    (*velSensitivity)(i,gradIndex) = vdot(i);

  return 0;
}

int 
Node::saveAccelSensitivity(const Vector &vdotdot, int gradIndex, int numGrads)
{
  // If the sensitivity matrices are not already created:
  if (accSensitivity == 0) {
    accSensitivity = new Matrix( numberDOF, numGrads );
  } 

  for (int i=0; i<numberDOF; i++ )
    (*accSensitivity)(i,gradIndex) = vdotdot(i);

  return 0;
}

double 
Node::getDispSensitivity(int dof, int gradIndex)
{
  if (dispSensitivity != 0)
    return (*dispSensitivity)(dof-1,gradIndex);
  else 
    return 0.0;
}

double 
Node::getVelSensitivity(int dof, int gradIndex)
{
  if (velSensitivity != 0)
    return (*velSensitivity)(dof-1,gradIndex);
  else
    return 0.0;
}

double 
Node::getAccSensitivity(int dof, int gradIndex)
{
  if (accSensitivity != 0)
    return (*accSensitivity)(dof-1,gradIndex);
  else
    return 0.0;
}
// AddingSensitivity:END /////////////////////////////////////////



const Vector &
Node::getReaction() {
  if (reaction == 0) {
    reaction = new Vector(numberDOF);
    if (reaction == 0) {
      opserr << "FATAL Node::getReaction() - out of memory\n";
      exit(-1);
    }
  }

  return *reaction;
}

int   
Node::addReactionForce(const Vector &add, double factor){

  // create rection vector if have not done so already
  if (reaction == 0) {
    reaction = new Vector(numberDOF);
    if (reaction == 0) {
      opserr << "WARNING Node::addReactionForce() - out of memory\n";
      return -1;
    }
  }

  // check vector of appropraie size
  if (add.Size() != numberDOF) {
    opserr << "WARNING Node::addReactionForce() - vector not of correct size\n";
    return -1;
  }

  if (factor == 1.0) 
    *reaction += add;
  else if (factor == -1.0)
    *reaction -= add;
  else
    *reaction = add * factor;

  return 0;
}

int   
Node::resetReactionForce(int flag){

  // create rection vector if have not done so already
  if (reaction == 0) {
    reaction = new Vector(numberDOF);
    if (reaction == 0) {
      opserr << "WARNING Node::addReactionForce() - out of memory\n";
      return -1;
    }
  }

  reaction->Zero();

  // add unbalance, the negative of applied forces hence the -=
  if (flag == 0) {
    *reaction -= this->getUnbalancedLoad();
  } if (flag == 1) {
    *reaction -= this->getUnbalancedLoadIncInertia();
  } else {
    if (mass != 0 && alphaM != 0) {
      if (alphaM != 0.0) {
	const Vector &theVel = this->getTrialVel(); // in case vel not created
	reaction->addMatrixVector(1.0, *mass, theVel, alphaM);
      }
    } 
  }
  return 0;
}

const Vector *
Node::getResponse(NodeResponseType responseType)
{
  const Vector *result = NULL;
  if (responseType == Disp) 
    result  = &(this->getDisp());
  else if (responseType == Vel) 
    return &(this->getVel());
  else if (responseType == Accel) 
    return &(this->getAccel());
  else if (responseType == IncrDisp) 
    return &(this->getIncrDisp());
  else if (responseType == IncrDeltaDisp) 
    return &(this->getIncrDeltaDisp());
  else if (responseType == Reaction) 
    return &(this->getReaction());
  else if (responseType == Unbalance) 
    return &(this->getUnbalancedLoad());
  else if (responseType == RayleighForces) {
    if (unbalLoadWithInertia == 0) {
      unbalLoadWithInertia = new Vector(this->getUnbalancedLoad());
    }
    if (alphaM != 0.0 && mass != 0) {
      const Vector &theVel = this->getTrialVel(); // in case vel not created
      unbalLoadWithInertia->addMatrixVector(0.0, *mass, theVel, -alphaM);
    } else
      unbalLoadWithInertia->Zero();

    return unbalLoadWithInertia;
  } else
    return NULL;

  return result;
}

void
Node::setCrds(double Crd1)
{
  if (Crd != 0 && Crd->Size() >= 1)
    (*Crd)(0) = Crd1;

  // Need to "setDomain" to make the change take effect. 
  Domain *theDomain = this->getDomain();
  ElementIter &theElements = theDomain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) {
    theElement->setDomain(theDomain);
  }
}

void
Node::setCrds(double Crd1, double Crd2)
{
  if (Crd != 0 && Crd->Size() >= 2) {
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;

    // Need to "setDomain" to make the change take effect. 
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      theElement->setDomain(theDomain);
    }
  }
}

void
Node::setCrds(double Crd1, double Crd2, double Crd3)
{
  if (Crd != 0 && Crd->Size() >= 3) {
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
    (*Crd)(2) = Crd3;

    // Need to "setDomain" to make the change take effect. 
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      theElement->setDomain(theDomain);
    }
  }
}

void
Node::setCrds(const Vector &newCrds) 
{
  if (Crd != 0 && Crd->Size() == newCrds.Size()) {
    (*Crd) = newCrds;

	return;

    // Need to "setDomain" to make the change take effect. 
    Domain *theDomain = this->getDomain();
    ElementIter &theElements = theDomain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      theElement->setDomain(theDomain);
    }
  }
}

int
Node::getDisplayCrds(Vector &res, double fact, int mode) 
{
  int ndm = Crd->Size();
  int resSize = res.Size();

  if (resSize < ndm)
    return -1;

  if (mode < 0) {
    int eigenMode = -mode;
    if ((theEigenvectors != 0) && ((*theEigenvectors).noCols() > eigenMode)) {
      if (displayLocation != 0)
	for (int i=0; i<ndm; i++)
	  res(i) = (*displayLocation)(i)+(*theEigenvectors)(i,eigenMode-1)*fact;
      else
	for (int i=0; i<ndm; i++)
	  res(i) = (*Crd)(i)+(*theEigenvectors)(i,eigenMode-1)*fact;
    }
  } else {    
  
    if (commitDisp != 0) {
      if (displayLocation != 0)
	for (int i=0; i<ndm; i++)
	  res(i) = (*displayLocation)(i)+(*commitDisp)(i)*fact;
      else
      for (int i=0; i<ndm; i++)
	res(i) = (*Crd)(i)+(*commitDisp)(i)*fact;
    } else {
      if (displayLocation != 0)
	for (int i=0; i<ndm; i++)
	  res(i) = (*displayLocation)(i);
      else
	for (int i=0; i<ndm; i++)
	  res(i) = (*Crd)(i);
    }
    
  }

  // zero rest
  for (int i=ndm; i<resSize; i++)
    res(i) = 0;

  return 0;
}

int
Node::setDisplayCrds(const Vector &theCrds) 
{
  if (theCrds.Size() != Crd->Size()) {
    return -1;
  }

  if (displayLocation == 0) {
    displayLocation = new Vector(theCrds);
  } else {
    *displayLocation = theCrds;
  }
  return 0;
}


//Add Pointer to NodalThermalAction id applicable------begin-----L.Jiang, [SIF]
NodalThermalAction*
Node::getNodalThermalActionPtr(void)
{
	return theNodalThermalActionPtr;
}
void
Node::setNodalThermalActionPtr(NodalThermalAction* theAction)
{
	theNodalThermalActionPtr = theAction;
}
//Add Pointer to NodalThermalAction id applicable-----end------L.Jiang, {SIF]

int
Node::setGlobalMatrices()
{
    if (index == -1) {
	for (int i=0; i<numMatrices; i++) {
	    if (theMatrices[i]->noRows() == numberDOF) {
		index = i;
		i = numMatrices;
	    }
	}
    }
    if (index == -1) {
	Matrix **nextMatrices = new Matrix *[numMatrices+1];
	if (nextMatrices == 0) {
	    opserr << "Element::getTheMatrix - out of memory\n";
	    exit(-1);
	}
	for (int j=0; j<numMatrices; j++)
	    nextMatrices[j] = theMatrices[j];
	Matrix *theMatrix = new Matrix(numberDOF, numberDOF);
	if (theMatrix == 0) {
	    opserr << "Element::getTheMatrix - out of memory\n";
	    exit(-1);
	}
	nextMatrices[numMatrices] = theMatrix;
	if (numMatrices != 0) 
	    delete [] theMatrices;
	index = numMatrices;
	numMatrices++;
	theMatrices = nextMatrices;
    }

    return 0;
}
