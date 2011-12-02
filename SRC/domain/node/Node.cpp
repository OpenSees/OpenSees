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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/Node.cpp,v $
                                                                        
                                                                        
// File: ~/domain/node/Node.C
//
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
#include <G3Globals.h>

void NodeItoa(int x, char *str);

// for FEM_Object Broker to use
Node::Node(int theClassTag)
:DomainComponent(0,theClassTag), 
 numberDOF(0), theDOF_GroupPtr(0), 
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    // for FEM_ObjectBroker, recvSelf() must be invoked on object
}    


Node::Node(int tag, int theClassTag)
:DomainComponent(tag,theClassTag), 
 numberDOF(0), theDOF_GroupPtr(0), 
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    // for subclasses - they must implement all the methods with
    // their own data structures.
}

Node::Node(int tag, int ndof, double Crd1)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    Crd = new Vector(1);
    (*Crd)(0) = Crd1;
}


//  Node(int tag, int ndof, double Crd1, double yCrd);
//	constructor for 2d nodes
Node::Node(int tag, int ndof, double Crd1, double Crd2)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    Crd = new Vector(2);
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
}


//  Node(int tag, int ndof, double Crd1, double Crd2, double zCrd);
//	constructor for 3d nodes

Node::Node(int tag, int ndof, double Crd1, double Crd2, double Crd3)
:DomainComponent(tag,NOD_TAG_Node), 
 numberDOF(ndof), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    Crd = new Vector(3);
    (*Crd)(0) = Crd1;
    (*Crd)(1) = Crd2;
    (*Crd)(2) = Crd3;    
}



Node::Node(const Node *otherNode)
:DomainComponent(otherNode->getTag(),NOD_TAG_Node), 
 numberDOF(otherNode->numberDOF), theDOF_GroupPtr(0),
 Crd(0), commitDisp(0), commitVel(0), commitAccel(0), 
 trialDisp(0), trialVel(0), trialAccel(0), unbalLoad(0), incrDisp(0),
 disp(0), vel(0), accel(0), dbTag1(0), dbTag2(0), dbTag3(0), dbTag4(0),
  R(0), mass(0), unbalLoadWithInertia(0), theEigenvectors(0)
{
    if (otherNode->Crd != 0) {
	Crd = new Vector(*(otherNode->Crd));
	if (Crd == 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for Crd\n";
	    exit(-1);
	}
    }

    if (otherNode->commitDisp != 0) {
	if (this->createDisp() < 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for displacement\n";
	    exit(-1);
	}
	for (int i=0; i<3*numberDOF; i++)
	    disp[i] = otherNode->disp[i];
    }    
    
    if (otherNode->commitVel != 0) {
	if (this->createVel() < 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for velocity\n";
	    exit(-1);
	}
	for (int i=0; i<2*numberDOF; i++)
	    vel[i] = otherNode->vel[i];
    }    
    
    if (otherNode->commitAccel != 0) {
	if (this->createAccel() < 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for acceleration\n";
	    exit(-1);
	}
	for (int i=0; i<2*numberDOF; i++)
	    accel[i] = otherNode->accel[i];
    }    

    if (otherNode->unbalLoad != 0){
	unbalLoad = new Vector(*(otherNode->unbalLoad));    
	if (unbalLoad == 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for Load\n";
	    exit(-1);
	}
    }    

    if (otherNode->mass != 0) {
	mass = new Matrix(*(otherNode->mass)) ;
	if (mass == 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for mass\n";
	    exit(-1);
	}
    }
    if (otherNode->R != 0) {
      R = new Matrix(*(otherNode->R));
	if (R == 0) {
	    cerr << " FATAL Node::Node(node *) - ran out of memory for R\n";
	    exit(-1);
	}
    }
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
	    cerr << "FATAL Node::getDisp() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getVel() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getAccel() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getTrialVel() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getTrialAccel() - ran out of memory\n";
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
	    cerr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getTrialDisp() -- ran out of memory\n";
	    exit(-1);
	}
    }    

    return *incrDeltaDisp;
}


int
Node::setTrialDisp(const Vector &newTrialDisp)
{
    // check vector arg is of correct size
    if (newTrialDisp.Size() != numberDOF) {
	    cerr << "WARNING Node::setTrialDisp() - incompatable sizes\n";
	    return -2;
	}    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialDisp(),
    // getDisp(), or incrTrialDisp()        
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    cerr << "FATAL Node::setTrialDisp() - ran out of memory\n";
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
	    cerr << "WARNING Node::setTrialVel() - incompatable sizes\n";
	    return -2;
    }    

    // construct memory and Vectors for trial and committed
    // accel on first call to this method, getTrialVEl(),
    // getVEl(), or incrTrialVel()        
    if (trialVel == 0) {
	if (this->createVel() < 0) {
	    cerr << "FATAL Node::setTrialVel() - ran out of memory\n";
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
	    cerr << "WARNING Node::setTrialAccel() - incompatable sizes\n";
	    return -2;
    }    

    // create a copy if no trial exists    
    if (trialAccel == 0) {
	if (this->createAccel() < 0) {
	    cerr << "FATAL Node::setTrialAccel() - ran out of memory\n";
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
	cerr << "WARNING Node::incrTrialDisp() - incompatable sizes\n";
	return -2;
    }    

    // create a copy if no trial exists andd add committed
    if (trialDisp == 0) {
	if (this->createDisp() < 0) {
	    cerr << "FATAL Node::incrTrialDisp() - ran out of memory\n";
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
	cerr << "WARNING Node::incrTrialVel() - incompatable sizes\n";
	return -2;
    }    

    // create Vectors and array if none exist and set trial
    if (trialVel == 0) {
	if (this->createVel() < 0) {
	    cerr << "FATAL Node::incrTrialVel - ran out of memory\n";
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
	cerr << "WARNING Node::incrTrialAccel() - incompatable sizes\n";
	return -2;
    }    

    // create a copy if no trial exists andd add committed    
    if (trialAccel == 0) {
	if (this->createAccel() < 0) {
	    cerr << "FATAL Node::incrTrialAccel() - ran out of memory\n";
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
	cerr << "Node::addunbalLoad - load to add of incorrect size ";
	cerr << add.Size() << " should be " <<  numberDOF << endl;	
	return -1;
    }

    // if no load yet create it and assign
    if (unbalLoad == 0) {
	unbalLoad = new Vector(add); 
	if (unbalLoad == 0) {
	    cerr << "FATAL Node::addunbalLoad - ran out of memory\n";
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
    cerr << "Node::addInertiaLoadToUnbalance - accelG not of correct dimension";
    return -1;
  }

  // if no load yet create it and assign
  if (unbalLoad == 0) {
      unbalLoad = new Vector(numberDOF); 
      if (unbalLoad == 0 || unbalLoad->Size() != numberDOF) {
	  cerr << "FATAL Node::addunbalLoad - ran out of memory\n";
	  exit(-1);
      }  
  }

  // form - fact * M*R*accelG and add it to the unbalanced load
  (*unbalLoad) -= ((*mass) * (*R) * accelG)*fact;

  return 0;
}



const Vector &
Node::getUnbalancedLoad(void) 
{
    // make sure it was created before we return it
    if (unbalLoad == 0) {
	unbalLoad = new Vector(numberDOF);
	if (unbalLoad == 0 || unbalLoad->Size() != numberDOF) {
	    cerr << "FATAL Node::getunbalLoad() -- ran out of memory\n";
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
	    cerr << "FATAL Node::getunbalLoadWithInertia -- ran out of memory\n";
	    exit(-1);
	}
    } else
      (*unbalLoadWithInertia) = this->getUnbalancedLoad();

    if (mass != 0) {
      const Vector &theAccel = this->getTrialAccel(); // in case accel not created
      unbalLoadWithInertia->addMatrixVector(1.0, *mass, theAccel, -1.0);
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

    // if we get here we are done
    return 0;
}


const Matrix &
Node::getMass(void) 
{
    // make sure it was created before we return it
    if (mass == 0) {
	mass = new Matrix(numberDOF,numberDOF);
	if (mass == 0 || mass->noRows() != numberDOF) {
	    cerr << "FATAL Node::GetMass() -- ran out of memory\n";
	    exit(-1);
	}
    }
    
    return *mass;
}



int
Node::setMass(const Matrix &newMass)
{
    // check right size
    if (newMass.noRows() != numberDOF || newMass.noCols() != numberDOF) {
	cerr << "Node::setMass - incompatable matrices\n";
	return -1;
    }	

    // create a matrix if no mass yet set
    if (mass == 0) {
	mass = new Matrix(newMass);
	if (mass == 0 || mass->noRows() != numberDOF) {
	    cerr << "FATAL Node::setMass - ran out of memory\n";
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
      cerr << "FATAL Node::setNumColR() - out of memory\n";
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
    cerr << "Node:setR() - R has not been initialised\n";
    return -1;
  }
  
  // ensure row, col in range (matrix assignment will catch this - extra work)
  if (row < 0 || row > numberDOF || col < 0 || col > R->noCols()) {
    cerr << "Node:setR() - row, col index out of range\n";
    return -1;
  }

  // do the assignment
  (*R)(row,col) = Value;
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
	    cerr << "Node::getunbalLoadWithInertia -- ran out of memory\n";
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
	cerr << "WARNING Node::getRV() - R and V of incompatable dimesions\n";
	cerr << "R: " << *R << "V: " << V;
	unbalLoadWithInertia->Zero();
	return *unbalLoadWithInertia;
    }    

    // determine the product
    (*unbalLoadWithInertia) = (*R) * V;
    return *unbalLoadWithInertia;
}


int 
Node::setNumEigenvectors(int numVectorsToStore)
{
  // ensure a positive number of vectors
  if (numVectorsToStore <= 0) {
      g3ErrorHandler->warning("Node::setNumEigenvectors() - %d < 0\n", 
			      numVectorsToStore);
      return -1;
  }    

  // if matrix not yet assigned or not of correct size delete old and create new
  if (theEigenvectors == 0 || theEigenvectors->noCols() != numVectorsToStore) {
    if (theEigenvectors != 0)
      delete theEigenvectors;


    theEigenvectors = new Matrix(numberDOF, numVectorsToStore);
    if (theEigenvectors == 0 || theEigenvectors->noCols() != numVectorsToStore) {
      g3ErrorHandler->warning("Node::setNumEigenvectors() - out of memory\n");
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
      g3ErrorHandler->warning("Node::setEigenvectors() - mode %d invalid\n", mode);
      return -1;
  }

  if (eigenVector.Size() != numberDOF) {
      g3ErrorHandler->warning("Node::setEigenvectors() - eigenvector of incorrect size\n");
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
  if (theEigenvectors == 0)
    g3ErrorHandler->fatal("Node::getEigenvectors() - eigenvectors have not been set\n");
  
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
      cerr << " Node::sendSelf() - failed to send ID data\n";
      return res;
    }

    res = theChannel.sendVector(dataTag, cTag, *Crd);
    if (res < 0) {
      cerr << " Node::sendSelf() - failed to send Vecor data\n";
      return res;
    }

    if (commitDisp != 0) {
	res = theChannel.sendVector(dbTag1, cTag, *commitDisp);	
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send Disp data\n";
	  return res;
	}
    }

    if (commitVel != 0) {
	res = theChannel.sendVector(dbTag2, cTag, *commitVel);		
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send Vel data\n";
	  return res;
	}
    }

    if (commitAccel != 0) {
	res = theChannel.sendVector(dbTag3, cTag, *commitAccel); 
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send Accel data\n";
	  return res;
	}
    }

    if (mass != 0) {
	res = theChannel.sendMatrix(dataTag, cTag, *mass);
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send Mass data\n";
	  return res;
	}
    }
    
    if (R != 0) {
	res = theChannel.sendMatrix(dataTag, cTag, *R);
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send R data\n";
	  return res;
	}
    }    
    
    if (unbalLoad  != 0) {
	res = theChannel.sendVector(dbTag4, cTag, *unbalLoad);	
	if (res < 0) {
	  cerr << " Node::sendSelf() - failed to send Load data\n";
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
      cerr << "Node::recvSelf() - failed to receive ID data\n";
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
    if (Crd == 0)
      Crd = new Vector(numberCrd);

    // check we did not run out of memory
    if (Crd == 0) {
      g3ErrorHandler->warning("Node::recvSelf() - out of memory creating Coordinate vector\n");
      return -1;
    }

    if (theChannel.recvVector(dataTag, cTag, *Crd) < 0) {
      g3ErrorHandler->warning("Node::recvSelf() - failed to receive the Coordinate vector\n");
      return -2;
    }
	
    if (data(2) == 0) {
      // create the disp vectors if node is a total blank
      if (commitDisp == 0)
	this->createDisp();

      // recv the committed disp
      if (theChannel.recvVector(dbTag1, cTag, *commitDisp) < 0) {
	g3ErrorHandler->warning("Node::recvSelf - failed to receive Disp data\n");
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
	g3ErrorHandler->warning("Node::recvSelf - failed to receive Velocity data\n");
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
	g3ErrorHandler->warning("Node::recvSelf - failed to receive Acceleration data\n");
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
	  g3ErrorHandler->warning("Node::recvData -- ran out of memory\n");
	  return -5;
	}
      }
      if (theChannel.recvMatrix(dataTag, cTag, *mass) < 0) {
	g3ErrorHandler->warning("Node::recvSelf() - failed to receive Mass data\n");
	return -6;
      }
    }            
    
    if (data(12) == 0) {
      // create a matrix for R
      int noCols = data(13);
      if (R == 0) {
	R = new Matrix(numberDOF, noCols);
	if (R == 0) {
	  g3ErrorHandler->warning("Node::recvData -- ran out of memory\n");
	  return -1;
	}
      }
      // now recv the R matrix
      if (theChannel.recvMatrix(dataTag, cTag, *R) < 0) {
	g3ErrorHandler->warning("Node::recvSelf() - failed to receive R data\n");
	return res;
      }
    }

    if (data(6) == 0) {
      // create a vector for the load
      if (unbalLoad == 0) {
	unbalLoad = new Vector(numberDOF);
	if (unbalLoad == 0) {
	  g3ErrorHandler->warning("Node::recvData -- ran out of memory\n");
	  return -10;
	}
      }
      if (theChannel.recvVector(dbTag4, cTag, *unbalLoad) < 0) {
	g3ErrorHandler->warning("Node::recvSelf() - failed to receive Load data\n");
	return res;
      }
    }        

    return 0;
}



void
Node::Print(ostream &s, int flag)
{
  if (flag == 0) { // print out everything
    s << "\n Node: " << this->getTag() << endl;
    if (mass != 0) s << "\tMass : " << *mass;
    s << "\tCoordinates  : " << *Crd;
    if (commitDisp != 0)         
	s << "\tcommitDisps: " << *commitDisp;
    if (commitVel != 0)     
	s << "\tVelocities   : " << *commitVel;
    if (commitAccel != 0)         
	s << "\tcommitAccels: " << *commitAccel;
    if (unbalLoad != 0)
      s << "\t unbalanced Load: " << *unbalLoad;

    if (theEigenvectors != 0)
	s << "\t Eigenvectors: \n" << *theEigenvectors;
    s << "\n"; 
  }
  else if (flag == 1) { // print out: nodeId displacements
    s << this->getTag() << "  " << *commitDisp;
  }
}
  
int
Node::displaySelf(Renderer &theRenderer, int displayMode, float fact)
{
    if (displayMode == 1) { // draw a text string containing tag
	const Vector &theDisp = this->getDisp();
	Vector position(*Crd);
	for (int i=0; i<Crd->Size(); i++) 
	    position(i) += theDisp(i)*fact;	
	
	char theText[20];
	NodeItoa(this->getTag(), theText);    
	return theRenderer.drawGText(position, theText, strlen(theText));
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
    g3ErrorHandler->warning("WARNING - Node::createDisp() %s %d\n",
			    "ran out of memory for array of size", 2*numberDOF);
    return -1;
  }
  for (int i=0; i<4*numberDOF; i++)
    disp[i] = 0.0;
    
  commitDisp = new Vector(&disp[numberDOF], numberDOF); 
  trialDisp = new Vector(disp, numberDOF);
  incrDisp = new Vector(&disp[2*numberDOF], numberDOF);
  incrDeltaDisp = new Vector(&disp[3*numberDOF], numberDOF);
  
  if (commitDisp == 0 || trialDisp == 0 || incrDisp == 0 || incrDeltaDisp == 0) {
    g3ErrorHandler->warning("WARNING - Node::createDisp() %s \n",
			    "ran out of memory creating Vectors(double *,int)");
    return -2;
  }
    
  return 0;
}


int
Node::createVel(void)
{
    vel = new double[2*numberDOF];
    
    if (vel == 0) {
      g3ErrorHandler->warning("WARNING - Node::createVel() %s %d\n",
			      "ran out of memory for array of size ",2*numberDOF);
	return -1;
    }
    for (int i=0; i<2*numberDOF; i++)
      vel[i] = 0.0;
    
    commitVel = new Vector(&vel[numberDOF], numberDOF); 
    trialVel = new Vector(vel, numberDOF);
    
    if (commitVel == 0 || trialVel == 0) {
      g3ErrorHandler->warning("WARNING - Node::createVel() %s\n",
			      "ran out of memory creating Vectors(double *,int) ");
      return -2;
    }
    
    return 0;
}

int
Node::createAccel(void)
{
    accel = new double[2*numberDOF];
    
    if (accel == 0) {
      g3ErrorHandler->warning("WARNING - Node::createAccel() %s %d\n",
				"ran out of memory for array of size ",2*numberDOF);
	return -1;
    }
    for (int i=0; i<2*numberDOF; i++)
	accel[i] = 0.0;
    
    commitAccel = new Vector(&accel[numberDOF], numberDOF);
    trialAccel = new Vector(accel, numberDOF);
    
    if (commitAccel == 0 || trialAccel == 0) {
	g3ErrorHandler->warning("WARNING - Node::createAccel() %s\n",
				"ran out of memory creating Vectors(double *,int)");
	return -2;
    }

    return 0;
}


char 
NodeItoc(int x)
{
 if (x == 1) return '1';
 if (x == 2) return '2';
 if (x == 3) return '3';
 if (x == 4) return '4';
 if (x == 5) return '5';
 if (x == 6) return '6';
 if (x == 7) return '7';
 if (x == 8) return '8';
 if (x == 9) return '9';
 return '0';
}

void
NodeItoa(int x, char *str)
{
  int y=x;
  while (y >= 10) 
    y = y/10;
  str[0] = NodeItoc(y);
  str[1] = '\0';
  if (x >= 10) {
    int z = x/10;
    z = x - 10*z;
    NodeItoa(z,&str[1]);
  }
}
