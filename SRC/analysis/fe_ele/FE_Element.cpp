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
                                                                        
// $Revision: 1.21 $
// $Date: 2010-02-04 01:21:04 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/fe_ele/FE_Element.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the code for implementing the methods
// of the FE_Element class interface.

#include <FE_Element.h>
#include <stdlib.h>

#include <Element.h>
#include <Domain.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Integrator.h>
#include <Subdomain.h>
#include <AnalysisModel.h>
#include <Matrix.h>
#include <Vector.h>

#define MAX_NUM_DOF 64

// static variables initialisation
Matrix FE_Element::errMatrix(1,1);
Vector FE_Element::errVector(1);
Matrix **FE_Element::theMatrices; // pointers to class wide matrices
Vector **FE_Element::theVectors;  // pointers to class widde vectors
int FE_Element::numFEs(0);           // number of objects

//  FE_Element(Element *, Integrator *theIntegrator);
//	construictor that take the corresponding model element.
FE_Element::FE_Element(int tag, Element *ele)
  :TaggedObject(tag),
   myDOF_Groups((ele->getExternalNodes()).Size()), myID(ele->getNumDOF()), 
   numDOF(ele->getNumDOF()), theModel(0), myEle(ele), 
   theResidual(0), theTangent(0), theIntegrator(0)
{
  if (numDOF <= 0) {
    opserr << "FE_Element::FE_Element(Element *) ";
    opserr << " element must have 1 dof " << *ele;
    exit(-1);
  }
    
    // get elements domain & check it is valid
    Domain *theDomain = ele->getDomain();
    if (theDomain == 0) {
	opserr << "FATAL FE_Element::FE_Element() - element has no domain "<< *ele;
	exit(-1);
    }

    // keep a pointer to all DOF_Groups
    int numGroups = ele->getNumExternalNodes();
    const ID &nodes = ele->getExternalNodes();

    for (int i=0; i<numGroups; i++) {
	Node *nodePtr =theDomain->getNode(nodes(i));
	if (nodePtr == 0) {
	    opserr << "FATAL FE_Element::FE_Element() - Node: ";
	    opserr <<  nodes(i) <<  "does not exist in the Domain\n";
	    opserr << *ele;
	    exit(-1);
	}
	
	DOF_Group *dofGrpPtr = nodePtr->getDOF_GroupPtr();
	if (dofGrpPtr != 0) 
	    myDOF_Groups(i) = dofGrpPtr->getTag();	
	else {
	    opserr << "FATAL FE_Element::FE_Element() - Node: ";
	    opserr <<  *nodePtr <<  " has no DOF_Group associated with it\n";
	    exit(-1);
	}
    }

    // if this is the first FE_Element we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numFEs == 0) {
	theMatrices = new Matrix *[MAX_NUM_DOF+1];
	theVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (theMatrices == 0 || theVectors == 0) {
	    opserr << "FE_Element::FE_Element(Element *) ";
	    opserr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    theMatrices[i] = 0;
	    theVectors[i] = 0;
	}
    }

    if (ele->isSubdomain() == false) {
	
	// if Elements are not subdomains, set up pointers to
	// objects to return tangent Matrix and residual Vector.

	if (numDOF <= MAX_NUM_DOF) {
	    // use class wide objects
	    if (theVectors[numDOF] == 0) {
		theVectors[numDOF] = new Vector(numDOF);
		theMatrices[numDOF] = new Matrix(numDOF,numDOF);
		theResidual = theVectors[numDOF];
		theTangent = theMatrices[numDOF];
		if (theResidual == 0 || theResidual->Size() != numDOF ||	
		    theTangent == 0 || theTangent->noCols() != numDOF)	{  
		    opserr << "FE_Element::FE_Element(Element *) ";
		    opserr << " ran out of memory for vector/Matrix of size :";
		    opserr << numDOF << endln;
		    exit(-1);
		}
	    } else {
		theResidual = theVectors[numDOF];
		theTangent = theMatrices[numDOF];
	    }
	} else {
	    // create matrices and vectors for each object instance
	    theResidual = new Vector(numDOF);
	    theTangent = new Matrix(numDOF, numDOF);
	    if (theResidual == 0 || theTangent ==0 ||
		theTangent ==0 || theTangent->noRows() ==0) {
	    
		opserr << "FE_Element::FE_Element(Element *) ";
		opserr << " ran out of memory for vector/Matrix of size :";
		opserr << numDOF << endln;
		exit(-1);
	    }
	} 
    } else {

	// as subdomains have own matrix for tangent and residual don't need
	// to set matrix and vector pointers to these objects
        theResidual = new Vector(numDOF);	
 	// invoke setFE_ElementPtr() method on Subdomain
	Subdomain *theSub = (Subdomain *)ele;
	theSub->setFE_ElementPtr(this);
    }
 
    // increment number of FE_Elements by 1
    numFEs++;
}


FE_Element::FE_Element(int tag, int numDOF_Group, int ndof)
  :TaggedObject(tag),
   myDOF_Groups(numDOF_Group), myID(ndof), numDOF(ndof), theModel(0),
   myEle(0), theResidual(0), theTangent(0), theIntegrator(0)
{
    // this is for a subtype, the subtype must set the myDOF_Groups ID array
    numFEs++;

    // if this is the first FE_Element we now
    // create the arrays used to store pointers to class wide
    // matrix and vector objects used to return tangent and residual
    if (numFEs == 0) {
	theMatrices = new Matrix *[MAX_NUM_DOF+1];
	theVectors  = new Vector *[MAX_NUM_DOF+1];
	
	if (theMatrices == 0 || theVectors == 0) {
	    opserr << "FE_Element::FE_Element(Element *) ";
	    opserr << " ran out of memory";	    
	}
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    theMatrices[i] = 0;
	    theVectors[i] = 0;
	}
    }
    
    // as subtypes have no access to the tangent or residual we don't set them
    // this way we can detect if subclass does not provide all methods it should
}



// ~FE_Element();    
//	destructor.
FE_Element::~FE_Element()
{

    // decrement number of FE_Elements
    numFEs--;

    // delete tangent and residual if created specially
    if (numDOF > MAX_NUM_DOF) {
	if (theTangent != 0) delete theTangent;
	if (theResidual != 0) delete theResidual;
    }

    // if this is the last FE_Element, clean up the
    // storage for the matrix and vector objects
    if (numFEs == 0) {
	for (int i=0; i<MAX_NUM_DOF; i++) {
	    if (theVectors[i] != 0)
		delete theVectors[i];
	    if (theMatrices[i] != 0)
		delete theMatrices[i];
	}	
	delete [] theMatrices;
	delete [] theVectors;
    }
}    


const ID &
FE_Element::getDOFtags(void) const 
{
    return myDOF_Groups;
}


// const ID &getID(void) const;
//	Method to return the current ID.

const ID &
FE_Element::getID(void) const
{
    return myID;
}

void 
FE_Element::setAnalysisModel(AnalysisModel &theAnalysisModel)
{
    theModel = &theAnalysisModel;
}

// void setID(int index, int value);
//	Method to set the corresponding index of the ID to value.

int
FE_Element::setID(void)
{
    int current = 0;
    
    if (theModel == 0) {
	opserr << "WARNING FE_Element::setID() - no AnalysisModel set\n";
	return -1;
    }
    
    int numGrps = myDOF_Groups.Size();
    for (int i=0; i<numGrps; i++) {
	int tag = myDOF_Groups(i);

	DOF_Group *dofPtr = theModel->getDOF_GroupPtr(tag);
	if (dofPtr == 0) {
	    opserr << "WARNING FE_Element::setID: 0 DOF_Group Pointer\n";
	    return -2;
	}
	    
	const ID &theDOFid = dofPtr->getID();

	for (int j=0; j<theDOFid.Size(); j++)  
	    if (current < numDOF)
		myID(current++) = theDOFid(j);
	    else {
		opserr << "WARNING FE_Element::setID() - numDOF and";
		opserr << " number of dof at the DOF_Groups\n";
		return -3;
	    }		
    }
    return 0;
}


const Matrix &
FE_Element::getTangent(Integrator *theNewIntegrator)
{
    theIntegrator = theNewIntegrator;
    
    if (myEle == 0) {
	opserr << "FATAL FE_Element::getTangent() - no Element *given ";
	opserr << "- subclasses must provide implementation - ";
	opserr << " - a 1x1 error matrix will be returned.\n";
	exit(-1);
    }

    if (myEle->isSubdomain() == false) {
      if (theNewIntegrator != 0)
	theNewIntegrator->formEleTangent(this);	    	    

      return *theTangent;
    } else {
      Subdomain *theSub = (Subdomain *)myEle;
      theSub->computeTang();	    
      return theSub->getTang();
    }
}

const Vector &
FE_Element::getResidual(Integrator *theNewIntegrator)
{
    theIntegrator = theNewIntegrator;

    if (theIntegrator == 0)
      return *theResidual;

    if (myEle == 0) {
	opserr << "FATAL FE_Element::getTangent() - no Element *given ";
	opserr << "- subclasses must provide implementation - ";
	opserr << " - an error Vector of order 1 will be returned.\n";
	exit(-1);
    }    

    if (myEle->isSubdomain() == false) {
      theNewIntegrator->formEleResidual(this);
      return *theResidual;
    } else {
      Subdomain *theSub = (Subdomain *)myEle;
      theSub->computeResidual();	    
      return theSub->getResistingForce();
    }
}



void  
FE_Element::zeroTangent(void)
{
    if (myEle != 0) {
	if (myEle->isSubdomain() == false)
	    theTangent->Zero();
	else {
	    opserr << "WARNING FE_Element::zeroTangent() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    
    }
}

void  
FE_Element::addKtToTang(double fact)
{
    if (myEle != 0) {
	
	// check for a quick return	
	if (fact == 0.0) 
	    return;
	else if (myEle->isSubdomain() == false)	    
	    theTangent->addMatrix(1.0, myEle->getTangentStiff(),fact);
	else {
	    opserr << "WARNING FE_Element::addKToTang() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    
    }
}

void  
FE_Element::addCtoTang(double fact)
{
    if (myEle != 0) {
	
	// check for a quick return	
	if (fact == 0.0) 
	  return;
	else if (myEle->isSubdomain() == false)	    	    
	  theTangent->addMatrix(1.0, myEle->getDamp(),fact);
	else {
	  opserr << "WARNING FE_Element::addCToTang() - ";
	  opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    	
    }
}
    
void  
FE_Element::addMtoTang(double fact)
{
    if (myEle != 0) {

	// check for a quick return	
	if (fact == 0.0) 
	  return;
	else if (myEle->isSubdomain() == false)	    	    
	  theTangent->addMatrix(1.0, myEle->getMass(),fact);
	else {
	  opserr << "WARNING FE_Element::addMToTang() - ";
	  opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    	
    }	
}    


void
FE_Element::addKiToTang(double fact)
{
  if (myEle != 0) {
    // check for a quick return	
    if (fact == 0.0) 
      return;
    else if (myEle->isSubdomain() == false)	    	    
      theTangent->addMatrix(1.0, myEle->getInitialStiff(), fact);
    else {
	opserr << "WARNING FE_Element::addKiToTang() - ";
	opserr << "- this should not be called on a Subdomain!\n";
    }    	    	    	
  }	
}

void
FE_Element::addKgToTang(double fact)
{
  if (myEle != 0) {
    // check for a quick return	
    if (fact == 0.0) 
      return;
    else if (myEle->isSubdomain() == false)	    	    
      theTangent->addMatrix(1.0, myEle->getGeometricTangentStiff(), fact);
    else {
	opserr << "WARNING FE_Element::addKgToTang() - ";
	opserr << "- this should not be called on a Subdomain!\n";
    }    	    	    	
  }	
}    

void
FE_Element::addKpToTang(double fact, int numP)
{
  if (myEle != 0) {
    // check for a quick return	
    if (fact == 0.0) 
      return;
    else if (myEle->isSubdomain() == false) {
      const Matrix *thePrevMat = myEle->getPreviousK(numP);
      if (thePrevMat != 0)
	theTangent->addMatrix(1.0, *thePrevMat, fact);
    } else {
      opserr << "WARNING FE_Element::addKpToTang() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }    	    	    	
  }	
}    

int
FE_Element::storePreviousK(int numP)
{
  int res = 0;
  if (myEle != 0) {
    res = myEle->storePreviousK(numP);
  } 
  
  return res;
}


void  
FE_Element::zeroResidual(void)
{
    if (myEle != 0) {
	if (myEle->isSubdomain() == false)
	    theResidual->Zero();
	else {
	    opserr << "WARNING FE_Element::zeroResidual() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    
    }
    else {
	opserr << "WARNING FE_Element::zeroResidual() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    
}

void  
FE_Element::addRtoResidual(double fact)
{
  if (myEle != 0) {
    // check for a quick return	
    if (fact == 0.0) 
      return;
    else if (myEle->isSubdomain() == false) {
      const Vector &eleResisting = myEle->getResistingForce();
      theResidual->addVector(1.0, eleResisting, -fact);
    }
    else {
      opserr << "WARNING FE_Element::addRtoResidual() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }    	    	    	
  }
  else {
    opserr << "WARNING FE_Element::addRtoResidual() - no Element *given ";
    opserr << "- subclasses must provide implementation\n";
    }    	        
}


void  
FE_Element::addRIncInertiaToResidual(double fact)
{
    if (myEle != 0) {
	// check for a quick return	
	if (fact == 0.0) 
	    return;
	else if (myEle->isSubdomain() == false) {
	  const Vector &eleResisting = myEle->getResistingForceIncInertia();
	  theResidual->addVector(1.0, eleResisting, -fact);
	}
	else {
	    opserr << "WARNING FE_Element::addRtoResidual() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    	
    }
    else {
	opserr << "WARNING FE_Element::addRtoResidual() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	        
}


const Vector &
FE_Element::getTangForce(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// zero out the force vector
	theResidual->Zero();

	// check for a quick return
	if (fact == 0.0) 
	    return *theResidual;

	// get the components we need out of the vector
	// and place in a temporary vector
	Vector tmp(numDOF);
	for (int i=0; i<numDOF; i++) {
	  int dof = myID(i);
	  if (dof >= 0)
	    tmp(i) = disp(myID(i));
	  else
	    tmp(i) = 0.0;
	}

	if (myEle->isSubdomain() == false) {
	    // form the tangent again and then add the force
	    theIntegrator->formEleTangent(this);
	    if (theResidual->addMatrixVector(1.0, *theTangent,tmp,fact) < 0) {
		opserr << "WARNING FE_Element::getTangForce() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }				
	}
	else {
	    Subdomain *theSub = (Subdomain *)myEle;
	    if (theResidual->addMatrixVector(1.0, theSub->getTang(),tmp,fact) < 0) {
		opserr << "WARNING FE_Element::getTangForce() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }						
	}
	return *theResidual;
    }
    else {
	opserr << "WARNING FE_Element::addTangForce() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
	return errVector;	
    }    	            
}



const Vector &
FE_Element::getK_Force(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// zero out the force vector
	theResidual->Zero();

	// check for a quick return
	if (fact == 0.0) 
	    return *theResidual;

	// get the components we need out of the vector
	// and place in a temporary vector
	Vector tmp(numDOF);
	for (int i=0; i<numDOF; i++) {
	  int dof = myID(i);
	  if (dof >= 0)
	    tmp(i) = disp(myID(i));
	  else
	    tmp(i) = 0.0;
	}

	if (theResidual->addMatrixVector(1.0, myEle->getTangentStiff(), tmp, fact) < 0){
	  opserr << "WARNING FE_Element::getKForce() - ";
	  opserr << "- addMatrixVector returned error\n";		 
	}		

	return *theResidual;
    }
    else {
	opserr << "WARNING FE_Element::getKForce() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
	return errVector;	
    }    	            
}


const Vector &
FE_Element::getKi_Force(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// zero out the force vector
	theResidual->Zero();

	// check for a quick return
	if (fact == 0.0) 
	    return *theResidual;

	// get the components we need out of the vector
	// and place in a temporary vector
	Vector tmp(numDOF);
	for (int i=0; i<numDOF; i++) {
	  int dof = myID(i);
	  if (dof >= 0)
	    tmp(i) = disp(myID(i));
	  else
	    tmp(i) = 0.0;
	}

	if (theResidual->addMatrixVector(1.0, myEle->getInitialStiff(), tmp, fact) < 0){
	  opserr << "WARNING FE_Element::getKForce() - ";
	  opserr << "- addMatrixVector returned error\n";		 
	}		

	return *theResidual;
    }
    else {
	opserr << "WARNING FE_Element::getKForce() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
	return errVector;	
    }    	            
}

const Vector &
FE_Element::getM_Force(const Vector &disp, double fact)
{

    if (myEle != 0) {    

	// zero out the force vector
	theResidual->Zero();

	// check for a quick return
	if (fact == 0.0) 
	    return *theResidual;

	// get the components we need out of the vector
	// and place in a temporary vector
	Vector tmp(numDOF);
	for (int i=0; i<numDOF; i++) {
	  int dof = myID(i);
	  if (dof >= 0)
	    tmp(i) = disp(myID(i));
	  else
	    tmp(i) = 0.0;
	}

	if (theResidual->addMatrixVector(1.0, myEle->getMass(), tmp, fact) < 0){
	  opserr << "WARNING FE_Element::getMForce() - ";
	  opserr << "- addMatrixVector returned error\n";		 
	}		


	return *theResidual;
    }
    else {
	opserr << "WARNING FE_Element::getMForce() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
	return errVector;	
    }    	            
}

const Vector &
FE_Element::getC_Force(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// zero out the force vector
	theResidual->Zero();

	// check for a quick return
	if (fact == 0.0) 
	    return *theResidual;

	// get the components we need out of the vector
	// and place in a temporary vector
	Vector tmp(numDOF);
	for (int i=0; i<numDOF; i++) {
	  int dof = myID(i);
	  if (dof >= 0)
	    tmp(i) = disp(myID(i));
	  else
	    tmp(i) = 0.0;
	}

	if (theResidual->addMatrixVector(1.0, myEle->getDamp(), tmp, fact) < 0){
	  opserr << "WARNING FE_Element::getDForce() - ";
	  opserr << "- addMatrixVector returned error\n";		 
	}		

	return *theResidual;
    }
    else {
	opserr << "WARNING FE_Element::getDForce() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
	return errVector;	
    }    	            
}



Integrator *
FE_Element::getLastIntegrator(void)
{
    return theIntegrator;
}

const Vector &
FE_Element::getLastResponse(void)
{
    if (myEle != 0) {
      if (theIntegrator != 0) {
	if (theIntegrator->getLastResponse(*theResidual,myID) < 0) {
	  opserr << "WARNING FE_Element::getLastResponse(void)";
	  opserr << " - the Integrator had problems with getLastResponse()\n";
	}
      }
      else {
	theResidual->Zero();
	opserr << "WARNING  FE_Element::getLastResponse()";
	opserr << " No Integrator yet passed\n";
      }
    
      Vector &result = *theResidual;
      return result;
    }
    else {
	opserr << "WARNING  FE_Element::getLastResponse()";
	opserr << " No Element passed in constructor\n";
	return errVector;
    }
}

void  
FE_Element::addM_Force(const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    // get the components we need out of the vector
	    // and place in a temporary vector
	    Vector tmp(numDOF);
	    for (int i=0; i<numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
		    tmp(i) = accel(loc);
		else
		    tmp(i) = 0.0;		
	    }	 
		
	    if (theResidual->addMatrixVector(1.0, myEle->getMass(), tmp, fact) < 0){
		opserr << "WARNING FE_Element::addM_Force() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addM_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addM_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}

void  
FE_Element::addD_Force(const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    // get the components we need out of the vector
	    // and place in a temporary vector
	    Vector tmp(numDOF);
	    for (int i=0; i<numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
		    tmp(i) = accel(loc);
		else
		    tmp(i) = 0.0;		
	    }	  
		
	    if (theResidual->addMatrixVector(1.0, myEle->getDamp(), tmp, fact) < 0){
		opserr << "WARNING FE_Element::addD_Force() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addD_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addD_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}

void  
FE_Element::addK_Force(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    // get the components we need out of the vector
	    // and place in a temporary vector
	    Vector tmp(numDOF);
	    for (int i=0; i<numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
		    tmp(i) = disp(loc);
		else
		    tmp(i) = 0.0;		
	    }	  
		
	    if (theResidual->addMatrixVector(1.0, myEle->getTangentStiff(), tmp, fact) < 0){
		opserr << "WARNING FE_Element::addK_Force() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addK_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addK_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}

void  
FE_Element::addKg_Force(const Vector &disp, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    // get the components we need out of the vector
	    // and place in a temporary vector
	    Vector tmp(numDOF);
	    for (int i=0; i<numDOF; i++) {
		int loc = myID(i);
		if (loc >= 0)
		    tmp(i) = disp(loc);
		else
		    tmp(i) = 0.0;		
	    }	  
		
	    if (theResidual->addMatrixVector(1.0, myEle->getGeometricTangentStiff(), tmp, fact) < 0){
		opserr << "WARNING FE_Element::addKg_Force() - ";
		opserr << "- addMatrixVector returned error\n";		 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addKg_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addKg_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}


void  
FE_Element::addLocalM_Force(const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    if (theResidual->addMatrixVector(1.0, myEle->getMass(),
					     accel, fact) < 0){

	      opserr << "WARNING FE_Element::addLocalM_Force() - ";
	      opserr << "- addMatrixVector returned error\n"; 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addLocalM_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addLocalM_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}

void  
FE_Element::addLocalD_Force(const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    if (theResidual->addMatrixVector(1.0, myEle->getDamp(),
					     accel, fact) < 0){

	      opserr << "WARNING FE_Element::addLocalD_Force() - ";
	      opserr << "- addMatrixVector returned error\n"; 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addLocalD_Force() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addLocalD_Force() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}


Element *
FE_Element::getElement(void)
{
  return myEle;
}


// AddingSensitivity:BEGIN /////////////////////////////////
void  
FE_Element::addResistingForceSensitivity(int gradNumber, double fact)
{
  theResidual->addVector(1.0, myEle->getResistingForceSensitivity(gradNumber), -fact);
}

void  
FE_Element::addM_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
  // Get the components we need out of the vector
  // and place in a temporary vector
  Vector tmp(numDOF);
  for (int i=0; i<numDOF; i++) {
    int loc = myID(i);
    if (loc >= 0) {
      tmp(i) = vect(loc);
    }
    else {
      tmp(i) = 0.0;
    }
  }
  if (theResidual->addMatrixVector(1.0, myEle->getMassSensitivity(gradNumber),tmp,fact) < 0) {
    opserr << "WARNING FE_Element::addM_ForceSensitivity() - ";
    opserr << "- addMatrixVector returned error\n";		 
  }
}

void  
FE_Element::addD_ForceSensitivity(int gradNumber, const Vector &vect, double fact)
{
  if (myEle != 0) {    
    
    // check for a quick return
    if (fact == 0.0) 
      return;
    if (myEle->isSubdomain() == false) {
      // get the components we need out of the vector
      // and place in a temporary vector
      Vector tmp(numDOF);
      for (int i=0; i<numDOF; i++) {
	int loc = myID(i);
	if (loc >= 0)
	  tmp(i) = vect(loc);
	else
	  tmp(i) = 0.0;		
      }	
      if (theResidual->addMatrixVector(1.0, myEle->getDampSensitivity(gradNumber), tmp, fact) < 0){
	opserr << "WARNING FE_Element::addD_ForceSensitivity() - ";
	opserr << "- addMatrixVector returned error\n";		 
      }		
    }
    else {
      opserr << "WARNING FE_Element::addD_ForceSensitivity() - ";
      opserr << "- this should not be called on a Subdomain!\n";
    }    	    	    				
  }
  else {
    opserr << "WARNING FE_Element::addD_ForceSensitivity() - no Element *given ";
    opserr << "- subclasses must provide implementation\n";
  }    	            
}

void  
FE_Element::addLocalD_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    if (theResidual->addMatrixVector(1.0, myEle->getDampSensitivity(gradNumber),
					     accel, fact) < 0){

	      opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
	      opserr << "- addMatrixVector returned error\n"; 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}

void  
FE_Element::addLocalM_ForceSensitivity(int gradNumber, const Vector &accel, double fact)
{
    if (myEle != 0) {    

	// check for a quick return
	if (fact == 0.0) 
	    return;
	if (myEle->isSubdomain() == false) {
	    if (theResidual->addMatrixVector(1.0, myEle->getMassSensitivity(gradNumber),
					     accel, fact) < 0){

	      opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
	      opserr << "- addMatrixVector returned error\n"; 
	    }		
	}
	else {
	    opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - ";
	    opserr << "- this should not be called on a Subdomain!\n";
	}    	    	    				
    }
    else {
	opserr << "WARNING FE_Element::addLocalD_ForceSensitivity() - no Element *given ";
	opserr << "- subclasses must provide implementation\n";
    }    	            
}





int  
FE_Element::commitSensitivity(int gradNum, int numGrads)
{
  myEle->commitSensitivity(gradNum, numGrads);
  
  return 0;
}

// AddingSensitivity:END ////////////////////////////////////


int  
FE_Element::updateElement(void)
{
  if (myEle != 0) {
    return myEle->update();
    opserr << "FE_Element::update()"; myEle->Print(opserr, 0);
  }

  
  return 0;
}
