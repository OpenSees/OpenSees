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
// $Date: 2008-03-10 18:25:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/node/NodalLoad.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the implementation of NodalLoad

#include <NodalLoad.h>
#include <stdlib.h>
#include <Domain.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>
#include <string>
#include <elementAPI.h>

// AddingSensitivity:BEGIN /////////////////////////////////////
Vector NodalLoad::gradientVector(1);
// AddingSensitivity:END ///////////////////////////////////////

NodalLoad::NodalLoad(int theClasTag)
:Load(0,theClasTag), 
 myNode(0), myNodePtr(0), load(0), konstant(false)
{
    // constructor for a FEM_ObjectBroker
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////
}

NodalLoad::NodalLoad(int tag, int node, int theClassTag)
:Load(tag,theClassTag), 
 myNode(node), myNodePtr(0), load(0), konstant(false)
{
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  // AddingSensitivity:END ///////////////////////////////////////
}

NodalLoad::NodalLoad(int tag, int node, const Vector &theLoad, bool isLoadConstant)
:Load(tag, LOAD_TAG_NodalLoad), 
 myNode(node), myNodePtr(0), load(0), konstant(isLoadConstant)
{
    load = new Vector(theLoad);    

    if (load == 0) {
	opserr << "FATAL NodalLoad::NodalLoad(int node, const Vector &theLoad) -";
	opserr << " ran out of memory for load on Node " << node << endln;
	exit(-1);
    }
    // AddingSensitivity:BEGIN /////////////////////////////////////
    parameterID = 0;
    // AddingSensitivity:END ///////////////////////////////////////
}


// ~NodalLoad()
// 	destructor

NodalLoad::~NodalLoad()
{
    if (load != 0)
	delete load;
}

void 
NodalLoad::setDomain(Domain *newDomain)
{
    // first get myNodePtr
  if (newDomain == 0)
    return;

  // invoke the ancestor class method
  this->DomainComponent::setDomain(newDomain);    

  /*
  if (newDomain != 0) {
    myNodePtr = newDomain->getNode(myNode);
    if (myNodePtr == 0) {
      opserr << *newDomain;
      opserr << "WARNING NodalLoad::setDomain() - No associated Node node " ;
      opserr << " for NodalLoad " << *this;
      //	opserr << *newDomain;

      return;
    }
  }
  */
}

int 
NodalLoad::getNodeTag(void) const
{
    return myNode;
}


void
NodalLoad::applyLoad(double loadFactor)
{
    if (myNodePtr == 0) {
      Domain *theDomain=this->getDomain();
      if ((theDomain == 0) || 
	  (myNodePtr = theDomain->getNode(myNode)) == 0) {
	opserr << "WARNING NodalLoad::applyLoad() - No associated Node node " ;
	opserr << " for NodalLoad " << *this;
	return;
      }
    }

    // add the load times the loadfactor to nodal unbalanced load
    if (konstant == false)
	myNodePtr->addUnbalancedLoad(*load,loadFactor);
    else
	myNodePtr->addUnbalancedLoad(*load,1.0);	
    
    //    opserr << "loadFactor: " << loadFactor << *myNodePtr;
}

void
NodalLoad::applyLoadSensitivity(double loadFactor)
{
    if (myNodePtr == 0) {
	Domain *theDomain=this->getDomain();
	if ((theDomain == 0) || 
	    (myNodePtr = theDomain->getNode(myNode)) == 0) {
	    opserr << "WARNING NodalLoad::applyLoadSensitivity() - No associated Node node " ;
	    opserr << " for NodalLoad " << *this;
	    return;
	}
    }

    // load sensitivity
    Vector loadsens(load->Size());
    if(parameterID == 0) return;
    if(parameterID > loadsens.Size()) return;
    loadsens(parameterID-1) = 1;

    // add the load times the loadfactor to nodal unbalanced load
    if (konstant == false)
	myNodePtr->addUnbalancedLoad(loadsens,loadFactor);
    else
	myNodePtr->addUnbalancedLoad(loadsens,1.0);	
    
    //    opserr << "loadFactor: " << loadFactor << *myNodePtr;
} 

int 
NodalLoad::sendSelf(int cTag, Channel &theChannel)
{
    int dataTag = this->getDbTag();
    ID data(5);
    data(0) = this->getTag();    
    data(1) = myNode;
    if (load != 0)
	data(2) = load->Size();
    else
	data(2) = 0;
    data(3) = konstant;
    data(4) = this->getLoadPatternTag();
    
    int result = theChannel.sendID(dataTag, cTag, data);
    if (result < 0) {
	opserr << "NodalLoad::sendSelf - failed to send data\n";
	return result;
    }

    if (load != 0){
	int result = theChannel.sendVector(dataTag, cTag, *load);
	if (result < 0) {
	    opserr << "NodalLoad::sendSelf - failed to Load data\n";
	    return result;
	}
    }    

    return 0;
}

int 
NodalLoad::recvSelf(int cTag, Channel &theChannel, 
		    FEM_ObjectBroker &theBroker)
{	
    int result;
    int dataTag = this->getDbTag();
    ID data(5);
    result = theChannel.recvID(dataTag, cTag, data);
    if (result < 0) {
      opserr << "NodalLoad::recvSelf() - failed to recv data\n";
      return result;
    }    
    this->setTag(data(0));
    myNode = data(1);
    int loadSize = data(2);
    konstant = (data(3) != 0);
    this->setLoadPatternTag(data(4));
    if (loadSize != 0) {
	load = new Vector(data(2));
	result = theChannel.recvVector(dataTag, cTag, *load);
	if (result < 0) {
	  opserr << "NodalLoad::recvSelf() - failed to recv load\n";
	  return result;
	}    
    }

    return 0;
}


void
NodalLoad::Print(OPS_Stream &s, int flag)
{
     s << "Nodal Load: " << myNode;
     if (load != 0)
	 s << " load : " << *load;
}


// AddingSensitivity:BEGIN /////////////////////////////////////
int
NodalLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"1") == 0) {
    param.setValue((*load)(0));
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"2") == 0) {
    param.setValue((*load)(1));
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"3") == 0) {
    param.setValue((*load)(2));
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"4") == 0) {
    param.setValue((*load)(3));
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"5") == 0) {
    param.setValue((*load)(4));
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"6") == 0) {
    param.setValue((*load)(5));
    return param.addObject(6, this);
  }

  return -1;
}

int
NodalLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case -1:
    return -1;
  case 1:
    (*load)(0) = info.theDouble;
    return 0;
  case 2:
    (*load)(1) = info.theDouble;
    return 0;
  case 3:
    (*load)(2) = info.theDouble;
    return 0;
  case 4:
    (*load)(3) = info.theDouble;
    return 0;
  case 5:
    (*load)(4) = info.theDouble;
    return 0;
  case 6:
    (*load)(5) = info.theDouble;
    return 0;

  default:
    return -1;
  }
}


int
NodalLoad::activateParameter(int passedParameterID)
{
  parameterID = passedParameterID;
  
  return 0;
}


const Vector &
NodalLoad::getExternalForceSensitivity(int gradNumber)
{
  gradientVector(0) = (double)parameterID;
  
  return gradientVector;
}


// AddingSensitivity:END //////////////////////////////////////

