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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/region/MeshRegion.cpp,v $
                                                                        
                                                                        
// Written: fmk 
//
// Description: This file contains the class definition for MeshRegion.
// A MeshRegion is a part of the domain which is defined by a set of
// Elements and Nodes (all the end nodes of the elements are in the region, 
// as are all elements whose end nodes are in the region)
//
// What: "@(#) MeshRegion.h, revA"


#include <MeshRegion.h>
#include <stdlib.h>
#include <string.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <ID.h>
#include <Domain.h>
#include <NodeRecorder.h>
#include <ElementRecorder.h>
#include <NodeIter.h>
#include <ElementIter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <SP_Constraint.h>
#include <elementAPI.h>

MeshRegion::MeshRegion(int tag) 
  :DomainComponent(tag, REGION_TAG_MeshRegion), 
  alphaM(0), betaK(0), betaK0(0), betaKc(0),
   theNodes(0), theElements(0), xEles(),
  currentGeoTag(0), lastGeoSendTag(-1), dbNod(0), dbEle(0)
{
    // does nothing
}

MeshRegion::MeshRegion(int tag, int cTag) 
  :DomainComponent(tag, cTag), alphaM(0), betaK(0), betaK0(0), betaKc(0),
  theNodes(0), theElements(0), xEles(),
  currentGeoTag(0), lastGeoSendTag(-1), dbNod(0), dbEle(0)
{
    // does nothing
}

MeshRegion::~MeshRegion() 
{
  if (theNodes != 0)
    delete theNodes;
  if (theElements != 0)
    delete theElements;
}

int 
MeshRegion::setNodes(const ID &theNods)
{
  // destroy the old lists
  if (theNodes != 0)
    delete theNodes;
  if (theElements != 0)
    delete theElements;

  // create new element & node lists
  Domain *theDomain = this->getDomain();
  if (theDomain == 0) {
    opserr << "MeshRegion::setNodes() - no domain yet set\n";
    return -1;
  }

  int numNodes = theNods.Size();
  theNodes = new ID(0, numNodes);
  theElements = new ID(0, numNodes);
  if (theNodes == 0 || theElements == 0) {
    opserr << "MeshRegion::setNodes() - ran out of memory\n";
    return -1;
  }

  // add nodes to the node list if in the domain
  int loc = 0;
  for (int i=0; i<numNodes; i++) {
    int nodeTag = theNods(i);
    Node *theNode = theDomain->getNode(nodeTag);
    if (theNode != 0) {
      if (theNodes->getLocation(nodeTag) < 0)
	(*theNodes)[loc++] = nodeTag;      
    }
  }

  // now loop over the ele to create the ele list
  // NOTE - ele added to list if all it's nodes are in the region
  loc = 0;

  ElementIter &theEles = theDomain->getElements();
  Element *theEle;

  // loop over all ele
  while ((theEle = theEles()) != 0) {

    int eleTag = theEle->getTag();
    
    // check to see if all external nodes in node list
    bool in = true;
    const ID &theEleNodes = theEle->getExternalNodes();
    int numNodes = theEleNodes.Size();

    for (int i=0; i<numNodes; i++) {
      int nodeTag = theEleNodes(i);
      if (theNodes->getLocation(nodeTag) < 0) {
	in = false;
	i = numNodes;
      }
    }
    
    // if they are all in the node list add the ele to ele list
    if (in == true) 
      (*theElements)[loc++] = eleTag;
  }

  return 0;
}

int
MeshRegion::setNodesOnly(const ID &theNods)
{
    // destroy the old node list
    if (theNodes != 0)
        delete theNodes;
    
    // create new node list
    Domain *theDomain = this->getDomain();
    if (theDomain == 0) {
        opserr << "MeshRegion::setNodesOnly() - no domain yet set\n";
        return -1;
    }
    
    int numNodes = theNods.Size();
    theNodes = new ID(0, numNodes);
    if (theNodes == 0) {
        opserr << "MeshRegion::setNodesOnly() - ran out of memory\n";
        return -1;
    }
    
    // add nodes to the node list if in the domain
    int loc = 0;
    for (int i = 0; i<numNodes; i++) {
        int nodeTag = theNods(i);
        Node *theNode = theDomain->getNode(nodeTag);
        if (theNode != 0) {
            if (theNodes->getLocation(nodeTag) < 0)
                (*theNodes)[loc++] = nodeTag;
        }
    }
    
    return 0;
}

int 
MeshRegion::setElements(const ID &theEles)
{
  // destroy the old lists
  if (theNodes != 0)
    delete theNodes;
  if (theElements != 0)
    delete theElements;

  // create new element & node lists
  int numEle = theEles.Size();

  theElements = new ID(0, numEle); // don't copy yet .. make sure ele in domain
  theNodes = new ID(0, numEle); // initial guess at size of ID
  if (theElements == 0 || theNodes == 0) {
    opserr << "MeshRegion::setElements() - ran out of memory\n";
    return -1;
  }

  // now loop over the elements in ele ID passed in to create the node & ele list
  // NOTE - only add those elements to the list that are in the domain
  // NOTE - node added to region if any element has it as an external node
  int locEle = 0;
  int locNode = 0;

  Domain *theDomain = this->getDomain();
  if (theDomain == 0) {
    opserr << "MeshRegion::setElements() - no domain yet set\n";
    return -1;
  }

  Element *theEle;
  for (int i=0; i<numEle; i++) {
    int eleTag = theEles(i);
    theEle = theDomain->getElement(eleTag);
    if (theEle != 0) {

      if (theElements->getLocation(eleTag) < 0)
	  (*theElements)[locEle++] = eleTag;

      const ID &theEleNodes = theEle->getExternalNodes();

      for (int i=0; i<theEleNodes.Size(); i++) {
	int nodeTag = theEleNodes(i);
	// add the node tag if not already there
	if (theNodes->getLocation(nodeTag) < 0)
	  (*theNodes)[locNode++] = nodeTag;
      }
    }
  }

  return 0;
}

int
MeshRegion::setElementsOnly(const ID &theEles)
{
    // destroy the old element list
    if (theElements != 0)
        delete theElements;

    // create new element list
    Domain *theDomain = this->getDomain();
    if (theDomain == 0) {
        opserr << "MeshRegion::setElementsOnly() - no domain yet set\n";
        return -1;
    }

    int numEles = theEles.Size();
    theElements = new ID(0, numEles);
    if (theElements == 0) {
        opserr << "MeshRegion::setElementsOnly() - ran out of memory\n";
        return -1;
    }

    // add elements to the ele list if in the domain
    int loc = 0;
    for (int i = 0; i<numEles; i++) {
        int eleTag = theEles(i);
        Element *theEle = theDomain->getElement(eleTag);
        if (theEle != 0) {
            if (theElements->getLocation(eleTag) < 0)
                (*theElements)[loc++] = eleTag;
        }
    }

    return 0;
}

const ID &
MeshRegion::getNodes(void)
{
  if (theNodes == 0)
    opserr << "FATAL::MeshRegion::getNodes(void) - no nodes yet set\n";
  
  return *theNodes;
}

const ID &
MeshRegion::getElements(void)
{
  if (theElements == 0)
    opserr << "FATAL::MeshRegion::getElements(void) - no elements yet set\n";
  
  return *theElements;
}

int
MeshRegion::setRayleighDampingFactors(double alpham, double betak, double betak0, double betakc)
{
  alphaM = alpham;
  betaK  = betak;
  betaK0 = betak0;
  betaKc = betakc;

  // now set the damping factors at the nodes & elements
  Domain *theDomain = this->getDomain();
  if (theDomain == 0) {
    opserr << "MeshRegion::setRayleighDampingFactors() - no domain yet set\n";
    return -1;
  }
  if (theElements != 0) {
    for (int i=0; i<theElements->Size(); i++) {
      int eleTag = (*theElements)(i);
      Element *theEle = theDomain->getElement(eleTag);
      if (theEle != 0)
	theEle->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
    }
  }

  if (theNodes != 0) {
    for (int i=0; i<theNodes->Size(); i++) {
      int nodTag = (*theNodes)(i);
      Node *theNode = theDomain->getNode(nodTag);
      if (theNode != 0)
	theNode->setRayleighDampingFactor(alphaM);
    }
  }

  return 0;
}

int
MeshRegion::setDamping(Damping *theDamping)
{
  // now set the damping factors at the nodes & elements
  Domain *theDomain = this->getDomain();
  if (theDomain == 0) {
    opserr << "MeshRegion::setDamping() - no domain yet set\n";
    return -1;
  }
  if (theElements != 0) {
    for (int i=0; i<theElements->Size(); i++) {
      int eleTag = (*theElements)(i);
      Element *theEle = theDomain->getElement(eleTag);
      if (theEle != 0 && theEle->setDamping(theDomain, theDamping)) {
        opserr << "MeshRegion::setDamping - failed to set damping for " << theEle->getClassType() << " Element #" << eleTag << endln;
      }
    }
  }

  return 0;
}

int 
MeshRegion::sendSelf(int commitTag, Channel &theChannel)
{
  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of LoadPattern
  ID lpData(6);
  lpData(0) = currentGeoTag;
  lpData(1) = this->getTag();

  int numEle= theElements->Size();
  int numNod = theNodes->Size();

  lpData(2) = numEle;
  lpData(3) = numNod;

  if (dbNod == 0) {
    dbNod = theChannel.getDbTag();
    dbEle = theChannel.getDbTag();
  } 

  lpData(4) = dbNod;
  lpData(5) = dbEle;

  if (theChannel.sendID(myDbTag, commitTag, lpData) < 0) {
   opserr << "MeshRegion::sendSelf - channel failed to send the initial ID\n";
    return -1;
  }

  // now check if data defining the objects in the LoadPAttern needs to be sent 
  // NOTE THIS APPROACH MAY NEED TO CHANGE FOR VERY LARGE PROBLEMS IF CHANNEL CANNOT
  // HANDLE VERY LARGE ID OBJECTS.
  if (lastGeoSendTag != currentGeoTag) {
    if (numNod != 0) 
      if (theChannel.sendID(dbNod, currentGeoTag, *theNodes) < 0) {
	opserr << "MeshRegion::sendSelf - channel failed to send the nodes\n";
	return -1;
      }
    if (numEle != 0) 
      if (theChannel.sendID(dbEle, currentGeoTag, *theElements) < 0) {
	opserr << "MeshRegion::sendSelf - channel failed to send the elements\n";
	return -1;
      }

    Vector dData(4);
    dData(0) = alphaM;
    dData(1) = betaK;
    dData(2) = betaK0;
    dData(3) = betaKc;

    if (theChannel.sendVector(dbEle, currentGeoTag, dData) < 0) {
     opserr << "MeshRegion::sendSelf - channel failed to send the elements\n";
      return -1;
    }  

    // set the lst send db tag so we don't have to send them again unless they change
    lastGeoSendTag = currentGeoTag;
  }    
  
  return 0;
}

int 
MeshRegion::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of MeshRegion
  ID lpData(6);

  if (theChannel.recvID(myDbTag, commitTag, lpData) < 0) {
   opserr << "MeshRegion::recvSelf - channel failed to recv the initial ID\n";
    return -1;
  }

  // only recv the nodes & ele if not current
  if (currentGeoTag != lpData(0)) {

    currentGeoTag = lpData(0);
    this->setTag(lpData(1));

    int numEle = lpData(2);
    int numNod = lpData(3);
    if (theNodes != 0) {
      delete theNodes;
      theNodes = 0;
    }
    if (theElements != 0) {
      delete theElements;      
      theElements = 0;
    }

    if (numEle != 0)
      theElements = new ID(numEle); 
    if (numNod != 0)
      theNodes = new ID(numNod); 

    if (numNod != 0) 
      if (theChannel.recvID(dbNod, currentGeoTag, *theNodes) < 0) {
	opserr << "MeshRegion::sendSelf - channel failed to recv the nodes\n";
	return -1;
      }
    if (numEle != 0) 
      if (theChannel.recvID(dbEle, currentGeoTag, *theElements) < 0) {
	opserr << "MeshRegion::sendSelf - channel failed to recv the elements\n";
	return -1;
      }

    Vector dData(4);

    if (theChannel.sendVector(dbEle, currentGeoTag, dData) < 0) {
     opserr << "MeshRegion::sendSelf - channel failed to send the elements\n";
      return -1;
    }  
    alphaM = dData(0);
    betaK  = dData(1);
    betaK0 = dData(2);
    betaKc = dData(3);
  }

  this->setRayleighDampingFactors(alphaM, betaK, betaK0, betaKc);
  return 0;
}

void 
MeshRegion::Print(OPS_Stream &s, int flag)
{
  s << "Region: " << this->getTag() << endln;
  if (theElements != 0)
    s << "Elements: "<< *theElements;
  if (theNodes != 0)
    s << "Nodes: "<< *theNodes;
  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0) {
    s << "rayleigh damping factors:: alphaM: " << alphaM << " betaK: ";
    s << betaK << " betaK0: " << betaK0 << endln;
  }
}
