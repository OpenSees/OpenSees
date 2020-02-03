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
                                                                        
// $Revision: 1.19 $
// $Date: 2010-06-09 17:43:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/pattern/LoadPattern.cpp,v $
                                                                        
// Written: fmk 07/99
// Revised:
//
// Purpose: This file contains the method definitions for class LoadPattern.
// LoadPattern is a container class.
//
// The LoadPattern interface:

#include <string.h>
#include <string>
#include <stdlib.h>

#include <LoadPattern.h>
#include <stdlib.h>
#include <ID.h>
#include <TimeSeries.h>
#include <NodalLoad.h>
#include <ElementalLoad.h>
#include <SP_Constraint.h>
#include <MapOfTaggedObjects.h>
#include <ElementalLoadIter.h>
#include <NodalLoadIter.h>
#include <SingleDomSP_Iter.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <GroundMotion.h>

#include <OPS_Globals.h>
#include <elementAPI.h>

void* OPS_LoadPattern()
{
    if(OPS_GetNumRemainingInputArgs() < 2) {
	opserr<<"insufficient number of args\n";
	return 0;
    }

    LoadPattern *thePattern = 0;

    // get tags
    int tags[2];
    int numData = 2;
    if(OPS_GetIntInput(&numData, &tags[0]) < 0) {
	opserr << "WARNING failed to get load pattern tag\n";
	return 0;
    }

    // get factor
    double fact = 1.0;
    if(OPS_GetNumRemainingInputArgs() > 1) {
	std::string type = OPS_GetString();
	if(type=="-fact" || type=="-factor") {
	    numData = 1;
	    if(OPS_GetDoubleInput(&numData,&fact) < 0) {
		opserr << "WARNING failed to get load pattern factor\n";
		return 0;
	    }
	}
    }

    // create pattern
    thePattern = new LoadPattern(tags[0], fact);
    TimeSeries *theSeries = OPS_getTimeSeries(tags[1]);

    // check 
    if(thePattern == 0 || theSeries == 0) {

	if(thePattern == 0) {
	    opserr << "WARNING - out of memory creating LoadPattern \n";
	} else {
	    opserr << "WARNING - problem creating TimeSeries for LoadPattern \n";
	}

	// clean up the memory and return an error
	if(thePattern != 0)
	    delete thePattern;
	return 0;
    }
    
    thePattern->setTimeSeries(theSeries);

    return thePattern;
}

LoadPattern::LoadPattern(int tag, int clasTag, double fact)
:DomainComponent(tag,clasTag),
 isConstant(1), loadFactor(0), scaleFactor(fact),
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0), lastChannel(0)
{
    // constructor for subclass
    theNodalLoads = new MapOfTaggedObjects();
    theElementalLoads = new MapOfTaggedObjects();
    theSPs = new MapOfTaggedObjects();

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }  
    // AddingSensitivity:BEGIN /////////////////////////////
    randomLoads = 0;
    dLambdadh = 0;
    // AddingSensitivity:END ///////////////////////////////
}


LoadPattern::LoadPattern()
:DomainComponent(0,PATTERN_TAG_LoadPattern),
 isConstant(1), loadFactor(0), scaleFactor(1.0),
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 dbSPs(0), dbNod(0), dbEle(0), 
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0), lastChannel(0)
{
    theNodalLoads = new MapOfTaggedObjects();
    theElementalLoads = new MapOfTaggedObjects();
    theSPs = new MapOfTaggedObjects();

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }
    // AddingSensitivity:BEGIN /////////////////////////////
    randomLoads = 0;
    dLambdadh = 0;
    // AddingSensitivity:END ///////////////////////////////
}


LoadPattern::LoadPattern(int tag, double fact)
:DomainComponent(tag,PATTERN_TAG_LoadPattern),
 isConstant(1), loadFactor(0.), scaleFactor(fact),
 theSeries(0), 
 currentGeoTag(0), lastGeoSendTag(-1),
 dbSPs(0), dbNod(0), dbEle(0), 
 theNodalLoads(0), theElementalLoads(0), theSPs(0),
 theNodIter(0), theEleIter(0), theSpIter(0), lastChannel(0)
{
    theNodalLoads = new MapOfTaggedObjects();
    theElementalLoads = new MapOfTaggedObjects();
    theSPs = new MapOfTaggedObjects();

    if (theNodalLoads == 0 || theElementalLoads == 0 || theSPs == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }    

    theEleIter = new ElementalLoadIter(theElementalLoads);    
    theNodIter = new NodalLoadIter(theNodalLoads);
    theSpIter = new SingleDomSP_Iter(theSPs);
    
    if (theEleIter == 0 || theNodIter == 0 || theSpIter == 0) {
	opserr << " LoadPattern::LoadPattern() - ran out of memory\n";
	exit(-1);
    }
    // AddingSensitivity:BEGIN /////////////////////////////
    randomLoads = 0;
    dLambdadh = 0;
    // AddingSensitivity:END ///////////////////////////////
}

    
// ~LoadPattern()
//	destructor

LoadPattern::~LoadPattern()
{
    if (theSeries != 0)
	delete theSeries;
    
    if (theNodalLoads != 0)
      delete theNodalLoads;

    if (theElementalLoads != 0)
      delete theElementalLoads;

    if (theSPs != 0)
      delete theSPs;

    if (theEleIter != 0)
      delete theEleIter;

    if (theNodIter != 0)
      delete theNodIter;

    if (theSpIter != 0)
      delete theSpIter;

    // AddingSensitivity:BEGIN /////////////////////////////
    if (randomLoads != 0)
      delete randomLoads;
    if (dLambdadh != 0)
      delete dLambdadh;
    // AddingSensitivity:END ///////////////////////////////
}


void
LoadPattern::setTimeSeries(TimeSeries *theTimeSeries)
{
    // invoke the destructor on the old TimeSeries
    if (theSeries != 0)
	delete theSeries;

    // set the pointer to the new series object
    theSeries = theTimeSeries;
}

void
LoadPattern::setDomain(Domain *theDomain)
{
    // if subclass does not implement .. check for 0 pointer
    if (theNodalLoads != 0) {
	NodalLoad *nodLoad;
	NodalLoadIter &theNodalIter = this->getNodalLoads();
	while ((nodLoad = theNodalIter()) != 0)
	    nodLoad->setDomain(theDomain);
    
	ElementalLoad *eleLoad;
	ElementalLoadIter &theElementalIter = this->getElementalLoads();
	while ((eleLoad = theElementalIter()) != 0)
	    eleLoad->setDomain(theDomain);    

	SP_Constraint *theSP;
	SP_ConstraintIter &theSpConstraints = this->getSPs();
	while ((theSP = theSpConstraints()) != 0)
	    theSP->setDomain(theDomain);
    }

    // now we set this load patterns domain
    this->DomainComponent::setDomain(theDomain);
}



bool
LoadPattern::addNodalLoad(NodalLoad *load)
{
    Domain *theDomain = this->getDomain();    

    bool result = theNodalLoads->addComponent(load);
    if (result == true) {
	if (theDomain != 0)
	    load->setDomain(theDomain);
	load->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else  
	opserr << "WARNING: LoadPattern::addNodalLoad() - load could not be added\n";

    return result;
}

bool
LoadPattern::addElementalLoad(ElementalLoad *load)
{
    Domain *theDomain = this->getDomain();
    
    bool result = theElementalLoads->addComponent(load);
    if (result == true) {
	if (theDomain != 0)
	    load->setDomain(theDomain);
	load->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else
	opserr << "WARNING: LoadPattern::addElementalLoad() - load could not be added\n";	
    
    return result;
}

bool
LoadPattern::addSP_Constraint(SP_Constraint *theSp)
{
    Domain *theDomain = this->getDomain();
    
    bool result = theSPs->addComponent(theSp);
    if (result == true) {
	if (theDomain != 0)
	    theSp->setDomain(theDomain);
	theSp->setLoadPatternTag(this->getTag());
	currentGeoTag++;
    } else
	opserr << "WARNING: LoadPattern::addSP_Constraint() - load could not be added\n";
    return result;
}


NodalLoadIter  &
LoadPattern::getNodalLoads(void)
{
    theNodIter->reset();
    return *theNodIter;
}
    
ElementalLoadIter &
LoadPattern::getElementalLoads(void)  
{
    theEleIter->reset();
    return *theEleIter;    
}

SP_ConstraintIter &
LoadPattern::getSPs(void)  
{
    theSpIter->reset();
    return *theSpIter;    
}

void
LoadPattern::clearAll(void)
{
    theElementalLoads->clearAll();
    theNodalLoads->clearAll();
    theSPs->clearAll();
    currentGeoTag++;
    lastChannel = 0;
    if (dLambdadh != 0) {
      dLambdadh->Zero();
    }
}

NodalLoad *
LoadPattern::removeNodalLoad(int tag)
{
    TaggedObject *obj = theNodalLoads->removeComponent(tag);
    if (obj == 0)
	return 0;
    NodalLoad *result = (NodalLoad *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;
}

ElementalLoad *
LoadPattern::removeElementalLoad(int tag)
{
    TaggedObject *obj = theElementalLoads->removeComponent(tag);
    if (obj == 0)
	return 0;

    ElementalLoad *result = (ElementalLoad *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;    
}    

SP_Constraint *
LoadPattern::removeSP_Constraint(int tag)
{
    TaggedObject *obj = theSPs->removeComponent(tag);
    if (obj == 0)
	return 0;
    SP_Constraint *result = (SP_Constraint *)obj;
    result->setDomain(0);
    currentGeoTag++;
    return result;    
}    


void
LoadPattern::applyLoad(double pseudoTime)
{
  // first determine the load factor
  if (theSeries != 0 && isConstant != 0) {
    loadFactor = theSeries->getFactor(pseudoTime);
    loadFactor *= scaleFactor;
  }

  NodalLoad *nodLoad;
  NodalLoadIter &theNodalIter = this->getNodalLoads();

  while ((nodLoad = theNodalIter()) != 0)
    nodLoad->applyLoad(loadFactor);
    
  ElementalLoad *eleLoad;
  ElementalLoadIter &theElementalIter = this->getElementalLoads();
  while ((eleLoad = theElementalIter()) != 0)
    eleLoad->applyLoad(loadFactor);

  SP_Constraint *sp;
  SP_ConstraintIter &theIter = this->getSPs();
  while ((sp = theIter()) != 0)
    sp->applyConstraint(loadFactor);
}

void
LoadPattern::setLoadConstant(void) 
{
  isConstant = 0;
}


void
LoadPattern::unsetLoadConstant(void) 
{
  isConstant = 1;
}

double
LoadPattern::getLoadFactor(void)
{
  if (theSeries != 0)
    return loadFactor;
  else
    return 0.0;
}

int
LoadPattern::sendSelf(int cTag, Channel &theChannel)
{
  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of LoadPattern
  int numNodLd, numEleLd, numSPs;
  ID lpData(11);

  numNodLd = theNodalLoads->getNumComponents();
  numEleLd = theElementalLoads->getNumComponents();
  numSPs = theSPs->getNumComponents();

  lpData(10) = this->getTag();
  lpData(0) = currentGeoTag;
  lpData(1) = numNodLd;
  lpData(2) = numEleLd;
  lpData(3) = numSPs;

  if (dbNod == 0) {
    dbNod = theChannel.getDbTag();
    dbEle = theChannel.getDbTag();
    dbSPs = theChannel.getDbTag();
  } 

  lpData(4) = dbNod;
  lpData(5) = dbEle;
  lpData(6) = dbSPs;

  lpData(7) = isConstant;

  if (theSeries != 0) {
    int dbtag = theSeries->getDbTag();
    int classtag = theSeries->getClassTag();
    if (dbtag == 0) {
      dbtag = theChannel.getDbTag();
      theSeries->setDbTag(dbtag);
    }
    lpData(8) = classtag;
    lpData(9) = dbtag;
  } else
    lpData(8) = -1;


  // see if we can save sending the vector containing just the load factor
  // will happen in parallel if sending the loadPattern .. not in database

  if (theChannel.sendID(myDbTag, cTag, lpData) < 0) {
    opserr << "LoadPattern::sendSelf - channel failed to send the initial ID\n";
    return -1;
  }    
  
 
    Vector data(2);
    data(0) = loadFactor;
    data(1) = scaleFactor;
    if (theChannel.sendVector(myDbTag, cTag, data) < 0) {
      opserr << "LoadPattern::sendSelf - channel failed to send the Vector\n";
      return -2;
  

  }

  if (theSeries != 0)
    if (theSeries->sendSelf(cTag, theChannel) < 0) {
      opserr << "LoadPattern::sendSelf - the TimeSeries failed to send\n";
      return -3;
    }

  // now check if data defining the objects in the LoadPAttern needs to be sent 
  // NOTE THIS APPROACH MAY NEED TO CHANGE FOR VERY LARGE PROBLEMS IF CHANNEL CANNOT
  // HANDLE VERY LARGE ID OBJECTS.

  /*
  if (theChannel.isDatastore() == 1) {
    static ID theLastSendTag(1);
    if (theChannel.recvID(myDbTag,0,theLastSendTag) == 0)
      lastGeoSendTag = theLastSendTag(0);
    else
      lastGeoSendTag = -1;
  }
  */

  if (lastChannel != theChannel.getTag() || lastGeoSendTag != currentGeoTag || theChannel.isDatastore() == 0) {

    lastChannel = theChannel.getTag();

    //
    // into an ID we are gonna place the class and db tags for each node so can rebuild
    // this ID we then send to the channel
    //

    // create the ID and get the node iter
    if (numNodLd != 0) {
      ID nodeData(numNodLd*2);
      NodalLoad *theNode;
      NodalLoadIter &theNodes = this->getNodalLoads();
      int loc =0;

      // loop over nodes in domain adding their classTag and dbTag to the ID
      while ((theNode = theNodes()) != 0) {
	nodeData(loc) = theNode->getClassTag();
	int dbTag = theNode->getDbTag();
	
	// if dbTag still 0 get one from Channel; 
	// if this tag != 0 set the dbTag in node
	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theNode->setDbTag(dbTag);
	}
	
	nodeData(loc+1) = dbTag;
	loc+=2;
      }    

      // now send the ID
      if (theChannel.sendID(dbNod, currentGeoTag, nodeData) < 0) {
	opserr << "LoadPattern::sendSelf - channel failed to send the NodalLoads ID\n";
	return -4;
      }
    }

    // we do the same for elemental loads as we did for nodal loads above .. see comments above!

    if (numEleLd != 0) {
      ID elementData(numEleLd*2);
      ElementalLoad *theEle;
      ElementalLoadIter &theElements = this->getElementalLoads();
      int loc = 0;
    
      while ((theEle = theElements()) != 0) {
	elementData(loc) = theEle->getClassTag();
	int dbTag = theEle->getDbTag();

	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theEle->setDbTag(dbTag);
	}
      
	elementData(loc+1) = dbTag;
	loc+=2;
      }

      // now send the ID
      if (theChannel.sendID(dbEle, currentGeoTag, elementData) < 0) {
	opserr << "Domain::send - channel failed to send the element ID\n";
	return -5;
      }
    }

    // we do the same for SP_Constraints as for NodalLoads above .. see comments above!
    
    if (numSPs != 0) {
      ID spData(numSPs*2);
      SP_Constraint *theSP;
      SP_ConstraintIter &theSPs = this->getSPs();
      int loc = 0;
    
      while ((theSP = theSPs()) != 0) {
	spData(loc) = theSP->getClassTag();
	int dbTag = theSP->getDbTag();

	if (dbTag == 0 && myDbTag != 0) {// go get a new tag and setDbTag in ele if this not 0 
	  dbTag = theChannel.getDbTag();
	  if (dbTag != 0)
	    theSP->setDbTag(dbTag);
	}
	
	spData(loc+1) = dbTag;
	loc+=2;
      }    

      if (theChannel.sendID(dbSPs, currentGeoTag, spData) < 0) {
	opserr << "LoadPAttern::sendSelf - channel failed sending SP_Constraint ID\n";
	return -6;
      }
    }

    // set the lst send db tag so we don't have to do all that again
    lastGeoSendTag = currentGeoTag;
    if (theChannel.isDatastore() == 1) {
      static ID theLastSendTag(1);
      theLastSendTag(0) = lastGeoSendTag;
      theChannel.sendID(myDbTag,0, theLastSendTag);
    }
  }

  // now we invoke sendSelf on all the NodalLoads, ElementalLoads and SP_Constraints
  // which have been added to the LoadCase
  NodalLoad *theNode;
  NodalLoadIter &theNodes = this->getNodalLoads();
  while ((theNode = theNodes()) != 0) {
    if (theNode->sendSelf(cTag, theChannel) < 0) {
      opserr << "LoadPattern::sendSelf - node with tag " << theNode->getTag() << " failed in sendSelf\n";
      return -7;
    }
  }

  ElementalLoad *theEle;
  ElementalLoadIter &theElements = this->getElementalLoads();
  while ((theEle = theElements()) != 0) {
    if (theEle->sendSelf(cTag, theChannel) < 0) {
      opserr << "LoadPattern::sendSelf - element with tag " << theEle->getTag() << " failed in sendSelf\n";
      return -8;
    }
  }

  SP_Constraint *theSP;
  SP_ConstraintIter &theSPs = this->getSPs();
  while ((theSP = theSPs()) != 0) {
    if (theSP->sendSelf(cTag, theChannel) < 0) {
      
      opserr << "LoadPattern::sendSelf - SP_Constraint: " << *theSP << " failed sendSelf\n";
      return -9;
    }
  }    

  // if we get here we are successful
  return 0;
}



int
LoadPattern::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{

  // get my current database tag
  // NOTE - dbTag equals 0 if not sending to a database OR has not yet been sent
  int myDbTag = this->getDbTag();

  // into an ID we place all info needed to determine state of LoadPattern
  int numNod, numEle, numSPs;
  ID lpData(11);

  if (theChannel.recvID(myDbTag, cTag, lpData) < 0) {
    opserr << "LoadPattern::recvSelf - channel failed to recv the initial ID\n";
    return -1;
  }

  isConstant = lpData(7);

  this->setTag(lpData(10));

 
    Vector data(2);
    if (theChannel.recvVector(myDbTag, cTag, data) < 0) {
      opserr << "LoadPattern::recvSelf - channel failed to recv the Vector\n";
      return -2;
    }
    loadFactor = data(0);
    scaleFactor = data(1);
  
  
  // read data about the time series
  if (lpData(8) != -1) {
    if (theSeries == 0) {
      theSeries = theBroker.getNewTimeSeries(lpData(8));
    } else if (theSeries->getClassTag() != lpData(8)) {
      delete theSeries;    
      theSeries = theBroker.getNewTimeSeries(lpData(8));
    }
    if (theSeries == 0) {
      opserr << "LoadPattern::recvSelf - failed to create TimeSeries\n";
      return -3;
    }
  
    theSeries->setDbTag(lpData(9));

    if (theSeries->recvSelf(cTag, theChannel, theBroker) < 0) {
      opserr << "LoadPattern::recvSelf - the TimeSeries failed to recv\n";
      return -3;
    }
  }

  /*
  if (theChannel.isDatastore() == 1) {
    static ID theLastSendTag(1);
    if (theChannel.recvID(myDbTag,0,theLastSendTag) == 0)
      lastGeoSendTag = theLastSendTag(0);
  }
  */

  if (lastChannel != theChannel.getTag() || currentGeoTag != lpData(0) || theChannel.isDatastore() == 0) {

    // clear out the all the components in the current load pattern
    this->clearAll();
    lastChannel = theChannel.getTag();
    currentGeoTag = lpData(0);

    numNod = lpData(1);
    numEle = lpData(2);
    numSPs = lpData(3);
    dbNod = lpData(4);
    dbEle = lpData(5);
    dbSPs = lpData(6);    

    // 
    // now we rebuild the nodal loads
    //
    
    // first get the information from the domainData about the nodes
    if (numNod != 0) {
      ID nodeData(2*numNod);

      // now receive the ID about the nodes, class tag and dbTags
      if (theChannel.recvID(dbNod, currentGeoTag, nodeData) < 0) {
	opserr << "LoadPAttern::recvSelf - channel failed to recv the NodalLoad ID\n";
	return -2;
      }

      // now for each NodalLoad we 1) get a new node of the correct type from the ObjectBroker
      // 2) ensure the node exists and set it's dbTag, 3) we invoke recvSelf on this new 
      // blank node and 4) add this node to the domain

      int loc = 0;

      for (int i=0; i<numNod; i++) {
	int classTag = nodeData(loc);
	int dbTag = nodeData(loc+1);
	
	NodalLoad *theNode = theBroker.getNewNodalLoad(classTag);

	if (theNode == 0) {
	  opserr << "LoadPattern::recv - cannot create NodalLoad with classTag " << classTag << endln;
	  return -2;
	}			
	
	theNode->setDbTag(dbTag);
	
	if (theNode->recvSelf(cTag, theChannel, theBroker) < 0) {
	  opserr << "LoadPattern::recvSelf - NodalLoad with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			

	if (this->addNodalLoad(theNode) == false) {
	  opserr << "LoadPattern::recvSelf - failed adding NodalLoad tagged " << theNode->getTag() << " into LP!\n";
	  return -3;
	}			
	  
	loc+=2;
      }   
    }

    // 
    // now we rebuild the ElementalLoads .. same as NodalLoads above .. see comments above
    //
    
    if (numEle != 0) {
      ID eleData(2*numEle);
      
      if (theChannel.recvID(dbEle, currentGeoTag, eleData) < 0) {
	opserr << "LoadPattern::recvSelf - channel failed to recv the EleLoad ID\n";
	return -2;
      }

      int loc = 0;
      for (int i=0; i<numEle; i++) {
	int classTag = eleData(loc);
	int dbTag = eleData(loc+1);
      
	ElementalLoad *theEle = theBroker.getNewElementalLoad(classTag);
	if (theEle == 0) {
	  opserr << "LoadPattern::recv - cannot create ElementalLoad with classTag " << classTag << endln;
	  return -2;
	}			

	theEle->setDbTag(dbTag);
	
	if (theEle->recvSelf(cTag, theChannel, theBroker) < 0) {
	  opserr << "LoadPattern::recvSelf - Ele with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			
	
	if (this->addElementalLoad(theEle) == false) {
	  opserr << "LoadPattern::recvSelf - could not add Ele with tag " << theEle->getTag() << " into LP!\n";
	  return -3;
	}			
	
	loc+=2;
      }
    }

    // 
    // now we rebuild the SP_Constraints .. same as nodes above .. see above if can't understand!!
    //
    
    if (numSPs != 0) {
      ID spData(2*numSPs);

      if (theChannel.recvID(dbSPs, currentGeoTag, spData) < 0) {
	opserr << "LoadPattern::recvSelf - channel failed to recv the SP_Constraints ID\n";
	return -2;
      }

      int loc = 0;
      for (int i=0; i<numSPs; i++) {
	int classTag = spData(loc);
	int dbTag = spData(loc+1);
      
	SP_Constraint *theSP = theBroker.getNewSP(classTag);
	if (theSP == 0) {
	  opserr << "LoadPattern::recv - cannot create SP_Constraint with classTag " << classTag << endln;
	  return -2;
	}			
	theSP->setDbTag(dbTag);
      
	if (theSP->recvSelf(cTag, theChannel, theBroker) < 0) {
	  opserr << "LoadPattern::recvSelf - SP_Constraint with dbTag " << dbTag << " failed in recvSelf\n";
	  return -2;
	}			
	
	if (this->addSP_Constraint(theSP) == false) {
	  opserr << "LoadPattern::recvSelf - could not add SP_Constraint with tag " << theSP->getTag()
		 << " into LP!\n";
				  
	  return -3;
	}			
	
	loc+=2;
      }
    }

    // now set the load pattern db count
    currentGeoTag = lpData(0);
    lastGeoSendTag  = currentGeoTag;

  } else {
    if (theSeries != 0)
      if (theSeries->recvSelf(cTag, theChannel, theBroker) < 0) {
	opserr << "LoadPattern::recvSelf - the TimeSeries failed to recv\n";
	return -3;
      }

    
    NodalLoad *theNode;
    NodalLoadIter &theNodes = this->getNodalLoads();
    while ((theNode = theNodes()) != 0) {
      if (theNode->recvSelf(cTag, theChannel, theBroker) < 0) {
	opserr << "LoadPattern::recvSelf - node with tag " << theNode->getTag() << " failed in recvSelf\n";
	return -7;
      }
    }

    ElementalLoad *theEle;
    ElementalLoadIter &theElements = this->getElementalLoads();
    while ((theEle = theElements()) != 0) {
      if (theEle->recvSelf(cTag, theChannel, theBroker) < 0) {
	opserr << "LoadPattern::recvSelf - element with tag " << theEle->getTag() << " failed in recvSelf\n";
	return -8;
      }
    }

    SP_Constraint *theSP;
    SP_ConstraintIter &theSPs = this->getSPs();
    while ((theSP = theSPs()) != 0) {
      if (theSP->recvSelf(cTag, theChannel, theBroker) < 0) {
	opserr << "LoadPattern::recvSelf - SP_Constraint tagged " << theSP->getTag() << "  failed recvSelf\n";
	return -9;
      }
    }    
  }

  // if we get here we are successful
  return 0;
}

void
LoadPattern::Print(OPS_Stream &s, int flag)
{
    s << "Load Pattern: " << this->getTag() << "\n";
    s << "  Scale Factor: " << scaleFactor << endln;
    if (theSeries != 0)
      theSeries->Print(s,flag);
    s << "  Nodal Loads: \n";
    theNodalLoads->Print(s,flag);
    s << "\n  Elemental Loads: \n";
    theElementalLoads->Print(s, flag);
    s << "\n  Single Point Constraints: \n";
    theSPs->Print(s, flag);
}


LoadPattern *
LoadPattern::getCopy(void)
{
  LoadPattern *theCopy = new LoadPattern(this->getTag());
  if (theCopy == 0) {
    opserr << "LoadPattern::getCopy() - ran out of memory\n";
    return theCopy; // in case fatal() does not exit
  }
  theCopy->loadFactor = loadFactor;
  theCopy->scaleFactor = scaleFactor;
  theCopy->isConstant = isConstant;
  theCopy->theSeries = theSeries;
  return theCopy;
}

int
LoadPattern::addMotion(GroundMotion &theMotion, int tag)
{
  opserr << "LoadPattern::addMotion() - cannot add GroundMotion - use MultiSupport Pattern instead\n";
  return -1;
}

GroundMotion *
LoadPattern::getMotion(int tag)
{
  return 0;
}


// AddingSensitivity:BEGIN ////////////////////////////////////
void
LoadPattern::applyLoadSensitivity(double pseudoTime)
{
    // P*dfactor/dh
    if (theSeries != 0 && isConstant != 0) {
	loadFactor = theSeries->getFactorSensitivity(pseudoTime);
	loadFactor *= scaleFactor;
    }
  
    NodalLoad *nodLoad;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    while ((nodLoad = theNodalIter()) != 0)
	nodLoad->applyLoad(loadFactor);

    // factor*dP/dh
    if (theSeries != 0 && isConstant != 0) {
	loadFactor = theSeries->getFactor(pseudoTime);
	loadFactor *= scaleFactor;
    }

    NodalLoadIter &theNodalIter2 = this->getNodalLoads();
    while ((nodLoad = theNodalIter2()) != 0)
	nodLoad->applyLoadSensitivity(loadFactor);
  
  // Don't include element loads and sp constraints for now
  /*
    ElementalLoad *eleLoad;
    ElementalLoadIter &theElementalIter = this->getElementalLoads();
    while ((eleLoad = theElementalIter()) != 0)
    eleLoad->applyLoad(loadFactor);
    
    SP_Constraint *sp;
    SP_ConstraintIter &theIter = this->getSPs();
    while ((sp = theIter()) != 0)
    sp->applyConstraint(loadFactor);
  */
}

int
LoadPattern::setParameter(const char **argv, int argc, Parameter &param)
{
    if (theSeries == 0) {
        opserr << "set/update/activate parameter is illegaly called in LoadPattern " << endln;
	return 0;
    }

    if (argc < 1)
      return -1;

    // Nodal load
    if (strstr(argv[0],"loadAtNode") != 0) {

      if (argc < 3)
	return -1;

        RVisRandomProcessDiscretizer = false;

        int nodeNumber = atoi(argv[1]);
        NodalLoad *thePossibleNodalLoad;
        NodalLoad *theNodalLoad = 0;
        NodalLoadIter &theNodalIter = this->getNodalLoads();

        while ((thePossibleNodalLoad = theNodalIter()) != 0) {
            if ( nodeNumber == thePossibleNodalLoad->getNodeTag() ) {
                theNodalLoad = thePossibleNodalLoad;
            }
        }

	if (theNodalLoad != 0)
	  return theNodalLoad->setParameter(&argv[2], argc-2, param);
	else
	  return -1;
    }

    else if (strstr(argv[0],"elementPointLoad") != 0 || strstr(argv[0],"elementLoad") != 0) {

      if (argc < 3)
	return -1;

      RVisRandomProcessDiscretizer = false;

      int eleNumber = atoi(argv[1]);
      ElementalLoad *theEleLoad = 0;
      ElementalLoadIter &theEleLoadIter = this->getElementalLoads();
      while ((theEleLoad = theEleLoadIter()) != 0) {
	int eleTag = theEleLoad->getElementTag();
	if (eleNumber == eleTag) {
	  return theEleLoad->setParameter(&argv[2], argc-2, param);
	}
      }

      return -1;
    }

    else if (strstr(argv[0],"randomProcessDiscretizer") != 0) {

      if (argc < 2)
	return -1;

        RVisRandomProcessDiscretizer = true;
        return theSeries->setParameter(&argv[1], argc-1, param);
    }

    // Unknown parameter
    else
      return -1;
}

int
LoadPattern::updateParameter(int parameterID, Information &info)
{
  if (theSeries == 0) {
    opserr << "set/update/activate parameter is illegaly called in LoadPattern " << endln;
  }
  
  opserr << "LoadPattern::updateParameter -- no parameters defined, this method should not be called" << endln;

  return 0;

  /*
  if (RVisRandomProcessDiscretizer) {
    return theSeries->updateParameter(parameterID,info);
  }
  else {
    NodalLoad *thePossibleNodalLoad = 0;
    NodalLoad *theNodalLoad = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    
    switch (parameterID) {
    case 1: case -1:  // Not implemented.
      return -1;
    default:
      if (parameterID > 1000  &&  parameterID < 2000)  {
	int nodeNumber = parameterID-1000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(1, info);
      }
      else if (parameterID > 2000  &&  parameterID < 3000)  {
	int nodeNumber = parameterID-2000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(2, info);
      }
      else if (parameterID > 3000  &&  parameterID < 4000)  {
	int nodeNumber = parameterID-3000;
	while ((thePossibleNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == thePossibleNodalLoad->getNodeTag() )  {
	    theNodalLoad = thePossibleNodalLoad;
	  }
	}
	return theNodalLoad->updateParameter(3, info);
            }
      else
	return -1;
    }
  }
  */
}





int
LoadPattern::activateParameter(int parameterID)
{
  if (theSeries == 0) {
    opserr << "set/update/activate parameter is illegaly called in LoadPattern " << endln;
  }
  
  opserr << "LoadPattern::activateParameter -- no parameters defined, this method should not be called" << endln;

  return 0;

  /*
  if (RVisRandomProcessDiscretizer) {
    return theSeries->activateParameter(parameterID);
  }
  else {
    
    
    // Don't set flag here in the load pattern itself.
    // (Assume there always may be random loads)
    
    NodalLoad *theNodalLoad = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    
    if (parameterID == 0) {
      
      // Go through all nodal loads and zero out gradientIdentifier
      // (Remember: the identifier is only zero if we are in
      // the process of zeroing out all sensitivity flags).
      while ((theNodalLoad = theNodalIter()) != 0)  {
	theNodalLoad->activateParameter(parameterID);
      }
      
    }
    else {
      
      // Find the right nodal load and set the flag
      if (parameterID > 1000  &&  parameterID < 2000)  {
	int nodeNumber = parameterID-1000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(1);
	  }
	}
      }
      else if (parameterID > 2000  &&  parameterID < 3000)  {
	int nodeNumber = parameterID-2000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(2);
	  }
	}
      }
      else if (parameterID > 3000  &&  parameterID < 4000)  {
	int nodeNumber = parameterID-3000;
	while ((theNodalLoad = theNodalIter()) != 0)  {
	  if ( nodeNumber == theNodalLoad->getNodeTag() )  {
	    theNodalLoad->activateParameter(3);
	  }
	}
      }
      else {
	opserr << "LoadPattern::gradient() -- error in identifier. " << endln;
      }
    }
  }
  return 0;
  */

}


const Vector &
LoadPattern::getExternalForceSensitivity(int gradNumber)
{

    // THIS METHOD IS CURRENTLY ONLY USED FOR THE STATIC CASE
    // IT SHOULD BE DELETED AND REPLACED BY THE DYNAMIC CASE

    // Initial declarations
    Vector tempRandomLoads(1);
    int sizeRandomLoads;

    // Start with a fresh return vector
    if (randomLoads == 0) {
        randomLoads = new Vector(1);
    }
    else {
        delete randomLoads;
        randomLoads = new Vector(1);
    }

    // Prepare the vector identifying which loads are random.
    NodalLoad *theNodalLoad = 0;
    NodalLoadIter &theNodalIter = this->getNodalLoads();
    int i;

    // Loop through the nodal loads to pick up possible contributions
    int nodeNumber;
    int dofNumber;
    while ((theNodalLoad = theNodalIter()) != 0)  {
        const Vector &gradientVector = theNodalLoad->getExternalForceSensitivity(gradNumber);
        if (gradientVector(0) != 0.0 ) {

            // Found a random load! Get nodeNumber and dofNumber
            nodeNumber = theNodalLoad->getNodeTag();
            dofNumber = (int)gradientVector(0);

            // Update the randomLoads vector
            sizeRandomLoads = randomLoads->Size();
            if (sizeRandomLoads == 1) {
                delete randomLoads;
                randomLoads = new Vector(2);
                (*randomLoads)(0) = (double)nodeNumber;
                (*randomLoads)(1) = (double)dofNumber;
            }
            else {
                tempRandomLoads = (*randomLoads);
                delete randomLoads;
                randomLoads = new Vector(sizeRandomLoads+2);
                for (i=0; i<sizeRandomLoads; i++) {
                    (*randomLoads)(i) = tempRandomLoads(i);
                }
                (*randomLoads)(sizeRandomLoads) = nodeNumber;
                (*randomLoads)(sizeRandomLoads+1) = dofNumber;
            }
        }
    }

    return (*randomLoads);
}

int
LoadPattern::saveLoadFactorSensitivity(double dlambdadh, int gradIndex, int numGrads)
{

  //opserr << "LoadPattern::savedlamdh " << gradIndex << ' ' << numGrads << endln;
  if (dLambdadh == 0) {
    dLambdadh = new Vector(numGrads);
  }

  if (dLambdadh == 0 || dLambdadh->Size() != numGrads) {
    if (dLambdadh != 0)
      delete dLambdadh;
    dLambdadh = new Vector(numGrads);
  }

  if (gradIndex >= 0 && gradIndex < numGrads) {
    (*dLambdadh)(gradIndex) = dlambdadh;
    return 0;
  }
  else {
    opserr << "LoadPattern::saveLoadFactorSensitivity -- gradIndex out of bounds" << endln;
    return -1;
  }
}

double
LoadPattern::getLoadFactorSensitivity(int gradIndex)
{
  if (dLambdadh != 0 && gradIndex >= 0 && gradIndex < dLambdadh->Size())
    return (*dLambdadh)(gradIndex);
  else
    return 0.0;
}

// AddingSensitivity:END //////////////////////////////////////
