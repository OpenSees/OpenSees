/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
                                                                        
// $Revision: 1.9 $
// $Date: 2008-08-26 15:43:43 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/ElementParameter.cpp,v $

#include <classTags.h>
#include <ElementParameter.h>
#include <DomainComponent.h>
#include <Domain.h>
#include <Element.h>
#include <Channel.h>
#include <Message.h>

ElementParameter::ElementParameter(int passedTag,
				   int eleTag,
				   const char **theArgv, int theArgc)
  :Parameter(passedTag, PARAMETER_TAG_ElementParameter),
   eleTags(1), numChannels(0), theChannels(0)
{
  argc = theArgc;
  argv = 0;
  argvSize = 0;
  
  if (argc != 0) {
    argv = new char *[argc];
    for (int i=0; i< argc; i++) {
      int lengthArgvI = strlen(theArgv[i]);
      argvSize += lengthArgvI + 1;
    }

    argv[0] = new char [argvSize];
    strcpy(argv[0], theArgv[0]);
    argvSize = strlen(theArgv[0]) + 1;

    for (int i=1; i< argc; i++) {
      int lengthArgvI = strlen(theArgv[i-1]);
      argv[i] = argv[i-1]+lengthArgvI+1;
      strcpy(argv[i], theArgv[i]);
      argvSize += lengthArgvI + 1;
    }
  }

  eleTags[0] = eleTag;
}


ElementParameter::ElementParameter()
  :Parameter(-1, PARAMETER_TAG_ElementParameter),
   numChannels(0), theChannels(0)
{
  argc = 0;
  argv = 0;
  argvSize = 0;
}

ElementParameter::~ElementParameter()
{
  if (argv != 0)
    delete [] argv[0]; // stored in 1 array

  delete [] argv;

  if (theChannels != 0)
    delete [] theChannels;
}


int 
ElementParameter::update(int newValue)
{
  return this->Parameter::update(newValue);
}

int 
ElementParameter::update(double newValue)
{
	/* THIS DOES NOT WORK, NOT ALL PARAMETERS ARE WORKING THE SAME WAY!
  if (numChannels != 0) {
    // due to stupid design of parameter class (aka addComponnet)
    // we have to always check if eleTags has changed
    static ID idData(1);
    if (numChannels > 0) {
      idData(0) = eleTags.Size();
      for (int i=0; i<numChannels; i++) {
	theChannels[i]->sendID(0,0,idData);
	theChannels[i]->sendID(0,0,eleTags);
      }
    } else {

      theChannels[0]->recvID(0,0,idData);      
      int numEle = idData(0);
      if (numEle == eleTags.Size())
	theChannels[0]->recvID(0,0,eleTags);
      else {
	eleTags.resize(numEle);
	theChannels[0]->recvID(0,0,eleTags);
	this->setDomain(theDomain);
      }
    }
  }
  */
  return this->Parameter::update(newValue);
}

int
ElementParameter::addComponent(int eleTag, const char **theArgv, int theArgc)
{
	opserr << "elementParameter::addComponent - hopefully not called\n";
  int numEle = eleTags.Size();
  eleTags[numEle] = eleTag;

  if (theDomain != 0) {
    Element *theEle = theDomain->getElement(eleTag);
    if (theEle != 0)
      return this->Parameter::addComponent(theEle, theArgv, theArgc);
  }

  if (argc != theArgc) {
    opserr << "ElementParameter::addComponent(int eleTag) " << eleTag << " argc passed differ from stored, won't work in SP\n";
    return 0;
  }
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i], theArgv[i]) != 0) {
      opserr << "ElementParameter::addComponent(int eleTag) " << eleTag << " argc passed differ from stored, won't work in SP\n";
    }
  }
  return 0;
}

void
ElementParameter::setDomain(Domain *aDomain)
{
  theDomain = aDomain;
  this->Parameter::clean();
 
  const char **theArgv = (const char **)argv;
  int numEle = eleTags.Size();
  for (int i=0; i<numEle; i++) {
    int eleTag = eleTags[i];
    Element *theEle = theDomain->getElement(eleTag);
    if (theEle != 0) {
      this->Parameter::addComponent(theEle, theArgv, argc);
	}
  }
}

int 
ElementParameter::sendSelf(int commitTag, Channel &theChannel)
{
  ID idData(4);
  idData(0) = this->getTag();
  idData(1) = eleTags.Size();
  idData(2) = argvSize;
  idData(3) = argc; 

  if (theChannel.sendID(0, commitTag, idData) < 0) {

  }

  if (theChannel.sendID(0, commitTag, eleTags) < 0) {

  }
  /*
  ID argData(argc-1);

  for (int i=0; i<argc-1; i++) {
    argData(i) = argv[i+1]-argv[i];
  }

  if (theChannel.sendID(0, commitTag, argData) < 0) {

  }
  */
  Message msgData(argv[0], argvSize);
  if (theChannel.sendMsg(0, commitTag, msgData) < 0) {

  }
  

  //
  // add the channel to theChanels array
  //

  Channel **newChannels = new Channel *[numChannels+1];
  for (int i=0; i<numChannels; i++)
    newChannels[i] = theChannels[i];
  
  newChannels[numChannels] = &theChannel;
  numChannels++;
  
  if (theChannels != 0)
    delete [] theChannels;
  theChannels = newChannels;
  return 0;
}

int 
ElementParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  ID idData(4);

  if (theChannel.recvID(0, commitTag, idData) < 0) {

  }

  this->setTag(idData(0));

  eleTags.resize(idData(1));

  if (theChannel.recvID(0, commitTag, eleTags) < 0) {

  }

  if (argv != 0) {
    if (argv[0] != 0)
      delete [] argv[0];

    delete [] argv;
  }

  argc = idData(3);
  argvSize = idData(2);
 
  /*
  ID argData(argc-1);
  if (theChannel.recvID(0, commitTag, argData) < 0) {

  }
  */

  argv = new char *[argc];
  argv[0] = new char[argvSize];
  Message msgData(argv[0], argvSize);
  if (theChannel.recvMsg(0, commitTag, msgData) < 0) {

  } 
  
  for (int i=0; i<argc-1; i++)
    argv[i+1] = argv[i] + strlen(argv[i])+1;

  if (theChannels != 0)
    delete [] theChannels;

  theChannels = new Channel *[1];
  theChannels[0] = &theChannel;
  numChannels = -1;
 
  return 0;
}

