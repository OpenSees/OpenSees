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
                                                                        
#include <classTags.h>
#include <ElementStateParameter.h>
#include <Node.h>
#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Channel.h>
#include <Message.h>

ElementStateParameter::ElementStateParameter(double value, 
					     const char **Argv, 
					     int Argc, 
					     int Flag, 
					     ID *theEle)
  :Parameter(0,PARAMETER_TAG_ElementStateParameter),
   currentValue(value),
   flag(Flag),
   argc(Argc), fromFree(1)
{
  if (theEle != 0)
    theEleIDs = new ID(*theEle);

  argv = new char *[argc];
  for (int i=0; i<argc; i++) {
    int length = strlen(Argv[i])+1;
    argv[i] = new char[length];
    strcpy(argv[i], Argv[i]);
  }
}


ElementStateParameter::ElementStateParameter()
  :Parameter(0,PARAMETER_TAG_ElementStateParameter),
   currentValue(0.0),
   theEleIDs(0), flag(0),
   argv(0), argc(0), fromFree(1)
{

}

ElementStateParameter::~ElementStateParameter()
{
  if (fromFree == 0) {
    if (argc != 0) {
      for (int i=0; i<argc; i++)
	delete argv[i];

      delete [] argv;
      
      if (theEleIDs != 0)
	delete theEleIDs;
    }
  }
}

void
ElementStateParameter::Print(OPS_Stream &s, int flag)  
{
  s << "ElementStateParameter, tag = " << this->getTag() << endln;
}

void
ElementStateParameter::setDomain(Domain *theDomain)
{
  Parameter *theParameter = new Parameter(0, 0, 0, 0);

  Element *theEle;
  ElementIter &theEles = theDomain->getElements();

  if (flag == 0) {

    //
    // setParameter on all ele in the domain
    //

    while (((theEle = theEles()) != 0)) {
      int theResult = theEle->setParameter((const char **)argv, argc, *theParameter);
      if (theResult != -1) {
	theParameter->update(currentValue);
	theParameter->clean();
      }
    }    
  } else if (flag == 1) {

    //
    // setParameter on all ele whose tags in theELeIDs
    //
    
    int numEle = theEleIDs->Size();
    for (int i=0; i<numEle; i++) {
      int eleTag = (*theEleIDs)(i);
      theEle = theDomain->getElement(eleTag);
      if (theEle != 0) {
	int theResult = theEle->setParameter((const char **)argv, argc, *theParameter);
	if (theResult != -1) {
	  theParameter->update(currentValue);
	  theParameter->clean();
	}	  
      }
    }
  } else {
    
    //
    // setParameter on all ele whose tags in range given by 2 tags in theEleIDs
    //
    
    int startEle = (*theEleIDs)(0);
    int endEle = (*theEleIDs)(1);
    while (((theEle = theEles()) != 0)) {
      int eleTag = theEle->getTag();
      if (eleTag >= startEle && eleTag <= endEle) {
	int theResult = theEle->setParameter((const char **)argv, argc, *theParameter);
	if (theResult != -1) {
	  theParameter->update(currentValue);
	  theParameter->clean();
	}
      }      
    }
  }    

  delete theParameter;

  return;
}
  



int 
ElementStateParameter::sendSelf(int commitTag, Channel &theChannel)
{
  static ID iData(3);
  iData(0) = flag;
  iData(1) = argc;
  if (theEleIDs != 0)
    iData(2) = theEleIDs->Size();
  else
    iData(2) = 0;

  theChannel.sendID(commitTag, 0, iData);

  static Vector dData(1);
  dData(0) = currentValue;
  theChannel.sendVector(commitTag, 0, dData);

  if (theEleIDs != 0)
    theChannel.sendID(commitTag, 0, *theEleIDs);

  ID argvData(argc);
  for (int j=0; j<argc; j++)
    argvData(j) = strlen(argv[j])+1;

  theChannel.sendID(commitTag, 0, argvData);      

  for (int j=0; j<argc; j++) {    
    Message theMessage((char *)argv[j], argvData(j));
    theChannel.sendMsg(commitTag, 0, theMessage);      
  }  

  return 0;
}

int 
ElementStateParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID iData(3);
  theChannel.recvID(commitTag, 0, iData);
  flag = iData(0);
  argc = iData(1);
  int numEle = iData(2);


  static Vector dData(1);
  theChannel.recvVector(commitTag, 0, dData);
  currentValue = dData(0);

  fromFree = 1;

  if (theEleIDs != 0) 
    delete theEleIDs;

  if (numEle != 0) {
    theEleIDs = new ID(numEle);
    theChannel.recvID(commitTag, 0, *theEleIDs);    
  } else
    theEleIDs = 0;

  ID argvData(argc);
  theChannel.recvID(commitTag, 0, argvData);      

  argv = new char *[argc];
  for (int j=0; j<argc; j++) {    
    int argLength = argvData[j];
    argv[j] = new char[argLength];
    if (argv[j] == 0) {
      opserr << "ElementRecorder::recvSelf() - out of memory\n";
      return -1;
    }
    Message theMessage((char *)argv[j], argLength);
    theChannel.recvMsg(commitTag, 0, theMessage);      
  }

  return 0;
}
