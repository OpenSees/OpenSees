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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-10-31 17:49:55 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/MatParameter.cpp,v $

// written: fmk

#include <classTags.h>
#include <MatParameter.h>
#include <DomainComponent.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Channel.h>
#include <Message.h>

MatParameter::MatParameter(int theTag, int materialTag, const char *materialParameterName)
:Parameter(theTag, PARAMETER_TAG_MatParameter),
 theParameterName(0), theMaterialTag(materialTag), theParameterID(-1)
{
  if (materialParameterName != 0) {
    theParameterName = new char[strlen(materialParameterName) +1];
    if (theParameterName == 0) {
      opserr << "MatParameter::MatParameter - out of memory for parameter: ";
      opserr << materialParameterName << endln;
    }
    strcpy(theParameterName, materialParameterName);
  }
}

MatParameter::MatParameter()
  :Parameter(), 
   theParameterName(0), theMaterialTag(0), theParameterID(-1)
{

}

MatParameter::~MatParameter()
{
  if (theParameterName != 0)
    delete [] theParameterName;
}

void
MatParameter::Print(OPS_Stream &s, int flag)  
{
  s << "MaterialParameter, tag = " << this->getTag() << endln;
}

void
MatParameter::setDomain(Domain *theDomain)  
{
  Element *theEle;
  ElementIter &theEles = theDomain->getElements();

  int theResult = -1;
  
  const char *theString[2]; 
  char materialIdTag[20];

  sprintf(materialIdTag,"%d",theMaterialTag);
  theString[0] = theParameterName;
  theString[1] = materialIdTag;

  // note because of the way this parameter is updated only need to find one in the domain
  while ((theEle = theEles()) != 0) {
      int result = theEle->setParameter(theString, 2, *this);
      if (result != -1)
	theResult = result;
  }
	
  if (theResult == -1) 
    opserr << "MatParameter::setDomain(Domain *theDomain) - NO RESULT\n";
}

int 
MatParameter::sendSelf(int commitTag, Channel &theChannel)
{
  static ID theData(3);
  theData[0] = this->getTag();
  theData[1] = theMaterialTag;
  if (theParameterName != 0)
    theData[2] = strlen(theParameterName);
  else
    theData[2] = 0;

  theChannel.sendID(commitTag, 0, theData);
  
  if (theParameterName != 0) {
    Message theMessage(theParameterName, strlen(theParameterName));
    theChannel.sendMsg(commitTag, 0, theMessage);
  }
  
  return 0;
}

int 
MatParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID theData(3);  
  theChannel.recvID(commitTag, 0, theData);
  this->setTag(theData[0]);
  theMaterialTag = theData[1];
  
  
  if (theData(2) != 0) {
    theParameterName = new char[theData(2)+1];
    theParameterName[theData(2)]='\0';
    Message theMessage(theParameterName, theData(2));
    theChannel.recvMsg(commitTag, 0, theMessage);	
	theParameterName[theData(2)+1]='\n';
  }
 
  return 0;
}

