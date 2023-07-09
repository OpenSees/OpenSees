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
                                                                        
// $Revision: 1.10 $
// $Date: 2009-08-25 23:26:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/component/MaterialStageParameter.cpp,v $

#include <classTags.h>
#include <MaterialStageParameter.h>
#include <DomainComponent.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Channel.h>

MaterialStageParameter::MaterialStageParameter(int theTag, int materialTag)
:Parameter(theTag, PARAMETER_TAG_MaterialStageParameter),
 theMaterialTag(materialTag)
{

}

MaterialStageParameter::MaterialStageParameter()
  :Parameter(), 
   theMaterialTag(0)
{

}

MaterialStageParameter::~MaterialStageParameter()
{

}

void
MaterialStageParameter::Print(OPS_Stream &s, int flag)  
{
  s << "MaterialStageParameter, tag = " << this->getTag() << endln;
}

void
MaterialStageParameter::setDomain(Domain *theDomain)  
{
  Element *theEle;
  ElementIter &theEles = theDomain->getElements();

  int theResult = -1;

  const char *theString[2];// = new const char*[2];
  char parameterName[21];
  char materialIdTag[10];
  sprintf(parameterName,"updateMaterialStage");
  sprintf(materialIdTag,"%d",theMaterialTag);
  theString[0] = parameterName;
  theString[1] = materialIdTag;

  // note because of the way this parameter is updated only need to find one in the domain
  while (((theEle = theEles()) != 0) && (theResult == -1)) {
    theResult = theEle->setParameter(theString, 2, *this);
  }

  if (theResult == -1)
    opserr << "WARNING: MaterialStageParameter::setDomain() - no effect with material tag " << theMaterialTag << endln;

  theResult = 0;

  return;
}

int 
MaterialStageParameter::sendSelf(int commitTag, Channel &theChannel)
{

  static ID theData(2);
  theData[0] = this->getTag();
  theData[1] = theMaterialTag;
  theChannel.sendID(commitTag, 0, theData);

  return 0;
}

int 
MaterialStageParameter::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static ID theData(2);  
  theChannel.recvID(commitTag, 0, theData);
  this->setTag(theData[0]);
  theMaterialTag = theData[1];
  return 0;
}

