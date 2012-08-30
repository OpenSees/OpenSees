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
#include <InitialStateParameter.h>
#include <DomainComponent.h>

#include <Domain.h>
#include <Element.h>
#include <ElementIter.h>
#include <Channel.h>

InitialStateParameter::InitialStateParameter(bool stateFlag)
  :Parameter(0, PARAMETER_TAG_InitialStateParameter), theDomain(0)
{
  if (stateFlag == true)
    flag = 1;
  else
    flag = 0;
}

InitialStateParameter::InitialStateParameter()
  :Parameter(0, PARAMETER_TAG_InitialStateParameter), 
   flag(0)
{

}

InitialStateParameter::~InitialStateParameter()
{

}

void
InitialStateParameter::Print(OPS_Stream &s, int flag)  
{

}

void
InitialStateParameter::setDomain(Domain *theDomain)  
{
  if (flag == 1)
    ops_InitialStateAnalysis = true;
  else
    ops_InitialStateAnalysis = false;    

  return;
}

int 
InitialStateParameter::sendSelf(int commitTag, Channel &theChannel)
{
  static ID theData(2);
  theData[0] = this->getTag();
  theData[1] = flag;
  theChannel.sendID(commitTag, 0, theData);

  return 0;
}

int 
InitialStateParameter::recvSelf(int commitTag, 
				Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
  static ID theData(2);  
  theChannel.recvID(commitTag, 0, theData);

  this->setTag(theData[0]);
  flag = theData[1];

  theData[0] = -1;
  theData[1] = -2;

  return 0;
}

