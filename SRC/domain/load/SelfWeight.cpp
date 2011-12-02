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
                                                                        
// Written: Chris McGann, U.Washington
//          02.2011
//
// Description: This file contains the implementation for the SelfWeight class.

#include <SelfWeight.h>
#include <Vector.h>

Vector SelfWeight::data(1);

SelfWeight::SelfWeight(int tag, int theElementTag)
  : ElementalLoad(tag, LOAD_TAG_SelfWeight, theElementTag)
{
}

SelfWeight::SelfWeight()
  : ElementalLoad(LOAD_TAG_SelfWeight)
{
}

SelfWeight::~SelfWeight()
{
}

const Vector &
SelfWeight::getData(int &type, double loadFactor)
{
  	type = LOAD_TAG_SelfWeight;

	return data;
}

int 
SelfWeight::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int 
SelfWeight::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
	return -1;
}

void 
SelfWeight::Print(OPS_Stream &s, int flag)
{
	s << "SelfWeight...";
	s << "  element acted on: " << eleTag << endln;;
}
