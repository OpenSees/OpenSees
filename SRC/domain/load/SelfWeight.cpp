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
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector SelfWeight::data(3);

SelfWeight::SelfWeight(int tag, double xf, double yf, double zf, int theElementTag)
  : ElementalLoad(tag, LOAD_TAG_SelfWeight, theElementTag),
  xFact(xf), yFact(yf), zFact(zf)
{
}

SelfWeight::SelfWeight()
  : ElementalLoad(LOAD_TAG_SelfWeight),
  xFact(0.0), yFact(0.0), zFact(0.0)
{
}

SelfWeight::~SelfWeight()
{
}

const Vector &
SelfWeight::getData(int &type, double loadFactor)
{
  	type = LOAD_TAG_SelfWeight;
    data(0) = xFact;
    data(1) = yFact;
    data(2) = zFact;

	return data;
}

int 
SelfWeight::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();

    static Vector vData(5);
    vData(0) = xFact;
    vData(1) = yFact;
    vData(2) = zFact;
    vData(3) = eleTag;
    vData(4) = this->getTag();

    int result = theChannel.sendVector(dbTag, commitTag, vData);
    if (result < 0) {
        opserr << "SelfWeight::sendSelf - failed to send data\n";
        return result;
    }

    return 0;
}

int 
SelfWeight::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();

    static Vector vData(5);

    int result = theChannel.recvVector(dbTag, commitTag, vData);
    if (result < 0) {
        opserr << "SelfWeight::recvSelf - failed to recv data\n";
        return result;
    }

    this->setTag(vData(4));
    xFact = vData(0);
    yFact = vData(1);
    zFact = vData(2);
    eleTag = (int)vData(3);

    return 0;
}

void 
SelfWeight::Print(OPS_Stream &s, int flag)
{
	s << "SelfWeight...";
	s << "  element acted on: " << eleTag << endln;;
	s << "  (xFact, yFact, zFact) = (" << xFact << ", " << yFact << ", " << zFact << ") " << endln;;
}
