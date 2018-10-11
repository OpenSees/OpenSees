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
// Description: This file contains the implementation for the LysmerVelocityLoader class.

#include <LysmerVelocityLoader.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector LysmerVelocityLoader::data(3);

LysmerVelocityLoader::LysmerVelocityLoader(int tag, int theElementTag, int dir_)
  : ElementalLoad(tag, LOAD_TAG_LysmerVelocityLoader, theElementTag), dir(dir_)
{
}

LysmerVelocityLoader::LysmerVelocityLoader()
  : ElementalLoad(LOAD_TAG_LysmerVelocityLoader), dir(0)
{
}

LysmerVelocityLoader::~LysmerVelocityLoader()
{
}

const Vector &
LysmerVelocityLoader::getData(int &type, double loadFactor)
{
	type = LOAD_TAG_LysmerVelocityLoader;

    data.Zero();
    data(dir -1) = loadFactor;

	return data;
}

int
LysmerVelocityLoader::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;
    static ID iddata(4);
    int dataTag = this->getDbTag();
    iddata(0) = this->getTag();
    iddata(1) = dataTag;
    iddata(2) = eleTag;
    iddata(3) = dir;

    res = theChannel.sendID(dataTag, commitTag, iddata);
    if (res < 0) {
        opserr << "WARNING LysmerVelocityLoader::sendSelf() - " << this->getTag() << " failed to send iddata\n";
        return res;
    }

    return res;
}

int
LysmerVelocityLoader::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;
    static ID iddata(4);
    int dataTag = this->getDbTag();

    res = theChannel.recvID(dataTag, commitTag, iddata);
    if (res < 0) {
        opserr << "WARNING LysmerVelocityLoader::recvSelf() - " << this->getTag() << " failed to receive iddata\n";
        return res;
    }
     this->setTag(iddata(0));
     eleTag = iddata(2);

    dir = iddata(3);

	return res;
}

void
LysmerVelocityLoader::Print(OPS_Stream &s, int flag)
{
	s << "LysmerVelocityLoader...";
	s << "  element acted on: " << eleTag << " dir = " << dir << endln;
}

