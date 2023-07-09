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
                                                                        
// Written: Felipe Elgueta and Jose Abell
//          01.2021
//
// Description: This file contains the implementation for IGAFollowerLoad, a load class for 
//               applying point forces inside load patterns for IGA analysis.


#include <IGAFollowerLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector IGAFollowerLoad::data(5);

IGAFollowerLoad::IGAFollowerLoad(int tag, double xi_, double eta_, double f1_, double f2_, double f3_, int patchTag)
  : ElementalLoad(tag, LOAD_TAG_IGAFollowerLoad, patchTag),
  xi(xi_),
  eta(eta_),
  f1(f1_),
  f2(f2_),
  f3(f3_)
{
    opserr << "IGAFollowerLoad::IGAFollowerLoad" 
        << "xi = " << xi << endln
        << "eta = " << eta << endln
        << "f1 = " << f1 << endln
        << "f2 = " << f2 << endln
        << "f3 = " << f3 << endln;
}

IGAFollowerLoad::IGAFollowerLoad()
  : ElementalLoad(LOAD_TAG_IGAFollowerLoad),
    xi(0),
  eta(0),
  f1(0),
  f2(0),
  f3(0)
{
}

IGAFollowerLoad::~IGAFollowerLoad()
{
}

const Vector &
IGAFollowerLoad::getData(int &type, double loadFactor)
{
  	type = LOAD_TAG_IGAFollowerLoad;
    data(0) = xi;
    data(1) = eta;
    data(2) = f1;
    data(3) = f2;
    data(4) = f3;

	return data;
}

int 
IGAFollowerLoad::sendSelf(int commitTag, Channel &theChannel)
{
	int dbTag = this->getDbTag();

    static Vector vData(7);
    vData(0) = xi;
    vData(1) = eta;
    vData(2) = f1;
    vData(3) = f2;
    vData(4) = f3;
    vData(5) = eleTag;
    vData(6) = this->getTag();

    int result = theChannel.sendVector(dbTag, commitTag, vData);
    if (result < 0) {
        opserr << "IGAFollowerLoad::sendSelf - failed to send data\n";
        return result;
    }

    return 0;
}

int 
IGAFollowerLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
	int dbTag = this->getDbTag();

    static Vector vData(7);

    int result = theChannel.recvVector(dbTag, commitTag, vData);
    if (result < 0) {
        opserr << "IGAFollowerLoad::recvSelf - failed to recv data\n";
        return result;
    }

    xi = vData(0);
    eta = vData(1);
    f1 = vData(2);
    f2 = vData(3);
    f3 = vData(4);
    eleTag = (int)vData(5);

    this->setTag(vData(6));

    return 0;
}

void 
IGAFollowerLoad::Print(OPS_Stream &s, int flag)
{
	s << "IGAFollowerLoad...";
	s << "  element acted on: " << eleTag << endln;;
	s << "  (xi, eta, f1, f2, f3) = (" 
        << xi << ", " 
        << eta << ", " 
        << f1 << ", " 
        << f2 << ", " 
        << f3 << ") " 
        << endln;
}
