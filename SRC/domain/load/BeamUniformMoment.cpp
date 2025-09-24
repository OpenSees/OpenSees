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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/BeamUniformMoment.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of BeamUniformMoment.

#include <BeamUniformMoment.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector BeamUniformMoment::data(3);

BeamUniformMoment::BeamUniformMoment(int tag, double x, double y, double z,
				     int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_BeamUniformMoment, theElementTag),
   mx(x), my(y), mz(z)
{

}

BeamUniformMoment::BeamUniformMoment()
  :ElementalLoad(LOAD_TAG_BeamUniformMoment),
   mx(0.0), my(0.0), mz(0.0)
{

}

BeamUniformMoment::~BeamUniformMoment()
{

}

const Vector &
BeamUniformMoment::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_BeamUniformMoment;
  data(0) = mx;
  data(1) = my;
  data(2) = mz;
  return data;
}


int 
BeamUniformMoment::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);
  vectData(4) = this->getTag();
  vectData(0) = mx;
  vectData(1) = my;
  vectData(2) = mz;
  vectData(3) = eleTag;

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "BeamUniformMoment::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
BeamUniformMoment::recvSelf(int commitTag, Channel &theChannel,
			    FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "BeamUniformMoment::recvSelf - failed to recv data\n";
    return result;
  }

  mx = vectData(0);;
  my = vectData(1);;
  mz = vectData(2);;
  eleTag = (int)vectData(3);
  this->setTag(vectData(4));

  return 0;
}

void 
BeamUniformMoment::Print(OPS_Stream &s, int flag)
{
  s << "BeamUniformMoment - Reference load: " << this->getTag() << endln;
  s << "  Torque (x): " << mx << endln;
  s << "  About y: " << my << endln;
  s << "  About z: " << mz << endln;
  s << "  Element: " << eleTag << endln;
}
