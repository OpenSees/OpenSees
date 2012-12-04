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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam3dPointLoad.

#include <Beam3dPointLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dPointLoad::data(4);

Beam3dPointLoad::Beam3dPointLoad(int tag, double py, double pz, double dist,
				 int theElementTag, double px)
  :ElementalLoad(tag, LOAD_TAG_Beam3dPointLoad, theElementTag),
   Py(py), Pz(pz), Px(px), x(dist)
{

}

Beam3dPointLoad::Beam3dPointLoad()
  :ElementalLoad(LOAD_TAG_Beam3dPointLoad),
   Py(0.0), Pz(0.0), Px(0.0), x(0.0)
{

}

Beam3dPointLoad::~Beam3dPointLoad()
{

}

const Vector &
Beam3dPointLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dPointLoad;
  data(0) = Py;
  data(1) = Pz;
  data(2) = Px;
  data(3) = x;
  return data;
}

int 
Beam3dPointLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(6);
  vectData(0) = Px;
  vectData(1) = Py;
  vectData(2) = Pz;
  vectData(3) = x;  
  vectData(4) = eleTag;
  vectData(5) = this->getTag();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPointLoad::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
Beam3dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(6);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPointLoad::recvSelf - failed to recv data\n";
    return result;

  }
  this->setTag(vectData(5));
  Px = vectData(0);;
  Py = vectData(1);;
  Pz = vectData(2);;
  x  = vectData(3);  
  eleTag = (int)vectData(4);

  return 0;
}

void 
Beam3dPointLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dPointLoad - Reference load" << endln;
  s << "  Transverse (y): " << Py << endln;
  s << "  Transverse (z): " << Pz << endln;
  s << "  Axial (x):      " << Px << endln;
  s << "  Element: " << eleTag << endln;;
}
