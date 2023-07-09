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
                                                                        
// $Revision: 1.8 $
// $Date: 2008-05-22 22:41:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam2dPointLoad.

#include <Beam2dPointLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector Beam2dPointLoad::data(3);

Beam2dPointLoad::Beam2dPointLoad(int tag, double Pt, double dist,
				 int theElementTag, double Pa)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPointLoad, theElementTag),
   Ptrans(Pt), Paxial(Pa), x(dist), parameterID(0)
{

}

Beam2dPointLoad::Beam2dPointLoad()
  :ElementalLoad(LOAD_TAG_Beam2dPointLoad),
   Ptrans(0.0), Paxial(0.0), x(0.0), parameterID(0)
{

}

Beam2dPointLoad::~Beam2dPointLoad()
{

}

const Vector &
Beam2dPointLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dPointLoad;
  data(0) = Ptrans;
  data(1) = Paxial;
  data(2) = x;
  return data;
}

int 
Beam2dPointLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);
  vectData(0) = Ptrans;
  vectData(1) = Paxial;
  vectData(2) = x;  
  vectData(3) = eleTag;
  vectData(4) = this->getTag();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPointLoad::sendSelf - failed to send data\n";
    return result;
  }
  
  return 0;
}

int 
Beam2dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPointLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag((int)vectData(4));
  Ptrans = vectData(0);
  Paxial = vectData(1);
  x      = vectData(2);  
  eleTag = (int)vectData(3);
  
  return 0;
}

void 
Beam2dPointLoad::Print(OPS_Stream &s, int flag)
{
    s << "Beam3dPointLoad - Reference load" << endln;
    s << "  Transverse: " << Ptrans << endln;
    s << "  Axial:      " << Paxial << endln;
    s << "  Relative Distance: " << x << endln;
    s << "  Element: " << eleTag << endln;;
}

int
Beam2dPointLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return 0;

  if (strcmp(argv[0],"Ptrans") == 0 || strcmp(argv[0],"P") == 0) {
    param.setValue(Ptrans);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"Paxial") == 0 || strcmp(argv[0],"N") == 0) {
    param.setValue(Paxial);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"x") == 0) {
    param.setValue(x);
    return param.addObject(3, this);
  }

  return 0;
}

int
Beam2dPointLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    Ptrans = info.theDouble;
    return 0;
  case 2:
    Paxial = info.theDouble;
    return 0;
  case 3:
    x = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
Beam2dPointLoad::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
Beam2dPointLoad::getSensitivityData(int gradNumber)
{
  data.Zero();

  switch(parameterID) {
  case 1:
    data(0) = 1.0;
    break;
  case 2:
    data(1) = 1.0;
    break;
  case 3:
    data(2) = 1.0;
    break;
  default:
    break;
  }

  return data;
}
