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
                                                                        
// $Revision$
// $Date$
// $Source$
                                                                        

#include <Beam3dPartialUniformLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector Beam3dPartialUniformLoad::data(5);

Beam3dPartialUniformLoad::Beam3dPartialUniformLoad(int tag, double wy, double wz, double wa,
						   double aL, double bL, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dPartialUniformLoad, theElementTag),
   wTransy(wy), wTransz(wz), wAxial(wa), aOverL(aL), bOverL(bL), parameterID(0)
{

}

Beam3dPartialUniformLoad::Beam3dPartialUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam3dPartialUniformLoad),
   wTransy(0.0), wTransz(0.0), wAxial(0.0), aOverL(0.0), bOverL(0.0), parameterID(0)
{

}

Beam3dPartialUniformLoad::~Beam3dPartialUniformLoad()
{

}

const Vector &
Beam3dPartialUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dPartialUniformLoad;
  data(0) = wTransy;
  data(1) = wTransz;
  data(2) = wAxial;
  data(3) = aOverL;
  data(4) = bOverL;
  return data;
}


int 
Beam3dPartialUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(7);
  vectData(0) = wTransy;
  vectData(1) = wTransz;  
  vectData(2) = wAxial;
  vectData(3) = eleTag;
  vectData(4) = this->getTag();
  vectData(5) = aOverL;
  vectData(6) = bOverL;

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPartialUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
Beam3dPartialUniformLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(7);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPartialUniformLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag(vectData(4));
  wTransy = vectData(0);
  wTransz = vectData(1);  
  wAxial = vectData(2);
  eleTag = vectData(3);
  aOverL = vectData(5);
  bOverL = vectData(6);

  return 0;
}

void 
Beam3dPartialUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dPartialUniformLoad - tag " << this->getTag() << endln;
  s << "  Transverse y: " << wTransy << endln;
  s << "  Transverse z: " << wTransz << endln;
  s << "  Axial:      " << wAxial << endln;
  s << "  Region:     " << aOverL << " to " << bOverL << endln;
  s << "  Element acted on: " << eleTag << endln;
}

int
Beam3dPartialUniformLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  if (strcmp(argv[0],"wTransy") == 0 || strcmp(argv[0],"wy") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"wTransz") == 0 || strcmp(argv[0],"wz") == 0)
    return param.addObject(5, this);  

  if (strcmp(argv[0],"wAxial") == 0 || strcmp(argv[0],"wx") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"aOverL") == 0 || strcmp(argv[0],"a") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"bOverL") == 0 || strcmp(argv[0],"b") == 0)
    return param.addObject(4, this);

  return -1;
}

int
Beam3dPartialUniformLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    wTransy = info.theDouble;
    return 0;
  case 5:
    wTransz = info.theDouble;
    return 0;    
  case 2:
    wAxial = info.theDouble;
    return 0;
  case 3:
    aOverL = info.theDouble;
    return 0;
  case 4:
    bOverL = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

int
Beam3dPartialUniformLoad::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
Beam3dPartialUniformLoad::getSensitivityData(int gradNumber)
{
  data.Zero();

  switch(parameterID) {
  case 1:
    data(0) = 1.0;
    break;
  case 5:
    data(1) = 1.0;
    break;    
  case 2:
    data(2) = 1.0;
    break;
  case 3:
    data(3) = 1.0;
    break;
  case 4:
    data(4) = 1.0;
    break;
  default:
    break;
  }

  return data;
}
