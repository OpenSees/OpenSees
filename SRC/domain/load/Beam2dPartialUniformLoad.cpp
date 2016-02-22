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
                                                                        

#include <Beam2dPartialUniformLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

Vector Beam2dPartialUniformLoad::data(4);

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad(int tag, double wt, double wa,
						   double aL, double bL, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPartialUniformLoad, theElementTag),
   wTrans(wt), wAxial(wa), aOverL(aL), bOverL(bL), parameterID(0)
{

}

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam2dPartialUniformLoad),
   wTrans(0.0), wAxial(0.0), aOverL(0.0), bOverL(0.0), parameterID(0)
{

}

Beam2dPartialUniformLoad::~Beam2dPartialUniformLoad()
{

}

const Vector &
Beam2dPartialUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dPartialUniformLoad;
  data(0) = wTrans;
  data(1) = wAxial;
  data(2) = aOverL;
  data(3) = bOverL;
  return data;
}


int 
Beam2dPartialUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(6);
  vectData(0) = wTrans;
  vectData(1) = wAxial;
  vectData(2) = eleTag;
  vectData(3) = this->getTag();
  vectData(4) = aOverL;
  vectData(5) = bOverL;

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPartialUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  return 0;
}

int 
Beam2dPartialUniformLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(4);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPartialUniformLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag(vectData(3));
  wTrans = vectData(0);
  wAxial = vectData(1);
  eleTag = vectData(2);
  aOverL = vectData(4);
  bOverL = vectData(5);

  return 0;
}

void 
Beam2dPartialUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dPartialUniformLoad - tag " << this->getTag() << endln;
  s << "  Transverse: " << wTrans << endln;
  s << "  Axial:      " << wAxial << endln;
  s << "  Region:     " << aOverL << " to " << bOverL << endln;
  s << "  Element acted on: " << eleTag << endln;
}

int
Beam2dPartialUniformLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  if (strcmp(argv[0],"wTrans") == 0 || strcmp(argv[0],"wy") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"wAxial") == 0 || strcmp(argv[0],"wx") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"aOverL") == 0 || strcmp(argv[0],"a") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"bOverL") == 0 || strcmp(argv[0],"b") == 0)
    return param.addObject(4, this);

  return -1;
}

int
Beam2dPartialUniformLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    wTrans = info.theDouble;
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
Beam2dPartialUniformLoad::activateParameter(int paramID)
{
  parameterID = paramID;

  return 0;
}

const Vector&
Beam2dPartialUniformLoad::getSensitivityData(int gradNumber)
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
  case 4:
    data(3) = 1.0;
    break;
  default:
    break;
  }

  return data;
}
