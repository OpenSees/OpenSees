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

Vector Beam2dPartialUniformLoad::data(6);

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad(int tag, double wt, double wa,
						   int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPartialUniformLoad, theElementTag),
   wTrans_a(wt), wTrans_b(wt), wAxial_a(wa), wAxial_b(wa), aOverL(0.0), bOverL(1.0), parameterID(0)
{

}

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad(int tag, double wt, double wa,
						   double aL, double bL, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPartialUniformLoad, theElementTag),
   wTrans_a(wt), wTrans_b(wt), wAxial_a(wa), wAxial_b(wa), aOverL(aL), bOverL(bL), parameterID(0)
{

}

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad(int tag, double wta, double wtb, double waa, double wab,
						   double aL, double bL, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPartialUniformLoad, theElementTag),
   wTrans_a(wta), wTrans_b(wtb), wAxial_a(waa), wAxial_b(wab), aOverL(aL), bOverL(bL), parameterID(0)
{

}

Beam2dPartialUniformLoad::Beam2dPartialUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam2dPartialUniformLoad),
   wTrans_a(0.0), wTrans_b(0.0), wAxial_a(0.0), wAxial_b(0.0), aOverL(0.0), bOverL(0.0), parameterID(0)
{

}

Beam2dPartialUniformLoad::~Beam2dPartialUniformLoad()
{

}

const Vector &
Beam2dPartialUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dPartialUniformLoad;
  data(0) = wTrans_a;
  data(1) = wTrans_b;
  data(2) = wAxial_a;
  data(3) = wAxial_b;
  data(4) = aOverL;
  data(5) = bOverL;
  return data;
}


int 
Beam2dPartialUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(8);
  vectData(0) = wTrans_a;
  vectData(1) = wTrans_b;  
  vectData(2) = wAxial_a;
  vectData(3) = wAxial_b;  
  vectData(4) = eleTag;
  vectData(5) = this->getTag();
  vectData(6) = aOverL;
  vectData(7) = bOverL;

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

  static Vector vectData(8);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPartialUniformLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag(vectData(5));
  wTrans_a = vectData(0);
  wTrans_b = vectData(1);
  wAxial_a = vectData(2);
  wAxial_b = vectData(3);
  eleTag = vectData(4);
  aOverL = vectData(6);
  bOverL = vectData(7);

  return 0;
}

void 
Beam2dPartialUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dPartialUniformLoad - tag " << this->getTag() << endln;
  s << "  Transverse: " << wTrans_a << ' ' << wTrans_b << endln;
  s << "  Axial:      " << wAxial_a << ' ' << wAxial_b << endln;
  s << "  Region:     " << aOverL << " to " << bOverL << endln;
  s << "  Element acted on: " << eleTag << endln;
}

int
Beam2dPartialUniformLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  if (strcmp(argv[0],"wTrans") == 0 || strcmp(argv[0],"wy") == 0) {
    param.setValue(wTrans_a);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"wTransA") == 0 || strcmp(argv[0],"wya") == 0) {
    param.setValue(wTrans_a);
    return param.addObject(5, this);
  }
  if (strcmp(argv[0],"wTransB") == 0 || strcmp(argv[0],"wyb") == 0) {
    param.setValue(wTrans_b);
    return param.addObject(6, this);
  }    
  if (strcmp(argv[0],"wAxial") == 0 || strcmp(argv[0],"wx") == 0) {
    param.setValue(wAxial_a);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"wAxialA") == 0 || strcmp(argv[0],"wxa") == 0) {
    param.setValue(wAxial_a);
    return param.addObject(7, this);
  }
  if (strcmp(argv[0],"wAxialB") == 0 || strcmp(argv[0],"wxb") == 0) {
    param.setValue(wAxial_b);
    return param.addObject(8, this);
  }    
  if (strcmp(argv[0],"aOverL") == 0 || strcmp(argv[0],"a") == 0) {
    param.setValue(aOverL);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"bOverL") == 0 || strcmp(argv[0],"b") == 0) {
    param.setValue(bOverL);
    return param.addObject(4, this);
  }

  return -1;
}

int
Beam2dPartialUniformLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    wTrans_a = info.theDouble;
    wTrans_b = info.theDouble;    
    return 0;
  case 2:
    wAxial_a = info.theDouble;
    wAxial_b = info.theDouble;    
    return 0;
  case 3:
    aOverL = info.theDouble;
    return 0;
  case 4:
    bOverL = info.theDouble;
    return 0;
  case 5:
    wTrans_a = info.theDouble;
    return 0;
  case 6:
    wTrans_b = info.theDouble;
    return 0;
  case 7:
    wAxial_a = info.theDouble;
    return 0;
  case 8:
    wAxial_b = info.theDouble;
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
    data(1) = 1.0;    
    break;
  case 2:
    data(2) = 1.0;
    data(3) = 1.0;    
    break;
  case 3:
    data(4) = 1.0;
    break;
  case 4:
    data(5) = 1.0;
    break;
  case 5:
    data(0) = 1.0;
    break;
  case 6:
    data(1) = 1.0;
    break;
  case 7:
    data(2) = 1.0;
    break;
  case 8:
    data(3) = 1.0;
    break;
  default:
    break;
  }

  return data;
}
