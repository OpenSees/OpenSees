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

Vector Beam3dPartialUniformLoad::data(8);

Beam3dPartialUniformLoad::Beam3dPartialUniformLoad(int tag, double wya, double wza, double waa,
						   double aL, double bL, double wyb, double wzb, double wab, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dPartialUniformLoad, theElementTag),
   wTransya(wya), wTransza(wza), wAxiala(waa), aOverL(aL), bOverL(bL), wTransyb(wyb), wTranszb(wzb), wAxialb(wab), parameterID(0)
{

}

Beam3dPartialUniformLoad::Beam3dPartialUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam3dPartialUniformLoad),
   wTransya(0.0), wTransza(0.0), wAxiala(0.0), aOverL(0.0), bOverL(0.0), wTransyb(0.0), wTranszb(0.0), wAxialb(0.0), parameterID(0)
{

}

Beam3dPartialUniformLoad::~Beam3dPartialUniformLoad()
{

}

const Vector &
Beam3dPartialUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dPartialUniformLoad;
  data(0) = wTransya;
  data(1) = wTransza;
  data(2) = wAxiala;
  data(3) = aOverL;
  data(4) = bOverL;
  data(5) = wTransyb;
  data(6) = wTranszb;
  data(7) = wAxialb;
  return data;
}


int 
Beam3dPartialUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();

  static Vector vectData(10);
  vectData(0) = wTransya;
  vectData(1) = wTransza;  
  vectData(2) = wAxiala;
  vectData(3) = eleTag;
  vectData(4) = this->getTag();
  vectData(5) = aOverL;
  vectData(6) = bOverL;
  vectData(7) = wTransyb;
  vectData(8) = wTranszb;
  vectData(9) = wAxialb;
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

  static Vector vectData(10);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPartialUniformLoad::recvSelf - failed to recv data\n";
    return result;
  }

  this->setTag(vectData(4));
  wTransya = vectData(0);
  wTransza = vectData(1);  
  wAxiala = vectData(2);
  eleTag = vectData(3);
  aOverL = vectData(5);
  bOverL = vectData(6);
  wTransyb = vectData(7);
  wTranszb = vectData(8);
  wAxialb = vectData(9);

  return 0;
}

void 
Beam3dPartialUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dPartialUniformLoad - tag " << this->getTag() << endln;
  s << "  Transverse y: at start:" << wTransya << " , at end:" << wTransyb << endln;
  s << "  Transverse z: at start:" << wTransza << " , at end:" << wTranszb << endln;
  s << "  Axial:        at start:" << wAxiala << " , at end:" << wAxialb << endln;
  s << "  Region:     " << aOverL << " to " << bOverL << endln;
  s << "  Element acted on: " << eleTag << endln;
}

int
Beam3dPartialUniformLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;
  
  if (strcmp(argv[0],"wTransya") == 0 || strcmp(argv[0],"wya") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"wTransza") == 0 || strcmp(argv[0],"wza") == 0)
    return param.addObject(5, this);  

  if (strcmp(argv[0],"wAxiala") == 0 || strcmp(argv[0],"wxa") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"aOverL") == 0 || strcmp(argv[0],"a") == 0)
    return param.addObject(3, this);

  if (strcmp(argv[0],"bOverL") == 0 || strcmp(argv[0],"b") == 0)
    return param.addObject(4, this);

  if (strcmp(argv[0], "wTransyb") == 0 || strcmp(argv[0], "wyb") == 0)
    return param.addObject(6, this);

  if (strcmp(argv[0], "wTranszb") == 0 || strcmp(argv[0], "wzb") == 0)
    return param.addObject(7, this);

  if (strcmp(argv[0], "wAxialb") == 0 || strcmp(argv[0], "wxb") == 0)
    return param.addObject(8, this);

  return -1;
}

int
Beam3dPartialUniformLoad::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    wTransya = info.theDouble;
    return 0;
  case 5:
    wTransza = info.theDouble;
    return 0;    
  case 2:
    wAxiala = info.theDouble;
    return 0;
  case 3:
    aOverL = info.theDouble;
    return 0;
  case 4:
    bOverL = info.theDouble;
    return 0;
  case 6:
    wTransyb = info.theDouble;
    return 0;
  case 7:
    wTranszb = info.theDouble;
    return 0;
  case 8:
    wAxialb = info.theDouble;
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
  case 6:
    data(5) = 1.0;
    break;
  case 7:
    data(6) = 1.0;
    break;
  case 8:
    data(7) = 1.0;
    break;
  default:
    break;
  }
  return data;
}
