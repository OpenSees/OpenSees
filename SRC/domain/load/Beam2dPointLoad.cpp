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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-09-05 23:08:00 $
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
				 const ID &theElementTags, double Pa)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPointLoad, theElementTags),
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
  const ID &theElements = this->getElementTags();

  static Vector vectData(4);
  vectData(0) = Ptrans;
  vectData(1) = Paxial;
  vectData(2) = x;  
  vectData(3) = theElements.Size();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPointLoad::sendSelf - failed to send data\n";
    return result;
  }

  result = theChannel.sendID(dbTag, commitTag, theElements);
  if (result < 0) {
    opserr << "Beam2dPointLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

int 
Beam2dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(4);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dPointLoad::sendSelf - failed to send data\n";
    return result;
  }

  Ptrans = vectData(0);
  Paxial = vectData(1);
  x      = vectData(2);  
  int numEle = (int)vectData(3);


  if (theElementTags == 0 || theElementTags->Size() != numEle) {
    if (theElementTags != 0)
      delete theElementTags;
    theElementTags = new ID(numEle);
    if (theElementTags == 0) {
      opserr << "Beam2dPointLoad::sendSelf - failed to create an ID\n";
      return -3;
    }
  }

  result = theChannel.recvID(dbTag, commitTag, *theElementTags);
  if (result < 0) {
    opserr << "Beam2dPointLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

void 
Beam2dPointLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dPointLoad - reference load : (" << Ptrans
    << ", " << Paxial << ") acting at : " << x << " relative to length\n";
  s << "  elements acted on: " << this->getElementTags();
}

int
Beam2dPointLoad::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 1)
    return -1;

  if (strcmp(argv[0],"Ptrans") == 0 || strcmp(argv[0],"P") == 0)
    return param.addObject(1, this);

  if (strcmp(argv[0],"Paxial") == 0 || strcmp(argv[0],"N") == 0)
    return param.addObject(2, this);

  if (strcmp(argv[0],"x") == 0)
    return param.addObject(3, this);

  return -1;
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
