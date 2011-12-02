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
                                                                        
// $Revision: 1.3 $
// $Date: 2003-02-14 23:00:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dUniformLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam3dUniformLoad.

#include <Beam3dUniformLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dUniformLoad::data(3);

Beam3dUniformLoad::Beam3dUniformLoad(int tag, double wY, double wZ, double wX,
				     const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam3dUniformLoad, theElementTags),
   wy(wY), wz(wZ), wx(wX)
{

}

Beam3dUniformLoad::Beam3dUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam3dUniformLoad),
   wy(0.0), wz(0.0), wx(0.0)
{

}

Beam3dUniformLoad::~Beam3dUniformLoad()
{

}

const Vector &
Beam3dUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam3dUniformLoad;
  data(0) = wy;
  data(1) = wz;
  data(2) = wx;
  return data;
}


int 
Beam3dUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  const ID &theElements = this->getElementTags();

  static Vector vectData(4);
  vectData(0) = wx;
  vectData(1) = wy;
  vectData(2) = wz;
  vectData(3) = theElements.Size();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  result = theChannel.sendID(dbTag, commitTag, theElements);
  if (result < 0) {
    opserr << "Beam3dUniformLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

int 
Beam3dUniformLoad::recvSelf(int commitTag, Channel &theChannel,
			    FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(4);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  wx = vectData(0);;
  wy = vectData(1);;
  wz = vectData(2);;
  int numEle = vectData(3);


  if (theElementTags == 0 || theElementTags->Size() != numEle) {
    if (theElementTags != 0)
      delete theElementTags;
    theElementTags = new ID(numEle);
    if (theElementTags == 0) {
      opserr << "Beam3dUniformLoad::sendSelf - failed to create an ID\n";
      return -3;
    }
  }

  result = theChannel.recvID(dbTag, commitTag, *theElementTags);
  if (result < 0) {
    opserr << "Beam3dUniformLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

void 
Beam3dUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dUniformLoad - Reference load" << endln;
  s << "  Transverse (y): " << wy << endln;
  s << "  Transverse (z): " << wz << endln;
  s << "  Axial (x):      " << wx << endln;
  s << "  Elements acted on: " << this->getElementTags();
}
