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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam3dPointLoad.

#include <Beam3dPointLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam3dPointLoad::data(4);

Beam3dPointLoad::Beam3dPointLoad(int tag, double py, double pz, double dist,
				 const ID &theElementTags, double px)
  :ElementalLoad(tag, LOAD_TAG_Beam3dPointLoad, theElementTags),
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
  const ID &theElements = this->getElementTags();

  static Vector vectData(5);
  vectData(0) = Px;
  vectData(1) = Py;
  vectData(2) = Pz;
  vectData(3) = x;  
  vectData(4) = theElements.Size();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPointLoad::sendSelf - failed to send data\n";
    return result;
  }

  result = theChannel.sendID(dbTag, commitTag, theElements);
  if (result < 0) {
    opserr << "Beam3dPointLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

int 
Beam3dPointLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(5);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam3dPointLoad::sendSelf - failed to send data\n";
    return result;
  }

  Px = vectData(0);;
  Py = vectData(1);;
  Pz = vectData(2);;
  x  = vectData(3);  
  int numEle = vectData(4);


  if (theElementTags == 0 || theElementTags->Size() != numEle) {
    if (theElementTags != 0)
      delete theElementTags;
    theElementTags = new ID(numEle);
    if (theElementTags == 0) {
      opserr << "Beam3dPointLoad::sendSelf - failed to create an ID\n";
      return -3;
    }
  }

  result = theChannel.recvID(dbTag, commitTag, *theElementTags);
  if (result < 0) {
    opserr << "Beam3dPointLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

void 
Beam3dPointLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam3dPointLoad - Reference load" << endln;
  s << "  Transverse (y): " << Py << endln;
  s << "  Transverse (z): " << Pz << endln;
  s << "  Axial (x):      " << Px << endln;
  s << "  Elements acted on: " << this->getElementTags();
}
