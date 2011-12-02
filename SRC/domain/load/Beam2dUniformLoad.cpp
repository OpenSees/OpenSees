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
                                                                        
// $Revision: 1.5 $
// $Date: 2003-02-14 23:00:57 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dUniformLoad.cpp,v $
                                                                        

// Written: fmk 
//
// Purpose: This file contains the class implementation of Beam2dUniformLoad.

#include <Beam2dUniformLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

Vector Beam2dUniformLoad::data(2);

Beam2dUniformLoad::Beam2dUniformLoad(int tag, double wt, double wa,
				     const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dUniformLoad, theElementTags),
   wTrans(wt), wAxial(wa)
{

}

Beam2dUniformLoad::Beam2dUniformLoad()
  :ElementalLoad(LOAD_TAG_Beam2dUniformLoad), wTrans(0.0), wAxial(0.0)
{

}

Beam2dUniformLoad::~Beam2dUniformLoad()
{

}

const Vector &
Beam2dUniformLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dUniformLoad;
  data(0) = wTrans;
  data(1) = wAxial;
  return data;
}


int 
Beam2dUniformLoad::sendSelf(int commitTag, Channel &theChannel)
{
  int dbTag = this->getDbTag();
  const ID &theElements = this->getElementTags();

  static Vector vectData(3);
  vectData(0) = wTrans;
  vectData(1) = wAxial;
  vectData(2) = theElements.Size();

  int result = theChannel.sendVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  result = theChannel.sendID(dbTag, commitTag, theElements);
  if (result < 0) {
    opserr << "Beam2dUniformLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

int 
Beam2dUniformLoad::recvSelf(int commitTag, Channel &theChannel,  FEM_ObjectBroker &theBroker)
{
  int dbTag = this->getDbTag();

  static Vector vectData(3);

  int result = theChannel.recvVector(dbTag, commitTag, vectData);
  if (result < 0) {
    opserr << "Beam2dUniformLoad::sendSelf - failed to send data\n";
    return result;
  }

  wTrans = vectData(0);;
  wAxial = vectData(1);;
  int numEle = vectData(2);

  if (theElementTags == 0 || theElementTags->Size() != numEle) {
    if (theElementTags != 0)
      delete theElementTags;
    theElementTags = new ID(numEle);
    if (theElementTags == 0) {
      opserr << "Beam2dUniformLoad::sendSelf - failed to create an ID\n";
      return -3;
    }
  }

  result = theChannel.recvID(dbTag, commitTag, *theElementTags);
  if (result < 0) {
    opserr << "Beam2dUniformLoad::sendSelf - failed to send element tags\n";
    return result;
  }
  
  return 0;
}

void 
Beam2dUniformLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dUniformLoad - Reference load" << endln;
  s << "  Transverse: " << wTrans << endln;
  s << "  Axial:      " << wAxial << endln;
  s << "  Elements acted on: " << this->getElementTags();
}
