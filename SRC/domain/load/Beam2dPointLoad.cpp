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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dPointLoad.cpp,v $
                                                                        
// Written: fmk 

// Purpose: This file contains the class implementation Beam2dPointLoad.

#include <Beam2dPointLoad.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>


Vector Beam2dPointLoad::data(3);

Beam2dPointLoad::Beam2dPointLoad(int tag, double Pt, double dist,
				 const ID &theElementTags, double Pa)
  :ElementalLoad(tag, LOAD_TAG_Beam2dPointLoad, theElementTags),
   Ptrans(Pt), Paxial(Pa), x(dist)
{

}

Beam2dPointLoad::Beam2dPointLoad()
  :ElementalLoad(LOAD_TAG_Beam2dPointLoad),
   Ptrans(0.0), Paxial(0.0), x(0.0)
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

  Ptrans = vectData(0);;
  Paxial = vectData(1);;
  x      = vectData(2);  
  int numEle = vectData(3);


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
