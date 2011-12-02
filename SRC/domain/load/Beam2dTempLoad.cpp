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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dTempLoad.cpp,v $
                                                                        
// Written: Scott R. Hamilton 15 July 2002

// Purpose: This file contains the class implementation Beam2dTempLoad 
// modeled after Beam2dPointLoad.cpp by fmk.

#include <Beam2dTempLoad.h>
#include <Vector.h>

Vector Beam2dTempLoad::data(4);

Beam2dTempLoad::Beam2dTempLoad(int tag, 
			       double temp1, double temp2, 
			       double temp3, double temp4, 
			       const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dTempLoad, theElementTags), 
  Ttop1(temp1),  Tbot1(temp2), Ttop2(temp3), Tbot2(temp4)
{

}
Beam2dTempLoad::Beam2dTempLoad(int tag, 
			       double temp1, 
			       const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dTempLoad, theElementTags), 
  Ttop1(temp1),  Tbot1(temp1), Ttop2(temp1), Tbot2(temp1)
{

}

Beam2dTempLoad::Beam2dTempLoad(int tag, 
			       double temp1, double temp2, 
			       const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dTempLoad, theElementTags), 
  Ttop1(temp1),  Tbot1(temp2), Ttop2(temp1), Tbot2(temp2)
{

}
Beam2dTempLoad::Beam2dTempLoad(int tag, const ID &theElementTags)
  :ElementalLoad(tag, LOAD_TAG_Beam2dTempLoad, theElementTags), 
  Ttop1(0.0), Tbot1(0.0), Ttop2(0.0), Tbot2(0.0)
{

}
Beam2dTempLoad::Beam2dTempLoad()
  :ElementalLoad(LOAD_TAG_Beam2dTempLoad), 
  Ttop1(0.0), Tbot1(0.0), Ttop2(0.0), Tbot2(0.0)
{

}
Beam2dTempLoad::~Beam2dTempLoad()
{

}

const Vector &
Beam2dTempLoad::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dTempLoad;
  data(0) = Ttop1;
  data(1) = Tbot1;
  data(2) = Ttop2;
  data(3) = Tbot2;
  return data;
}

int 
Beam2dTempLoad::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam2dTempLoad::recvSelf(int commitTag, Channel &theChannel,  
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
Beam2dTempLoad::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dTempLoad - reference load : " << Ttop1 << " change in temp at top of node 1 : " << Tbot1 << " change in temp at bottom of node 1\n";
  s <<  Ttop2 << " change in temp at top of node 2 : " << Tbot2 << " change in temp at bottom of node 2\n";
  s << "  elements acted on: " << this->getElementTags();
}
