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
                                                                                                                                               
// $Revision: 1.1 $
// $Date: 2011-07-18 10:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dThermalAction.cpp,v $


//Modified by Jian Zhang, [Univeristy of Edinburgh]
//Modified by Panagiotis Kotsovinos, [Univeristy of Edinburgh]


// Description: This file contains the class implementation for Beam2dThermalAction.
// Beam2dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Beam2dThermalAction.h>
#include <Vector.h>

Vector Beam2dThermalAction::data(18);

Beam2dThermalAction::Beam2dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2,
					 double t3, double locY3, double t4, double locY4,
					 double t5, double locY5, double t6, double locY6,
					 double t7, double locY7, double t8, double locY8,
					 double t9, double locY9, 
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), 
   T1(t1),LocY1(locY1),T2(t2),LocY2(locY2),T3(t3),LocY3(locY3),T4(t4),LocY4(locY4),
   T5(t5),LocY5(locY5),T6(t6),LocY6(locY6),T7(t7),LocY7(locY7),T8(t8),LocY8(locY8),
   T9(t9),LocY9(locY9)
{
  
}


Beam2dThermalAction::Beam2dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2,
					 double t3, double locY3, double t4, double locY4,
					 double t5, double locY5, int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), 
   T1(t1),LocY1(locY1),T3(t2),LocY3(locY2),T5(t3),LocY5(locY3),T7(t4),LocY7(locY4),
   T9(t5),LocY9(locY5)
{
  T2 = (T1 + T3)/2;
  LocY2 = (LocY1 + LocY3)/2;
  T4 = (T3 + T5)/2;
  LocY4 = (LocY3 + LocY5)/2;
  T6 = (T5 + T7)/2;
  LocY6 = (LocY5 + LocY7)/2;
  T8 = (T7 + T9)/2;
  LocY8 = (LocY7 + LocY9)/2;
}


Beam2dThermalAction::Beam2dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2, 
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), 
   T1(t1),LocY1(locY1),T9(t2),LocY9(locY2)
   
{
  T2 = T1 - 1*(T1 - T9)/8;
  T3 = T1 - 2*(T1 - T9)/8;
  T4 = T1 - 3*(T1 - T9)/8;
  T5 = T1 - 4*(T1 - T9)/8;
  T6 = T1 - 5*(T1 - T9)/8;
  T7 = T1 - 6*(T1 - T9)/8;
  T8 = T1 - 7*(T1 - T9)/8;
  
  LocY2 = LocY1 + 1*(LocY9 - LocY1)/8;
  LocY3 = LocY1 + 2*(LocY9 - LocY1)/8;
  LocY4 = LocY1 + 3*(LocY9 - LocY1)/8;
  LocY5 = LocY1 + 4*(LocY9 - LocY1)/8;
  LocY6 = LocY1 + 5*(LocY9 - LocY1)/8;
  LocY7 = LocY1 + 6*(LocY9 - LocY1)/8;
  LocY8 = LocY1 + 7*(LocY9 - LocY1)/8;
}

Beam2dThermalAction::~Beam2dThermalAction()
{
  indicator=0;
}

const Vector &
Beam2dThermalAction::getData(int &type, double loadFactor)
{
  type = LOAD_TAG_Beam2dThermalAction;
  data(0) = T1;
  data(1) = LocY1;
  data(2) = T2;
  data(3) = LocY2;
  data(4) = T3;
  data(5) = LocY3;
  data(6) = T4;
  data(7) = LocY4;
  data(8) = T5;
  data(9) = LocY5;
  data(10) = T6;
  data(11) = LocY6;
  data(12) = T7;
  data(13) = LocY7;
  data(14) = T8;
  data(15) = LocY8;
  data(16) = T9;
  data(17) = LocY9;

  return data;
}


int 
Beam2dThermalAction::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam2dThermalAction::recvSelf(int commitTag, Channel &theChannel,  
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

// do it later
void 
Beam2dThermalAction::Print(OPS_Stream &s, int flag)
{
  s << "Beam2dThermalAction - reference load : " << T1 <<" change  temp of bot\n";
  s <<  T2 << " change  temp at top\n";
  s << "  element acted on: " << eleTag << endln;
}

