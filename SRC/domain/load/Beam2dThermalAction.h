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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dThermalAction.h,v $

//Modified by Jian Zhang, [Univeristy of Edinburgh]
//Modified by Panagiotis Kotsovinos, [Univeristy of Edinburgh]


// Description: This file contains the class definition for Beam2dThermalAction.
// Beam2dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#ifndef Beam2dThermalAction_h
#define Beam2dThermalAction_h


#include <ElementalLoad.h>

class Beam2dThermalAction : public ElementalLoad
{
  public:
  // Constructors based on 9, 5 or 2 temperature points 
  // t-temperature; locY-coordinate through the depth of section
  Beam2dThermalAction(int tag,
		      double t1, double locY1, double t2, double locY2,
		      double t3, double locY3, double t4, double locY4,
		      double t5, double locY5, double t6, double locY6,
		      double t7, double locY7, double t8, double locY8,
		      double t9, double locY9, 
		      int theElementTag);
  
  Beam2dThermalAction(int tag,
		      double t1, double locY1, double t2, double locY2,
		      double t3, double locY3, double t4, double locY4,
		      double t5, double locY5, int theElementTag);
  
  Beam2dThermalAction(int tag,
		      double t1, double locY1, double t2, double locY2,
		      int theElementTag);
  
  Beam2dThermalAction();    
  
  ~Beam2dThermalAction();
  
  const Vector &getData(int &type, double loadFactor);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double T1; //Temperature
  double LocY1; // Location through the depth of section
  double T2;
  double LocY2;
  double T3;
  double LocY3;
  double T4;
  double LocY4;
  double T5;
  double LocY5;
  double T6;
  double LocY6;
  double T7;
  double LocY7;
  double T8;
  double LocY8;
  double T9;
  double LocY9;
  static Vector data; // data for temperature and locations

  //--Adding a factor vector for FireLoadPattern [-BEGIN-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
  static Vector factors;
  int indicator; //indicator if fireloadpattern was called
  double Factor2;
  double Factor3;
  double Factor4;
  double Factor5;
  double Factor6;
  double Factor7;
  double Factor8;
  double Factor9;
  //--Adding a factor vector for FireLoadPattern [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
 };


#endif

