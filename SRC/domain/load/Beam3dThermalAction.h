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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dThermalAction.h,v $


// Created: Jian Jiang UoE, July,2011
// //Modified by Liming Jiang [http://openseesforfire.github.io]


// Description: This file contains the class definition for Beam3dThermalAction.
// Beam3dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#ifndef Beam3dThermalAction_h
#define Beam3dThermalAction_h


#include <ElementalLoad.h>
#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>


class Beam3dThermalAction : public ElementalLoad
{
  public:
	  // GR added
	  Beam3dThermalAction(int tag, double indata[], int theElementTag);

  // Constructors based on 9, 5 or 2 temperature points
  // t-temperature; locY-coordinate through the depth of section
  Beam3dThermalAction(int tag,
                double t1, double locY1, double t2, double locY2,
                double t3, double locY3, double t4, double locY4,
                double t5, double locY5, double t6, double t7, double locZ1,
                double t8, double t9, double locZ2, double t10, double t11, double locZ3,
                double t12, double t13, double locZ4, double t14,double t15, double locZ5,
		        int theElementTag);

 Beam3dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2,
					 double t3, double locY3, double t4, double locY4,
					 double t5, double locY5, double t6, double locY6,
					 double t7, double locY7, double t8, double locY8,
					 double t9, double locY9, 
					 int theElementTag);

  Beam3dThermalAction(int tag,
                double locY1, double locY2, double locZ1, double Z2,
                TimeSeries* theSeries,
		        int theElementTag);
  // for receiving 9 data points temperature definition
 Beam3dThermalAction(int tag,
               const Vector& locs, TimeSeries* theSeries,
		        int theElementTag);

  Beam3dThermalAction(int tag, int theElementTag);

  Beam3dThermalAction();

  ~Beam3dThermalAction();


  const Vector &getData(int &type, double loadFactor);
  virtual void applyLoad(const Vector &loadFactors);
  virtual void applyLoad(double loadFactor);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double Temp[25]; //Initial Temperature for using plain patterns
  double TempApp[25]; // Temperature applied
  double Loc[10]; // 5 Locsthrough the depth of section+ 5 locs through the width
  static Vector data; // data for temperature and locations
  int ThermalActionType;

  //--The BeamThermalAction are modified by Liming and having a new structure for applying the fire action
 int indicator; //indicator if fireloadpattern was called
  Vector Factors;
  TimeSeries* theSeries;
};

#endif

