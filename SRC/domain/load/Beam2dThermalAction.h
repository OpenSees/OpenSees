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

//Modified by Jian Zhang, [University of Edinburgh]
//Modified by Panagiotis Kotsovinos, [University of Edinburgh]
//Modified by Liming Jiang [http://openseesforfire.github.io]

// Description: This file contains the class definition for Beam2dThermalAction.
// Beam2dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#ifndef Beam2dThermalAction_h
#define Beam2dThermalAction_h

class TimeSeries;

#include <ElementalLoad.h>
#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>


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
					 double locY1, double locY2,
					 TimeSeries* theSeries, int theElementTag
					 );
  
  Beam2dThermalAction(int tag, 
					 const Vector& locs,
					 TimeSeries* theSeries, int theElementTag
					 );
  Beam2dThermalAction(int tag, int theElementTag);
  
  Beam2dThermalAction();    
  
  ~Beam2dThermalAction();
  
  const Vector &getData(int &type, double loadFactor);
  virtual void applyLoad(const Vector &loadFactors); 
  virtual void applyLoad(double loadFactor); 
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double Temp[9]; //Initial Temperature 
  double TempApp[9]; // Temperature applied
  double Loc[9]; // Location through the depth of section
  static Vector data; // data for temperature and locations

  int ThermalActionType;


  //--Adding a factor vector for FireLoadPattern [-BEGIN-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
 int indicator; //indicator if fireloadpattern was called
  Vector Factors;
  TimeSeries* theSeries;
  //--Adding a factor vector for FireLoadPattern [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
 };


#endif

