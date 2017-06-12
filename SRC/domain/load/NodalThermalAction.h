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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/NodalThermalAction.h,v $


// Description: This file contains the class definition for NodalThermalAction.
// NodalThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


//Created by Liming Jiang

#ifndef NodalThermalAction_h
#define NodalThermalAction_h

class TimeSeries;
class Vector;
#include <NodalLoad.h>
#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>
#include <Vector.h>

class NodalThermalAction : public NodalLoad
{
  public:
  // Constructors based on 9, 5 or 2 temperature points 
  // t-temperature; locY-coordinate through the depth of section
  NodalThermalAction(int tag,int theNodeTag,
		      const Vector& locy, TimeSeries* theSeries, Vector* crds=0
		      );
  NodalThermalAction(int tag,int theNodeTag,
		      double locY1, double locY2, double locZ1, double locZ2,
		      TimeSeries* theSeries, Vector* crds=0
		      );
  
  NodalThermalAction(int tag, int theNodeTag,
		      double t1, double locY1, double t2, double locY2, Vector* crds =0
		       );
  
  NodalThermalAction();    
  
  ~NodalThermalAction();
      
  const Vector &getData(int& type);
  const Vector &getCrds(void);
  virtual void applyLoad(const Vector &loadFactors); 
  virtual void applyLoad(double loadFactor); 

  int getThermalActionType(void);
  //int sendSelf(int commitTag, Channel &theChannel);  
  //int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double Temp[15]; //Temperature
  double TempApp[15]; //Temperature applied
  double Loc[10]; // Location through the depth of section
  Vector data; // data for temperature and locations

  //--Adding a factor vector for FireLoadPattern [-BEGIN-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
  //static Vector factors;
  int indicator; //indicator if fireloadpattern was called
  int ThermalActionType; //1:Default, for 2D Beam or Shell elements, 2: for 3D beam element
  Vector Factors;
  Vector Crds;
  TimeSeries* theSeries;
 
  //--Adding a factor vector for FireLoadPattern [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
 };


#endif

