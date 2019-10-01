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

//Added by Liming Jiang, [University of Edinburgh]
//This class is a wrapper processing the thermal action definition


#ifndef ThermalActionWrapper_h
#define ThermalActionWrapper_h

class TimeSeries;

#include <ElementalLoad.h>
#include <TimeSeries.h>
#include <PathTimeSeriesThermal.h>
#include <NodalThermalAction.h>
#include <Matrix.h>

class NodalThermalAction;
class Matrix;

class ThermalActionWrapper : public ElementalLoad
{
  public:
 
 //To set up different interpolation scheme
  
  ThermalActionWrapper(int tag,int EleTag,
		       NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2
			   
		      );

  ThermalActionWrapper(int tag, int EleTag ,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2,
					  NodalThermalAction* theNodalTA3
					 );
  
  ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2,
					 NodalThermalAction* theNodalTA3 ,  NodalThermalAction* theNodalTA4
					
					 );

   ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2,
					 NodalThermalAction* theNodalTA3 ,  NodalThermalAction* theNodalTA4,
					  NodalThermalAction* theNodalTA5 
					 );
  ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2,
					 NodalThermalAction* theNodalTA3 ,  NodalThermalAction* theNodalTA4,
					  NodalThermalAction* theNodalTA5 , NodalThermalAction* theNodalTA6 
					
					 );
  
  ThermalActionWrapper();    
  
  ~ThermalActionWrapper();
  
  const Vector &getData(int &type, double loadFactor);
   
  const Vector &getIntData(const Vector& theIntCrds);

  int setRatios(const Vector& theRatio);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  Matrix NodalLocs; // Location through the depth of section
  NodalThermalAction** theNodalTA;
  Vector theRatios;   //store the ratios
  int ThermalActionType;
  int NumData;
  int ndm;
  double ConstLoc;
  double Transpoint;
  Vector IntData;
  //--Adding a factor vector for FireLoadPattern [-END-]: by L.J&P.K(university of Edinburgh)-07-MAY-2012-///
 };


#endif

