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


//Modified by Jian Zhang, [University of Edinburgh]
//Modified by Panagiotis Kotsovinos, [University of Edinburgh]
//Modified by Liming Jiang [http://openseesforfire.github.io]

// Description: This file contains the class implementation for Beam2dThermalAction.
// Beam2dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Beam2dThermalAction.h>
#include <Vector.h>
#include <Element.h>
Vector Beam2dThermalAction::data(18);

Beam2dThermalAction::Beam2dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2,
					 double t3, double locY3, double t4, double locY4,
					 double t5, double locY5, double t6, double locY6,
					 double t7, double locY7, double t8, double locY8,
					 double t9, double locY9, 
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag), 
   ThermalActionType(LOAD_TAG_Beam2dThermalAction),theSeries(0)
{
  Temp[0]=t1; Temp[1] = t2; Temp[2] = t3; Temp[3] = t4; Temp[4] = t5;
  Temp[5]=t6; Temp[6] = t7; Temp[7] = t8; Temp[8] = t9; 
  Loc[0]=locY1; Loc[1] = locY2; Loc[2] = locY3; Loc[3] = locY4; Loc[4] = locY5;
  Loc[5]=locY6; Loc[6] = locY7; Loc[7] = locY8; Loc[8] = locY9;

  Factors.Zero();
  indicator=1; //without path timeseries defined;
}


Beam2dThermalAction::Beam2dThermalAction(int tag,
					 double locY1, double locY2,
					 TimeSeries* theSeries,int theElementTag
					 )
:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag),theSeries(theSeries),
ThermalActionType(LOAD_TAG_Beam2dThermalAction)
{
  Loc[0]=locY1;
  Loc[8]=locY2;
  
  for(int i= 1; i<8;i++){
		Loc[i]=Loc[0]-i*(Loc[0]-Loc[8])/8;
	}


  for(int i=0 ;i<9;i++) {
		Temp[i]=0;
		TempApp[i]=0;
	}
   Factors.Zero();
   indicator=2 ;// Independent timeseries were created;

}

Beam2dThermalAction::Beam2dThermalAction(int tag,
					 const Vector& locs,
					 TimeSeries* theSeries,int theElementTag
					 )
:ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag),theSeries(theSeries),
ThermalActionType(LOAD_TAG_Beam2dThermalAction)
{
	if (locs.Size()!=9){
		opserr<<" WARNING::Beam2DThermalAction constructor failed to get 9 loc values"<<endln;

	}
  
  for(int i= 0; i<9;i++){
		Loc[i]=locs(i);
	}


  for(int i=0 ;i<9;i++) {
		Temp[i]=0;
		TempApp[i]=0;
	}
   Factors.Zero();
   indicator=2 ;// Independent timeseries were created;

}

Beam2dThermalAction::Beam2dThermalAction(int tag,  
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam2dThermalAction, theElementTag),ThermalActionType(LOAD_TAG_NodalThermalAction),theSeries(0)  
{
	 for(int i=0 ;i<9;i++) {
		Temp[i]=0;
		TempApp[i]=0;
		Loc[i]=0;
	}
	Factors.Zero();
        indicator=3 ;// USing Nodal Thermal Action;
}

Beam2dThermalAction::~Beam2dThermalAction()
{
  indicator=0;
    //if(theSeries!=0)
	  //delete theSeries;
  
  theSeries=0;
}

const Vector &
Beam2dThermalAction::getData(int &type, double loadFactor)
{
  type = ThermalActionType;

  for(int i=0; i<9;i++) {
		data(2*i) = TempApp[i];
		data(2*i+1)= Loc[i];
	}
  Factors.Zero();
	return data;
}

void
Beam2dThermalAction::applyLoad(const Vector &factors)
{
	for(int i=0; i<9 ;i++) {
		TempApp[i]= Temp[i]*factors(i);
	}
	if (theElement != 0)
    theElement->addLoad(this, factors(0));
}

void
Beam2dThermalAction::applyLoad(double loadfactor)
{
	// first determine the load factor
	if (indicator==2) {
		for(int i=0;i<9;i++) {
		  Factors=((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		  //PathTimeSeriesThermal returns absolute temperature;
		  TempApp[i]=Factors(i);
		}
	}else if (indicator ==1){
		for(int i=0;i<9;i++) {
		  TempApp[i]=Temp[i]*loadfactor;
		}
	}

	if (theElement != 0)
    theElement->addLoad(this, loadfactor);
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
  s << "Beam2dThermalAction - reference load : " << Temp[0] <<" change  temp of bot\n";
  s <<  Temp[8] << " change  temp at top\n";
  s << "  element acted on: " << eleTag << endln;
}

