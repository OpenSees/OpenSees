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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam3dThermalAction.cpp,v $


// Written: Jian Jiang, UoE
// //Modified by Liming Jiang [http://openseesforfire.github.io]


// Description: This file contains the class implementation for Beam3dThermalAction.
// Beam3dThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <Beam3dThermalAction.h>
#include <Vector.h>
#include <Element.h>
Vector Beam3dThermalAction::data(35);

// It allows for linear interpolation in a 5x5 grid (5 locs in Z and Y)
Beam3dThermalAction::Beam3dThermalAction(int tag,
	double indata[],
	int theElementTag)
	:ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag),
	ThermalActionType(LOAD_TAG_Beam3dThermalAction), theSeries(0)
{
	for (int i = 0; i < 5; i++) {
		Loc[i] = indata[i]; Loc[i + 5] = indata[i + 5];
	}
	for (int i = 0; i < 25; i++) {
		Temp[i] = indata[i + 10];
	}

	Factors.Zero();
	indicator = 6; // without path timeseries defined;
}

//Basically there are 5 datapoints respectively in the top flange , the web , and the bottom flange . 
// And 5 loc data for defining the zones along y direction, and another 5 for z direction.
Beam3dThermalAction::Beam3dThermalAction(int tag,
                         double t1, double locY1, double t2, double locY2,
                         double t3, double locY3, double t4, double locY4,
                         double t5, double locY5, double t6, double t7, double locZ1,
                         double t8, double t9, double locZ2, double t10, double t11, double locZ3,
                         double t12, double t13, double locZ4, double t14, double t15,double locZ5,
			             int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag),
  ThermalActionType(LOAD_TAG_Beam3dThermalAction), theSeries(0)
{
  Temp[0]=t1; Temp[1] = t2; Temp[2] = t3; Temp[3] = t4; Temp[4] = t5;
  Temp[5]=t6; Temp[6] = t8; Temp[7] = t10; Temp[8] = t12; Temp[9] = t14;
  Temp[10]=t7; Temp[11] = t9; Temp[12] = t11; Temp[13] = t13; Temp[14] = t15;

  Loc[0]=locY1; Loc[1] = locY2; Loc[2] = locY3; Loc[3] = locY4; Loc[4] = locY5;
  Loc[5]=locZ1; Loc[6] = locZ2; Loc[7] = locZ3; Loc[8] = locZ4; Loc[9] = locZ5;

  Factors.Zero();
  indicator=1; //without path timeseries defined;
}

Beam3dThermalAction::Beam3dThermalAction(int tag, 
					 double t1, double locY1, double t2, double locY2,
					 double t3, double locY3, double t4, double locY4,
					 double t5, double locY5, double t6, double locY6,
					 double t7, double locY7, double t8, double locY8,
					 double t9, double locY9, 
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag), 
   ThermalActionType(LOAD_TAG_Beam3dThermalAction), theSeries(0)
{
  Temp[0]=t1; Temp[1] = t2; Temp[2] = t3; Temp[3] = t4; Temp[4] = t5;
  Temp[5]=t6; Temp[6] = t7; Temp[7] = t8; Temp[8] = t9; 
  Loc[0]=locY1; Loc[1] = locY2; Loc[2] = locY3; Loc[3] = locY4; Loc[4] = locY5;
  Loc[5]=locY6; Loc[6] = locY7; Loc[7] = locY8; Loc[8] = locY9;

  Factors.Zero();
  indicator=5; //without path timeseries defined;
}

Beam3dThermalAction::Beam3dThermalAction(int tag,
                         double locY1, double locY2, double locZ1, double locZ2,
                         TimeSeries* theSeries, 
                         int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag),theSeries(theSeries),
  ThermalActionType(LOAD_TAG_Beam3dThermalAction)
{
    Loc[0]=locY1;Loc[4]=locY2; Loc[5]=locZ1 ;Loc[9] = locZ2;

    for(int i= 1; i<4;i++){
		Loc[i]=Loc[0]+i*(Loc[4]-Loc[0])/4;    //locs through the depth
                Loc[5+i]=Loc[5]+i*(Loc[9]-Loc[5])/4;  //locs through the width
	}
	Factors.Zero();
    for(int i=0 ;i<15;i++) {
		Temp[i]=0;   //Here the original temp is set as 1, which will be factorized by the value obtained from 
		TempApp[i]=0;
	}

	indicator=2 ;// Using PathTimeSeriesThermal for elemental thermal action
	
}


//Using 9 data points for defining temperature distribution
Beam3dThermalAction::Beam3dThermalAction(int tag,
                         const Vector& locs,
                         TimeSeries* theSeries, 
                         int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag),theSeries(theSeries),
  ThermalActionType(LOAD_TAG_Beam3dThermalAction)
{
	if (locs.Size()!=9){
		opserr<<" WARNING::Beam3DThermalAction constructor failed to get 9 loc values"<<endln;

	}
   
   for(int i= 0; i<9;i++){
		Loc[i]=locs(i);
	}
	
	Factors.Zero();
    for(int i=0 ;i<15;i++) {
		Temp[i]=0;   //Here the original temp is set as 1, which will be factorized by the value obtained from 
		TempApp[i]=0;
	}

	indicator=4 ;// Using PathTimeSeriesThermal and 9 data points for elemental thermal action
	
}

Beam3dThermalAction::Beam3dThermalAction(int tag,  
					 int theElementTag)
  :ElementalLoad(tag, LOAD_TAG_Beam3dThermalAction, theElementTag),ThermalActionType(LOAD_TAG_NodalThermalAction), theSeries(0)
{
	 Factors.Zero();
	 for(int i=0 ;i<15;i++) {
		Temp[i]=0;
		TempApp[i]=0;
	}
	 indicator=3 ;// USing Nodal Thermal Action;

}




Beam3dThermalAction::~Beam3dThermalAction()
{
  indicator=0;
  if(theSeries!=0)
    theSeries=0;
}

const Vector &
Beam3dThermalAction::getData(int &type, double loadFactor)
{
 type = ThermalActionType;
  //data(0) = T1;
  //data(1) = LocY1;
  //....
  //data(8) = T5;
 // data(9) = LocY5;
  //data(10) = T6;
  //data(11) = T7;
  //data(13) = LocZ1;
  //...
  //data(22) = T14;
  //data(23) = T15;
  //data(24) = LocZ5;
 if (indicator ==4){
	 //using 9 data points
	 data.resize(18);
	  for(int i=0; i<9;i++) {
		data(2*i) = TempApp[i];
		data(2*i+1)= Loc[i];
	}
	  
 }
 else if (indicator == 6) {
	 data.resize(35);
	 for (int i = 0; i < 5; i++) {
		 data(i) = Loc[i];            //5 locs through y
		 data(i + 5) = Loc[i + 5];    //5 locs through z
	 }
	 for (int j = 0; j < 25; j++) {
		 data(j + 10) = TempApp[j];   //25 temps
	 }

 }
 else{
	 data.resize(25);
  for(int i=0; i<5;i++) {
          data(2*i) = TempApp[i];            //5 temps through y
          data(2*i+1)= Loc[i];               //5 locs through y
          data(3*i+10) = TempApp[i+5];       //5 temps through Z in bottom flange
          data(3*i+11)= TempApp[i+10];       //5 temps through Z in top flange
          data(3*i+12)= Loc[i+5];            //5 locs through Z
	}
 
 }

  Factors.Zero();
  return data;
}


void 
Beam3dThermalAction::applyLoad(const Vector &factors)
{
	if (indicator ==4||indicator ==5){
	   for(int i=0; i<9 ;i++) {
		 TempApp[i]= Temp[i]*factors(i);
	   }
	}
	else{
	   for(int i=0; i<25 ;i++) {
		  TempApp[i]= Temp[i]*factors(i);
	   }
	}


	if (theElement != 0)
           theElement->addLoad(this, factors(0));
}

void 
Beam3dThermalAction::applyLoad(double loadfactor)
{
	// type 2& 4 are defined with 9 data points
	//type 1& 5 are defined with the given temperature
	if (indicator==2) {
		//Looking for loadfactors from timeseries;
		Factors=((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		for(int i=0;i<15;i++) {
		   //PathTimeSeriesThermal returns absolute temperature;
		  TempApp[i]=Factors(i);
		}
	}else if(indicator ==1 || indicator == 6) {
		for(int i=0;i<25;i++) {
		  TempApp[i]=Temp[i]*loadfactor;
		}
	}
	else if(indicator ==4) {
		 Factors=((PathTimeSeriesThermal*)theSeries)->getFactors(loadfactor);
		for(int i=0;i<9;i++) {
		  //PathTimeSeriesThermal returns absolute temperature;
		  TempApp[i]=Factors(i);
		}
	}
	else if(indicator ==5) {
		for(int i=0;i<9;i++) {
		  TempApp[i]=Temp[i]*loadfactor;
		}
	}

	if (theElement != 0)
    theElement->addLoad(this, loadfactor);
	//This will be called to send load factor to element
}
int 
Beam3dThermalAction::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
Beam3dThermalAction::recvSelf(int commitTag, Channel &theChannel,  
			 FEM_ObjectBroker &theBroker)
{
  return -1;
}

// do it later
void 
Beam3dThermalAction::Print(OPS_Stream &s, int flag)
{
	if (indicator == 4 || indicator == 5) {
		s << "Beam3dThermalAction - reference load : " << TempApp[0] << " at bot\n";
		s << TempApp[8] << " at top\n";
		s << "  element acted on: " << eleTag << endln;
	}
	else if (indicator == 6) {
		s << "Beam3dThermalAction - reference load : " << TempApp[0] << " at bot\n";
		s << TempApp[24] << " at top\n";
		s << "  element acted on: " << eleTag << endln;
	}
	else {
		s << "Beam3dThermalAction - reference load : " << TempApp[0] << " at bot\n";
		s << TempApp[5] << " at top\n";
		s << "  element acted on: " << eleTag << endln;
	}
  
}

