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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/NodalThermalAction.cpp,v $


//Written by Liming Jiang [http://openseesforfire.github.io]


// Description: This file contains the class implementation for NodalThermalAction.
// NodalThermalAction is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <NodalThermalAction.h>
#include <Vector.h>


NodalThermalAction::NodalThermalAction(int tag, int theNodeTag,
					 double t1, double locY1, double t2, double locY2, Vector* crds)
  :NodalLoad(tag, theNodeTag,LOAD_TAG_NodalThermalAction),data(18),Crds(0),
  ThermalActionType(1),theSeries(0)
{
	Temp[0] = t1;
	Temp[8] = t2;
	Loc[0]=locY1;
	Loc[8]=locY2;
	for(int i= 1; i<8;i++){
		Temp[i]=Temp[0]-i*(Temp[0]-Temp[8])/8;
		Loc[i]=Loc[0]-i*(Loc[0]-Loc[8])/8;
	}
	Factors.Zero();
	

	if(crds!=0)
		Crds = (*crds);

 
}

NodalThermalAction::NodalThermalAction(int tag, int theNodeTag,
					 const Vector& locy,
					 TimeSeries* theSeries,Vector* crds 
					 )
  :NodalLoad(tag, theNodeTag,  LOAD_TAG_NodalThermalAction),theSeries(theSeries),Crds(0),
  data(18),ThermalActionType(1)  
{
	for(int i=0 ;i<15;i++) {
		Temp[i]=1;
		TempApp[i]=0;
	}
	for(int i =0;i<9;i++)
		Loc[i]=locy(i);

	Factors.Zero();
	
	if(crds!=0)
		Crds = (*crds);
}

NodalThermalAction::NodalThermalAction(int tag, int theNodeTag,
					 double locY1, double locY2,double locZ1,double locZ2,
					 TimeSeries* theSeries, Vector* crds
					 )
  :NodalLoad(tag, theNodeTag,  LOAD_TAG_NodalThermalAction),theSeries(theSeries),Crds(0),
  data(25),ThermalActionType(2)  
{
	Loc[0]=locY1;Loc[4]=locY2; Loc[5]=locZ1 ;Loc[9] = locZ2;
	for(int i= 1; i<4;i++){
		Loc[i]=Loc[0]+i*(Loc[4]-Loc[0])/4;    //locs through the depth
        Loc[5+i]=Loc[5]+i*(Loc[9]-Loc[5])/4;  //locs through the width
	}

    for(int i=0 ;i<15;i++) {
		Temp[i]=1;   //Here the original temp is set as 1, which will be factorized by the value obtained from 
		TempApp[i]=0;
	}
	Factors.Zero();

	if(crds!=0)
		Crds = (*crds);

	
}

NodalThermalAction::~NodalThermalAction()
{
  indicator=0;
  if(theSeries!=0)
	  delete theSeries;
  
  theSeries=0;
}

const Vector &
NodalThermalAction::getData(int &type)
{
	type=LOAD_TAG_NodalThermalAction;
	if(ThermalActionType==1){
	 for(int i=0; i<9;i++) {
		data(2*i) = TempApp[i];
		data(2*i+1)= Loc[i];
		
	 }
	}
	else if(ThermalActionType==2){
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
	for(int i=0; i<5;i++) {
          data(2*i) = TempApp[i];            //5 temps through y
          data(2*i+1)= Loc[i];               //5 locs through y
          data(3*i+10) = TempApp[i+5];       //5 temps through Z in bottom flange
          data(3*i+11)= TempApp[i+10];       //5 temps through Z in top flange
          data(3*i+12)= Loc[i+5];            //5 locs through Z
		}
	}
	else {
		opserr<<"NodalThermalAction::getData, ThermalActionType tag "
			  <<ThermalActionType<<"is invalid"<<endln;
	}
	Factors.Zero();
	return data;
	//data to store Temperature and location of 9 datapoints
}


void 
NodalThermalAction::applyLoad(const Vector &factors) 
{
	opserr<<"NodalThermalAction::applyLoad(Vector& factors) should not be called)"<<endln;
}

void 
NodalThermalAction::applyLoad(double time)
{
	// first determine the load factor
	if (theSeries!= 0) {
		Factors=((PathTimeSeriesThermal*)theSeries)->getFactors(time);
		//opserr<<"factors"<<Factors<<endln;
		for(int i=0;i<15;i++) {
		  TempApp[i]=Factors(i);
		  if((ThermalActionType==1)&&(i==8))
			break;
		  else if((ThermalActionType==2)&&(i==14))
			break;
		}
	}else{
		for(int i=0;i<15;i++) {
		  TempApp[i]=Temp[i]*time;
		  if((ThermalActionType==1)&&(i==8))
			break;
		  else if((ThermalActionType==2)&&(i==14))
			break;
		}
	}
}

int 
NodalThermalAction::getThermalActionType()
{
  return ThermalActionType;
}

const Vector&
NodalThermalAction::getCrds()
{
  return Crds;
}
/*
int 
NodalThermalAction::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
NodalThermalAction::recvSelf(int commitTag, Channel &theChannel,  
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}
*/

// do it later
void 
NodalThermalAction::Print(OPS_Stream &s, int flag)
{
  s << "NodalThermalAction: " << this->getNodeTag() << endln;
}

