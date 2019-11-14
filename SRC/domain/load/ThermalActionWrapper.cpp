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
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/ThermalActionWrapper.cpp,v $


//Added by Liming Jiang, [University of Edinburgh]



// Description: This file contains the class implementation for ThermalActionWrapper.
// ThermalActionWrapper is a thermal field class created to store the temperature
// distribution through the depth of section defined by temperature and location.


#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>
#include <ThermalActionWrapper.h>
#include <Vector.h>
#include <Element.h>
#include <math.h>


//Vector ThermalActionWrapper::IntData;

//two nodes
ThermalActionWrapper::ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2)
  :ElementalLoad(tag, LOAD_TAG_ThermalActionWrapper, EleTag),ThermalActionType(LOAD_TAG_ThermalActionWrapper),
  theRatios(0),NumData(0),IntData(0),ConstLoc(0),Transpoint(0)
{
	 theNodalTA = new NodalThermalAction*[2];
	 theNodalTA[0] = theNodalTA1;
	 theNodalTA[1] = theNodalTA2;

	 ndm = (theNodalTA[0]->getCrds()).Size();
	 
	 NodalLocs.Zero();
	 NodalLocs.resize(2,ndm);

	 for (int i=0; i<2; i++){
		 for (int j=0;j<ndm; j++){
			NodalLocs(i,j)=(theNodalTA[i]->getCrds())(j);
		}
	 }

	 if((theNodalTA[0]->getThermalActionType())!= (theNodalTA[1]->getThermalActionType())){
		 opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is not consistent"<<endln;
	 }
	 else{
		 if(theNodalTA[0]->getThermalActionType()==1){
			NumData=9; //This is for 2D beam
		 }
		 else if( theNodalTA[0]->getThermalActionType()==2  ){
			NumData = 15; //This is for 3D I-section Beam
		 }
		 else {
            opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is unable to be identified"<<endln;
		 }
	  }

	
	 
	 
}

//three nodes

ThermalActionWrapper::ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2,  
					  NodalThermalAction* theNodalTA3)
  :ElementalLoad(tag, LOAD_TAG_ThermalActionWrapper, EleTag),ThermalActionType(LOAD_TAG_ThermalActionWrapper) ,
  theRatios(0),NumData(0),IntData(0),ConstLoc(0),Transpoint(0)
{
	 theNodalTA = new NodalThermalAction*[3];
	 theNodalTA[0] = theNodalTA1;
	 theNodalTA[1] = theNodalTA2;
	 theNodalTA[2] = theNodalTA3;

	ndm = (theNodalTA[0]->getCrds()).Size();
	 
	 NodalLocs.Zero();
	 NodalLocs.resize(3,ndm);

	 for (int i=0; i<3; i++){
		 for (int j=0;j<ndm; j++){
			NodalLocs(i,j)=(theNodalTA[i]->getCrds())(j);
		}
	 }

	 if((theNodalTA[0]->getThermalActionType())!= (theNodalTA[2]->getThermalActionType())){
		 opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is not consistent"<<endln;
	 }
	 else{
		 if(theNodalTA[0]->getThermalActionType()==1){
			NumData=9; //This is for 2D beam
		 }
		 else if( theNodalTA[0]->getThermalActionType()==2  ){
			NumData = 15; //This is for 3D I-section Beam
		 }
		 else {
            opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is unable to be identified"<<endln;
		 }
	  }
	
	 
}

//four nodes

ThermalActionWrapper::ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2, 
					  NodalThermalAction* theNodalTA3 , NodalThermalAction* theNodalTA4 )
  :ElementalLoad(tag, LOAD_TAG_ThermalActionWrapper, EleTag),ThermalActionType(LOAD_TAG_ThermalActionWrapper), 
  theRatios(0),NumData(0),ConstLoc(0),Transpoint(0)
{
	 theNodalTA = new NodalThermalAction*[5];
	 theNodalTA[0] = theNodalTA1;
	 theNodalTA[1] = theNodalTA2;
	 theNodalTA[2] = theNodalTA3;
	 theNodalTA[3] = theNodalTA4;

	 ndm = (theNodalTA[0]->getCrds()).Size();
	 
	 NodalLocs.Zero();
	 NodalLocs.resize(4,ndm);

	 for (int i=0; i<4; i++){
		 for (int j=0;j<ndm; j++){
			NodalLocs(i,j)=(theNodalTA[i]->getCrds())(j);
		}
	 }

	 if((theNodalTA[0]->getThermalActionType())!= (theNodalTA[3]->getThermalActionType())){
		 opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is not consistent"<<endln;
	 }
	 else{
		 if(theNodalTA[0]->getThermalActionType()==1){
			NumData=9; //This is for 2D beam
		 }
		 else if( theNodalTA[0]->getThermalActionType()==2  ){
			NumData = 15; //This is for 3D I-section Beam
		 }
		 else {
            opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is unable to be identified"<<endln;
		 }
	  }
		 
}

//five nodes

ThermalActionWrapper::ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2, 
					  NodalThermalAction* theNodalTA3 , NodalThermalAction* theNodalTA4 ,
					  NodalThermalAction* theNodalTA5)
  :ElementalLoad(tag, LOAD_TAG_ThermalActionWrapper, EleTag),ThermalActionType(LOAD_TAG_ThermalActionWrapper), 
  theRatios(0),NumData(0),ConstLoc(0),Transpoint(0)
{
	 theNodalTA = new NodalThermalAction*[5];
	 theNodalTA[0] = theNodalTA1;
	 theNodalTA[1] = theNodalTA2;
	 theNodalTA[2] = theNodalTA3;
	 theNodalTA[3] = theNodalTA4;
	 theNodalTA[4] = theNodalTA5;

	 ndm = (theNodalTA[0]->getCrds()).Size();
	 
	 NodalLocs.Zero();
	 NodalLocs.resize(5,ndm);

	 for (int i=0; i<5; i++){
		 for (int j=0;j<ndm; j++){
			NodalLocs(i,j)=(theNodalTA[i]->getCrds())(j);
		}
	 }

	 if((theNodalTA[0]->getThermalActionType())!= (theNodalTA[4]->getThermalActionType())){
		 opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is not consistent"<<endln;
	 }
	 else{
		 if(theNodalTA[0]->getThermalActionType()==1){
			NumData=9; //This is for 2D beam
		 }
		 else if( theNodalTA[0]->getThermalActionType()==2  ){
			NumData = 15; //This is for 3D I-section Beam
		 }
		 else {
            opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is unable to be identified"<<endln;
		 }
	  }
		 
}

//six nodes

ThermalActionWrapper::ThermalActionWrapper(int tag, int EleTag,
					  NodalThermalAction* theNodalTA1 ,  NodalThermalAction* theNodalTA2, 
					  NodalThermalAction* theNodalTA3 , NodalThermalAction* theNodalTA4 ,
					  NodalThermalAction* theNodalTA5, NodalThermalAction* theNodalTA6)
  :ElementalLoad(tag, LOAD_TAG_ThermalActionWrapper, EleTag),ThermalActionType(LOAD_TAG_ThermalActionWrapper), 
  theRatios(0),NumData(0),ConstLoc(0),Transpoint(0)
{
	 theNodalTA = new NodalThermalAction*[6];
	 theNodalTA[0] = theNodalTA1;
	 theNodalTA[1] = theNodalTA2;
	 theNodalTA[2] = theNodalTA3;
	 theNodalTA[3] = theNodalTA4;
	 theNodalTA[4] = theNodalTA5;
	 theNodalTA[5] = theNodalTA6;

	 ndm = (theNodalTA[0]->getCrds()).Size();
	 
	 NodalLocs.Zero();
	 NodalLocs.resize(6,ndm);

	 for (int i=0; i<6; i++){
		 for (int j=0;j<ndm; j++){
			NodalLocs(i,j)=(theNodalTA[i]->getCrds())(j);
		}
	 }

	 if((theNodalTA[0]->getThermalActionType())!= (theNodalTA[5]->getThermalActionType())){
		 opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is not consistent"<<endln;
	 }
	 else{
		 if(theNodalTA[0]->getThermalActionType()==1){
			NumData=9; //This is for 2D beam
		 }
		 else if( theNodalTA[0]->getThermalActionType()==2  ){
			NumData = 15; //This is for 3D I-section Beam
		 }
		 else {
            opserr<<"WARNING::ThermalActionWrapper: Thermal action type derived from NodalThermalAction is unable to be identified"<<endln;
		 }
	  }
		 
}


ThermalActionWrapper::~ThermalActionWrapper()
{
  
}

int 
ThermalActionWrapper::setRatios(const Vector& theRatio){
	//Clean the ratios
	if(theRatio!=0){
		theRatios.Zero();
		theRatios = theRatio;
	}
	else{
		opserr<<"WARNIGN::ThermalActionWrapper received invalid ratios"<<endln;
		return -1;
	}
  
  //check the num of interpolation points
  if (theRatios.Size()!=NodalLocs.noRows()) {
    opserr<<"WARNIGN::ThermalActionWrapper received an incompatible ratio"<<endln;
    return -2;
  }
  
  int NumNodTA = NodalLocs.noRows();
  
  if(theRatios(NumNodTA-1)>2){
    opserr<<"WARNING::ThermalActionWrapper received a ration vector ends up with "<<theRatios(NumNodTA-1)<< " , which should be 1.0 or 2.0"<<endln;
    return -2;
  }
  else if(theRatios(0)<0){
	  opserr<<"WARNING::ThermalActionWrapper received a ration vector ends up with "<<theRatios(NumNodTA-1)<< " , which should be 0 or greater"<<endln;
    return -2;
  }
  else {
    if(theRatios(0)>0)
		ConstLoc = theRatios(0);

	for(int i=1; i<NumNodTA-1;i++){
		if(theRatios(i)<-1e-6){
			if(Transpoint !=0)
				opserr<<"WARNING::ThermalActionWrapper over defined the transition point"<<endln;
			else
				Transpoint= (-1)*theRatios(i);
		}

	}
    //for (int i=1; i<NumNodTA-1; i++){
     // for (int j=0;j<ndm; j++){
     //   NodalLocs(i,j)= NodalLocs(0,j)+ theRatios(i)*(NodalLocs(NumNodTA-1,j)-NodalLocs(0,j));
     // }
    //}
    //end of updating locs
    
  }

	return 0;
}

const Vector &
ThermalActionWrapper::getData(int &type, double loadFactor)
{
  type = ThermalActionType;

  Vector data = 0;

	return data;
}






const Vector&
ThermalActionWrapper::getIntData(const Vector& locs)
{
 int NumNodalTA = NodalLocs.noRows();
 
 int ndm = NodalLocs.noCols();
 //Vector IntData;
 IntData.Zero();
  double ratio[6];
  double r1, r2, r3, r4, r5;
  double constStart=0;
  double constEnd =1.0;
  
#ifdef _DEBUG
	 //opserr<<locs<<endln;
#endif

  double sQdistInt =0;double sQr=0 ;double sQr1=0; 
  for(int i=0;i<ndm;i++){
		sQdistInt += (locs(i)- NodalLocs(0,i))*(locs(i)- NodalLocs(0,i));
		sQr += (NodalLocs(NumNodalTA-1,i)-NodalLocs(0,i))*(NodalLocs(NumNodalTA-1,i)-NodalLocs(0,i));
		//sQr1 += (NodalLocs(1,i)-NodalLocs(0,i))*(NodalLocs(1,i)-NodalLocs(0,i));
   }
   
   double distInt = sqrt(sQdistInt);
   double  r = sqrt(sQr);
   constStart = ConstLoc*r;
   constEnd = r;
   if(Transpoint>1e-6)
	   Transpoint= Transpoint*r;

 if(NumNodalTA==2){
	 // obtain the ratio

	if(theRatios(1)-1<0)
		constEnd= theRatios(1)*r;
   
   //ratio[0]= (locs(crdi)-NodalLocs(1,crdi))/(NodalLocs(0,crdi)-NodalLocs(1,crdi));
   //ratio[1]= (locs(crdi)-NodalLocs(0,crdi))/(NodalLocs(1,crdi)-NodalLocs(0,crdi));

   ratio[0]= (distInt-r)/(constStart-r);
   ratio[1]= (distInt-constStart)/(r-constStart);
   
   int type;
   
   Vector data0 = theNodalTA[0]->getData( type);
   Vector data1 = theNodalTA[1]->getData( type);

#ifdef _DEBUG
 		//opserr<<"NodalT0: "<<data0<<endln<<"NodalT1: "<<data1;
#endif
   
   if(NumData==9){
     IntData.resize(18);
     
     for(int i =0; i<9;i++){
       if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6){
         opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
       }
       else{
		   if(distInt<constStart){
			IntData(2*i+1)=data0(2*i+1);
			IntData(2*i)=data0(2*i);  // reversed ratio
			}
		   else if(distInt>constEnd){
			IntData(2*i+1)=data1(2*i+1);
			IntData(2*i)=data1(2*i);  // reversed ratio
			}
		   else{
			IntData(2*i+1)=data0(2*i+1);
			IntData(2*i)=data0(2*i)*ratio[0]+data1(2*i)*ratio[1];  // reversed ratio
		   }
         
       }
     }
     
   }
   //end of numData==9;
   else if(NumData ==15)
   {
     IntData.resize(25);
     for(int i =0; i<5;i++){
       if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6&&fabs(data0(3*i+12)-data1(3*i+12))>1e-6){
         opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
       }
       else{
		   if(distInt<constStart){
			IntData(2*i) = data0(2*i);          //5 temps through y
			IntData(2*i+1)= data0(2*i+1);               //5 locs through y
			IntData(3*i+10) =data0(3*i+10);       //5 temps through Z in bottom flange
			IntData(3*i+11)= data0(3*i+10);       //5 temps through Z in top flange
			IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
			}
		   else if(distInt>constEnd){
			IntData(2*i) = data1(2*i);          //5 temps through y
			IntData(2*i+1)= data1(2*i+1);               //5 locs through y
			IntData(3*i+10) =data1(3*i+10);       //5 temps through Z in bottom flange
			IntData(3*i+11)= data1(3*i+10);       //5 temps through Z in top flange
			IntData(3*i+12)= data1(3*i+12);            //5 locs through Z
			}
		   else{
			IntData(2*i) = data0(2*i)*ratio[0]+data1(2*i)*ratio[1];            //5 temps through y
			IntData(2*i+1)= data0(2*i+1);               //5 locs through y
			IntData(3*i+10) =data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1];       //5 temps through Z in bottom flange
			IntData(3*i+11)= data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1];       //5 temps through Z in top flange
			IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
		   }
         
       }
     }
   }
   //end of numData==18;

 }
 else if(NumNodalTA==3){
   //for 3-point interpolation
	 // obtain the ratio

    #ifdef _DEBUG
 		//opserr<<"NodalLocs: "<<NodalLocs<<endln;
	#endif
  r1 = theRatios(1)*r;
   ratio[0]= (distInt-r1)*(distInt-r)/(constStart-r1)/(constStart-r);
   ratio[1]= (distInt-constStart)*(distInt-r)/(r1-constStart)/(r1-r);
   ratio[2] = (distInt-constStart)*(distInt-r1)/(r-constStart)/(r-r1);
   // ratio[2]= (locs(crdi)-NodalLocs(0,crdi))*(locs(crdi)-NodalLocs(1,crdi))/(NodalLocs(2,crdi)-NodalLocs(0,crdi))/(NodalLocs(2,crdi)-NodalLocs(1,crdi));
		
	   int type;
   
	   Vector data0 = theNodalTA[0]->getData( type);
     Vector data1 = theNodalTA[1]->getData( type);
	   Vector data2 = theNodalTA[2]->getData( type);
	   #ifdef _DEBUG
 		//opserr<<"NodalT0: "<<data0<<endln<<"NodalT2: "<<data1;
		#endif

	 if(NumData==9){
		 IntData.resize(18);
		
		for(int i =0; i<9;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data0(2*i);  // reversed ratio
				}
				else{
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2];  // reversed ratio
				}
		   }
		 }

	 }
	 //end of numData==9;
	 else if(NumData ==15)
	 {
		 IntData.resize(25);	
		for(int i =0; i<5;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6&&fabs(data0(3*i+12)-data1(3*i+12))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i) = data0(2*i);            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10);       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+10);       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
				else{
				IntData(2*i) = data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2];            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1]+data2(3*i+10)*ratio[2];       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1]+data2(3*i+10)*ratio[2];       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
		   }
		 }
	 }
	 //end of numData==18;


 }
 else if(NumNodalTA==4){
    r1 = theRatios(1)*r; r2 = theRatios(2)*r;
   ratio[0]= (distInt-r1)*(distInt-r)*(distInt-r2)/(constStart-r1)/(constStart-r)/(constStart-r2);
   ratio[1]= (distInt-constStart)*(distInt-r)*(distInt-r2)/(r1-constStart)/(r1-r)/(r1-r2);
   ratio[2] = (distInt-constStart)*(distInt-r)*(distInt-r1)/(r2-constStart)/(r2-r)/(r2-r1);
   ratio[3] = (distInt-constStart)*(distInt-r1)*(distInt-r2)/(r-constStart)/(r-r1)/(r-r2);
  
   // ratio[2]= (locs(crdi)-NodalLocs(0,crdi))*(locs(crdi)-NodalLocs(1,crdi))/(NodalLocs(2,crdi)-NodalLocs(0,crdi))/(NodalLocs(2,crdi)-NodalLocs(1,crdi));
		
	   int type;
   
	   Vector data0 = theNodalTA[0]->getData( type);
     Vector data1 = theNodalTA[1]->getData( type);
	   Vector data2 = theNodalTA[2]->getData( type);
	   Vector data3 = theNodalTA[3]->getData( type);
	   #ifdef _DEBUG
 		//opserr<<"NodalT0: "<<data0<<endln<<"NodalT2: "<<data1;
		#endif

	 if(NumData==9){
		 IntData.resize(18);
		
		for(int i =0; i<9;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data0(2*i);  // reversed ratio
				}
				else{
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2]+data3(2*i)*ratio[3];  // reversed ratio
				}
		   }
		 }

	 }
	 //end of numData==9;
	 else if(NumData ==15)
	 {
		 IntData.resize(25);	
		for(int i =0; i<5;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6&&fabs(data0(3*i+12)-data1(3*i+12))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i) = data0(2*i);            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10);       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+10);       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
				else{
				IntData(2*i) = data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2]+data3(2*i)*ratio[3];            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1]+data2(3*i+10)*ratio[2]+data3(3*i+10)*ratio[3];       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+11)*ratio[0]+data1(3*i+11)*ratio[1]+data2(3*i+11)*ratio[2]+data3(3*i+11)*ratio[3];       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
		   }
		 }
	 }
	 //end of numData==18;
 }
 //for num nodes =4

 else if(NumNodalTA==5){
    r1 = theRatios(1)*r; r3 = theRatios(3)*r; 
   if(distInt>constStart&&distInt<Transpoint){
	ratio[0]= (distInt-r1)*(distInt-Transpoint)/(constStart-r1)/(constStart-Transpoint);
	ratio[1]= (distInt-constStart)*(distInt-Transpoint)/(r1-constStart)/(r1-Transpoint);
	ratio[2] = (distInt-constStart)*(distInt-r1)/(Transpoint-constStart)/(Transpoint-r1);
   }
   else if(distInt>Transpoint&&distInt<r){
	ratio[0]= (distInt-r3)*(distInt-r)/(Transpoint-r3)/(Transpoint-r);
	ratio[1]= (distInt-Transpoint)*(distInt-r)/(r3-Transpoint)/(r3-r);
	ratio[2] = (distInt-Transpoint)*(distInt-r3)/(r-Transpoint)/(r-r3);
   }
   else{
	opserr<<"WARNIGN ThermalActionWrapper received a int Value "<<distInt<<" out of range"<<endln;
	   }
  
   // ratio[2]= (locs(crdi)-NodalLocs(0,crdi))*(locs(crdi)-NodalLocs(1,crdi))/(NodalLocs(2,crdi)-NodalLocs(0,crdi))/(NodalLocs(2,crdi)-NodalLocs(1,crdi));
		
	   int type;
   
	   Vector data0 = theNodalTA[0]->getData( type);
     Vector data1 = theNodalTA[1]->getData( type);
	   Vector data2 = theNodalTA[2]->getData( type);
	   Vector data3 = theNodalTA[3]->getData( type);
	    Vector data4 = theNodalTA[4]->getData( type);
	   #ifdef _DEBUG
 		//opserr<<"NodalT0: "<<data0<<endln<<"NodalT2: "<<data1;
		#endif

	 if(NumData==9){
		 IntData.resize(18);
		
		for(int i =0; i<9;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data0(2*i);  // reversed ratio
				}
				else if(distInt<Transpoint){
				IntData(2*i+1)=data3(2*i+1);
				IntData(2*i)=data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2];  // reversed ratio
				}
				else{
				IntData(2*i+1)=data0(2*i+1);
				IntData(2*i)=data2(2*i)*ratio[0]+data3(2*i)*ratio[1]+data4(2*i)*ratio[2];  // reversed ratio
				}
		   }
		 }

	 }
	 //end of numData==9;
	 else if(NumData ==15)
	 {
		 IntData.resize(25);	
		for(int i =0; i<5;i++){
			if(fabs(data0(2*i+1)-data1(2*i+1))>1e-6&&fabs(data0(3*i+12)-data1(3*i+12))>1e-6){
				opserr<<"Warning:The NodalThermalAction in dispBeamColumn2dThermalNUT "<<this->getTag()
			      << "incompatiable loc input for datapoint "<< i << endln;
			 }
			else{
				if(distInt<constStart){
				IntData(2*i) = data0(2*i);            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10);       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+10);       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
				else if(distInt<Transpoint){
				IntData(2*i) = data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data1(2*i)*ratio[2];            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1]+data2(3*i+10)*ratio[2]+data3(3*i+10)*ratio[3];       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+11)*ratio[0]+data1(3*i+11)*ratio[1]+data2(3*i+11)*ratio[2]+data3(3*i+11)*ratio[3];       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
				else
				{
				IntData(2*i) = data0(2*i)*ratio[0]+data1(2*i)*ratio[1]+data2(2*i)*ratio[2]+data3(2*i)*ratio[3];            //5 temps through y
				IntData(2*i+1)= data0(2*i+1);               //5 locs through y
				IntData(3*i+10) =data0(3*i+10)*ratio[0]+data1(3*i+10)*ratio[1]+data2(3*i+10)*ratio[2]+data3(3*i+10)*ratio[3];       //5 temps through Z in bottom flange
				IntData(3*i+11)= data0(3*i+11)*ratio[0]+data1(3*i+11)*ratio[1]+data2(3*i+11)*ratio[2]+data3(3*i+11)*ratio[3];       //5 temps through Z in top flange
				IntData(3*i+12)= data0(3*i+12);            //5 locs through Z
				}
		   }
		 }
	 }
	 //end of numData==18;
 }
 //for num nodes =4

 else{

 }
#ifdef _DEBUG
 //opserr<<"ThermalActionWrapper:Interpolated TempData: "<<IntData<<endln;
#endif
	return IntData;
}




int 
ThermalActionWrapper::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int 
ThermalActionWrapper::recvSelf(int commitTag, Channel &theChannel,  
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

// do it later
void 
ThermalActionWrapper::Print(OPS_Stream &s, int flag)
{
  s << "ThermalActionWrapper"<< endln;
}

