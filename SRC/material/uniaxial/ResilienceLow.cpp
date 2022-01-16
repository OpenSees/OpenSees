// Written: Q. Gu  and Z.Zeng
// Created: 2014
//Reference:Seismic Performance of SRC Columns with High Ratio of Encased Steel written by YIN Xiaowei,supervised by Prof. LU Xilin,Tongji University 2012 
// Description: This file contains the class implementation for 
// ResilienceLow.
#include <ResilienceLow.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
void * OPS_ADD_RUNTIME_VPV(OPS_ResilienceLow) {

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int    iData[1];
  double dData[5];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterialtag" << endln;
    return 0;
  }
  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 5) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceLow " << iData[0] << "  PY DPmax Pmax Ke Kd" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceLow " << iData[0] << "  PY DPmax Pmax Ke Kd" << endln;
    return 0;
  }

  

  // Parsing was successful, allocate the material
  theMaterial = new ResilienceLow(iData[0], dData[0], dData[1], 
			    dData[2], dData[3], dData[4] 
			    );

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ResilienceLow\n";
    return 0;
  }
     
  return theMaterial;
}
ResilienceLow::ResilienceLow(int tag, double pPY, double pDPmax, double pPmax, double pKe,double pKd)
  :UniaxialMaterial(tag,MAT_TAG_ResilienceLow),
   strain(0.0), stress(0.0), tangent(0.0),
   Cstrain(0.0), Cstress(0.0), Ctangent(0.0)
{


	this->PY = pPY;
	this->DPmax = pDPmax;
	this->Pmax = pPmax;
	this->Ke= pKe;
	this->Kd= pKd;

	mode = 1;
	Cmode =1;
	strainRFMode11=-DPmax;
	stressRFMode11=-Pmax;
	strainRFMode6=DPmax;
	stressRFMode6=Pmax;

	return;

}

ResilienceLow::~ResilienceLow()
{

	return; 
}

int 
ResilienceLow::setTrialStrain(double pStrain, double strainRate) //// determine trial stress and tangent!!!!!!!!!!!!!!!!
{
	  strain = pStrain;
	  stress = Cstress;
	  tangent =Ctangent;
	  mode = Cmode; 
	  Flag=CFlag;
	  Di=CDi;
	  DY=PY/Ke;
	 mode= this->determineState();
	 tangent=this->getTangent();
	return 0;
}
double ResilienceLow::getStress(void){
	return stress;
};
double ResilienceLow::getTangent(){
	switch (mode){
		case 1:{tangent=Ke;}// mode1����1,3
			break; 
		case 2:{tangent=(Pmax-PY)/(DPmax-DY);}	// ����1,2
			break;
		case 3:{tangent=(stressRFMode2+PY)/(strainRFMode2+DY);}	// ����3,5  mode3
			break;
		case 4:{tangent=(Pmax-PY)/(DPmax-DY);}	// mode4������3,4
			break;
		case 5:
			{tangent=(stressRFMode4-PY)/(strainRFMode4-DY);}		// mode5,����4,5
			break;
		case 6:
			{	if(stress>=0.55*Pmax){        
               tangent=-Kd;}
            else
			{  stress=0;}
				}// 6 ����6,7
			break;
		case 7:{tangent=Kui;}// mode 7 ����7,8
			break;
		case 8:{tangent=Kri;}
			break;
		case 9:{tangent=Kui;}
			break;
		case 10:{tangent=Kri;}
			break;
		case 11:
			{	if(stress<=-0.55*Pmax)        
			{ tangent=-Kd;}
            else
			{stress=0;}
				}
			break;
	}// switch}
	return tangent;
};
 double ResilienceLow::getInitialTangent(void)
{
  return Ke;
}

 /*double ResiliencelMaterial::getTangent(void)
{

 double savedStrain = strain;
 double stress1 = this->getStress();
 strain += 1e-7;
 this->setTrialStrain(strain, 0);
 double stress2 = this->getStress();
 tangent = (stress2-stress1)/1e-7;

 this->setTrialStrain(savedStrain, 0);

  return tangent;
}*/

double  ResilienceLow::getStrain(void)
{
  return strain;
}
int 
ResilienceLow::commitState(void)
{
  Cstrain  = strain;
  Cstress  = stress;
  Ctangent = tangent;
  Cmode = mode; 
  CFlag=Flag;
  CKui=Kui;
  CKri=Kri;
  CDi=Di;
  CstrainRFMode2 = strainRFMode2;
  CstressRFMode2 = stressRFMode2;
  CstrainRFMode4 = strainRFMode4;
  CstressRFMode4 = stressRFMode4;
  CstrainRFMode6 = strainRFMode6;
  CstressRFMode6 = stressRFMode6; 
  CstrainRFMode8 = strainRFMode8;
  CstressRFMode8 = stressRFMode8;
  CstrainRFMode10 = strainRFMode10;
  CstressRFMode10 = stressRFMode10;
  CstrainRFMode11 = strainRFMode11;
  CstressRFMode11 = stressRFMode11;

  return 0;
}
int ResilienceLow::revertToLastCommit(void){return 0;}
int ResilienceLow::revertToStart(void){return 0;}

UniaxialMaterial *ResilienceLow::getCopy(void){
  ResilienceLow *theCopy = new ResilienceLow(this->getTag(), PY, DPmax, Pmax,Ke, Kd);
  return theCopy;
}

int ResilienceLow::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int ResilienceLow::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void ResilienceLow::Print(OPS_Stream &s, int flag)
{
  s << "ResilienceLow : " << this->getTag();
  return;
}

int ResilienceLow::determineState(){
   double dd = strain - Cstrain; 
   if (fabs(dd)<1.0e-14) {
      stress = Cstress;
	  tangent = Ctangent;
       return mode;
   }
   switch (mode){
   // ��һ�μ���

// mode1����1,3
	case 1:
		if(dd>=0){
			if(strain>DY){        
                mode=2;
				determineState( );}
            else
                stress=Ke*strain;
		}
		else {
			if(strain<-DY){
                mode=4;
				determineState( );}
            else
                stress = Ke*strain;
         }
		break; 
        
// ����1,2
	case 2:
		if(dd>=0){
			if(strain>DPmax){
                mode=6;
				determineState( );}
            else
                stress=(Pmax-PY)/(DPmax-DY)*(strain-DY)+PY;
		}
		else{
            strainRFMode2=Cstrain;
            stressRFMode2=Cstress;
            mode=3;
            determineState( );
		}
		break;
       
// ����3,5  mode3
	case 3:
		if(dd>=0)
			{
				if(strain>strainRFMode2)
				{
                mode=2;
				determineState( );
				}
				else
                stress=(stressRFMode2+PY)/(strainRFMode2+DY)*(strain-strainRFMode2)+stressRFMode2;
			}
        else
		{
			if(strain<-DY)
			{
                mode=4;
				determineState( );
			}
            else
                stress=(stressRFMode2+PY)/(strainRFMode2+DY)*(strain-strainRFMode2)+stressRFMode2;
		}
        break;
        
// mode4������3,4
	case 4:
		if(dd>=0)
		{
            strainRFMode4=Cstrain;
            stressRFMode4=Cstress;
            mode=5;
			determineState( );
		}           
        else
		{
			if(strain<-DPmax)
			{
                mode=11;
                determineState( ); 
			}
            else
                stress=(Pmax-PY)/(DPmax-DY)*(strain+DPmax)-Pmax;
		}
		break;
        
        // mode5,����4,5
	case 5:
        if(dd>=0)
		{
            if(strain>DY)
			{mode=2;
             determineState( );}
            else
                stress=(stressRFMode4-PY)/(strainRFMode4-DY)*(strain-DY)+PY;
		}
        else
		{
			if(strain<strainRFMode4)
			{mode=4;
			determineState( );}
            else
                stress=(stressRFMode4-PY)/(strainRFMode4-DY)*(strain-DY)+PY;
		}
        break;
        // 6 ����6,7
	case 6:
		if(dd>=0)
		{stress=-Kd*(strain-DPmax)+Pmax;
			if(strain>Di)
			{Di=strain;}
			 if(stress<0.55*Pmax)
			 { stress=0.55*Pmax;}
		}
		else
		{
            strainRFMode6=Cstrain;
            stressRFMode6=Cstress;
			Flag=6;
            mode=7;
			determineState( );
		}
        break;
        
        // mode 7 ����7,8

	case 7:
		if(strain>0.45*Pmax/Kd+DPmax)
			{ Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
			Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
        else
			{Kri=(0.5*Pmax-Kd*(Di-DPmax))/(Di-0.5*Pmax/Ke);
			Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
		if(Kri<0)
		{   Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
            Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di); }
        if(Kui>Ke||Kui<0)
		{Kui=Ke;}

		if(Flag==6)
			{stress=Kui*(strain-strainRFMode6)+stressRFMode6;
				if(dd>=0){
					if(strain>strainRFMode6){
                         mode=6;
						 determineState( );}}
				else 
				{if(stress<Kri*(strain+0.5*Pmax/Ke)-0.5*Pmax)
				{ mode=8;
                     determineState( ); }
					}}
      else 
		 { stress=Kui*(strain-strainRFMode8)+stressRFMode8;
			if(dd>=0){
                  if(stress>Kri*(strain-0.5*Pmax/Ke)+0.5*Pmax)
				  { mode=10;
				  determineState( );}
				  stressRFMode6=-Kd*(strainRFMode6-DPmax)+Pmax;
				if(stressRFMode6<0.55*Pmax)
					stressRFMode6=0.55*Pmax;
				  if(strain>strainRFMode6||stress>stressRFMode6){
                         mode=6;
                         determineState( );
				  }}
              else
			  { if(strain<strainRFMode8){
                     mode=8;
                     determineState( );
			}}}
 break;


	case 8:
			if(dd>=0)
			{		strainRFMode8=Cstrain;
					stressRFMode8=Cstress;
					Flag=8;
					mode=7;
					determineState(); }
				else
				{   if(strain<-0.45*Pmax/Kd-DPmax){
					  Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
					  Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
				else{
					  Kri=(0.5*Pmax-Kd*(Di-DPmax))/(Di-0.5*Pmax/Ke);
					  Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);
				}
				if(Kri<0){
					Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
					Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);
				}
					stress=Kri*(strain+0.5*Pmax/Ke)-0.5*Pmax;
			 if(strain<-DPmax&&stress<-Kd*(strain+DPmax)-Pmax)
			 {mode=11;
			 determineState( );}}

	break;

case 9:
		if(strain<-0.45*Pmax/Kd-DPmax){
              Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
			  Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
		else{
              Kri=(0.5*Pmax-Kd*(Di-DPmax))/(Di-0.5*Pmax/Ke);
              Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
		if(Kri<0){
            Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
            Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
            if(Kui>Ke||Kui<0)
			{  Kui=Ke;}

			if(Flag==11){
                stress=Kui*(strain-strainRFMode11)+stressRFMode11;
				if(dd<0){
					if(strain<strainRFMode11){
                        mode=11;
                        determineState( );
					}}
                else
				{if(stress>Kri*(strain-0.5*Pmax/Ke)+0.5*Pmax){
                        mode=10;
                        determineState( );
					}}}
            else
			{stress=Kui*(strain-strainRFMode10)+stressRFMode10;
             if(dd>=0){
					if(strain>strainRFMode10){
                        mode=10;
                        determineState( );
					}}
				else{
                    if(stress<Kri*(strain+0.5*Pmax/Ke)-0.5*Pmax)
					{  mode=8;
                          determineState( );}
					stressRFMode11=-Kd*(strainRFMode11+DPmax)-Pmax;
					if(stressRFMode11>-0.55*Pmax)
					stressRFMode11=-0.55*Pmax;
					if(strain<strainRFMode11||stress<stressRFMode11){
                        mode=11;
						determineState( );}}}
break;

case 10:
		if(dd>=0)
		{if(strain>0.45*Pmax/Kd+DPmax){
              Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
			  Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);}
            else
			{ Kri=(0.5*Pmax-Kd*(Di-DPmax))/(Di-0.5*Pmax/Ke);
              Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);
		}
		if(Kri<0){
            Kri=0.05*Pmax/(Di-0.5*Pmax/Ke);
            Kui=(1.25*Pmax-Kd*(Di-DPmax))/(0.5*Pmax/Ke-0.25*Pmax/Kri+Di);
		}
            stress=Kri*(strain-0.5*Pmax/Ke)+0.5*Pmax;
            
            if(strain>DPmax&&stress>-Kd*(strain-DPmax)+Pmax)
			{mode=6;
                 determineState();          
			}}
        else
		{	strainRFMode10=Cstrain;
            stressRFMode10=Cstress;
            Flag=10;
            mode=9;
			determineState();}

break;

case 11:
        if(dd>=0)
		{
            strainRFMode11=Cstrain;
            stressRFMode11=Cstress;
            Flag=11;
			 mode=9;
			determineState();
		}
        else
		{
           stress=-Kd*(strain+DPmax)-Pmax;
		   if(fabs(strain)>Di){
                Di=fabs(strain);
		   }
		   if(stress>-0.55*Pmax){
                stress=-0.55*Pmax;
                }}
		break;
  } 
  return mode;
};
