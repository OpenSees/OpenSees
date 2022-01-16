 /* ********************************************************************
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center	      **
**                                                                    **
** (C) Copyright by Quan Gu and Yongdou Liu @ Xiamen University	      **
**								      **
** ****************************************************************** */
// Created: 2014

// Referenced from "Seismic Performance of SRC Columns with High Ratio of Encased Steel ".Prof. LU Xilin.    and  YIN Xiaowei 
// Description: This file contains the class implementation for ResilienceMaterialHR

#include <ResilienceMaterialHR.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <Information.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <OPS_Globals.h>
void * OPS_ADD_RUNTIME_VPV(OPS_ResilienceMaterialHR)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  int    iData[1];
  double dData[7];
  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ResilienceMaterialHR tag" << endln;
    return 0;
  }
  numData = OPS_GetNumRemainingInputArgs();
  if (numData != 7) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceMaterialHR " << iData[0] << " DY PY DPmax Pmax Ke Kd coefficient" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args, want: uniaxialMaterial ResilienceMaterialHR " << iData[0] << " DY PY DPmax Pmax Ke Kd coefficient" << endln;
    return 0;
  }

  

  // Parsing was successful, allocate the material
  theMaterial = new ResilienceMaterialHR(iData[0], dData[0], dData[1], 
			    dData[2], dData[3], dData[4], 
			    dData[5], dData[6]);

  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ResilienceMaterialHR\n";
    return 0;
  }

  return theMaterial;
}

ResilienceMaterialHR::ResilienceMaterialHR(int tag, double  pDY, double pPY, double pDPmax, double pPmax, double pKe,double pKd,double pK)
  :UniaxialMaterial(tag,MAT_TAG_ResilienceMaterialHR),
   strain(0.0), stress(0.0), tangent(0.0),
   Cstrain(0.0), Cstress(0.0), Ctangent(0.0)
{
  this->DY = pDY;
  this->PY = pPY;
  this->DPmax = pDPmax;
  this->Pmax = pPmax;
  this->Ke= pKe;
  mode = 0;
  Cmode =0;
  strainRFMode13=-DPmax;
  stressRFMode13=-Pmax;
  strainRFMode6=DPmax;
  stressRFMode6=Pmax;
  coefficient = pK; 
  Kd = pKd;
  Ku=Ke;  
  Kr=coefficient*Ke;
}

ResilienceMaterialHR::~ResilienceMaterialHR()
{
  return; 
}

//------main part of the ResilenseMaterial which is aimed at  determining stress and tangent------------------------------------
int 
ResilienceMaterialHR::setTrialStrain(double pStrain, double strainRate) 
{
  strain = pStrain;
  stress = Cstress;
  tangent =Ctangent;
  mode = Cmode; 
  mode= this->determineState();
  tangent=this->getTangent();
  return 0;
}

double 
ResilienceMaterialHR::getStress(void){
  return stress;
};

double 
ResilienceMaterialHR::getTangent(){
  switch (mode){
  case 1:{tangent=Ke;}
    break; 
  case 2:{tangent=(Pmax-PY)/(DPmax-DY);}
    break;
  case 3:{tangent=(stressRFMode2+PY)/(strainRFMode2+DY);}
    break;
  case 4:{tangent=(Pmax-PY)/(DPmax-DY);}	
    break;
  case 5:
    {tangent=(stressRFMode4-PY)/(strainRFMode4-DY);}
    break;
  case 6:{tangent=-Kd;}
    break;
  case 7:{tangent=Ku;}
    break;
  case 8:{tangent=Kr;}
    break;
  case 9:
    {tangent=(-0.85*Pmax-stressRFMode13)/(strainP9-strainRFMode13);}
    break;
  case 10:{tangent=Ku;}
    break;
  case 11:{tangent=Kr;}
    break;
  case 12:{tangent=(0.85*Pmax-stressRFMode6)/(strainP12-strainRFMode6);}
    break;
  case 13:{tangent=-Kd;}
    break;
  }// switch}
  return tangent;
};

double 
ResilienceMaterialHR::getInitialTangent(void)
{
  return Ke;
}

double 
ResilienceMaterialHR::getStrain(void)
{
  return strain;
}

int 
ResilienceMaterialHR::commitState(void)
{
  Cstrain  = strain;
  Cstress  = stress;
  Ctangent = tangent;
  Cmode = mode; 
  CstrainP8 = strainP8;
  CstrainP9 = strainP9;
  CstrainP11 = strainP11;
  CstrainP12 = strainP12;
  CstrainP10 = strainP10;
  CstrainP7 = strainP7;
  CstrainRFMode9 = strainRFMode9;
  CstressRFMode9 = stressRFMode9;
  CstrainRFMode12 = strainRFMode12;
  CstressRFMode12 = stressRFMode12;
  CstrainRFMode2 = strainRFMode2;
  CstressRFMode2 = stressRFMode2;
  CstrainRFMode4 = strainRFMode4;
  CstressRFMode4 = stressRFMode4;
  CstrainRFMode6 = strainRFMode6;
  CstressRFMode6 = stressRFMode6;
  CstressP7 = stressP7;
  CstressP10 = stressP10;
  CstrainRFMode13 = strainRFMode13;
  CstressRFMode13 = stressRFMode13;
  return 0;
}

int 
ResilienceMaterialHR::revertToLastCommit(void)
{
  strain = Cstrain;
  stress = Cstress;
  tangent =Ctangent;
  mode = Cmode; 
  return 0;
}

int 
ResilienceMaterialHR::revertToStart(void)
{
  return -1;
}


UniaxialMaterial *
ResilienceMaterialHR::getCopy(void)
{
  ResilienceMaterialHR *theCopy = new ResilienceMaterialHR(this->getTag(),  DY, PY, DPmax, Pmax, Ke, Kd,coefficient);
  return theCopy;
}

int 
ResilienceMaterialHR::sendSelf(int cTag, Channel &theChannel)
{
  opserr << "ResilienceMaterialHR::sendSelf() - does not send self\n";
  return -1;
}

int 
ResilienceMaterialHR::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
ResilienceMaterialHR::Print(OPS_Stream &s, int flag)
{
  s << "ResilienceMaterialHR : " << this->getTag();
  return;
}

//-----------the core of the ResilienceMaterialHR which compute the stress under definite strian and strain history-------------
int 
ResilienceMaterialHR::determineState()
{
  double dd = strain - Cstrain; 
  if (fabs(dd)<1.0e-14) {
    stress = Cstress;
    tangent = Ctangent;
    return mode;
  }
  switch (mode){
  case 0:
    if(fabs(strain)<=DY) {
      mode=1;
      stress=Ke*strain;}
    else if(strain>DY&&strain<=DPmax){
      mode=2;
      stress=(Pmax-PY)/(DPmax-DY)*(strain-DY)+PY;}
    else if(strain>DPmax){
      mode=6;
      stress =-Kd*(strain-DPmax)+Pmax;}
    else if(strain < -DY && strain >= -DPmax){
      mode=4;
      stress=(Pmax-PY)/(DPmax-DY)*(strain+DPmax)-Pmax;}
    else{
      mode=13;
      stress = -Kd*(strain+DPmax)-Pmax;
    }
    break; 
    
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
    

  case 2:
    if(dd>0){
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
    

  case 3:
    if(dd>0){
      if(strain>strainRFMode2){
	mode=2;
	determineState( );
      }
      else
	stress=(stressRFMode2+PY)/(strainRFMode2+DY)*(strain-strainRFMode2)+stressRFMode2;
    }
    else{
      if(strain<-DY){
	mode=4;
	determineState( );
      }
      else
	stress=(stressRFMode2+PY)/(strainRFMode2+DY)*(strain-strainRFMode2)+stressRFMode2;
    }
    break;
    

  case 4:
    if(dd>0){
      strainRFMode4=Cstrain;
      stressRFMode4=Cstress;
      mode=5;
      determineState( );
    }           
    else
      {
	if(strain<-DPmax){
	  mode=13;
	  determineState( ); 
	}
	else
	  stress=(Pmax-PY)/(DPmax-DY)*(strain+DPmax)-Pmax;
      }
    break;
    
  case 5:
    if(dd>=0){
      if(strain>DY){
	mode=2;
	determineState( );
      }
      else
	stress=(stressRFMode4-PY)/(strainRFMode4-DY)*(strain-DY)+PY;
    }
    else{
      if(strain<strainRFMode4){
	mode=4;
	determineState( );
      }
      else
	stress=(stressRFMode4-PY)/(strainRFMode4-DY)*(strain-DY)+PY;
    }
    break;

  case 6:
    if(dd>=0){
      stress=-Kd*(strain-DPmax)+Pmax;
    }
    else
      {
	strainRFMode6=Cstrain;
	stressRFMode6=Cstress;
	strainP7=strainRFMode6;
	stressP7=stressRFMode6;
	mode=7;
	determineState( );
      }
    break;
    
  case 7:
    if(dd>0){
      if(strainP7>=strainRFMode6){
	if(strain>strainRFMode6){
	  mode=6;
	  determineState( );
	}
	else
	  stress=Ku*(strain-strainP7)+stressP7;
      }
      else
	{
	  if(strain>strainRFMode12){
	    mode=12;
	    determineState();
	  }
	  else
	    stress=Ku*(strain-strainP7)+stressP7;
	}                      
    }  //if dd>0
    else {
      strainP8=-stressP7/Ku+strainP7;
      if(strain>strainP8){
	stress=Ku*(strain-strainP7)+stressP7;
      }
      else{
	mode=8;
	determineState();
      }
    }
    break;
    
    
  case 8:
    if(dd>0){
      if(strain>strainP8){
	mode=7;
	determineState();
      }
      else
	stress=Kr*(strain-strainP8);
    }
    else{
      strainP9=-0.85*Pmax/Kr+strainP8;
      if(strain>strainP9){
	stress=Kr*(strain-strainP8);
      }
      else{
	mode=9;
	determineState();
      }
    }
    break;
    
  case 9:
    if(dd>0){
      strainRFMode9=Cstrain;
      stressRFMode9=Cstress;
      stressP10=stressRFMode9;
      strainP10=strainRFMode9;
      mode=10;
      determineState();
    }
    else{
      if(strain>strainRFMode13){
	stress=(-0.85*Pmax-stressRFMode13)/(strainP9-strainRFMode13)*(strain-strainRFMode13)+stressRFMode13;
      }
      else{
	mode=13;
	determineState();
      }
    }
    break;
    
  case 10:
    if(dd>0){
      strainP11=-stressP10/Ku+strainP10;
      if(strain>strainP11){
	mode=11;
	determineState();
      }
      else
	stress=Ku*(strain-strainP10)+stressP10;
    }
    else{
      if(strainP10<=strainRFMode13){
	if(strain<strainRFMode13){
	  mode=13;
	  determineState();
	}
	else
	  stress=Ku*(strain-strainP10)+stressP10;
      }
      else{
	if(strain<strainRFMode9){
	  mode=9;
	  determineState();
	}
	else
	  stress=Ku*(strain-strainP10)+stressP10;
      }
    }
    break;
    
  case 11:
    if(dd>0){
      strainP12=0.85*Pmax/Kr+strainP11;
      if(strain>strainP12){
	mode=12;
	determineState();
      }
      else
	stress=Kr*(strain-strainP11);
    }
    else{
      if(strain>strainP11){
	stress=Kr*(strain-strainP11);
      }
      else{
	mode=10;
	determineState();
      }
    }
    break;
    
  case 12:
    if(dd>0){
      if(strain>strainRFMode6){
	mode=6;
	determineState();
      }
      else
	stress=(0.85*Pmax-stressRFMode6)/(strainP12-strainRFMode6)*(strain-strainRFMode6)+stressRFMode6;
    }
    else{
      strainRFMode12=Cstrain;
      stressRFMode12=Cstress;
      strainP7=strainRFMode12;
      stressP7=stressRFMode12;
      mode=7;
      determineState();
    }
    break;
    
  case 13:
    if(dd>0){
      strainRFMode13=Cstrain;
      stressRFMode13=Cstress;
      strainP10=strainRFMode13;
      stressP10=stressRFMode13;
      mode=10;
      determineState();
    }
    else{
      stress=-Kd*(strain+DPmax)-Pmax;
    }
    break;
  } // switch
  
  return mode;
};
