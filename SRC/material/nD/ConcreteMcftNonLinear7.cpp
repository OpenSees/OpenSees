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
                                                                        
// $Revision: 1.19 $
// $Date: 2007/11/30 23:34:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01.cpp,v $
                                                                        
// Written: Murat 
// Created: 03/08
// 
// Description: This file contains the class implementation for 
// Concrete with . 
//


#include <NDMaterial.h>
#include <ConcreteMcftNonLinear7.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h> 
#include <Parameter.h>
#include <Information.h>
#include <float.h>

#include <elementAPI.h>

void *OPS_ConcreteMcftNonlinear7()
{
  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 9) {
    opserr << "ERROR not enough input args: nDMaterial ConcreteMcftNonlinear7 tag? fcu? ecu? Ec? fcr? Esv? fyv? alphaV? RoV?" << endln;
    return 0;
  }

  int tag;
  double dData[8];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "ERROR nDMaterial ConcreteMcftNonlinear7 - unable to read matTag" << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ERROR nDMaterial ConcreteMcftNonlinear7 - unable to read inputs" << endln;    
    return 0;
  }
  
  NDMaterial *theMaterial = new ConcreteMcftNonLinear7(tag,
						       dData[0], dData[1], dData[2], dData[3],
						       dData[4], dData[5], dData[6], dData[7]);
  if (theMaterial == 0) {
    opserr << "ERROR - could not create nDMaterial ConcreteMcftNonlinear7 with tag " << tag << endln;
    return 0;
  }
  
  return theMaterial;
}

// constructor: 
ConcreteMcftNonLinear7 ::ConcreteMcftNonLinear7  
(int tag, double fcui, double ecui, double Eci, double fcri, double Esvi, double fyvi, double alphaVi, double RoVi)
   :NDMaterial(tag, ND_TAG_ConcreteMcftNonLinear7 ), 
   fcu(fcui), ecu(ecui), Ec(Eci), fcr(fcri), Esv(Esvi), fyv(fyvi), alphaV(alphaVi), RoV(RoVi), epsf(2), Dri(2,2),
	Dr(2,2), sigf(2), fx(0.0), fxy(0.0), fy(0.0), DrP(2,2), e1(0.0), e2(0.0), sig1(0.0), sig2(0.0), theta(0.0), 
	exmin(0.0), exmax(0.0), ex(0.0),ey(0.0), exy(0.0), exymin(0.0), exymax(0.0),
	loadpath(0.0), loadpathP(0.0), crackLabelP(0.0), sigfsens(2),sigfsensP(2)

	

{



exP = 0.0;
eyP = 0.0;
exyP = 0.0;
fxP = 0.0; 
fxyP =0.0;
e1P = 0.0;
e2P = 0.0;

e1maxP = 0.0;
e2minP = 0.0;
e1max = 0.0;
e2min = 0.0;
fe1max = 0.0;
fe2min = 0.0;

exminP  =0.0;
exmaxP  =0.0;
eyminP  =0.0;
eymaxP  =0.0;
exyminP =0.0;
exymaxP =0.0;
exminc  =0.0;
exmaxc  =0.0;
eyminc  =0.0;
eymaxc  =0.0;
exyminc =0.0;
exymaxc =0.0;

Dr(0,0) = Ec;
Dr(0,1) = 0.0;
Dr(1,0) = 0.0;
Dr(1,1) = Ec/2;

 this->revertToLastCommit();
//opserr << " check1 " << endln;
//Sensitivity
parameterID = 0;
//SHVs = new Matrix(12,3);
 SHVs = 0;

//   SHVsP = 0;
//opserr << " check2 " << endln;
}


ConcreteMcftNonLinear7 ::ConcreteMcftNonLinear7 (void)
  :NDMaterial(0, ND_TAG_ConcreteMcftNonLinear7), epsf(2), Dri(2,2),
   Dr(2,2), sigf(2),sigfsens(2),
   e1(0.0), e2(0.0), sig1(0.0), sig2(0.0), theta(0.0), fx(0.0), fxy(0.0), fy(0.0), 
   exmin(0.0), exmax(0.0), ex(0.0),ey(0.0), exy(0.0), exymin(0.0), exymax(0.0), 
   loadpath(0.0), loadpathP(0.0), 	crackLabelP(0.0), sigfsensP(2)
 {
   exP = 0.0;
   eyP=0.0;
   exyP = 0.0;
   fxP = 0.0; 
   fxyP =0.0;
   e1P = 0.0;
   e2P = 0.0;
   e1maxP = 0.0;
   e2minP = 0.0;
   
   e1max = 0.0;
   e2min = 0.0;
   fe1max = 0.0;
   fe2min = 0.0;
   
   exminP = 0.0;
   exmaxP = 0.0;
   eyminP=0.0;
   eymaxP =0.0;
   exyminP = 0.0;
   exymaxP = 0.0;
   Dr(0,0) = Ec;
   Dr(0,1) = 0.0;
   Dr(1,0) = 0.0;
   Dr(1,1) = 0.5*Ec;
   exminc  =0.0;
   exmaxc  =0.0;
   eyminc  =0.0;
   eymaxc  =0.0;
   exyminc =0.0;
   exymaxc =0.0;
   
   DrP(0,0) = 0.0;
   DrP(0,1) = 0.0;
   DrP(1,0) = 0.0;
   DrP(1,1) = 0.0;
   
   
   this->revertToLastCommit();
   //Sensitivity
   parameterID = 0; 
   //SHVs = new Matrix(12,3);
   SHVs = 0;
}

ConcreteMcftNonLinear7 ::~ConcreteMcftNonLinear7 (void)
{
  if (SHVs != 0) 
    delete SHVs;	
}


NDMaterial*
ConcreteMcftNonLinear7 ::getCopy (void)
{
  ConcreteMcftNonLinear7  *theCopy;
  theCopy = new ConcreteMcftNonLinear7  (this->getTag(), fcu, ecu, Ec, fcr, Esv, fyv, alphaV, RoV);
  return theCopy;
}

NDMaterial*
ConcreteMcftNonLinear7 ::getCopy (const char *type)
{
  return this->getCopy();
}




int
ConcreteMcftNonLinear7 ::setTrialStrain (const Vector &strain)
{
  epsf = strain;
  return 0;
}

int
ConcreteMcftNonLinear7 ::setTrialStrain (const Vector &strain, const Vector &rate)
{
  epsf = strain;
  return 0;
}

int
ConcreteMcftNonLinear7 ::setTrialStrainIncr (const Vector &strain)
{
  epsf += strain;
  ////////opserr << " check7 " << endln;
  return 0;
}

int
ConcreteMcftNonLinear7 ::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsf += strain;
  ////////opserr << " check8 " << endln;
  return 0;
}


const Vector&
ConcreteMcftNonLinear7 ::getStrain (void)
{
  ////////opserr << " check9 " << endln;
  return epsf;
}


const Matrix&
ConcreteMcftNonLinear7 ::getTangent (void)
{
  //opserr << " check10 " << endln;
  //opserr << Dr <<endln; 
  return Dr;
}

const Matrix&
ConcreteMcftNonLinear7 ::getInitialTangent (void)
{
  Dri.Zero();
  Dri(0,0) = Ec;
  Dri(0,1) = 0.0;
  Dri(1,0) = 0.0;
  Dri(1,1) = Ec/2;
  
  ////////opserr << " check11 " << endln;
  return Dri;
}

const Vector&
ConcreteMcftNonLinear7 ::getStress (void)
{
  ////opserr << "GetStress ----------------------------" << endln;
  ////getchar();
  exmin = exminP;
  exmax = exmaxP;
  eymin = eyminP;
  eymax = eymaxP;
  exymin = exyminP;
  exymax = exymaxP;
  
  nE  = Ec/(Ec-fcu/ecu);
  
  //fiber strains
  ex  = epsf(0);
  exy = epsf(1);
  
  
  
  if(ex == 0.0 && exy == 0.0) {
    
    loadpath = 1.0;
    sigf(0) = 0.000;
    sigf(1) = 0.000;
    Dr(0,0) = Ec;
    Dr(0,1) = 1.0;
    Dr(1,0) = 1.0;
    Dr(1,1) = Ec/2;
    loadpath = 1.0;
    
    //////opserr << "loadpath = " <<  loadpath << endln;			
    
    //this->revertToStart();
    return sigf;
  } else if( (ex!=0.0 && exy!=0.0) && (ex==exP && exy==exyP)) {
    loadpath = 2.0;
    sigf(0) = fxP;
    sigf(1) = fxyP;
    Dr(0,0) = DrP(0,0);
    Dr(0,1) = DrP(0,1);
    Dr(1,0) = DrP(1,0);
    Dr(1,1) = DrP(1,1);
    
    crackLabel = crackLabelP;
    //////opserr << "loadpath = " <<  loadpath << endln;					
    //	 this->revertToLastCommit();
    return sigf;
  }  else if ((ex>0.10) || (ex<-0.10)){
    loadpath = 3.0;
    //////opserr << "loadpath = " <<  loadpath << endln;			
    Dr(0,0) =  0.0001;
    Dr(0,1) =  0.0;
    Dr(1,0) =  0.0;
    Dr(1,1) =  0.0001;
    
    if(ex<0) {
      sigf(0) = -0.0000001;
    } else {
      sigf(0) = 0.00000001;
    }
    if(exy<0) {
      sigf(1) = -0.0000001;
    } else {
      sigf(1) = 0.00000001;
    }
    
    return sigf;
  } else {
    
    Loadf();
    
    return sigf;
  }	 	
  
}

void
ConcreteMcftNonLinear7 ::Loadf(void)
{
	
  double exdelta = ex - exP;
  double exydelta = exy - exyP;
  //----------------------------------------------------------------------------------------------------------
  //compatibility equations depend on the sign of exy strain because e1 is always defined as Tension
  //where e2 is always defined as Compression
  
  sig1= 0;
  sig2 = 0;
  
  if ( (exy != 0) && (fabs(exy/ex) > 0.01) ) {
    loadpath = 4.1;
    //Update the principle axis Angle LOOP
    this->ForwardAngleSearch();
    
    Dr(0,0) = this->tangentstifness00(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV,e1max,e2min,fe1max,e1max,fe2min,e2min);
    Dr(0,1) = this->tangentstifness01(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV,e1max,e2min,fe1max,e1max,fe2min,e2min);
    Dr(1,0) = this->tangentstifness10(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV,e1max,e2min,fe1max,e1max,fe2min,e2min);
    Dr(1,1) = this->tangentstifness11(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV,e1max,e2min,fe1max,e1max,fe2min,e2min);
    
    sigf(0) = fx;
    sigf(1) = fxy;
  } else {
    
    if (ex<0) {
      //uniaxial compression
      loadpath = 4.2;
      e2 = ex;
      ey = 0.0;
      e1 = 0.0;
      fy = 0.0;
      FinalAnglex = 0.001;
      
      
      if (exP < 0) {
	e2P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*0.0) + exyP/2*sin(2*0.0);
	e1P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*0.0) - exyP/2*sin(2*0.0);
      } else {
	e1P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*0.0) + exyP/2*sin(2*0.0);
	e2P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*0.0) - exyP/2*sin(2*0.0);
      }
      
      //maximum and minimum strain envelopes on the same proncipal direction
      e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
      e2min = (exmin - eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);
      
      e2min = exmin;
      
      if(e2min<0 ) {
	fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));
      } else {
	fe2min = 0.0;
	
      }
      
      
      if (e2 <= e2min) {
	sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));
	
      } else {
	sig2 = fe2min + fe2min/e2min*(e2-e2min);
	
      }
      
      fx = sig2;
      fxy = Ec/2*exy;
      
      //recorded variables
      crackLabel = 0;	
      Strain1 = e1;
      Strain2 = e2;
      Sigma1 = sig1;
      Sigma2 = sig2;
      epsy = ey;
      
      
      sigf(0) = fx;
      sigf(1) = fxy;
      if (e2 <= e2min) {
	Dr(0,0) = -((e2*pow(Ec,2)*pow(e2/ecu,-1 + Ec/(Ec - fcu/ecu))*fcu)/(pow(ecu,2)*pow(Ec - fcu/ecu,2)*
									   pow(-1 + pow(e2/ecu,Ec/(Ec - fcu/ecu)) + Ec/(Ec - fcu/ecu),2))) + (Ec*fcu)/
	  (ecu*(Ec - fcu/ecu)*(-1 + pow(e2/ecu,Ec/(Ec - fcu/ecu)) + Ec/(Ec - fcu/ecu)));
	
	
      } else {
	Dr(0,0) = fe2min/e2min;
      }
      
      Dr(0,1) = 0.0;
      Dr(1,0) = 0.0;
      Dr(1,1) = Ec/2;
      
      


				

    } else if (ex>0) {
      //uniaxial tension
      loadpath = 4.3;
      e1 = ex;
      e2 = 0.0;
      ey = 0.0;
      fy = 0.0;
      FinalAnglex = 89.999;
      
      if (exP < 0) {
	e2P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*0.0) + exyP/2*sin(2*0.0);
	e1P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*0.0) - exyP/2*sin(2*0.0);
      } else {
	e1P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*0.0) + exyP/2*sin(2*0.0);
	e2P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*0.0) - exyP/2*sin(2*0.0);
      }
      
      //maximum and minimum strain envelopes on the same proncipal direction
      e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
      e2min = (exmin + eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);
      
      e1max = exmax;
      
      
      if (e1max >0 && e1max <= fcr/Ec) {
	fe1max = Ec*e1max;
	
      }  else if (e1max <0 ) {
	fe1max= 0.0;
	
      } else {
	fe1max = fcr / (1+sqrt(500.0*(e1max)));
	
	
      }
      
      
      //stress-strain
      //principal tension
      if (e1 >= e1max) {
	if (e1 <= fcr/Ec) {
	  sig1 = Ec*e1;
	} else {
	  sig1 = fcr / (1+sqrt(500.0*(e1)));
	}
      } else {
	sig1 = fe1max + fe1max/e1max*(e1-e1max);
	
      }
      
      fx = sig1;
      fy =  0.0;
      fxy = Ec/2*exy;
      
      if(e1 >= fcr/Ec) { 
	crackLabel = 1;
      } else {
	crackLabel = 0;
      }
      Strain1 = e1;
      Strain2 = e2;
      Sigma1 = sig1;
      Sigma2 = sig2;
      epsy = ey;
      
      
      sigf(0) = fx;
      sigf(1) = fxy;
      
      if (e1 >= e1max) {
	if (e1 <= fcr/Ec) {
	  Dr(0,0) = Ec;					
	} else {
	  Dr(0,0) = (-5*sqrt(5.0)*fcr)/(pow(1 + 10*sqrt(5.0)*sqrt(e1),2)*sqrt(e1));
	}
      } else {
	if(e1 <= fcr/Ec) {
	  Dr(0,0) = Ec;
	} else {
	  Dr(0,0) = fe1max/e1max;
	}	
      }
      Dr(0,1)= 0.0;
      Dr(1,0)= 0.0;
      Dr(1,1)= Ec/2;
      
      
      
      
      
    } else { 
      ////////////opserr << "C2-C" <<endln;
      loadpath = 4.4;
      e1 = fabs(exy);
      e2 = -fabs(exy);
      ey = 0.0;
      
      this->StressEnvelope(e1, e2, e1P, e2P, e1max, e2min);
      
      
      fx = 0.000;
      fy = 0.000;
      
      if (exy>0) {
	fxy = (sig1-sig2)/2;
      } else {
	fxy = -(sig1-sig2)/2;
      }
      
      FinalAnglex = 0.78;
      if(e1 >= fcr/Ec) { 
	crackLabel = 1;
      } else {
	crackLabel = 0;
      }
      Strain1 = e1;
      Strain2 = e2;
      Sigma1 = sig1;
      Sigma2 = sig2;
      epsy = ey;
      
      
      sigf(0) = fx;
      sigf(1) = fxy;
      
      Dr(0,0) = Ec;
      Dr(0,1) = 0.0;
      Dr(1,0) = 0.0;
      Dr(1,1) = (sig1-sig2)/2/exy;
    }
    
  }
  
  
  //Update max-min strains
  if (ex>0 && ex>exmax)
    exmax = ex;
  if (ex<0 && ex<exmin)
    exmin = ex;
  if (ey>0 && ey>eymax)
    eymax = ey;
  if (ey<0 && ey<eymin)
    eymin = ey;
  if (exy>0 && exy>exymax)
    exymax = exy;
  if (exy<0 && exy<exymin)
    exymin = exy;
  if (fabs(exymin)>exymax)
    exymax= fabs(exymin);
  
  
  exminc  = exmin;
  exmaxc  = exmax;
  eyminc  = eymin;
  eymaxc  = eymax;
  exyminc = exymin;
  exymaxc = exymax;
  
}


void 
ConcreteMcftNonLinear7 ::ForwardAngleSearch(void)
{
  //opserr << "ForwardAngleSearch" << endln;
  InitCrackAngle = 0.000001;
  int counter1 = 0;
  double pi = 3.141592654;
  double degreetorad = 3.141592654*2/360.0;
  nE  = Ec/(Ec-fcu/ecu);
  
  //angle search algorithm boundaries
  double bound1 = InitCrackAngle*degreetorad;
  double bound2 = (90-InitCrackAngle)*degreetorad;
  double ResiStress = 1;
  double ResiStressCheck = 1000.0 ;
  double toleranceEpsy = 0.000001;
  double stepsize;
  double nstep = 90;
  
  //stepsize and increments in angle
  stepsize = (bound2-bound1)/nstep;
  theta = bound1+ counter1*stepsize;
  
  //final values
  double e1f, e2f, eyf, fxf, fyf, fxyf, thetaf;
  int countStep = 1; 
  
  
  while((fabs(ResiStress) > toleranceEpsy) ) {
    // strain compatibility: current step for trial theta
    if ( exy > 0) {
      //CASE 1 formulation 
      e2 = ex - exy * tan(theta) / 2;				
    } else if (exy < 0 ) {
      //CASE 2 formulation
      e2 = ex + exy * tan(theta) / 2;
    }
    if (e2<0) {
      double tan2 = tan(theta)*tan(theta);
      e1 = (ex - e2 + ex*tan2)/ tan2;
      ey = e1 + e2 - ex;
      //previous step strains:
      if (exP < 0) {
	e2P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*theta) + exyP/2*sin(2*theta);
	e1P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*theta) - exyP/2*sin(2*theta);
      } else {
	e1P = (exP + eyP)/2 + (exP - eyP)/2*cos(2*theta) + exyP/2*sin(2*theta);
	e2P = (exP + eyP)/2 - (exP - eyP)/2*cos(2*theta) - exyP/2*sin(2*theta);
      }
      if (e1P <0)
	e1P =0.0;
      if (e2P>0)
	e2P = 0.0;
      
      //maximum and minimum strain envelopes on the same proncipal direction
      e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*theta) + exymax/2*sin(2*theta);
      e2min = (exmin + eymin)/2 - (exmin - eymin)/2 * cos(2*theta) + exymax/2*sin(2*theta);
      
      this->StressEnvelope(e1, e2, e1P, e2P, e1max, e2min);
      //Final Stress in Concrete Fiber(in terms of principal stresses) 
      //fy will be checked if its in equilbrm.
      
      if (exy<0) {
	fxy = -1*(sig1 - sig2)/2 * sin(2*theta);
	fx = sig2 - fxy*tan(theta);
	fy = sig1 + fxy*tan(theta);		
      } else if (exy>0) {
	fxy = (sig1 - sig2)/2 * sin(2*theta);
	fx = sig2 + fxy*tan(theta);
	fy = sig1 - fxy*tan(theta);
      }
      ResiStress = Esv * RoV * ey + fy; 
      
      //CheckPoint 1
      if (countStep  > 2 && ResiStress * ResiStressCheck < 0 ) {
	bound1 = theta - stepsize;
	bound2 = theta + stepsize;
	nstep = 10;
	stepsize = (bound2-bound1)/nstep;
	counter1 = 0;
      }
      
      //CheckPoint 2
      if ( ResiStress < toleranceEpsy ) {
	FinalAnglex = theta;
	Strain1 = e1;
	Strain2 = e2;
	Sigma1 = sig1;
	Sigma2 = sig2;
	epsy = ey;
      }
      
      //CheckPoint 3
      if( (countStep > 2) && (ResiStressCheck < 0 && ResiStress < 0) && (ResiStress < ResiStressCheck)) {
	e1 = e1f;
	e2 = e2f;
	ey = eyf;
	fx = fxf;
	fy = fyf;
	fxy= fxyf;
	theta = thetaf;
	
	FinalAnglex = theta;
	if(e1 >= fcr/Ec) { 
	  crackLabel = 1;
	} else {
	  crackLabel = 0;
	}
	Strain1 = e1;
	Strain2 = e2;
	Sigma1 = sig1;
	Sigma2 = sig2;
	epsy = ey;
	sigf(0) = fx;
	sigf(1) = fxy;
	break;
      }
      
      //CheckPoint 4
      if (countStep == 90) {
	break;
      }
      counter1++;		//count theta stepsize
      
      if(fabs(ResiStress) > toleranceEpsy)
	theta = bound1 + counter1*stepsize;
      
      ResiStressCheck = ResiStress;	
      countStep++;	//count computation steps
      
			// record previous step parameters
      e1f = e1;
      e2f = e2;
      eyf = ey;
      fxf = fx;
      fyf = fy;
      fxyf= fxy;
      thetaf = theta;
      
      
    } else if (e2>0) {
      
      counter1++;		//count theta stepsize
      theta = bound1 + counter1*stepsize;
      
      
      ResiStressCheck = ResiStress;	
      countStep++;	
      
    }
  }
}



void 
ConcreteMcftNonLinear7 ::StressEnvelope(double e1, double e2, double e1P, double e2P, double e1max, double e2min)
{
  // stress-strain and sensitivity 
  if(e1max>0 ) {
    
    if (e1max <= fcr/Ec) {
      fe1max = Ec*e1max;
    } else {
      fe1max = fcr / (1+sqrt(500.0*(e1max)));
    }
  } else { 
    fe1max=0.0;
  }
  
  if(e2min<0 ) {
    fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));
  } else {
    fe2min = 0.0;
  }
  //stress-strain
  //principal tension
  if(e1>0) {
    if (e1 >= e1max) {
      if (e1 <= fcr/Ec) {
	sig1 = Ec*e1;				
	loadpath = 4.11;
	
      } else {
	sig1 = fcr / (1+sqrt(500.0*(e1)));				
	loadpath = 4.12;
      }
    } else {
      sig1 = fe1max + fe1max/e1max*(e1-e1max);
      loadpath = 4.14;
    }
  } else {
    sig1 = Ec*e1;		
    loadpath = 4.15;
  }
  //principal compression
  if (e2 <= e2min) {
    sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));		
    loadpath = 4.16;
  } else {
    sig2 = fe2min + fe2min/e2min*(e2-e2min);
    loadpath = 4.17;
  }
  
}

int
ConcreteMcftNonLinear7 ::commitState (void)
{
  
  exP   = epsf(0);
  exyP   = epsf(1);
  eyP   = ey;
  fxP    = sigf(0);
  fxyP   = sigf(1);
  exminP = exminc;
  exmaxP = exmaxc;
  eyminP = eyminc;
  eymaxP = eymaxc;	
  exyminP=exyminc;
  exymaxP =exymaxc;
  loadpathP = loadpath;
  DrP(0,0) = Dr(0,0);
  DrP(0,1) = Dr(0,1);
  DrP(1,0) = Dr(1,0);
  DrP(1,1) = Dr(1,1);
  sigfsensP(0) = sigfsens(0);
  sigfsensP(1) = sigfsens(1);
  //SHVsP = SHVs;
  
  return 0;
}

int
ConcreteMcftNonLinear7 ::revertToLastCommit (void)
{
  
  epsf(0) = exP ; 
  epsf(1) = exyP;
  ey  = eyP;
  sigf(0) = fxyP;
  sigf(1)	= fxyP   ;
  exmin  = exminP;
  exmax  = exmaxP;
  eymin  = eyminP;
  eymax  = eymaxP;	
  exymin = exyminP;
  exymax = exymaxP;
  loadpath = loadpathP;
  Dr(0,0) = DrP(0,0);
  Dr(0,1) = DrP(0,1);
  Dr(1,0) = DrP(1,0);
  Dr(1,1) = DrP(1,1);
  sigfsens(0) = sigfsensP(0);
  sigfsens(1) = sigfsensP(1);
  //SHVs = SHVsP;
  
  return 0;
}

int
ConcreteMcftNonLinear7 ::revertToStart (void)
{
  epsf.Zero();
  sigf.Zero();
  exP   = 0.0;
  exyP   = 0.0;
  eyP   = 0.0;
  fxP    = 0.0;
  fxyP   = 0.0;
  
  exminP = 0.0;
  exmaxP = 0.0;
  eyminP = 0.0;
  eymaxP = 00.0;	
  exyminP=0.0;
  exymaxP =0.0;
  loadpathP = 0.0;
  DrP(0,0) = Ec;
  DrP(0,1) = 0.0;
  DrP(1,0) = 0.0;
  DrP(1,1) = Ec/2;
  
  if (SHVs !=0) {SHVs->Zero();}
  parameterID=0;
  return 0;
}


const char*
ConcreteMcftNonLinear7 ::getType (void) const
{
  //////////opserr << " check16 " << endln;
  return "BeamFiber2d";
}

int 
ConcreteMcftNonLinear7 ::getOrder(void) const
{
  //////////opserr << " check17 " << endln;
  return 2;
}

void
ConcreteMcftNonLinear7 ::Print (OPS_Stream &s, int flag)
{
  //////////opserr << " check18 " << endln;
  return;
}

int
ConcreteMcftNonLinear7 ::sendSelf(int commitTag, Channel &theChannel)
{
  //////////opserr << " check19 " << endln;
  return 0;
}

int
ConcreteMcftNonLinear7 ::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  //////////opserr << " check20 " << endln;
  return 0;
}


Response *
ConcreteMcftNonLinear7 ::setResponse (const char **argv, int argc, 
				      OPS_Stream &s)
{

  Response *theRes = NDMaterial::setResponse(argv, argc, s);
  if (theRes != 0)
    return theRes;
  else {
    if (strcmp(argv[0], "crackAngle") == 0) {
      return new MaterialResponse(this, 10, Vector(5));
    }
    else if (strcmp(argv[0], "fiberStress") == 0) {
      return new MaterialResponse(this, 11, Vector(8));
    }
    else
      return 0;
  }
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
ConcreteMcftNonLinear7 ::setParameter(const char **argv, int argc, Parameter &param)
{
  if (strcmp(argv[0],"fcu") == 0) {// Compressive strength
    //////////opserr << " check22 " << endln;
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"fcr") == 0) {// Vertical Reinforcement Ratio
    //////////opserr << " check22 " << endln;
    return param.addObject(2, this);
  } 
  else if (strcmp(argv[0],"Ec") == 0) {// Vertical Reinforcement Ratio
    //////////opserr << " check22 " << endln;
    return param.addObject(3, this);
  }  
  
  return -1;
}                         

int
ConcreteMcftNonLinear7 ::updateParameter(int parameterID, Information &info)
{
  switch (parameterID) {
  case 1:
    this->fcu = info.theDouble;
    break;
  case 2:
    this->fcr = info.theDouble;
    break;
  case 3:
    this->Ec = info.theDouble;
    break;
  default:
    break;
  }
  return 0;
}

int
ConcreteMcftNonLinear7 ::activateParameter(int passedParameterID)
{  
  parameterID = passedParameterID;
	return 0;
}

const Vector&
ConcreteMcftNonLinear7 ::getStressSensitivity(int gradNumber, bool conditional)
{
	//  sensitivity parameter sensitivity initiate
    	double fcusens = 0.0;
		double fcrsens = 0.0;
		double Ecsens  = 0.0;
	//initiate the gradients
		double fxsens =0.0; // stressSensitivity
		double fxysens=0.0;
		double fysens =0.0;
		double exsens =0.0; // strainSensitivity
		double exysens=0.0;
		double eysens=0.0;

		double nEsens = 0.0; //material parameter nE
		double sig1sens   = 0.0;
		double sig2sens   = 0.0;
		double e1sens     = 0.0;
		double e2sens     = 0.0;
		double e1maxsens  = 0.0;
		double e2minsens  = 0.0;
		double fe1maxsens = 0.0;
		double fe2minsens = 0.0;

	// committed Sensitivity params
		double fxsensC    =0.0; 
		double fxysensC   =0.0;
		double fysensC    =0.0;
		double exsensC    =0.0; 
		double exysensC   =0.0;
		double eysensC    =0.0;
		double exminsensC = 0.0;
		double exmaxsensC = 0.0;
		double eyminsensC = 0.0;
		double eymaxsensC = 0.0;
		double exyminsensC = 0.0;
		double exymaxsensC = 0.0;
   
  

if(SHVs != 0) {
 fxsensC     = (*SHVs)(0,  gradNumber);
 fxysensC    = (*SHVs)(1,  gradNumber);
 fysensC     = (*SHVs)(2,  gradNumber);
 exsensC     = (*SHVs)(3,  gradNumber);
 exysensC    = (*SHVs)(4,  gradNumber);
 eysensC     = (*SHVs)(5,  gradNumber);
 exminsensC  = (*SHVs)(6,  gradNumber);
 exmaxsensC  = (*SHVs)(7,  gradNumber);
 eyminsensC  = (*SHVs)(8,  gradNumber);
 eymaxsensC  = (*SHVs)(9,  gradNumber);
 exyminsensC = (*SHVs)(10, gradNumber);
 exymaxsensC = (*SHVs)(11, gradNumber);
}

//check sensitivity parameter
	if (  parameterID == 1 ) {
		  fcusens = 1.0;
	} else if ( parameterID == 2 ) {
		  fcrsens = 1.0;
    } else if ( parameterID == 3 ) {
	      Ecsens  = 1.0;
    } else {
	  sigfsens(0) = 0.0;
	  sigfsens(1) = 0.0;
	  return sigfsens;
    }
  
	//assign current variables	
	//assign committed strain envelope
	 ex  = epsf(0);
	 exy = epsf(1);
	 exmin = exminP;
	 exmax = exmaxP;
	 eymin = eyminP;
	 eymax = eymaxP;
	 exymin = exyminP;
	 exymax = exymaxP;

	nE  = Ec/(Ec-fcu/ecu);

		if (ex ==0 && exy == 0){
			sigfsens.Zero();
			return sigfsens;
		}


	//gradient required for each ParamID 
	nEsens  = Ecsens/(Ec-fcu/ecu)-Ec*(Ecsens-fcusens/ecu)/(Ec-fcu/ecu)/(Ec-fcu/ecu);


	// Forward Angle Search to get sens of principal axis stresses
	if ( (exy != 0) && (fabs(exy/ex) > 0.01) ) {
	  loadpath = 4.1;

	 
	//Update the principle axis Angle LOOP
	//this->ForwardAngleSearch();

	InitCrackAngle = 0.000001;
	int counter1 = 0;
    double pi = 3.141592654;
	double degreetorad = 3.141592654*2/360.0;
	nE  = Ec/(Ec-fcu/ecu);

	//angle search algorithm boundaries
	double bound1 = InitCrackAngle*degreetorad;
	double bound2 = (90-InitCrackAngle)*degreetorad;
	double ResiStress = 1;
	double ResiStressCheck = 1000.0 ;
	double toleranceEpsy = 0.000001;
	double stepsize;
	double nstep = 90;

	//stepsize and increments in angle
	stepsize = (bound2-bound1)/nstep;
	theta = bound1+ counter1*stepsize;

	//final values
	double e1f, e2f, eyf, fxf, fyf, fxyf, thetaf;

	int countStep = 1; 
	

	while((fabs(ResiStress) > toleranceEpsy) ) {

	    // derive strain compatibility: current step for trial theta
			if ( exy > 0) {
				//CASE 1 formulation 
				e2 = ex - exy * tan(theta) / 2;		
				e2sens = exsens - exysens*tan(theta)/2;
				//////////opserr << "e2sens (1583)= " <<  e2sens << endln;
			} else if (exy < 0 ) {
				//CASE 2 formulation
				e2 = ex + exy * tan(theta) / 2;
				e2sens = exsens+exysens*tan(theta)/2;
				//////////opserr << "e2sens (1589)= " <<  e2sens << endln;
			}

		if (e2<0) {
			double tan2 = tan(theta)*tan(theta);
			e1 = (ex - e2 + ex*tan2)/ tan2;
			e1sens = (exsens - e2sens + exsens*tan2)/ tan2;
			//////////opserr << "e1sens  (1597)= " <<  e1sens << endln;
			ey = e1 + e2 - ex;
			eysens = e1sens + e2sens - exsens;
			//////////opserr << "eysens (1600)= " <<  eysens << endln;
	
			//derive maximum and minimum strain envelopes on the same proncipal direction
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*theta) + exymax/2*sin(2*theta);
			e2min = (exmin + eymin)/2 - (exmin - eymin)/2 * cos(2*theta) + exymax/2*sin(2*theta);
	
			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*theta) + exymaxsensC/2*sin(2*theta);
			e2minsens = (exminsensC + eyminsensC)/2 - (exminsensC - eyminsensC)/2 * cos(2*theta) + exymaxsensC/2*sin(2*theta);
	
		if(e1max>0 ) {
		
		    if (e1max <= fcr/Ec) {
				fe1max = Ec*e1max;
				fe1maxsens = Ecsens*e1max+Ec*e1maxsens;
			} else {
			fe1max = fcr / (1+sqrt(500.0*(e1max)));
			fe1maxsens = fcrsens / (1+sqrt(500.0*e1max))- fcr/pow((1+sqrt(500.0*e1max)),2)*(1/2*pow(500*e1max,-0.5))*(500*e1maxsens);
			}
		} else { 
			fe1max=0.0;
			fe1maxsens=0.0;
		}

	if(e2min<0 ) {
			fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));
			fe2minsens = (e2minsens/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE))+(e2min/ecu)*fcusens*nE/(nE-1.0+pow(e2min/ecu,nE))+
					(e2min/ecu)*fcu*nEsens/(nE-1.0+pow(e2min/ecu,nE))+ 
					-(e2min/ecu)*fcu*nE/pow((nE-1.0+pow(e2min/ecu,nE)),2)*(nEsens+
					pow(e2min/ecu,nE)*(nEsens* log(e2min/ecu)/log(2.718281828459)+nE*(e2minsens/e2min)));

	} else {
			fe2min = 0.0;	
			fe2minsens = 0.0;
	}

//stress-strain
	//principal tension
	if(e1>0) {
		if (e1 >= e1max) {
			if (e1 <= fcr/Ec) {
				sig1 = Ec*e1;
				sig1sens = Ecsens*e1+Ec*e1sens;
			} else {
			sig1 = fcr / (1+sqrt(500.0*(e1)));		
			sig1sens = fcrsens / (1+sqrt(500.0*(e1)))- fcr/pow(1+sqrt(500.0*e1),2)*(1/2*pow(500*e1,-0.5))*(500*e1sens);
			}				
		} else {
			sig1 = fe1max + fe1max/e1max*(e1-e1max);
				sig1sens = fe1maxsens*e1/e1max + fe1max*(e1sens/e1max-e1*e1maxsens/e1max/e1max);
		}
	} else {
		sig1 = Ec*e1;
		sig1sens = Ecsens*e1+Ec*e1sens;
	}

	//principal compression
	if (e2 <= e2min) {
		sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));
		sig2sens = (e2sens/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE))+(e2/ecu)*fcusens*nE/(nE-1.0+pow(e2/ecu,nE))+
					(e2/ecu)*fcu*nEsens/(nE-1.0+pow(e2/ecu,nE))+ 
					-(e2/ecu)*fcu*nE/pow((nE-1.0+pow(e2/ecu,nE)),2)*(nEsens+
					pow(e2/ecu,nE)*(nEsens* log(e2/ecu)/log(2.718281828459)+nE*(e2sens/e2)));
	} else {
		sig2 = fe2min + fe2min/e2min*(e2-e2min);
		sig2sens = fe2minsens*e2/e2min + fe2min*(e2sens/e2min-e2*e2minsens/e2min/e2min); 
   	}

		//Derive Final Stress in Concrete Fiber(in terms of principal stresses) 
		//fy will be checked if its in equilbrm.
		
			if (exy<0) {
				fxy = -1*(sig1 - sig2)/2 * sin(2*theta);
				fx = sig2 - fxy*tan(theta);
				fy = sig1 + fxy*tan(theta);		

				fxysens= -1*(sig1sens - sig2sens)/2*sin(2*theta);
				fxsens = sig2sens - fxysens*tan(theta);
				fysens = sig1sens + fxysens*tan(theta); 
				
			} else if (exy>0) {
				fxy = (sig1 - sig2)/2 * sin(2*theta);
				fx = sig2 + fxy*tan(theta);
				fy = sig1 - fxy*tan(theta);
		
				fxysens= (sig1sens - sig2sens)/2*sin(2*theta);
				fxsens = sig2sens + fxysens*tan(theta);
				fysens = sig1sens - fxysens*tan(theta); 
				
			}
				ResiStress = Esv * RoV * ey + fy; 

				
//CheckPoint 1

			if (countStep  > 2 && ResiStress * ResiStressCheck < 0 ) {
					bound1 = theta - stepsize;
					bound2 = theta + stepsize;
					nstep = 10;
					stepsize = (bound2-bound1)/nstep;
					counter1 = 0;
			}

//CheckPoint 2

			if ( ResiStress < toleranceEpsy ) {
					FinalAnglex = theta;
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;
			}

//CheckPoint 3

			if( (countStep > 2) && (ResiStressCheck < 0 && ResiStress < 0) && (ResiStress < ResiStressCheck)) {
				
					e1 = e1f;
					e2 = e2f;
					ey = eyf;
					fx = fxf;
					fy = fyf;
					fxy= fxyf;
					theta = thetaf;
					
					FinalAnglex = theta;
					if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;

					sigfsens(0) = fxsens;
					sigfsens(1) = fxysens;

			
			
					break;
 
			}

//CheckPoint 4
			if (countStep == 90) {
				break;
			}

			counter1++;		//count theta stepsize
			
			if(fabs(ResiStress) > toleranceEpsy)
			theta = bound1 + counter1*stepsize;
		
			ResiStressCheck = ResiStress;	
			countStep++;	//count computation steps
			
			// record previous step parameters
			e1f = e1;
			e2f = e2;
			eyf = ey;
			fxf = fx;
			fyf = fy;
			fxyf= fxy;
			thetaf = theta;
	

} else if (e2>0) {

			counter1++;		//count theta stepsize
			theta = bound1 + counter1*stepsize;


			ResiStressCheck = ResiStress;	
			countStep++;	
			
		}
}
			sigfsens(0) = fxsens;
			sigfsens(1) = fxysens;
	   
  } else {

	  if (ex<0) {
//derive uniaxial compression
				
	loadpath = 4.2;
				e2 = ex;
				ey = 0.0;
				e1 = 0.0;
				fy = 0.0;
				FinalAnglex = 0.001;
	
				e2sens = exsens;
				exysens = 0.0;
				eysens = 0.0;
				e1sens = 0.0;
				fysens = 0.0;
		
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
			e2min = (exmin - eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);
			
			e2min = exmin;

					
			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*0.0) + exymaxsensC/2*sin(2*0.0);
			e2minsens = (exminsensC - eyminsensC)/2 + (exminsensC - eyminsensC)/2 * cos(2*0.0) - exymaxsensC/2*sin(2*0.0);
			
			e2minsens = exminsensC;

		
			if(e2min<0 ) {
				
				fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));
			fe2minsens = (e2minsens/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE))+(e2min/ecu)*fcusens*nE/(nE-1.0+pow(e2min/ecu,nE))+
						(e2min/ecu)*fcu*nEsens/(nE-1.0+pow(e2min/ecu,nE))+ 
						-(e2min/ecu)*fcu*nE/pow((nE-1.0+pow(e2min/ecu,nE)),2)*(nEsens+
						pow(e2min/ecu,nE)*(nEsens* log(e2min/ecu)/log(2.718281828459)+nE*(e2minsens/e2min)));
					//////////opserr << "fe2minsens (1866)= " <<  fe2minsens << endln;
			} else {
				fe2min = 0.0;
				fe2minsens = 0.0;
					//////////opserr << "fe2minsens (1871)= " <<  fe2minsens << endln;
			}

	
			if (e2 <= e2min) {
				sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));
				sig2sens=	(e2sens/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE))+(e2/ecu)*fcusens*nE/(nE-1.0+pow(e2/ecu,nE))+
							(e2/ecu)*fcu*nEsens/(nE-1.0+pow(e2/ecu,nE))+ 
							-(e2/ecu)*fcu*nE/pow((nE-1.0+pow(e2/ecu,nE)),2)*(nEsens+
							pow(e2/ecu,nE)*(nEsens* log(e2/ecu)/log(2.718281828459)+nE*(e2sens/e2)));
					
						//////////opserr << "sig2sens  (1878)= " <<  sig2sens << endln;
								
				fxysens  =1/2*Ecsens*exy + 1/2*Ec*exysens;
						//////////opserr << "fxysens (1883)= " <<  fxysens << endln;

				
			} else {
				sig2 = fe2min + fe2min/e2min*(e2-e2min);
				sig2sens = fe2minsens*e2/e2min + fe2min*(e2sens/e2min-e2*e2minsens/e2min/e2min); 
				//////////opserr << "sig2sens (1888)= " <<  sig2sens << endln;	
				fxysens  =1/2*Ecsens*exy + 1/2*Ec*exysens;
				//////////opserr << "fxsens (1892)= " <<  fxsens << endln;
			}
				
				    fxsens = sig2sens;
					fxysens = 1/2*Ecsens*exy + 1/2*Ec*exysens;
			
				sigfsens(0) = fxsens;
				sigfsens(1) = fxysens;
	} else if (ex>0) {

//derive uniaxial tension
		  //opserr << "UNIAX positive ex"<< endln;
		  
				loadpath = 4.3;
				e1 = ex;
				e2 = 0.0;
				ey = 0.0;
				fy = 0.0;
			
				e1sens = exsens;
				e2sens = 0.0;
				eysens = 0.0;
				fysens = 0.0;
				exysens = 0.0;
				FinalAnglex = 89.999;
//derive maximum and minimum strain envelopes on the same pronipal direction
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
			e2min = (exmin + eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);

			e1max = exmax;
			

			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*0.0) + exymaxsensC/2*sin(2*0.0);
			e2minsens = (exminsensC + eyminsensC)/2 + (exminsensC - eyminsensC)/2 * cos(2*0.0) - exymaxsensC/2*sin(2*0.0);
		
			e1maxsens = exmaxsensC ;
	

			if (e1max >0 && e1max <= fcr/Ec) {
				fe1max = Ec*e1max;
				fe1maxsens = Ecsens*e1max+Ec*e1maxsens;
				//////////opserr << "fe1maxsens (1932)= " <<  fe1maxsens << endln;
				
				
			}  else if (e1max <=0 ) {
				fe1max= 0.0;
				fe1maxsens  = 0.0;
					//////////opserr << "fe1maxsens (1936)= " <<  fe1maxsens << endln;
				
			} else {
				fe1max = fcr / (1+sqrt(500.0*(e1max)));
				fe1maxsens = fcrsens/(1+sqrt(500.0*e1max))+
							-fcr/(1+sqrt(500.0*e1max))/(1+sqrt(500.0*e1max))*(0.5/pow(500*e1max,0.5)*500*e1maxsens);
		
			}
	

			//stress-strain
			//principal tension
			if (e1 >= e1max) {
						if (e1 <= fcr/Ec) {
								sig1 = Ec*e1;
								sig1sens = Ecsens*e1+Ec*e1sens;
						//////////opserr << "fxsens (1951)= " <<  fxsens << endln;

						} else {
								sig1 = fcr / (1+sqrt(500.0*(e1)));
								sig1sens = fcrsens / (1+sqrt(500.0*(e1)))- fcr/pow(1+sqrt(500.0*e1),2)*(1/2*pow(500*e1,-0.5))*(500*e1sens);	
						//////////opserr << "fxsens (1956)= " <<  fxsens << endln;
						}
			} else {
						sig1 = fe1max + fe1max/e1max*(e1-e1max);
						sig1sens = fe1maxsens*e1/e1max + fe1max*(e1sens/e1max-e1*e1maxsens/e1max/e1max);
						//////////opserr << "fxsens (1967)= " <<  fxsens << endln;
			}				
				
				fxsens = sig1sens;		
				fxysens  =1/2*Ecsens*exy + Ec/2*exysens;	
				sigfsens(0) = fxsens;
				sigfsens(1) = fxysens;
	
			} else {
	//opserr << "C2-C" <<endln;
	loadpath = 4.4;
				e1 = 0.0;
				e2 = 0.0;
				ey = 0.0;
				fx = 0.000;
				fy = 0.000;
				fxy =0.000;
				sig1=fx;
				sig2=fy;
				
	
				e1sens = 0.0;
				e2sens = 0.0;
				eysens = 0.0;
				fxsens = 0.000;
				fysens = 0.000;
				fxysens =0.000;
				sig1sens=fxsens;
				sig2sens=fysens;
				FinalAnglex = 0.001;
				
				if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
			
				sigfsens(0) = fxsens;
				sigfsens(1) = fxysens;
				//////opserr << " line 1941" << sigfsens << endln; 

	}
			
	}


//Update max-min strains
	if (ex>0 && ex>exmax){
	//	exmax = ex;
		exmaxsensC = exsens;
	}
	if (ex<0 && ex<exmin){
	//	exmin = ex;
		exminsensC = exsens;
	}
	if (ey>0 && ey>eymax){
	//	eymax = ey;
		eymaxsensC = eysens;
	}
	if (ey<0 && ey<eymin) {
	//	eymin = ey;
		eyminsensC = eysens;
	}
	if (exy>0 && exy>exymax) {
	//	exymax = exy;
		exymaxsensC = exysens;
	}
	if (exy<0 && exy<exymin) {
	//	exymin = exy;
		exyminsensC = exysens;
	}
	if (fabs(exymin)>exymax){
	//	exymax= fabs(exy);
		exymaxsensC= fabs(exysens);
	}



	return sigfsens;

}



const Matrix&
ConcreteMcftNonLinear7 ::getInitialTangentSensitivity(int gradNumber)
{
  static Matrix dDridh(2,2);
  
  dDridh.Zero();
  
  double sens00 = 0.0;
  
  if (  parameterID == 1 ) {
    sens00 = 1.0;
  }
  else if ( parameterID == 2 ) {
    sens00 = 1.0;
  }
  else if ( parameterID == 3 ) {
    sens00 = 1.0;
  } else {
    sens00 =0.0;
  }
  dDridh(0,0) = sens00;
  dDridh(0,1) = 0.0;
  dDridh(1,0) = 0.0;
  dDridh(1,1) = sens00/2;
  
  return dDridh;
}

int
ConcreteMcftNonLinear7 ::commitSensitivity(const Vector &strainGradient, int gradNumber, int numGrads)
{

	//  sensitivity parameter sensitivity initiate
		double fcusens = 0.0;
		double fcrsens = 0.0;
		double Ecsens = 0.0;
    
	//initiate the gradients
		double fxsens =0.0; // stressSensitivity
		double fxysens=0.0;
		double fysens =0.0;

		double nEsens = 0.0; //material parameter nE
	//initiate  principal gradients
		double sig1sens   = 0.0;
		double sig2sens	  = 0.0;
		double e1sens     = 0.0;
		double e2sens     = 0.0;
		double e1maxsens  = 0.0;
		double e2minsens  = 0.0;
		double fe1maxsens = 0.0;
		double fe2minsens = 0.0;

	// initiate committed Sensitivity
		double fxsensC     = 0.0; 
		double fxysensC    = 0.0;
		double fysensC     = 0.0;
		double exsensC     = 0.0; 
		double exysensC    = 0.0;
		double eysensC     = 0.0;
		double exminsensC  = 0.0;
		double exmaxsensC  = 0.0;
		double eyminsensC  = 0.0;
		double eymaxsensC  = 0.0;
		double exyminsensC = 0.0;
		double exymaxsensC = 0.0;
  
//call committed stress-strain sensitivities
	 if(SHVs == 0) {
		SHVs = new Matrix(12,numGrads);
	 }else {

		fxsensC     = (*SHVs)(0, gradNumber);
		fxysensC    = (*SHVs)(1, gradNumber);
		fysensC     = (*SHVs)(2, gradNumber);
		exsensC     = (*SHVs)(3, gradNumber);
		exysensC    = (*SHVs)(4, gradNumber);
		eysensC     = (*SHVs)(5, gradNumber);
		exminsensC  = (*SHVs)(6, gradNumber);
		exmaxsensC  = (*SHVs)(7, gradNumber);
		eyminsensC  = (*SHVs)(8, gradNumber);
		eymaxsensC  = (*SHVs)(9, gradNumber);
		exyminsensC = (*SHVs)(10, gradNumber);
		exymaxsensC = (*SHVs)(11, gradNumber);
	 }




	if (  parameterID == 1 ) {
			fcusens = 1.0;
	 } else if ( parameterID == 2 ) {
	        fcrsens = 1.0;
	 } else if ( parameterID == 3 ) {
			Ecsens = 1.0;
     } else {
	        return 0;
	}

	 // strainSensitivity
		double exsens  = strainGradient(0); 
		double exysens = strainGradient(1);
		double eysens  =0.0;
	//Start sensitivity calculations
	 ex  = epsf(0);
	 exy = epsf(1);
	 exmin = exminP;
	 exmax = exmaxP;
	 eymin = eyminP;
	 eymax = eymaxP;
	 exymin = exyminP;
	 exymax = exymaxP;

	
	nE  = Ec/(Ec-fcu/ecu);
	
	//gradient required for each ParamID 
	nEsens  = Ecsens/(Ec-fcu/ecu)-Ec*(Ecsens-fcusens/ecu)/(Ec-fcu/ecu)/(Ec-fcu/ecu);

//checkpoint A	
//Stress Cases
	if(ex == 0.0 && exy == 0.0) {
			
			  (*SHVs)(0, gradNumber) = 0.0;
			  (*SHVs)(1, gradNumber) = 0.0; 
			  (*SHVs)(2, gradNumber) = 0.0;
			  (*SHVs)(3, gradNumber) = 0.0;     
			  (*SHVs)(4, gradNumber) =  0.0;
			  (*SHVs)(5, gradNumber) = 0.0;
			  (*SHVs)(6, gradNumber) = 0.0;
			  (*SHVs)(7, gradNumber) =  0.0;
			  (*SHVs)(8, gradNumber) = 0.0;
			  (*SHVs)(9, gradNumber) = 0.0;
			  (*SHVs)(10, gradNumber) = 0.0;
			  (*SHVs)(11, gradNumber) = 0.0;

	} else if( (ex!=0.0 && exy!=0.0) && (ex==exP && exy==exyP)) {
	
			    (*SHVs)(0, gradNumber) = fxsensC;
		        (*SHVs)(1, gradNumber) = fxysensC; 
			    (*SHVs)(2, gradNumber) = fysensC;
		        (*SHVs)(3, gradNumber) = exsensC;     
                (*SHVs)(4, gradNumber) = exysensC;
				(*SHVs)(5, gradNumber) = eysensC;
		    	(*SHVs)(6, gradNumber) = exminsensC;
			    (*SHVs)(7, gradNumber) = exmaxsensC;
				(*SHVs)(8, gradNumber) = eyminsensC;
				(*SHVs)(9, gradNumber) = eymaxsensC;
				(*SHVs)(10, gradNumber) = exyminsensC;
				(*SHVs)(11, gradNumber) = exymaxsensC;

	}  else if ((ex>0.10) || (ex<-0.10)){

			  (*SHVs)(0, gradNumber) = 0.0;
			  (*SHVs)(1, gradNumber) = 0.0; 
			  (*SHVs)(2, gradNumber) = 0.0;
			  (*SHVs)(3, gradNumber) = 0.0;     
			  (*SHVs)(4, gradNumber) =  0.0;
			  (*SHVs)(5, gradNumber) = 0.0;
			  (*SHVs)(6, gradNumber) = 0.0;
			  (*SHVs)(7, gradNumber) =  0.0;
			  (*SHVs)(8, gradNumber) = 0.0;
			  (*SHVs)(9, gradNumber) = 0.0;
			  (*SHVs)(10, gradNumber) = 0.0;
			  (*SHVs)(11, gradNumber) = 0.0;
		
	
	} else {

	// Forward Angle Search to get sens of principal axis stresses
	if ( (exy != 0) && (fabs(exy/ex) > 0.01) ) {
	  loadpath = 4.1;

	 // ////opserr << " loadpath sens" << loadpath << endln; 
			//Update the principle axis Angle LOOP
		//this->ForwardAngleSearch();

	InitCrackAngle = 0.000001;
	int counter1 = 0;
    double pi = 3.141592654;
	double degreetorad = 3.141592654*2/360.0;
	nE  = Ec/(Ec-fcu/ecu);

	//angle search algorithm boundaries
	double bound1 = InitCrackAngle*degreetorad;
	double bound2 = (90-InitCrackAngle)*degreetorad;
	double ResiStress = 1;
	double ResiStressCheck = 1000.0 ;
	double toleranceEpsy = 0.000001;
	double stepsize;
	double nstep = 90;

	//stepsize and increments in angle
	stepsize = (bound2-bound1)/nstep;
	theta = bound1+ counter1*stepsize;

	//final values
	double e1f, e2f, eyf, fxf, fyf, fxyf, thetaf;

	int countStep = 1; 
	

	while((fabs(ResiStress) > toleranceEpsy) ) {

	    // derive strain compatibility: current step for trial theta
			if ( exy > 0) {
				//CASE 1 formulation 
				e2 = ex - exy * tan(theta) / 2;		
				e2sens = exsens - exysens*tan(theta)/2;
			} else if (exy < 0 ) {
				//CASE 2 formulation
				e2 = ex + exy * tan(theta) / 2;
				e2sens = exsens+exysens*tan(theta)/2;
			}

		if (e2<0) {
			double tan2 = tan(theta)*tan(theta);
			e1 = (ex - e2 + ex*tan2)/ tan2;
			e1sens = (exsens - e2sens + exsens*tan2)/ tan2;
			ey = e1 + e2 - ex;
			eysens = e1sens + e2sens - exsens;
	
			//derive maximum and minimum strain envelopes on the same proncipal direction
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*theta) + exymax/2*sin(2*theta);
			e2min = (exmin + eymin)/2 - (exmin - eymin)/2 * cos(2*theta) + exymax/2*sin(2*theta);
	
			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*theta) + exymaxsensC/2*sin(2*theta);
			e2minsens = (exminsensC + eyminsensC)/2 - (exminsensC - eyminsensC)/2 * cos(2*theta) + exymaxsensC/2*sin(2*theta);
	
		if(e1max>0 ) {
		
		    if (e1max <= fcr/Ec) {
				fe1max = Ec*e1max;
				fe1maxsens = Ecsens*e1max+Ec*e1maxsens;
				
			} else {
				fe1max = fcr / (1+sqrt(500.0*(e1max)));
				fe1maxsens = fcrsens / (1+sqrt(500.0*e1max))- fcr/pow((1+sqrt(500.0*e1max)),2)*(1/2*pow(500*e1max,-0.5))*(500*e1maxsens);
			}
	} else { 
			fe1max=0.0;
			fe1maxsens=0.0;
	}

	if(e2min<0 ) {
				
			fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));

			fe2minsens = (e2minsens/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE))+(e2min/ecu)*fcusens*nE/(nE-1.0+pow(e2min/ecu,nE))+
					(e2min/ecu)*fcu*nEsens/(nE-1.0+pow(e2min/ecu,nE))+ 
					-(e2min/ecu)*fcu*nE/pow((nE-1.0+pow(e2min/ecu,nE)),2)*(nEsens+
					pow(e2min/ecu,nE)*(nEsens* log(e2min/ecu)/log(2.718281828459)+nE*(e2minsens/e2min)));
	} else {
			fe2min = 0.0;	
			fe2minsens = 0.0;
	}

		

//stress-strain
	//principal tension
	if(e1>0) {
		if (e1 >= e1max) {
			if (e1 <= fcr/Ec) {
				sig1 = Ec*e1;
				sig1sens = Ecsens*e1+Ec*e1sens;
			} else {
			sig1 = fcr / (1+sqrt(500.0*(e1)));		
			sig1sens = fcrsens / (1+sqrt(500.0*(e1)))- fcr/pow(1+sqrt(500.0*e1),2)*(1/2*pow(500*e1,-0.5))*(500*e1sens);
			}
				
		} else {
			sig1 = fe1max + fe1max/e1max*(e1-e1max);
				sig1sens = fe1maxsens*e1/e1max + fe1max*(e1sens/e1max-e1*e1maxsens/e1max/e1max);
		}
	} else {
		sig1 = Ec*e1;
		sig1sens = Ecsens*e1+Ec*e1sens;
	}


	//principal compression
	if (e2 <= e2min) {
		sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));
		sig2sens = (e2sens/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE))+(e2/ecu)*fcusens*nE/(nE-1.0+pow(e2/ecu,nE))+
					(e2/ecu)*fcu*nEsens/(nE-1.0+pow(e2/ecu,nE))+ 
					-(e2/ecu)*fcu*nE/pow((nE-1.0+pow(e2/ecu,nE)),2)*(nEsens+
					pow(e2/ecu,nE)*(nEsens* log(e2/ecu)/log(2.718281828459)+nE*(e2sens/e2)));
	
	} else {
		sig2 = fe2min + fe2min/e2min*(e2-e2min);
		sig2sens = fe2minsens*e2/e2min + fe2min*(e2sens/e2min-e2*e2minsens/e2min/e2min); 
   	}

	//Derive Final Stress in Concrete Fiber(in terms of principal stresses) 
	//fy will be checked if its in equilbrm.
		
			if (exy<0) {
				fxy = -1*(sig1 - sig2)/2 * sin(2*theta);
				fx = sig2 - fxy*tan(theta);
				fy = sig1 + fxy*tan(theta);		

				fxysens= -1*(sig1sens - sig2sens)/2*sin(2*theta);
				fxsens = sig2sens - fxysens*tan(theta);
				fysens = sig1sens + fxysens*tan(theta); 

			} else if (exy>0) {
				fxy = (sig1 - sig2)/2 * sin(2*theta);
				fx = sig2 + fxy*tan(theta);
				fy = sig1 - fxy*tan(theta);
		
				fxysens= (sig1sens - sig2sens)/2*sin(2*theta);
				fxsens = sig2sens + fxysens*tan(theta);
				fysens = sig1sens - fxysens*tan(theta); 
			}
				ResiStress = Esv * RoV * ey + fy; 
			
				
//CheckPoint 1

			if (countStep  > 2 && ResiStress * ResiStressCheck < 0 ) {
					bound1 = theta - stepsize;
					bound2 = theta + stepsize;
					nstep = 10;
					stepsize = (bound2-bound1)/nstep;
					counter1 = 0;
			}

//CheckPoint 2

			if ( ResiStress < toleranceEpsy ) {
					FinalAnglex = theta;
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;
			}

//CheckPoint 3

			if( (countStep > 2) && (ResiStressCheck < 0 && ResiStress < 0) && (ResiStress < ResiStressCheck)) {
				
					e1 = e1f;
					e2 = e2f;
					ey = eyf;
					fx = fxf;
					fy = fyf;
					fxy= fxyf;
					theta = thetaf;
					
					FinalAnglex = theta;
					if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;

					break;
 
			}

//CheckPoint 4
			if (countStep == 90) {
				break;
			}

			counter1++;		//count theta stepsize
			
			if(fabs(ResiStress) > toleranceEpsy)
			theta = bound1 + counter1*stepsize;
		
			ResiStressCheck = ResiStress;	
			countStep++;	//count computation steps
			
			// record previous step parameters
			e1f = e1;
			e2f = e2;
			eyf = ey;
			fxf = fx;
			fyf = fy;
			fxyf= fxy;
			thetaf = theta;
	

} else if (e2>0) {

			counter1++;		//count theta stepsize
			theta = bound1 + counter1*stepsize;


			ResiStressCheck = ResiStress;	
			countStep++;	
			
		}
}
	
  } else {

	  if (ex<0) {
//derive uniaxial compression
				
	loadpath = 4.2;
				e2 = ex;
				ey = 0.0;
				e1 = 0.0;
				fy = 0.0;
				FinalAnglex = 0.001;
	
				e2sens = exsens;
				exysens = 0.0;
				eysens = 0.0;
				e1sens = 0.0;
				fysens = 0.0;

		//derive maximum and minimum strain envelopes on the same proncipal direction
	
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
			e2min = (exmin - eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);
			
			e2min = exmin;

					
			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*0.0) + exymaxsensC/2*sin(2*0.0);
			e2minsens = (exminsensC - eyminsensC)/2 + (exminsensC - eyminsensC)/2 * cos(2*0.0) - exymaxsensC/2*sin(2*0.0);
			
			e2minsens = exminsensC;

		if(e2min<0 ) {
				
			fe2min = (e2min/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE));
			fe2minsens = (e2minsens/ecu)*fcu*nE/(nE-1.0+pow(e2min/ecu,nE))+(e2min/ecu)*fcusens*nE/(nE-1.0+pow(e2min/ecu,nE))+
						(e2min/ecu)*fcu*nEsens/(nE-1.0+pow(e2min/ecu,nE))+ 
						-(e2min/ecu)*fcu*nE/pow((nE-1.0+pow(e2min/ecu,nE)),2)*(nEsens+
						pow(e2min/ecu,nE)*(nEsens* log(e2min/ecu)/log(2.718281828459)+nE*(e2minsens/e2min)));

		} else {
				fe2min = 0.0;
				fe2minsens = 0.0;
		}

	
			if (e2 <= e2min) {
				
				sig2 = (e2/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE));
				sig2sens=	(e2sens/ecu)*fcu*nE/(nE-1.0+pow(e2/ecu,nE))+(e2/ecu)*fcusens*nE/(nE-1.0+pow(e2/ecu,nE))+
							(e2/ecu)*fcu*nEsens/(nE-1.0+pow(e2/ecu,nE))+ 
							-(e2/ecu)*fcu*nE/pow((nE-1.0+pow(e2/ecu,nE)),2)*(nEsens+
							pow(e2/ecu,nE)*(nEsens* log(e2/ecu)/log(2.718281828459)+nE*(e2sens/e2)));
					
						//////////opserr << "sig2sens  (1878)= " <<  sig2sens << endln;
								
				fxysens  =1/2*Ecsens*exy + 1/2*Ec*exysens;
						//////////opserr << "fxysens (1883)= " <<  fxysens << endln;

				
			} else {
				sig2 = fe2min + fe2min/e2min*(e2-e2min);
				sig2sens = fe2minsens*e2/e2min + fe2min*(e2sens/e2min-e2*e2minsens/e2min/e2min); 
				fxysens  =1/2*Ecsens*exy + 1/2*Ec*exysens;
			}
			    fxsens = sig2sens;
				fxysens = 1/2*Ecsens*exy + 1/2*Ec*exysens;

	  } else if (ex>0) {
//derive uniaxial tension

	loadpath = 4.3;
				e1 = ex;
				e2 = 0.0;
				ey = 0.0;
				fy = 0.0;
			
				e1sens = exsens;
				e2sens = 0.0;
				eysens = 0.0;
				fysens = 0.0;
				exysens = 0.0;
				FinalAnglex = 89.999;
			
			

		//derive maximum and minimum strain envelopes on the same pronipal direction
			e1max = (exmax + eymax)/2 + (exmax - eymax)/2 * cos(2*0.0) + exymax/2*sin(2*0.0);
			e2min = (exmin + eymin)/2 + (exmin - eymin)/2 * cos(2*0.0) - exymax/2*sin(2*0.0);

			e1max = exmax;
			

			e1maxsens = (exmaxsensC + eymaxsensC)/2 + (exmaxsensC - eymaxsensC)/2 * cos(2*0.0) + exymaxsensC/2*sin(2*0.0);
			e2minsens = (exminsensC + eyminsensC)/2 + (exminsensC - eyminsensC)/2 * cos(2*0.0) - exymaxsensC/2*sin(2*0.0);
		
			e1maxsens = exmaxsensC ;
		
			if (e1max >0 && e1max <= fcr/Ec) {
				fe1max = Ec*e1max;
				fe1maxsens = Ecsens*e1max+Ec*e1maxsens;
				//////////opserr << "fe1maxsens (1932)= " <<  fe1maxsens << endln;
				
				
			}  else if (e1max <=0 ) {
				fe1max= 0.0;
				fe1maxsens  = 0.0;
					//////////opserr << "fe1maxsens (1936)= " <<  fe1maxsens << endln;
				
			} else {
				fe1max = fcr / (1+sqrt(500.0*(e1max)));
				fe1maxsens = fcrsens/(1+sqrt(500.0*e1max))+
							-fcr/(1+sqrt(500.0*e1max))/(1+sqrt(500.0*e1max))*(0.5/pow(500*e1max,0.5)*500*e1maxsens);
		
			}
	

			//stress-strain
			//principal tension
			if (e1 >= e1max) {
						if (e1 <= fcr/Ec) {
								sig1 = Ec*e1;
								sig1sens = Ecsens*e1+Ec*e1sens;
				
						} else {
								sig1 = fcr / (1+sqrt(500.0*(e1)));
								sig1sens = fcrsens / (1+sqrt(500.0*(e1)))- fcr/pow(1+sqrt(500.0*e1),2)*(1/2*pow(500*e1,-0.5))*(500*e1sens);
						}
			} else {
					if(e1 <= fcr/Ec) {
						sig1 = Ec*e1;
						sig1sens = Ecsens*e1 + Ec*e1sens;
					} else {
						sig1 = fe1max + fe1max/e1max*(e1-e1max);
						sig1sens = fe1maxsens*e1/e1max + fe1max*(e1sens/e1max-e1*e1maxsens/e1max/e1max);
					}

			}				
				
				fxsens = sig1sens;		
				fxysens  =1/2*Ecsens*exy + Ec/2*exysens;	
		} else {
/////////opserr << "C2-C" <<endln;
	loadpath = 4.4;
				e1 = 0.0;
				e2 = 0.0;
				ey = 0.0;
				fx = 0.000;
				fy = 0.000;
				fxy =0.000;
				sig1=fx;
				sig2=fy;
				
	
				e1sens = 0.0;
				e2sens = 0.0;
				eysens = 0.0;
				fxsens = 0.000;
				fysens = 0.000;
				fxysens =0.000;
				sig1sens=fxsens;
				sig2sens=fysens;
				FinalAnglex = 0.001;
				
				if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
		}	
}


//Update max-min strains
if (ex>0 && ex>exmax){
		exmax = ex;
		exmaxsensC = exsens;
}
if (ex<0 && ex<exmin){
		exmin = ex;
		exminsensC = exsens;
}
if (ey>0 && ey>eymax){
		eymax = ey;
		eymaxsensC = eysens;
}
if (ey<0 && ey<eymin){
		eymin = ey;
		eyminsensC = eysens;
}
if (exy>0 && exy>exymax){
		exymax = exy;
		exymaxsensC = exysens;
}
if (exy<0 && exy<exymin){
		exymin = exy;
		exyminsensC = exysens;
}
if (fabs(exymin)>exymax){
		exymax= fabs(exy);
		exymaxsensC= fabs(exysens);
}

  (*SHVs)(0, gradNumber)  = fxsensC;
  (*SHVs)(1, gradNumber)  = fxysensC; 
  (*SHVs)(2, gradNumber)  = fysensC;
  (*SHVs)(3, gradNumber)  = exsensC;     
  (*SHVs)(4, gradNumber)  = exysensC;
  (*SHVs)(5, gradNumber)  = eysensC;
  (*SHVs)(6, gradNumber)  = exminsensC;
  (*SHVs)(7, gradNumber)  = exmaxsensC;
  (*SHVs)(8, gradNumber)  = eyminsensC;
  (*SHVs)(9, gradNumber)  = eymaxsensC;
  (*SHVs)(10, gradNumber) = exyminsensC;
  (*SHVs)(11, gradNumber) = exymaxsensC;
}
	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////

int 
ConcreteMcftNonLinear7 ::getResponse (int responseID, Information &matInformation)
{
double epsx = ex;
double epsxy = exy;
double epsy = ey;
double sigx = fx;
double sigxy = fxy;
double angl = FinalAnglex;
double cL = crackLabel;

	static Vector crackInfo(6);
	if (responseID == 10) {
		crackInfo(0) = epsx;
		crackInfo(1) = epsxy;
		crackInfo(2) = sigx;
		crackInfo(3) = sigxy;
		crackInfo(4) = angl;
		crackInfo(5) = epsy;
		matInformation.setVector(crackInfo);
	} 
	static Vector prinStress(8);
	if (responseID == 11) {
		prinStress(0) = e1;
		prinStress(1) = e2;
		prinStress(2) = Sigma1;
		prinStress(3) = Sigma2;
		prinStress(4) = Dr(0,0);
		prinStress(5) = Dr(0,1);
		prinStress(6) = Dr(1,0);
		prinStress(7) = Dr(1,1);
		matInformation.setVector(prinStress);
	}
	return 0;
}

double 
ConcreteMcftNonLinear7 ::tangentstifness00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min)
{
	
  double cott = 1/tan(theta);
  double sect = 1/cos(theta);
  double csct = 1/sin(theta);
  double d00;
  
  if (exy>0) {
    if(e1>e1max || e1<0){
      if(e1<=fcr/Ec) {
	if(e2<=e2min) {
	  //AD1-d00
	  d00 = -((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
		  (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
	    (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	    (sin(2*theta)*tan(theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				      (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
				      (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	    ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
						       (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
						       (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
	     ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
	      (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
	      (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
					2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
					(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
					(2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
					(Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
	      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
					 (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
	      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
				       (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
	    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
	     2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
	     Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	     (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
				       2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
				       (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				       (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				       (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	     (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
					(Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	     cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
				      (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
	} 
	else {
	  //AF1-d00
	  d00 = fe2min/e2min + ((Ec - fe2min/e2min)*sin(2*theta)*tan(theta))/2. - ((Ec + Esv*RoV - ((Ec - fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
										   (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (pow(sect,2)*sin(2*theta)*
															    (-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
										    cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
													     Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
										    (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
													      Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
													      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
	    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
	     2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
	     (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
					Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
	     cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
	     Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	     (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
				       Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
				       2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
	}
      } 
      else {
	if(e2<=e2min) {
	  //BD1-d00
	  d00= -((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
		 (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
	    (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	    (sin(2*theta)*tan(theta)*((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				      (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
				      (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
				      (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
							 pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
	    ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					   pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
	      (sin(2*theta)*tan(theta)*((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
					(pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
					(Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
					(5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
							   pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
	     ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
	      (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
	      (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
					(2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
					(Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
					(5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
							  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
					(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					 pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
	      (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
					   (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
					 fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
	      cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
					 (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
				       fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
	    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	     (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
			       2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
	     (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
	      pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
	     (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				       (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				       (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
				       (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
							 2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
				       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
	     (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
					  (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
					fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
	     cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
					(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
				      fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
	} 
	else {
	  //BF1-d00
	  d00=fe2min/e2min + (sin(2*theta)*tan(theta)*(-(fe2min/e2min) - (5*sqrt(5.0)*fcr)/
						       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
							pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
	    ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					   pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
	      (sin(2*theta)*tan(theta)*(-(fe2min/e2min) - (5*sqrt(5.0)*fcr)/
					(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					 pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
	     (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (sin(2*theta)*tan(theta)*
						      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (5*sqrt(5.0)*fcr*
											      (pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
											       2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
						       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
							pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
	      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
					 fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
	      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
				       fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
	    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	     (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
			       2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
	     (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
	      pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
	     (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
				       (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
							 2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
				       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
					pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
	     (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
					fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
	     cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
				      fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
	}
      }
    } 
    else {		
      if(e2<=e2min) {
	//CD1-d00
	d00 = -((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
		(pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
	  (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	  (sin(2*theta)*tan(theta)*(fe1max/e1max + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				    (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
				    (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	  ((fe1max/e1max + Esv*RoV - (sin(2*theta)*tan(theta)*(fe1max/e1max + 
							       (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
							       (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
							       (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
	   ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
	    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
	    (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
	    (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
				       (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				       (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
	    cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
				     (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				     (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
	    (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				      (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				      (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				      (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
					       2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
	  ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		    2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
	   Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		    2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	   (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
				      (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				      (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
	   cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
				    (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				    (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
	   (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				     (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				     (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
				     (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
					      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
      } 
      else {
	//CF1-d00
	d00 = fe2min/e2min + ((fe1max/e1max - fe2min/e2min)*sin(2*theta)*tan(theta))/2. - 
	  ((fe1max/e1max + Esv*RoV - ((fe1max/e1max - fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
	   (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (pow(sect,2)*sin(2*theta)*
						    (fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
						     (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
	    cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
				     (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
	    (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
				      (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
					       2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
	  ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		    2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
	   Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		    2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
	   (pow(sect,2)*sin(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
				      (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
	   cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
				    (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
	   (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
				     (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
					      2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
      }		
    }
    
  } 
  else {
    if(e1>e1max||e1<0){
      if(e1<=fcr/Ec) {
	if(e2<=e2min) {
	  //AD2-d00
	  d00 = -((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
		  (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
	    (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
	    (sin(2*theta)*tan(theta)*(-Ec - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				      (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				      (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	    ((Ec + Esv*RoV + (sin(2*theta)*tan(theta)*(-Ec - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
						       (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
						       (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
	     (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
	      (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
	      (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
	      (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
					2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
					(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
					(2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
					(Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	      (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
					 (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
	      cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
				       (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
	    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
	     2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
	     Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
		      2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
	     (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
				       2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
				       (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
				       (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
				       (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
	     (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
					(Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
	     cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
				      (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
	} else {
	  //AF2-d00
	  d00 = fe2min/e2min - ((-Ec + fe2min/e2min)*sin(2*theta)*tan(theta))/2. - 
	    ((Ec + Esv*RoV + ((-Ec + fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (pow(sect,2)*sin(2*theta)*
           (fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
        cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD2-d00
				d00 = -((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
      (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
   (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
   (sin(2*theta)*tan(theta)*(-((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
           (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
        (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(-((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
                (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF2-d00
					d00= fe2min/e2min - (sin(2*theta)*tan(theta)*(fe2min/e2min + (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(fe2min/e2min + (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (sin(2*theta)*tan(theta)*
           ((exy*fe2min*pow(sect,2))/(2.*e2min) + (5*sqrt(5.0)*fcr*
                (pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))); 
				}
			}
			} else {
			if(e2<=e2min) {
				//CD2-d00
				d00 = -((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
      (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
   (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
   (sin(2*theta)*tan(theta)*(-(fe1max/e1max) - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((fe1max/e1max + Esv*RoV + (sin(2*theta)*tan(theta)*(-(fe1max/e1max) - 
             (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
        cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF2-d00
				d00 = fe2min/e2min - ((-(fe1max/e1max) + fe2min/e2min)*sin(2*theta)*tan(theta))/2. - 
   ((fe1max/e1max + Esv*RoV + ((-(fe1max/e1max) + fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (pow(sect,2)*sin(2*theta)*
           (-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
             (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
        cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
	}


return d00;
}

double 
ConcreteMcftNonLinear7 ::tangentstifness01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d01;

	if (exy>0) {
		if(e1>e1max || e1<0){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD1-d01
					d01 = (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
   (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
				} else {
					//AF1-d01
					d01 = -(fe2min*tan(theta))/(2.*e2min) + (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (pow(sect,2)*sin(2*theta)*
           (-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
        cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.) ;
				}
			} else {
			if(e2<=e2min) {
					//BD1-d01
				d01 = (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
   (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
                (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF1-d01
					d01 = -(fe2min*tan(theta))/(2.*e2min) + (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (sin(2*theta)*tan(theta)*
           ((exy*fe2min*pow(sect,2))/(2.*e2min) - (5*sqrt(5.0)*fcr*
                (pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				}
			}
			}  else {
			if(e2<=e2min) {
				//CD1-d01
				d01 = (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
   (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   (((fe1max*cott)/(2.*e1max) + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - 
             (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
        (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
        cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF1-d01
				d01 = -(fe2min*tan(theta))/(2.*e2min) + (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   (((fe1max*cott)/(2.*e1max) + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (-(exy*fe2min*pow(sect,2))/(2.*e2min) + (pow(sect,2)*sin(2*theta)*
           (fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
             (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
        cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
		
	} else {
		if(e1>e1max || e1<0){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD2-d01
					d01 = -(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
   (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
             (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
        cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
         (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
				} else {
					//AF2-d01
					d01 = (fe2min*tan(theta))/(2.*e2min) - (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) + (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (pow(sect,2)*sin(2*theta)*
           (fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
        cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD2-d01
				d01 = -(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
   (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF2-d01
					d01 = (fe2min*tan(theta))/(2.*e2min) - (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (sin(2*theta)*tan(theta)*
           ((exy*fe2min*pow(sect,2))/(2.*e2min) + (5*sqrt(5.0)*fcr*
                (pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				}
			}
			}  else {
			if(e2<=e2min) {
				//CD2-d01
				d01 = -(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
    (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
   (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
   (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((-(fe1max*cott)/(2.*e1max) + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - 
             (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
        cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF2-d01
				d01 = (fe2min*tan(theta))/(2.*e2min) - (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   ((-(fe1max*cott)/(2.*e1max) + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      ((exy*fe2min*pow(sect,2))/(2.*e2min) - (pow(sect,2)*sin(2*theta)*
           (-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
             (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
        cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
        (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
	}

return d01;
}



double 
ConcreteMcftNonLinear7 ::tangentstifness10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d10;

		if (exy>0) {
		if(e1>e1max || e1<0){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD1-d10
					d10 = (sin(2*theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))) ;
				} else {
					//AF1-d10
					d10 = ((Ec - fe2min/e2min)*sin(2*theta))/2. - ((Ec + Esv*RoV - ((Ec - fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      (cos(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD1-d10
				d10 = (sin(2*theta)*((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF1-d10
					d10 = (sin(2*theta)*(-(fe2min/e2min) - (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fe2min/e2min) - (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				}
			}
			} else if(e1<fcr/Ec){
			if(e2<=e2min){
				//AD1-d10
				d10 = (sin(2*theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*(Ec + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))) ;
				
			} else {
				//AF1-d10
				d10 = ((Ec - fe2min/e2min)*sin(2*theta))/2. - ((Ec + Esv*RoV - ((Ec - fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      (cos(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
			}
		} else {
			if(e2<=e2min) {
				//CD1-d10
				d10 = (sin(2*theta)*(fe1max/e1max + (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((fe1max/e1max + Esv*RoV - (sin(2*theta)*tan(theta)*(fe1max/e1max + 
             (pow(Ec,2)*fcu*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) - 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (cos(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF1-d10
				d10 = ((fe1max/e1max - fe2min/e2min)*sin(2*theta))/2. - ((fe1max/e1max + Esv*RoV - ((fe1max/e1max - fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      (cos(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.) ;
			}		
		}
		
	} else {
		if(e1>e1max || e1<0){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD2-d10
					d10 = (sin(2*theta)*(-Ec - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((Ec + Esv*RoV + (sin(2*theta)*tan(theta)*(-Ec - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
         (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
				} else {
					//AF2-d10
					d10 = ((-Ec + fe2min/e2min)*sin(2*theta))/2. - ((Ec + Esv*RoV + ((-Ec + fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      (cos(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD2-d10
				d10 = (sin(2*theta)*(-((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
           (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
        (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(-((pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
                (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2))) + 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF2-d10
					d10 = (sin(2*theta)*(fe2min/e2min + (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(fe2min/e2min + (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				}
			}
			} else {
			if(e2<=e2min) {
				//CD2-d10
				d10 = (sin(2*theta)*(-(fe1max/e1max) - (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((fe1max/e1max + Esv*RoV + (sin(2*theta)*tan(theta)*(-(fe1max/e1max) - 
             (pow(Ec,2)*fcu*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu)/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (cos(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF2-d10
				d10 = ((-(fe1max/e1max) + fe2min/e2min)*sin(2*theta))/2. - ((fe1max/e1max + Esv*RoV + ((-(fe1max/e1max) + fe2min/e2min)*sin(2*theta)*tan(theta))/2.)*
      (cos(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
	}
return d10;
}


double 
ConcreteMcftNonLinear7 ::tangentstifness11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV, double e1P, double e2P, double fe1max, double e1max, double fe2min, double e2min)
{
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d11;

	if (exy>0) {
		if(e1>e1max || e1<0 ){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD1-d11
					d11 = (sin(2*theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (Ec*fcu*(ex - (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
				} else {
					//AF1-d11
					d11 = (sin(2*theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (cos(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD1-d11
				d11 = (sin(2*theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
              (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF1-d11
					d11 = (sin(2*theta)*((fe2min*tan(theta))/(2.*e2min) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				}
			}
			} else {
			if(e2<=e2min) {
				//CD1-d11
				d11 = (sin(2*theta)*((fe1max*cott)/(2.*e1max) - (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   (((fe1max*cott)/(2.*e1max) + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - 
             (pow(Ec,2)*fcu*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (cos(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - (Ec*fcu*(ex - (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex - (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF1-d11
				d11 = (sin(2*theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   (((fe1max*cott)/(2.*e1max) + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (cos(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (pow(sect,2)*sin(2*theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
           (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. - 
      cos(2*theta)*tan(theta)*(fe1max - fe2min - (fe2min*(-e2P + ex - (exy*tan(theta))/2.))/e2min + 
         (fe1max*(-e1P + pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) - 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (fe1max*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
		
	} else {
		if(e1>e1max || e1<0){
			if(e1<=fcr/Ec) {
				if(e2<=e2min) {
					//AD2-d11
					d11 = (sin(2*theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      ((sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
        cos(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta))) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      (pow(sect,2)*sin(2*theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
           (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. + 
      cos(2*theta)*tan(theta)*(-(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
         (Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))));
				} else {
					//AF2-d11
					d11 = (sin(2*theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) + (sin(2*theta)*tan(theta)*((Ec*cott)/2. + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (cos(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) + 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/2.);
				}
			} else {
			if(e2<=e2min) {
					//BD2-d11
				d11 = (sin(2*theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/(ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*((Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));
				} else {
					//BF2-d11
					d11 = (sin(2*theta)*((fe2min*tan(theta))/(2.*e2min) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
   ((Esv*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
        (sin(2*theta)*tan(theta)*((fe2min*tan(theta))/(2.*e2min) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) + 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
      (pow(sect,2)*sin(2*theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
      cos(2*theta)*tan(theta)*(fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))); 
				}
			}
			} else {
			if(e2<=e2min) {
				//CD2-d11
				d11 = (sin(2*theta)*((fe1max*cott)/(2.*e1max) - (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
         (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
        (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2. - 
   ((-(fe1max*cott)/(2.*e1max) + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) - 
             (pow(Ec,2)*fcu*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*fcu*tan(theta))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu))))))/2.)*
      (cos(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
              (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
             (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
            (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + (Ec*fcu*(ex + (exy*tan(theta))/2.))/
          (ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*(-(pow(Ec,2)*exy*fcu*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + Ec/(Ec - fcu/ecu)))/
            (2.*pow(ecu,2)*pow(Ec - fcu/ecu,2)*pow(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)),2)) + 
           (Ec*exy*fcu*pow(sect,2))/(2.*ecu*(Ec - fcu/ecu)*(-1 + Ec/(Ec - fcu/ecu) + pow((ex + (exy*tan(theta))/2.)/ecu,Ec/(Ec - fcu/ecu)))) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			} else {
				//CF2-d11
				d11 = (sin(2*theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2. - 
   ((-(fe1max*cott)/(2.*e1max) + Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (sin(2*theta)*tan(theta)*((fe1max*cott)/(2.*e1max) + (fe2min*tan(theta))/(2.*e2min)))/2.)*
      (cos(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
        (sin(2*theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
             (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.))/
    ((fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max + 
      Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) + 
      (pow(sect,2)*sin(2*theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
           (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2. + 
      cos(2*theta)*tan(theta)*(-fe1max + fe2min + (fe2min*(-e2P + ex + (exy*tan(theta))/2.))/e2min - 
         (fe1max*(-e1P + pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max) + 
      (sin(2*theta)*tan(theta)*((exy*fe2min*pow(sect,2))/(2.*e2min) - 
           (fe1max*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/e1max))/2.);
			}		
		}
	}

return d11;
}
