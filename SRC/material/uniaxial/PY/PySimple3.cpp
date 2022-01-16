/* *********************************************************************
**    Module:	PySimple3.cpp 
**
**    Purpose:	Provide a p-y spring for OpenSees to better capture
**              small-strain nonlinear and viscoelastic behavior.
**
**    Developed by Benjamin Turner and Scott J. Brandenberg
**    (C) Copyright 2013 and 2017, All Rights Reserved.
**
** ****************************************************************** */

// $Revision: 1.0 to include viscoelasticity and improve numerical stability
// $Date: 	2017/02/22
// $By:		BJT
// $Source: /OpenSees/SRC/material/uniaxial/PY/PySimple3.cpp

// Written: Rev. 0, SJB
// Created: Rev. 0, Dec 2013
// tested and checked: BJT 2015
//
// Description: This file contains the class implementation for PySimple3

#include <stdlib.h>
#include <math.h>

#include "PySimple3.h"
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <elementAPI.h>

//this allows the global time variable to be accessed
#include <cstdlib>

// Controls on internal iteration between spring components
const int PYmaxIterations = 200;
const double PYtolerance = 1.0e-12;

void * OPS_ADD_RUNTIME_VPV(OPS_PySimple3)
{
  int numdata = OPS_GetNumRemainingInputArgs();
  if (numdata < 5) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial PySimple3 tag? pult? pyield? ke? C? dashpot? " << endln;
    return 0;
  }
    
  int idata[1];
  numdata = 1;
  if (OPS_GetIntInput(&numdata, idata) < 0) {
    opserr << "WARNING invalid int inputs\n";
    return 0;
  }
    
  double ddata[5] = {0,0,0,0,0};
  numdata = OPS_GetNumRemainingInputArgs();
  if (numdata > 5) numdata = 5;
  if (OPS_GetDoubleInput(&numdata, ddata) < 0) {
    opserr << "WARNING invalid double inputs\n";
    return 0;
  }
    
  UniaxialMaterial *theMaterial = 0;
  theMaterial = new PySimple3(idata[0], MAT_TAG_PySimple3, ddata[0], ddata[1],
			      ddata[2], ddata[3], ddata[4]);
    
  return theMaterial;
}



/////////////////////////////////////////////////////////////////////
//	Constructor with data

PySimple3::PySimple3(int tag, int classtag, double p_ult, double p_yield,
		     double k_e, double C_, double dashpot_)
  :UniaxialMaterial(tag,classtag),
   pult(p_ult), pyield(p_yield), ke(k_e), C(C_), dashpot(dashpot_)
{
  // Initialize PySimple variables and history variables
  //
  this->revertToStart();
  initialTangent = Ttangent;
}

/////////////////////////////////////////////////////////////////////
//	Default constructor

PySimple3::PySimple3()
  :UniaxialMaterial(0,0),
   pult(0.0), pyield(0.0), ke(0.0), C(0.0), dashpot(0.0)
{
}

/////////////////////////////////////////////////////////////////////
//	Default destructor
PySimple3::~PySimple3()
{
  // Does nothing
}

/////////////////////////////////////////////////////////////////////
//	Residual function for plastic component
double
PySimple3::getResidual(double ke, double Cp, double Tp, double dy, double pu, double C, double Tpin, double dashpot, double tstepCurrent, double dyELast, double CyeTotal, double tstepLast, double Pveguess, double bump)
{
  signdy = sign(dy);
  double signVEguess = sign(Pveguess-Cp);
  if(pu>Tp*signdy)
    {	
      if(tstep != 0.0)
	{		
	  return    (C*ke*(dy-((Tp-Cp+(dashpot*(dyELast/tstepLast)-bump))/(ke+(dashpot/tstepCurrent)))))+(Tp-Cp)+(Tpin-pu*signVEguess)*(log(pu-Cp*signVEguess)-log(pu-Tp*signVEguess));
	}
      else{
	// for static analysis with tstep = 0 the previous form would have blown up...
	return ((Tp-Cp)*(1.0-1.0/C)+((Tpin-pu*signdy)*(log(pu-Tp*signdy)-log(pu-Cp*signdy)))/C - ke*dy);	
      }
    }	
  else
    return 0.0;
}

/////////////////////////////////////////////////////////////////////
//	Sign function
int
PySimple3::sign(double val)
{
  if (val > 0) return 1;
  if (val < 0) return -1;
  return 0;
}

/////////////////////////////////////////////////////////////////////
// Main function
int 
PySimple3::setTrialStrain (double y, double yRate)
{
  // Set trial values for displacement and load in the material
  int Mtag = this->getTag();
  if (fabs(y-Cy) <= 0.0000000001)
    //	Avoid rounding errors
    {
      y = Cy;
    }
  dy = y - Cy;	// trial displacement increment
  yLast = Cy-dyLast;
  sysTimeStep = ops_Dt;
  tstep = sysTimeStep;
  tstep1 = tstep;		//these will only differ from the EQ record timestep when sub-stepping is used for first yielding
  tstep2 = CtstepLast;
  if(dyELast == 0.0) {
    tstep2 = tstep;		//avoid blowup computation of previous dashpot force
    CtstepLast = tstep;
  }	
  TyRate = dy/tstep;
  Ty = y;
  signdy = sign(dy);	//note sign(dy) gives same value as sign(yRate)
  signdyLast = sign(dyLast);
  TpinF = CpinF;
  TpinB = CpinB;
  double CpinFB4check = CpinF;
  double CpinBB4check = CpinB;
  TpinUse = CpinLast;
  Tyin = Cyin;
  TLastYieldDir = CLastYieldDir;

  if(yRate == 0.0)
    {
      P1 = Cp+ke*dy;	//this is the viscoelastic force predictor which assumes all the displacement is elastic. If yielding occurs, a revised version of this equation will be used which considers only the viscoeslatic portion of the deformation.
    }
  else
    {
      P1 = Cp+ke*dy+((dashpot/tstep1)*dy)-((dashpot/tstep2)*dyELast);	
    }
  if (signdy!=0 && signdy != TLastYieldDir && sign(P1-Cp)!=signdy) {
    P1 = Cp+(0.000001*signdy);	//this ensures that when the displacement reverses, the load does not continue in the previous direction. This would only occur in unusual situations where deceleration yielding or perhaps a load reversal had also occurred during the previous step.	
  }	
  Tp = P1;
  // 	Ttangent
  if(yRate == 0.0)
    {
      Ttangent = ke;
    }
  else if (dy != 0.0 && dashpot == 0.0) {
    Ttangent = ke;	
  }
  else if (dy != 0.0) {
    Ttangent = ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*dy)));
  }
  else if (dy == 0.0 && dashpot != 0.0 && sign(Tp-Cp) != 0.0){
    // for case where dy is zero, if previous dashpot force was nonzero, the tangent will either be pos. or neg. infinity, since the force will
    // change without y changing. This is handled as a special case by the code below, but for now the tangent is saved as ke
    Ttangent = Ctangent;
  }
  else {
    Ttangent = ke;
  }
  Tpalpha = Cpalpha;
  f = fabs(Tp-Cpalpha)-pyield;	//yield function
  if (f <= 0.0)
    {
      TyeTotal =  (Tp-Cp+ke*CyeTotal+dashpot*CyeTotal/tstep1+dashpot*dyELast/tstep2)/(ke+dashpot/tstep1);
      TdyE    = TyeTotal-CyeTotal;
      return 0;	//no yielding
    }
	
  // Return mapping if yielding occurs
  if(f>0)
    {
      double CpUse = Cp;
      int signdyLast = sign(dyLast);
      int signP = sign(Tp-Cp);
      double dyELastUse = dyELast;
      double CyeTotalUse = CyeTotal;
      double bump = 0.0;
      /////If the current P is very close to Pult and displacement increment is in same dir. as last step,
      // // or the displacement increment is zero but the viscoelastic guess still lands above Pult,
      // // there can be convergence issues because the residual for the correct value of next P gets very
      // // small. Instead just assign Pult minus tol. as next P.
      if(  (fabs(Cp) > 0.99*pult && fabs(P1) > fabs(pult) && signdy == signdyLast) || (fabs(Cp) > 0.99*pult && fabs(P1) > fabs(pult) && signdy == 0)  )
	{
	  Tp = (pult-PYtolerance)*signdy;
	  if (signP == 1) {TpinUse = TpinF;}
	  else if (signP == -1) {TpinUse = TpinB;}
	  else {TpinUse = CpinLast;}				
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if (dy == 0.0 && dashpot == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  Tpalpha = Tp - pyield*signdy;
	  // Compute viscoelastic displacement increment
	  TLastYieldDir = signP;
	  return 0;	
	}

      // Test for deceleration yielding
      //Special case... if displacement increment is in same dir. as last converged step but velocity is much less
      // than the velocity in the previous step (or zero), the decreased force in the damper component could cause the
      // viscoelastic force to decrease relative to Cp (this can only occur if the decrease in dashpot force is large
      // relative to the increase in force in the elastic component, so with a relatively large dashpot coeff. and large 
      // deceleration). This is somewhat counterituitive, because you would expect continuing to displace in the same
      // dir. would keep increasing the load, but a large loss of damper force could cause the cumulative load to drop.
      // If this drop is large enough, yielding could occur in the opposite direction of displacement, so the previous 
      // if statements for re-assigning Tpin and Tyin won't work because they rely on sign(dy). Instead, need to use
      // the term sign(P1-Cp), where P1 is the viscoelastic predictor guess.
      //	If sign(P1-Cp) and signdy have opposite signs, this condition has occurred...
	
      if(signP != signdy && signP != 0 && dy == 0.0 && signP != CLastYieldDir)
	// special case for dy = zero
	{
	  //compute new Tpin based on signP instead of signdy...
	  TpinUse = Cpalpha + pyield*signP;
	  Tyin = Cy;	
	  //Compute the change in viscoelastic displacement that occurs when the decrease in dashpot force causes P
	  //	to drop from the current P to the yield surface. First compute time this takes by proportioning forces:
	  TLastYieldDir = signP;
	  pn1_a = TpinUse;									//lower bracket is current committed value of p
	  pn1_b = (1.0-0.5*PYtolerance)*signP*pult;			//upper bracket is pult, but with opposite sign
	  //...then don't change anything, just send in the regular values...
	  CpUse = TpinUse;
	  tstep1 = tstep;
	  tstep2 = tstep;
	  dyELastUse = dyELast;
	  CyeTotalUse = CyeTotal;
	  bump=Cp-TpinUse;
	}
      else if(signP != signdy && signP != 0 && signP != CLastYieldDir)
	//	Declaration yielding for a non-zero dy case, still flip the search direction:
	{
	  TpinUse = Cpalpha - pyield*signdy;
	  Tyin = Cy + (TpinUse-Cp)/(-1.0*((Tp-Cp)/dy));
	  CpUse = TpinUse;
	  TLastYieldDir = (-1)*signdy;
	  pn1_a = TpinUse;
	  pn1_b = (1.0-0.5*PYtolerance)*(-1.0)*signdy*pult;	//upper bracket is pult, but with opposite sign
	  bump = Cp-TpinUse;
	}
      else
	{
	  // if no deceleration yielding, only re-assign Tpin if previous yield was not in the current direction
	  if(TLastYieldDir != signdy && signdy!= 0)
	    {
	      TpinUse = Cpalpha + pyield*signdy;
	      Tyin = Cy + (TpinUse-Cp)/((Tp-Cp)/dy);
	    }
	  pn1_a = Cp;											//lower bracket is current committed value of p
	  if(signP == 0)
	    {
	      pn1_b = (1.0-0.5*PYtolerance)*CLastYieldDir*pult;			//upper bracket is pult
	      TLastYieldDir = CLastYieldDir;
	    }
	  else{
	    pn1_b = (1.0-0.5*PYtolerance)*signP*pult;			//upper bracket is pult
	    TLastYieldDir = signP;
	  }
	}
	
      // Bumping routine and update backstress
      //If declaration yielding occurred, bumped values have already been computed and this loop will not be entered.
      if((pn1_a < TpinUse && pn1_b > TpinUse) || (pn1_a > TpinUse && pn1_b < TpinUse)) // test to see if brackets are on opposite sides of yield surface
	{
	  pn1_a = TpinUse;
	  CpUse = TpinUse;
	  bump = TpinUse-Cp;	
	}	
      //Determine if backstress should be updated
      double TpinUseB4check = TpinUse;
      if (signP == 1 && TpinUse < CpinF) {
	TpinF = TpinUse;
      }
      else if (signP ==1 && TpinUse >= CpinF && CLastYieldDir != 0){
	TpinUse = CpinF;
	TpinF = CpinF;
      }
      else if (signP ==1 && TpinUse >= CpinF && CLastYieldDir == 0){
	TpinF = TpinUse;	
      }
      if (signP == -1 && TpinUse > CpinB) {
	TpinB = TpinUse;
      }
      else if (signP ==-1 && TpinUse <= CpinB && CLastYieldDir != 0){
	TpinUse = CpinB;
	TpinB = CpinB;
      }
      else if (signP ==-1 && TpinUse <= CpinB && CLastYieldDir == 0){
	TpinB = TpinUse;
      }		
      if (signP == 0){
	TpinUse = CpinLast;
      }

      //Compute residuals 	
      R1 = getResidual(ke,CpUse,pn1_a,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);
      R2 = getResidual(ke,CpUse,pn1_b,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);
		
      // Iterate to find P next			
      //first check if either of the residuals computed at the starting brackets satisfies the tolerance			
      if(fabs(R1)<PYtolerance*pult)
	{
	  Tp = pn1_a;
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  //check if declaration cause a decrease in force
	  signP = sign(Tp-Cp);
	  if(signP != signdy)
	    {
	      signPalphaNew = (-1)*signdy;
	    }
	  else 
	    {
	      signPalphaNew = signdy;
	    }
	  Tpalpha = Tp - pyield*signPalphaNew;
	  return 0;
	}
      if(fabs(R2)<PYtolerance*pult)
	{
	  Tp = pn1_b;
	  TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	  TdyE    = TyeTotal-CyeTotalUse;
	  //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	  if(dy==0.0 && TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else if (dy==0.0 && dashpot != 0.0){
	    Ttangent = Ctangent;
	  }
	  else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
	    Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	  }
	  else{
	    Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	  }
	  //check if declaration cause in decrease in force
	  signP = sign(Tp-Cp);
	  if(signP != signdy)
	    {
	      signPalphaNew = (-1)*signdy;
	    }
	  else 
	    {
	      signPalphaNew = signdy;
	    }
	  Tpalpha = Tp - pyield*signPalphaNew;
	  return 0;
	}
			
      // If R1 and R2 weren't within tolerance, check to see if they have the same sign...
      if(sign(R1)==sign(R2))
	{
	  // 	Check for a condition in which, for a very large trial displacement step (keep in mind the solver may send a trial displacement step that is significantly larger than the eventual converged disp. step), the state goes from being either elastic or mostly elastic all the way to the yield surface. If this jump happens from far away from the yield surface, the residual for the upper bound guess can be very large and won't satisfy tolerance, even though it is essentially the correct guess. Test for this condition by nudging the upper bound guess even closer to pult and seeing if the residual continues to decrease.
	  pn1_b = (1.0-0.0000000001*PYtolerance)*TLastYieldDir*pult;
	  R2 = getResidual(ke,CpUse,pn1_b,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);
	  if(fabs(R2)<PYtolerance*pult)
	    {
	      Tp = (1.0-0.5*PYtolerance)*TLastYieldDir*pult;	//keep consistent with other cases to avoid tiny fluctuations near the yield surface
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if declaration cause in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return 0;
	    }	
	}

      //start Ridder's loop, initialize variables
      double pn1_3 = 0.0;
      double pn1_4 = 0.0;
      double R3 = 0.0;
      double R4 = 0.0;
      double S = 0.0;
      int i=0;
      while((fabs(R1)>PYtolerance*pult) && (fabs(R2)>PYtolerance*pult))
	{
	  pn1_3 = (pn1_a+pn1_b)/2.0;
	  R3    = getResidual(ke,CpUse,pn1_3,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);
	  if(fabs(R3)<PYtolerance*pult)
	    {
	      Tp = pn1_3;
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if declaration causes in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return 0;
	    }
	  S = (sqrt(R3*R3 - R1*R2));
	  pn1_4 = pn1_3+(pn1_3-pn1_a)*((sign(R1-R2)*R3)/(sqrt(R3*R3 - R1*R2)));
	  R4    = getResidual(ke,CpUse,pn1_4,dy,pult,C,TpinUse,dashpot,tstep1,dyELastUse,CyeTotalUse,tstep2,P1,bump);	
	  if(fabs(R4)<PYtolerance*pult)
	    {
	      Tp = pn1_4;
	      TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
	      TdyE    = TyeTotal-CyeTotalUse;
	      //for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
	      if(dy==0.0 && TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else if (dy==0.0 && dashpot != 0.0){
		Ttangent = Ctangent;
	      }
	      else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
	      }
	      else{
		Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
	      }
	      //check if declaration cause in decrease in force
	      signP = sign(Tp-Cp);
	      if(signP != signdy)
		{
		  signPalphaNew = (-1)*signdy;
		}
	      else 
		{
		  signPalphaNew = signdy;
		}
	      Tpalpha = Tp - pyield*signPalphaNew;
	      return 0;
	    }
	  if(sign(R3) != sign(R4))
	    {
	      pn1_a = pn1_3;
	      pn1_b = pn1_4;
	      R1 = R3;
	      R2 = R4;
	    }
	  else
	    {
	      if(sign(R1) != sign(R4))
		{
		  pn1_a = pn1_a;
		  pn1_b = pn1_4;
		  R1 = R1;
		  R2 = R4;
		}
	      else
		{
		  if(sign(R2) != sign(R4))
		    {
		      pn1_a = pn1_4;
		      pn1_b = pn1_b;
		      R1 = R4;
		      R2 = R2;
		    }
		  else
		    {
		      //printf ("none of the Ridder's values had opposite signs-- error! \n");*/
		    }
		}
	    }
	  i++;
	  if(i==PYmaxIterations){
					
	    if(fabs(pn1_a) > 0.995*pult && fabs(pn1_b) > 0.995*pult)
	      {
		Tp = (pn1_a+pn1_b)/2.0;
		if(fabs(Tp) >= pult){
		  Tp = (pult-PYtolerance)*signdy;
		}
		TyeTotal = (Tp-CpUse+ke*CyeTotalUse+dashpot*CyeTotalUse/tstep1+dashpot*dyELastUse/tstep2)/(ke+dashpot/tstep1);
		TdyE    = TyeTotal-CyeTotalUse;
		//for tangent, consider entire increment (i.e. not just the post-bump behavior if first yield occurred)
		if(dy==0.0 && TdyE == 0.0){
		  Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
		}
		else if (dy==0.0 && dashpot != 0.0){
		  Ttangent = Ctangent;
		}
		else if ((dy == 0.0 && dashpot == 0.0) || TdyE == 0.0){
		  Ttangent = 	1.0/((1.0/(ke)) + (1.0/(C*ke*((Tp-pult*sign(signdy))/(TpinUse-Tp)))));
		}
		else{
		  Ttangent = 	1.0/((1.0/(ke+dashpot*((1.0/tstep1)-(dyELast/(tstep2*TdyE))))) + (1.0/(C*ke*((Tp-pult*sign(Tp-Cp))/(TpinUse-Tp)))));
		}
		//check if declaration cause in decrease in force
		signP = sign(Tp-Cp);
		if(signP != signdy)
		  {
		    signPalphaNew = (-1)*signdy;
		  }
		else 
		  {
		    signPalphaNew = signdy;
		  }
		Tpalpha = Tp - pyield*signPalphaNew;
		// Compute viscoelastic displacement increment
		return 0;	
	      }
	    else {
	      opserr << "Ridder's method for material tag " << this->getTag() << " failed to find a working value for P in " << PYmaxIterations << " iterations." << endln;
	      Tp = Cp;
	      Ttangent = Ctangent;
	      Tpalpha = Cpalpha;
	      TdyE = dy;
	      return 0;
	    }
	  }
	}	
			
			
    }	//end return mapping algorithm
  return 0;
}

/////////////////////////////////////////////////////////////////////
double 
PySimple3::getStress(void)
{
  return this->Tp;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple3::getTangent(void)
{
  return this->Ttangent;
}


/////////////////////////////////////////////////////////////////////
double 
PySimple3::getInitialTangent(void)
{
  return this->initialTangent;	
}

/////////////////////////////////////////////////////////////////////
double 
PySimple3::getDampTangent(void)
{
  //	The damping tangent is dp/d(yRate), i.e. the change in force per change in velocity. 
  // 	The damping tangent is just the dashpot coeff if load remains in elastic region. 
  // 	If yielding occurs, proportion by the displacement in the viscoelastic component.
	
  double DampTangent;
	
  if (tstep == 0.0) {
    DampTangent = 0.0;	//for static analysis, no dependence on velocity
  }
  else if (dashpot == 0.0) {
    DampTangent = 0.0;
  }
  else if (dy == 0.0 && TdyE == 0.0) {
    DampTangent = (Tp-Cp);
  }
  else if (dy == 0.0) {
    DampTangent = (Tp-Cp);
  }
  else {
    DampTangent	=	dashpot*(TdyE/dy);
  }
  return DampTangent;
}
/////////////////////////////////////////////////////////////////////
double 
PySimple3::getStrain(void)
{
  return this->Ty;
}

/////////////////////////////////////////////////////////////////////
double 
PySimple3::getStrainRate(void)
{
  return this->TyRate;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple3::commitState(void)
{
  // Commit trial state variables
  Cp = Tp;
  dyLast = Ty-Cy;
  Cy = Ty;
  CpinF = TpinF;
  CpinB = TpinB;
  CpinLast = TpinUse;
  Cyin = Tyin;
  Cpalpha = Tpalpha;
  Ctangent = Ttangent;
  dyELast = TdyE;
  CyeTotal = TyeTotal;
  CtstepLast = tstep1;
  CLastYieldDir = TLastYieldDir;
  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple3::revertToLastCommit(void)
{
  Tp = Cp;
  Ty = Cy;
  TpinF = CpinF;
  TpinB = CpinB;
  TpinUse = CpinLast;
  Tyin = Cyin;
  Tpalpha = Cpalpha;
  Ttangent = Ctangent;
  dyLast = dyLast;
  dyELast = dyELast;
  TyeTotal = CyeTotal;
  tstep1 = CtstepLast;
  TLastYieldDir = CLastYieldDir;
  return 0;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple3::revertToStart(void)
{
  dy = 0.0;
  CLastYieldDir = 0;
  TLastYieldDir = 0;
  signdy = 0;
  Cp = 0.0;
  Cy = 0.0;
  Cyin = 0.0;
  CpinF = 0.0;
  CpinB = 0.0;
  CpinLast = 0.0;
  Ty = 0.0;
  Tyin = 0.0;
  Tp = 0.0;
  TyRate = 0.0;
  Ttangent = ke;
  TpinF = 0.0;
  TpinB = 0.0;
  TpinUse = 0.0;
  Tpalpha = 0.0;
  dyLast = 0.0;
  signdyLast = 0.0;
  dyELast = 0.0;
  lam = 0;		//lambda, plastic multiplier
  lamLB = 0;
  lamUB = 0;		//lower- and upper-bound guesses for lambda
  ypDot = 0;		//plastic displacement rate
  P1 = 0;			//Force (dP1/dt) in viscoelastic spring component
  P2 = 0;			//Force (dP2/dt) in plastic spring component
  tstep = 0;		//time step, computed from yRate
  yLast = 0;		//y from the timestep before the currently committed timestep, i.e. y(i)=y, y(i-1)=Cy, y(i-2)=yLast
  ypRate = 0;		//plastic deformation rate
  yeRate = 0;		//elastic deformation rate
  TdyP = 0;		//plastic deformation increment
  TdyE = 0;		//elastic deformation increment
  signdy = 0;
  residualLam = 0;	//residual when iterating to determine correct value of plastic multiplier
  Rlam1 = 0;		//residual for lambda using lower- and upper-bound guesses
  Rlam2 = 0;
  sysTimeStep = 0;	//system time
  bumped = 0;
  pP1 = 0;
  dyELastUse = 0;
  signBump = 0.0;
  P1veGuess = 0.0;
  CyeTotal = 0.0;
  TyeTotal = 0.0;
  CtstepLast = 0.0;
  tstep1 = 0.0;
  tstep2 = 0.0;
	
  this->commitState();
  return 0;
}

/////////////////////////////////////////////////////////////////////
UniaxialMaterial *
PySimple3::getCopy(void)
{
  PySimple3 *theCopy;			// pointer to a PySimple3 class
  theCopy = new PySimple3();	// new instance of this class
  *theCopy= *this;			// theCopy (dereferenced) = this (dereferenced pointer)
  return theCopy;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple3::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(18);
  
  data(0) = this->getTag();
  data(1) = pult;
  data(2) = pyield;
  data(3) = ke;
  data(4) = C;
  data(5) = dashpot;
  
  data(6) = Cp;
  data(7) = Cy;
  data(8) = CpinF;
  data(9) = Cpalpha;
  data(10) = Cyin;
  data(11) = Ctangent;
  data(12) = dyLast;
  data(13) = dyELast;
  data(14) = CyeTotal;
  data(15) = CtstepLast;
  data(16) = CpinB;
  data(17) = CpinLast;


  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "PySimple3::sendSelf() - failed to send data\n";

  return res;
}

/////////////////////////////////////////////////////////////////////
int 
PySimple3::recvSelf(int cTag, Channel &theChannel, 
		    FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(18);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
    opserr << "PySimple3::recvSelf() - failed to receive data\n";
    this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    pult	 = data(1);
    pyield   = data(2);
    ke       = data(3);
    C     	 = data(4);
    dashpot	 = data(5);
    Cp       = data(6);
    Cy		 = data(7);
    CpinF	 = data(8);
    Cpalpha  = data(9);
    Cyin	 = data(10);
    Ctangent = data(11);
    dyLast	 = data(12);
    dyELast  = data(13);
    CyeTotal  = data(14);
    CtstepLast  = data(15);
    CpinB		= data(16);
    CpinLast	= data(17);

    // set the trial quantities
    this->revertToLastCommit();
  }

  return res;
}

/////////////////////////////////////////////////////////////////////
void 
PySimple3::Print(OPS_Stream &s, int flag)
{
  s << "PySimple3, tag: " << this->getTag() << endln;
  s << "  pult: " << pult << endln;
  s << "  pyield: " << pyield << endln;
  s << "  ke: " << ke << endln;
  s << "  C: " << C << endln;
  s << "  dashpot: " << dashpot << endln;
  s << "  disp: " << Cy << endln;
  s << "  force: " << Cp << endln;
  s << "  CpinF: " << CpinF << endln;
  s << "  Cy: " << Cy << endln;
  s << "  Cpalpha: " << Cpalpha << endln;
  s << "  Cyin: " << Cyin << endln;
  s << "  Ctangent: " << Ctangent << endln;
  s << "  dyLast: " << dyLast << endln;
  s << "  dyELast: " << dyLast << endln;
  s << "  CyeTotal: " << CyeTotal << endln;
  s << "  CpinB: " << CpinB << endln;
  s << "  CpinLast: " << CpinLast << endln;
	
}
