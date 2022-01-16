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

// $Revision: 0 $
// $Date: May 2015 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/BilinearOilDamper.h,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: May 2015
// Updated: May 2015
// Revision: A

// Description: This file contains the class interface for 
// Oil Damper Model Relationship of the form  before relief valve ==> F = K*u_s = C*V_d  after relief valve ==> F= K*u_s = Fr + p*C*(V_d-Fr/C)
//
// References: 
// Akcelyan, S., and Lignos, D.G. (2015), “Adaptive Numerical Method Algorithms for Nonlinear Viscous and Bilinear Oil Damper Models Under Random Vibrations”, ASCE Journal of Engineering Mechanics, (under review)
// Kasai, K., Takahashi, O., and Sekiguchi, Y. (2004). "JSSI manual for building passive control technology part-10 time-history analysis model for nonlinear oil dampers." Proc., The 13th World Conference on Earthquake Engineering, Vancouver, B.C., Canada.

// Variables:
// $K: Elastic stiffness of linear spring (to model the axial flexibility of a viscous damper (brace and damper portion)
// $C: Viscous damping coefficient of the damper
// $Fr: Relief load
// $p: post-relief viscous damping coefficient ratio, (p=C2/C)
// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2 = Finite differences)
// $RelTol: Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol: Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-10)
// $MaxHalf: Maximum number of sub-step iterations within an integration step (default value 15)

#include <BilinearOilDamper.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>
#include <OPS_Globals.h>
#include <elementAPI.h>
#include <Parameter.h>

static int numBilinearOilDamperMaterials = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_BilinearOilDamper)
{
  if (numBilinearOilDamperMaterials == 0) {
    numBilinearOilDamperMaterials++;
    opserr << "BilinearOilDamper Model by Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[9];
  int numData = 1;
        // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  BilinearOilDamper tag" << endln;
    return 0;
  }
  // Check if the input variables 
  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 2 && numData != 4 && numData != 5 && numData != 9 ) {
    opserr << "Invalid #args, want: uniaxialMaterial BilinearOilDamper " 
	   << iData[0] <<  " K? C? <Fr? p?> <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    return 0;
  }

  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args want: uniaxialMaterial BilinearOilDamper " 
	   << iData[0] << " K? C? <Fr? p?> <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    
    return 0;   
  }

    if (numData == 2) {
		//Default variables
      dData[2] = 1.0;
      dData[3] = 1.0;
      dData[4] = 0.0;
      dData[5] = 1;
      dData[6] = 0.000001;
      dData[7] = 0.0000000001;
      dData[8] = 15;
    }
    if (numData == 4){
      //Default variables
      dData[4] = 0.0;
      dData[5] = 1;
      dData[6] = 0.000001;
      dData[7] = 0.0000000001;
      dData[8] = 15;
    }
    if (numData == 5){
      //Default variables
      dData[5] = 1;
      dData[6] = 0.000001;
      dData[7] = 0.0000000001;
      dData[8] = 15;
    }
    
    // Parsing was successful, allocate the material with zero index
    theMaterial = new BilinearOilDamper(iData[0], 
					dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8]);
    
    if (theMaterial == 0) {
      opserr << "WARNING could not create uniaxialMaterial of type BilinearOilDamper Material\n";
      return 0;
    }
    
    return theMaterial;
}


BilinearOilDamper::BilinearOilDamper(int tag, double k, double c, double fr, double pp, double lgap, double nm, double reltol, double abstol, double maxhalf)
 :UniaxialMaterial(tag,MAT_TAG_BilinearOilDamper), K(k), C(c), Fr(fr), p(pp), LGap(lgap), NM(nm), RelTol(reltol), AbsTol(abstol), MaxHalf(maxhalf)       
{
  if (p < 0.0) {
    opserr << "BilinearOilDamper::BilinearOilDamper -- p < 0.0, setting to 0.0\n";
    p = 0.0;
  }
  
  //initialize variables
  this->revertToStart();
}

BilinearOilDamper::BilinearOilDamper()
 :UniaxialMaterial(0,MAT_TAG_BilinearOilDamper),
  K(0.0), C(0.0), Fr(0.0), p(0.0), LGap(0.0),  NM(0.0), RelTol(0.0), AbsTol(0.0), MaxHalf(0.0)         
{
  this->revertToStart();
}

BilinearOilDamper::~BilinearOilDamper()
{
  // does nothing
}

int 
BilinearOilDamper::setTrialStrain(double strain, double strainRate)
{
  //all variables to the last commit
  this->revertToLastCommit();
  
  
  // Determine the strain rate and acceleration 
  double dStrain = (strain - Tstrain);
  double Vel, fd0, acc, vel1, vel0;
  if (fabs(strainRate) == 0.0) { //static analysis
    Vel = 0.0;
    acc = 0.0;
    
  } else { 
    
    Vel = strainRate;
    acc = (Vel - TVel)/ops_Dt;
  }
  double smin = pow(0.5,MaxHalf);
  double s = 1.0;
  double stot = 0.0;
  double it = 0.0;
  
  // Method1 : Dormand Pronce Method (ODE Solver)
  if (NM == 1.0) { 
    double h, yt, eps, error;
    vel0 = TVel;  // Velocity of the previous step.
    fd0 = Tstress; 
    
    while (it < 1.0) { //iteration
      h = s * ops_Dt; // Time step 
      vel1 = vel0 + acc * h; // Velocity at the time step h
    
      DormandPrince(vel0, vel1, fd0, h, yt, eps, error);
      
      
      // Error check: Adaptive Step Size
      if ((eps <= RelTol) || (s == smin) || (fabs(error) <= AbsTol)) {
	vel0 = vel1;
	fd0 = yt;
	stot = stot+s;
      } else {
	if (s > smin) {
	  s=0.5*s; // step gets smaller -now try this step again.
	} else {
	  s=smin;
	}
      }
      
      if (stot == 1.0) { // The total internal stepsize reached dt
	it=1.0;
      }
 }
    
    if (p == 0) {
      if (fabs(fd0) > Fr) {
	fd0 = sgn(fd0)*Fr;
      }
    }
    
  }
  
  // Method2 : Numerical Integration
  if (NM == 2.0) { 
    double h, fdk1, fdk2, eps, error;
    
    
    while (it < 1.0) { //iteration
      for ( int k = 1; k < 3 ; k = k + 1) {
	if (k == 1) {
	  h = s * ops_Dt; // Time step 
	} else {
	  h = s/(s+1) * ops_Dt; // Time step 
	}
	
	vel0 = TVel;  // Velocity of the previous step.
	fd0 = Tstress;
	for ( int j = 1; j < (ops_Dt/h + 1) ; j = j + 1) {
	  vel1 = vel0 + acc * h; // Velocity at the time step h
	  double fd1 = (K*vel1*h+fd0)/(1+K*h/C);
	  
	  if (fd1 > Fr){ // if the force exceeds the relief force
	    if (p == 0.) {
	      fd1 = Fr;
	    } else {
	      fd1 = (K*vel1*h+Fr*(1-p)*K*h/(p*C)+fd0)/(1+K*h/(p*C));
	    }
	  }
	  
	  if (fd1 < -Fr) { // if the force exceeds the relief force
	    if (p == 0.) {
	      fd1 = -Fr;
	    } else {
	      fd1 = (K*vel1*h-Fr*(1-p)*K*h/(p*C)+fd0)/(1+K*h/(p*C));
	    }
	  }
	  
	  fd0 = fd1;
	  vel0 = vel1;
	}
	if (k == 1) {
	  fdk1 = fd0; // Solution for h = s * ops_Dt; 
	} else {
	  fdk2 = fd0; // Solution for  h = s/(s+1) * ops_Dt;
	} 
	
	error = fdk2 - fdk1;
	eps = fabs(error/fdk2);
	
      }
      // Error check: Adaptive Step Size
      if ((eps <= RelTol) || (s == smin) || (fabs(error) <= AbsTol)) {
	it = 1;
      } else {
	if (s > smin) {
	  s= 0.5*s; // step gets smaller -now try this step again.
	} else {
	  s= smin;
	}
      }
    }
  }
  
  // End of Methods
  
  
  // Effect of gap start 
  
  if (LGap > 0.) {
    
    double dStrain = (strain - Tstrain);
    
    if ((fd0 > 0) && (Tstress < 0)) {  //from negative to positive
      
      Tpugr = Tstrain + dStrain * fabs(fd0)/fabs(fd0 - Tstress);  // approximate displacement for gap initiation
      Tnugr = 0.;
      
      if (fabs(strain-Tpugr) < LGap) {
	fd0 = 0.;
      }
    }  
    
    if ((fd0 < 0) && (Tstress > 0)) {  //from positive to negative
      
      Tnugr = Tstrain + dStrain * fabs(fd0)/fabs(fd0 - Tstress);  // approximate displacement for gap initiation
      Tpugr = 0.;
      
      if (fabs(strain-Tnugr) < LGap) {
	fd0 = 0.;
      }
    }
    
    // After gap initiation
    
    if  ((fabs(Tpugr) > 0.) && (Tstress == 0)) {   //from negative to positive
      
      if ((strain > Tpugr) && ((strain-Tpugr) < LGap)) {
	fd0 = 0.;
      }
    }
    
    
    
    if  ((fabs(Tnugr) > 0.) && (Tstress == 0)) {   //from positive to negative
      
      if ((strain < Tnugr) && ((strain-Tnugr) > -LGap)) {
	fd0 = 0.;
      }
    }
    
  }
  // Effect of gap end 
  
  
  
  Tstress = fd0; // Stress 
  TVel = Vel;
  Tstrain = strain;
  
  return 0;
}

double BilinearOilDamper::getStress(void)
{
  return  Tstress;
}

double BilinearOilDamper::getTangent(void)
{
  // Why is this return 0.0? -- MHS
  return 0.0;
}

double BilinearOilDamper::getInitialTangent(void)
{
  //return 0.0;
  return K; // MHS
}

double BilinearOilDamper::getDampTangent(void)
{
  //return 0.0;
  return C; // Needs to return something, something is better than nothing -- MHS
}


double 
BilinearOilDamper::getStrain(void)
{
  return Tstrain;
}

double 
BilinearOilDamper::getStrainRate(void)
{
  return TVel;
}

int 
BilinearOilDamper::commitState(void)
{
  //commit trial  variables
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  CVel = TVel;
  
  return 0;
}

int 
BilinearOilDamper::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  TVel = CVel;
  Cpugr = Tpugr;
  Cnugr = Tnugr;
  
  return 0;
}

int 
BilinearOilDamper::revertToStart(void)
{
  // Initialize state variables
  Tstrain=0.0;
  Tstress=0.0;
  Ttangent = 0.0;
  TVel = 0.0;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Ctangent = 0.0;
  CVel = 0.0;
  Tpugr = Cpugr;
  Tnugr = Cnugr;
  return 0;
}

double
BilinearOilDamper::sgn(double dVariable){ 
  if (dVariable<0.0){
    return -1.0;
  }else{
    return 1.0;
  }
}

int 
BilinearOilDamper::DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error){
  double k1, k2, k3, k4, k5, k6, k7;
  
  k1 = f(vel0, y0) * h;
  
  k2 = f((vel1 - vel0)*(1./5.) + vel0, y0 + (1./5.)*k1) * h;

  k3 = f((vel1 - vel0)*(3./10.) + vel0, y0 + (3./40.)*k1 + (9./40.)*k2) * h;

  k4 = f((vel1 - vel0)*(4./5.) + vel0, y0 + (44./45.)*k1 + (-56./15.)*k2 + (32./9.)*k3) * h;

  k5 = f((vel1 - vel0)*(8./9.) + vel0, y0 + (19372.0/6561.0)*k1 + (-25360.0/2187.0)*k2 + (64448.0/6561.0)*k3 + (-212.0/729.0)*k4) * h;

  k6 = f((vel1 - vel0)*(1.) + vel0, y0 + (9017.0/3168.0)*k1 + (-355.0/33.0)*k2 + (46732.0/5247.0)*k3 + (49.0/176.0)*k4 + (-5103.0/18656.0)*k5) * h;

  yt = y0 + (35./384.)*k1 + (500./1113.)*k3 + (125./192.)*k4 + (-2187./6784.)*k5 + (11./84.)*k6;

  k7 = f((vel1 - vel0)*(1.) + vel0, yt) * h;

  error = (71./57600.)*k1 + (-71./16695.)*k3 + (71./1920.)*k4 + (-17253./339200.)*k5 + (22./525.)*k6 + (-1./40.)*k7;

  eps = fabs(error/ yt);


 return 0;
}

double
BilinearOilDamper::f(double v, double fd){

	 if ((fabs(fd) < Fr) || (p == 0)) {

		return ( v - (fd/C) )*K;

	 } else {

		return ( v - ((sgn(fd)*(p-1.0)*Fr+fd)/(p*C)))*K;
 
	}
}

UniaxialMaterial *
BilinearOilDamper::getCopy(void)
{
    BilinearOilDamper *theCopy = new BilinearOilDamper(this->getTag(), K, C, Fr, p, LGap, NM, RelTol, AbsTol, MaxHalf);
    // Converged state variables
        theCopy->Cstrain = Cstrain;
        theCopy->Cstress = Cstress;
        theCopy->Ctangent = Ctangent;
        theCopy->CVel = CVel;
		theCopy->Cpugr = Cpugr;
		theCopy->Cnugr = Cnugr;

        // Trial state variables
		theCopy->Tstrain = Tstrain;
		theCopy->Tstress = Tstress; 
        theCopy->Ttangent = Ttangent;
        theCopy->TVel = TVel;
		theCopy->Tpugr = Tpugr;
		theCopy->Tnugr = Tnugr;

    return theCopy;
}

int 
BilinearOilDamper::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(16);
  data(0) = this->getTag();

  // Material properties
  data(1) = K;
  data(2) = C;
  data(3) = Fr;
  data(4) = p;
  data(5) = LGap;
  data(6) = NM;
  data(7) = RelTol;
  data(8) = AbsTol;
  data(9) = MaxHalf;
  
  // State variables from last converged state
  data(10) = Cstrain;
  data(11) = Cstress;
  data(12) = Ctangent;
  data(13) = CVel;
  data(14) = Cpugr;
  data(15) = Cnugr;
        
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "BilinearOilDamper::sendSelf() - failed to send data\n";

  return res;
}

int 
BilinearOilDamper::recvSelf(int cTag, Channel &theChannel, 
                               FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(16);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "BilinearOilDamper::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
        
        // Material properties
        K = data(1);
        C = data(2);
		Fr = data(3);
		p = data(4);
		LGap = data(5);
	    NM = data(6);
        RelTol = data(7);
		AbsTol = data(8);
		MaxHalf = data(9);
        
        // State variables from last converged state 
        Cstrain = data(10);
        Cstress = data(11);
        Ctangent = data(12);
        CVel = data(13);
        Cpugr = data(14);  
		Cnugr = data(15); 


          
  }
    
  return res;
}

void 
BilinearOilDamper::Print(OPS_Stream &s, int flag)
{
    s << "BilinearOilDamper tag: " << this->getTag() << endln;
    s << "  K: " << K << endln; 
    s << "  C: " << C << endln;
    s << "  Fr: " << Fr << endln;
	s << "  p: " << p << endln;
    s << "  LGap: " << LGap << endln; 
	s << "  NM: " << NM << endln; 
    s << "  RelTol: " << RelTol << endln;
	s << "  AbsTol: " << AbsTol << endln;
    s << "  MaxHalf: " << MaxHalf << endln;    
}


int
BilinearOilDamper::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0 || strcmp(argv[0],"K") == 0) {
    param.setValue(K);
    return param.addObject(1, this);
  }
  if (strcmp(argv[0],"eta") == 0 || strcmp(argv[0],"C") == 0) {
    param.setValue(C);
    return param.addObject(2, this);
  }
  if (strcmp(argv[0],"Fr") == 0) {
    param.setValue(Fr);
    return param.addObject(3, this);
  }
  if (strcmp(argv[0],"p") == 0) {
    param.setValue(p);
    return param.addObject(4, this);
  }
  if (strcmp(argv[0],"Lgap") == 0 || strcmp(argv[0],"LGap") == 0) {
    param.setValue(LGap);
    return param.addObject(5, this);
  }  
  return -1;
}


int 
BilinearOilDamper::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    K = info.theDouble;
    return 0;
  case 2:
    C = info.theDouble;
    return 0;
  case 3:
    Fr = info.theDouble;
    return 0;
  case 4:
    p = info.theDouble;
    return 0;
  case 5:
    LGap = info.theDouble;
    return 0;
  default:
    return -1;
  }
}
