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
                                                                        
// $Revision: C $
// $Date: May 2015 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousDamper.cpp,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: January 2013
// Updated: May 2015
// Revision: C
//
// Description: This file contains the class interface for 
// Viscous Damper Model Relationship of the form F = K*u_s = C*pow(V_d,alpha)
// Reference: 
// Akcelyan, S., and Lignos, D.G. (2015), “Adaptive Numerical Method Algorithms for Nonlinear Viscous and Bilinear Oil Damper Models Under Random Vibrations”, ASCE Journal of Engineering Mechanics, (under review)
// Kasai K, Oohara K. (2001). “Algorithm and Computer Code To Simulate Response of Nonlinear Viscous Damper”. Proceedings Passively Controlled Structure Symposium 2001, Yokohama, Japan.

// Variables:
// $K: Elastic stiffness of linear spring (to model the axial flexibility of a viscous damper (brace and damper portion)
// $C: Viscous damping coefficient of the damper
// $Alpha: Viscous damper exponent
// $LGap: gap length to simulate the gap length due to the pin tolerance
// $NM:	Employed adaptive numerical algorithm (default value NM = 1; 1 = Dormand-Prince54, 2=6th order Adams-Bashforth-Moulton, 3=modified Rosenbrock Triple)
// $RelTol:	Tolerance for absolute relative error control of the adaptive iterative algorithm (default value 10^-6)
// $AbsTol:	Tolerance for absolute error control of adaptive iterative algorithm (default value 10^-6)
// $MaxHalf: Maximum number of sub-step iterations within an integration step, h=dt*(0.5)^MaxHalf (default value 15)


#include <math.h>
#include <elementAPI.h>
#include <ViscousDamper.h>
#include <Vector.h>
#include <Channel.h>
#include <OPS_Globals.h>
#include <Parameter.h>

static int numViscousDamperMaterials = 0;

void *
OPS_ViscousDamper(void)

{
  if (numViscousDamperMaterials == 0) {
    numViscousDamperMaterials++;
    opserr << "ViscousDamper Model by Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University\n";
  }
  
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  double dData[8];
  int numData = 1;
        // Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ViscousDamper tag" << endln;
    return 0;
  }
  // Check if we have 3 or 4 or 8 input variables
  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 3 && numData != 4 && numData != 8) {
    opserr << "Invalid #args, want: uniaxialMaterial ViscousDamper " << iData[0] <<  " K? C? Alpha? <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    return 0;
  }
  
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid #args want: uniaxialMaterial ViscousDamper " << iData[0] <<  " K? C? Alpha? <LGap?> <NM? RelTol? AbsTol? MaxHalf?>" << endln;
    
    return 0;   
  }
  
  if (numData == 3) {
    // Default variables
    dData[3] = 0.0;
    dData[4] = 1;
    dData[5] = 0.000001;
    dData[6] = 0.0000000001;
    dData[7] = 15;
  }
  
  if (numData == 4) {
    // Default variables
    dData[4] = 1;
    dData[5] = 0.000001;
    dData[6] = 0.0000000001;
    dData[7] = 15;
  }
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new ViscousDamper(iData[0], 
                                  dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ViscousDamper Material\n";
    return 0;
  }
  
  return theMaterial;
}


ViscousDamper::ViscousDamper(int tag, double k, double c, double a, double lgap, double nm, double reltol, double abstol, double maxhalf)
:UniaxialMaterial(tag,MAT_TAG_ViscousDamper), K(k), C(c), Alpha(a), LGap(lgap), NM(nm), RelTol(reltol), AbsTol(abstol), MaxHalf(maxhalf)    
{
    if (Alpha < 0.0) {
      opserr << "ViscousDamper::ViscousDamper -- Alpha < 0.0, setting to 1.0\n";
      Alpha = 1.0;
    }
        
    //initialize variables
    this->revertToStart();
}

ViscousDamper::ViscousDamper()
:UniaxialMaterial(0,MAT_TAG_ViscousDamper),
 K(0.0), C(0.0), Alpha(0.0), LGap(0.0),  NM(0.0), RelTol(0.0), AbsTol(0.0), MaxHalf(0.0)     
{
        this->revertToStart();
}

ViscousDamper::~ViscousDamper()
{
  // does nothing
}

int 
ViscousDamper::setTrialStrain(double strain, double strainRate)
{
  //all variables to the last commit
  this->revertToLastCommit();
  
  // Determine the strain rate and acceleration 
  
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
  fd0 = Tstress; 

  double h, yt, eps, error;
  vel0 = TVel;  // Velocity of the previous step.


  while (it < 1.0) { //iteration
    h = s * ops_Dt; // Time step 
    vel1 = vel0 + acc * h; // Velocity at the time step h
    
    
    // Selection of Numerical Method to solve the ODE
    if (NM == 1.0) {
      DormandPrince(vel0, vel1, fd0, h, yt, eps, error);
    }
    if (NM == 2.0) {
      ABM6(vel0, vel1, fd0, h, yt, eps, error);
    }
    if (NM == 3.0) {
      ROS(vel0, vel1, fd0, h, yt, eps, error);
    }
    
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
    
    // After gap inititon
    
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

double ViscousDamper::getStress(void)
{
  return  Tstress;
}

double ViscousDamper::getTangent(void)
{
  return 0;
}

double ViscousDamper::getInitialTangent(void)
{
  return 0;
}

double ViscousDamper::getDampTangent(void)
{
  
  return 0; 
}


double 
ViscousDamper::getStrain(void)
{
    return Tstrain;
}

double 
ViscousDamper::getStrainRate(void)
{
  return TVel;
}

int 
ViscousDamper::commitState(void)
{
  //commit trial  variables
  Cstrain = Tstrain;
  Cstress = Tstress;
  Ctangent = Ttangent;
  CVel = TVel;
  Cpugr = Tpugr;
  Cnugr = Tnugr; 
  
  return 0;
}

int 
ViscousDamper::revertToLastCommit(void)
{
  Tstrain = Cstrain;
  Tstress = Cstress;
  Ttangent = Ctangent;
  TVel = CVel;
  Tpugr = Cpugr;
  Tnugr = Cnugr;
  
  return 0;
}

int 
ViscousDamper::revertToStart(void)
{
  // Initialize state variables
  Tstrain=0.0;
  Tstress=0.0;
  Ttangent = 0.0;
  TVel = 0.0;
  Tpugr = 0.0;
  Tnugr = 0.0;
  
  Cstrain=0.0;
  Cstress = 0.0;
  Ctangent = 0.0;
  CVel = 0.0;
  Cpugr = 0.0;
  Cnugr = 0.0;
  
  return 0;
}

double
ViscousDamper::sgn(double dVariable){ 
  if (dVariable<0.0){
    return -1.0;
  }else{
    return 1.0;
  }
}

int 
ViscousDamper::DormandPrince(double vel0, double vel1, double y0, double h, double& yt, double& eps, double& error){
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

int 
ViscousDamper::ABM6(double vel0, double vel1, double y0, double h, double& y6, double& eps, double& error){
 double f0, f1, f2, f3, f4, f5, f6, y11, y2, y3, y4, y5, yp6;

    h = h/6.0;

    f0 = f((vel1 - vel0)*(0./6.) + vel0, y0);

    y11 = y0 + h*f0; //predictor

    f1 = f((vel1 - vel0)*(1./6.) + vel0, y11);

    y11 = y0 + h*f1; //corrector

    y2 = y11 + 0.5*h*(3.*f1 - 1.*f0); //predictor

    f2 = f((vel1 - vel0)*(2./6.) + vel0, y2);

    y2 = y11 + 0.5*h*(f2 + f1);  //corrector

    y3 = y2 + h/12.*(23.*f2 - 16.*f1 + 5.*f0); //predictor

    f3 = f((vel1 - vel0)*(3./6.) + vel0, y3);

    y3 = y2 + h/12.*(5.*f3 + 8.*f2 - 1.*f1); //corrector

    y4 = y3 + h/24.*(55.*f3 - 59.*f2 + 37.*f1 - 9.*f0); //predictor

    f4 = f((vel1 - vel0)*(4./6.) + vel0, y4);

    y4 = y3 + h/24.*(9.*f4 + 19.*f3 -5.*f2 + f1); // corrector

    y5 = y4 + h/720.*(1901.*f4 - 2774.*f3 + 2616.*f2 - 1274.*f1 + 251.*f0); // predictor

    f5 = f((vel1 - vel0)*(5./6.) + vel0, y5);

    y5 = y4 + h/720.*(251.*f5 + 646.*f4 -264.*f3 + 106.*f2 -19.*f1); // corrector

    yp6 = y5 + h/1440.*(4277.*f5 -7923.*f4 +9982.*f3 -7298.*f2 + 2877.*f1 -475.*f0); // predictor

    f6 = f((vel1 - vel0)*(6./6.) + vel0, yp6);

    y6 = y5 + (h/1440.)*(475.*f6 +1427.*f5  -798.*f4 + 482.*f3 -173.*f2 + 27.*f1); // corrector

	error = (yp6-y6);

    eps = fabs(error/y6);


 return 0;
}

int
ViscousDamper::ROS(double vel0, double vel1, double y0, double h, double& y2, double& eps, double& error){
  double k1, k2, k3, J, T, d, e32, W, f0, f1, f2, y3;
    J = -K / (Alpha*C);
    T = K;
    d = 1. / (2. + sqrt(2.));
    e32 = 6. + sqrt(2.);
    W = 1. - h*d*J;
    f0 = f(vel0, y0);
    k1 = (f0 + h*d*T)/W;
    f1 = f((vel1 - vel0)*(0.5) + vel0, y0 + (0.5)*k1*h);
    k2 = (f1 - k1)/W + k1;
    y2  = y0 + h*k2;
    f2 = f(vel1, y2);
    k3 = 1./W *(f2 - e32*(k2 - f1) - 2.*(k1 - f0) + h*d*T);
    error = h/6.*(k1 -2.*k2 + k3);
    y3 = y2 + error;
    eps = fabs(error/(y3));

 return 0;
}


double
ViscousDamper::f(double v, double fd){
    return (v - sgn(fd) * pow(fabs(fd)/C,1.0/Alpha))*K;
}

UniaxialMaterial *
ViscousDamper::getCopy(void)
{
    ViscousDamper *theCopy = new ViscousDamper(this->getTag(), K, C, Alpha, LGap, NM, RelTol, AbsTol, MaxHalf);
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
ViscousDamper::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(15);
  data(0) = this->getTag();

  // Material properties
  data(1) = K;
  data(2) = C;
  data(3) = Alpha;
  data(4) = LGap;
  data(5) = NM;
  data(6) = RelTol;
  data(7) = AbsTol;
  data(8) = MaxHalf;
  
  // State variables from last converged state
  data(9) = Cstrain;
  data(10) = Cstress;
  data(11) = Ctangent;
  data(12) = CVel;
  data(13) = Cpugr;
  data(14) = Cnugr;

        
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ViscousDamper::sendSelf() - failed to send data\n";

  return res;
}

int 
ViscousDamper::recvSelf(int cTag, Channel &theChannel, 
                               FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(15);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "ViscousDamper::recvSelf() - failed to receive data\n";
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
        
        // Material properties
        K = data(1);
        C = data(2);
		Alpha = data(3);
		LGap = data(4);
	    NM = data(5);
        RelTol = data(6);
		AbsTol = data(7);
		MaxHalf = data(8);
        
        // State variables from last converged state 
        Cstrain = data(9);
        Cstress = data(10);
        Ctangent = data(11);
        CVel = data(12);
		Cpugr = data(13);
		Cnugr = data(14);
          
  }
    
  return res;
}

int
ViscousDamper::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"E") == 0 || strcmp(argv[0],"K") == 0) {
    param.setValue(K);
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"eta") == 0 || strcmp(argv[0],"C") == 0) {
    param.setValue(C);
    return param.addObject(4, this);
  }
  return -1;
}

int 
ViscousDamper::updateParameter(int parameterID, Information &info)
{
  switch(parameterID) {
  case 1:
    K = info.theDouble;
    return 0;
  case 4:
    C = info.theDouble;
    return 0;
  default:
    return -1;
  }
}

void 
ViscousDamper::Print(OPS_Stream &s, int flag)
{
    s << "ViscousDamper tag: " << this->getTag() << endln;
    s << "  K: " << K << endln; 
    s << "  C: " << C << endln;
    s << "  Alpha: " << Alpha << endln;
	s << "  LGap: " << LGap << endln; 
	s << "  NM: " << NM << endln; 
    s << "  RelTol: " << RelTol << endln;
	s << "  AbsTol: " << AbsTol << endln;
    s << "  MaxHalf: " << MaxHalf << endln;
        
}
