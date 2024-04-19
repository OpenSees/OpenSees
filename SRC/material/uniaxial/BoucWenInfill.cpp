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

// Written by Stefano Sirotti (stefano.sirotti@unimore.it)
// Created on January 2022
//
// Description: 
// This file contains the class definition for 
// BoucWenInfill material. BoucWenInfill material provides 
// a hysteretic uniaxial material described by the Bouc-Wen law 
// with stiffness, strength degradation and pinching. 
// Particularly suitable to simulate hysteresis of infill panels.
//
// Reference: 
// Sirotti, S., Pelliciari, M., Di Trapani, F.,
// Briseghella, B., Carlo Marano, G., Nuti, C., & Tarantino, A. M. (2021).
// Development and validation of new Bouc–Wen data-driven hysteresis model 
// for masonry infilled RC frames. 
// Journal of Engineering Mechanics, 147(11), 04021092.
//

#include <elementAPI.h>
#include <BoucWenInfill.h>

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <Parameter.h>


void *
OPS_BoucWenInfill(void)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 16) {
	opserr << "WARNING: Insufficient arguments\n";
	opserr << "Want: uniaxialMaterial BoucWenInfill tag? mass? alpha? beta0? eta0?" << endln
		<< "n? k? xy? deltak? deltaf? psi? Zs? As? epsp? tol? maxNumIter?" << endln;
	return 0;
    }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData1[1];
  double dData[14];
  int    iData2[1];
  int numData;

  numData = 1;
  if (OPS_GetIntInput(&numData, iData1) != 0) {
    opserr << "WARNING invalid uniaxialMaterial BoucWenInfill tag" << endln;
    return 0;
  }

  numData = 14;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "WARNING invalid Double Values\n";
    return 0;	
  }

  numData = 1;
  if (OPS_GetIntInput(&numData, iData2) != 0) {
    opserr << "WARNING invalid maxNumIter" << endln;
    return 0;
  }

  theMaterial = new BoucWenInfill(iData1[0], dData[0], dData[1], dData[2],
  dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], dData[9], 
  dData[10], dData[11], dData[12], dData[13], iData2[0]);       

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type BoucWenInfill\n";
    return 0;
  }

  return theMaterial;
}


BoucWenInfill::BoucWenInfill(int tag,
	double p_mass,
	double p_alpha,
	double p_beta0,
	double p_eta0,
	double p_n,
	double p_k,
	double p_xy,
	double p_deltak,
	double p_deltaf,
	double p_psi,
	double p_Zs,
	double p_As,
	double p_epsp,
	double ptolerance,
	int pMaxNumIter)
 :UniaxialMaterial(tag, MAT_TAG_BoucWenInfill),
  xmax(0.0), xmaxp(0.0), mass(p_mass), alpha(p_alpha), beta0(p_beta0), eta0(p_eta0), n(p_n),
  k(p_k), xy(p_xy), deltak(p_deltak), deltaf(p_deltaf), psi(p_psi), 
  Zs(p_Zs), As(p_As), epsp(p_epsp), tolerance(ptolerance), maxNumIter(pMaxNumIter)
{
  // Initialize variables
  this->revertToStart();
}

BoucWenInfill::BoucWenInfill()
 :UniaxialMaterial(0, MAT_TAG_BoucWenInfill),
  xmax(0.0), xmaxp(0.0), mass(0.0), alpha(0.0), beta0(0.0), eta0(0.0), n(0.0),
  k(0.0), xy(0.0), deltak(0.0), deltaf(0.0), psi(0.0), 
  Zs(0.0), As(0.0), epsp(0.0), tolerance(0.0), maxNumIter(0)
{

}


BoucWenInfill::~BoucWenInfill()
{
  // does nothing
}

double 
BoucWenInfill::signum(double value)
{
  if (value > 0.0) {
    return 1.0;
  }
  else {
    return -1.0;
  }
}


int 
BoucWenInfill::setTrialStrain (double strain, double strainRate)
{
	// Set trial strain and compute strain increment
	Tstrain = strain;
	double dStrain = Tstrain - Cstrain;

	// Initial declarations (make sure not to declare class variables here!)
	double TDamIndex, Tpk, TA, Tbetak, Tbetaf, TPow, Psi, Phi, Tznew, Tzold, F;
	double f_z, geps_z, Te_z, TDamIndex_z, Tpk_z, TA_z, Tbetak_z, Tbetaf_z, TPow_z, Phi_z, F_z;
	double geps_x, a_x, Te_x, TDamIndex_x, Tpk_x, TA_x, Tbetak_x, Tbetaf_x, Phi_x, F_x;
	double geps_z_x, Te_z_x, TDamIndex_z_x, Tpk_z_x, TA_z_x, Tbetak_z_x, Tbetaf_z_x, Phi_z_x, F_z_x;
	double f, a, geps;
	double numerF_z, denomF_z;
	double numerF_z_x, denomF_z_x;

	// Newton-Raphson scheme to solve for z_{i+1} = z1
	int count = 0;
	double startPoint = 0.01;
	Tz = startPoint;
	Tzold = startPoint;
	Tznew = 1.0;  
	while ( ( fabs(Tzold-Tznew) > tolerance ) && count < maxNumIter) {

 		// Trial energy and updating xmax
		Te = Ce + (1.0-alpha)*k/mass*dStrain*Tz;
		if ((fabs(Tstrain)) >= (fabs(xmax))) {
			xmax = fabs(Tstrain);
		}
		else {
			xmax = xmaxp;
		}

		//Evaluate Slip-Lock functions
		f = exp(-pow(Tz/Zs, 2));
		a = As*fabs(xmax)/xy;
		geps = 1.0 - exp((-0.5)*pow(Te/epsp,2));

		// Damage index and degrading functions 
		TDamIndex = Te*mass/(k*pow(xy,2)) + fabs(xmax)/xy;
		Tpk = exp(-psi*TDamIndex);
		TA = exp(-deltak*TDamIndex*Tpk);
		Tbetak = beta0*TA;
		Tbetaf = exp(n*deltaf*TDamIndex);
		Psi = (eta0 + signum(dStrain*Tz));
		TPow = pow(fabs(Tz),(n));
		Phi = TA - TPow*Tbetak*Tbetaf*Psi;

		// Root function F
		F = Tz - Cz - dStrain*Phi/(1 + a*f*geps*Phi);

		// Evaluate function derivative F_z ( _z = derivative with respect to z_n+1 )
		Te_z = (1.0-alpha)*k/mass*dStrain;
		f_z = -2*Tz/pow(Zs,2)*f;
	    geps_z = exp((-0.5)*pow(Te/epsp,2))*Te*Te_z/pow(epsp,2);
		TDamIndex_z = Te_z*mass/(k*pow(xy,2));
		Tpk_z = Tpk*(-psi*TDamIndex_z);
		TA_z = TA*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z);
		Tbetak_z = beta0*TA_z;
		Tbetaf_z = Tbetaf*(n*deltaf*TDamIndex_z);
		if (Tz == 0.0) {
			TPow = 0.0;
			TPow_z = 0.0;
		}
		else {
			TPow = pow(fabs(Tz),(n));
			TPow_z = n*pow(fabs(Tz),(n-1))*signum(Tz);
		}
		Phi_z = TA_z - (TPow_z*Tbetak*Tbetaf + TPow*Tbetak_z*Tbetaf + TPow*Tbetak*Tbetaf_z)*Psi;
		F_z = 1.0 - dStrain*(Phi_z*(1 + a*f*geps*Phi) - Phi*(a*f_z*geps*Phi + a*f*geps_z*Phi + a*f*geps*Phi_z))/pow(1+a*f*geps*Phi,2);
		numerF_z = Phi_z - a*f_z*geps*pow(Phi,2) - a*f*geps_z*pow(Phi,2);
		denomF_z = pow(1 + a*f*geps*Phi,2);
		
		// Issue warning if derivative is zero
		if ( fabs(F_z)<1.0e-10 ) {
			opserr << "WARNING: BoucWenInfill::setTrialStrain() -- zero derivative " << endln
				<< " in Newton-Raphson scheme" << endln;
		}

		// Take a Newton step
		Tznew = Tz - F/F_z;


		// Update the root (but keep the old for convergence check)
		Tzold = Tz; 
		Tz = Tznew;

		// Update counter
		count++;


		// Issue warning if we didn't converge
		if (count == maxNumIter) {

			opserr << "WARNING: BoucWenInfill::setTrialStrain() -- did not" << endln
				<< " find the root z_{i+1}, after " << maxNumIter << " iterations" << endln
				<< " and norm: " << fabs(Tzold-Tznew) << endln;
		}

		// Compute stress
		
		Tstress = alpha*k*Tstrain + (1-alpha)*k*Tz;


		// Compute deterioration parameters
	
		Te = Ce + (1-alpha)*k/mass*dStrain*Tz;
		TDamIndex = Te*mass/(k*pow(xy,2)) + fabs(xmax)/xy;
		Tpk = exp(-psi*TDamIndex);
		TA = exp(-deltak*TDamIndex*Tpk);
		Tbetak = beta0*TA;
		Tbetaf = exp(n*deltaf*TDamIndex);
		if (Tz == 0.0) {
			TPow = 0.0;
			TPow_z = 0.0;
		}
		else {
			TPow = pow(fabs(Tz),(n));
			TPow_z = n*pow(fabs(Tz),(n-1))*signum(Tz);
		}
        Psi = eta0 + signum(dStrain*Tz);
		Phi =  TA - TPow*Tbetak*Tbetaf*Psi;
		
		// Compute tangent ( _x = derivative with respect to x_n+1 )
		
		// Compute the derivative of F with respect to x_n+1
		if (Tz != 0.0) {
			Te_x = (1-alpha)*k/mass*Tz;
			if (xmax == Tstrain) {
				TDamIndex_x = Te_x*mass/(k*pow(xy,2)) + 1/xy;
				a_x = As/xy;
			}
			else {
				TDamIndex_x = Te_x*mass/(k*pow(xy,2));
				a_x = 0;
			}
			geps_x = exp((-0.5)*pow(Te/epsp,2))*Te*Te_x/pow(epsp,2);
			Tpk_x = Tpk*(-psi*TDamIndex_x);
			TA_x = TA*(-deltak*TDamIndex_x*Tpk - deltak*TDamIndex*Tpk_x);
			Tbetak_x = beta0*TA_x;
			Tbetaf_x = Tbetaf*(n*deltaf*TDamIndex_x);
			Phi_x = TA_x - (Tbetak_x*Tbetaf + Tbetak*Tbetaf_x)*TPow*Psi;
			F_x = - Phi/(1 + a*f*geps*Phi) - dStrain*(Phi_x*(1 + a*f*geps*Phi) - Phi*(a_x*f*geps*Phi + a*f*geps_x*Phi + a*f*geps*Phi_x))/pow(1 + a*f*geps*Phi,2);
			
		// Compute the derivative of F_z with respect to x_n+1
			Te_z_x = (1-alpha)*k/mass;
			geps_z_x = - geps_x*Te*Te_z/pow(epsp,2) + exp((-0.5)*pow(Te/epsp,2))*Te_x*Te_z/pow(epsp,2) + exp((-0.5)*pow(Te/epsp,2))*Te*Te_z_x/pow(epsp,2);
			TDamIndex_z_x = Te_z_x*mass/(k*pow(xy,2));
			Tpk_z_x = Tpk_x*(-psi*TDamIndex_z) + Tpk*(-psi*TDamIndex_z_x);
			TA_z_x = TA_x*(-deltak*TDamIndex_z*Tpk - deltak*TDamIndex*Tpk_z) - TA*(deltak*TDamIndex_z_x*Tpk + deltak*TDamIndex_z*Tpk_x + deltak*TDamIndex_x*Tpk_z + deltak*TDamIndex*Tpk_z_x);
			Tbetak_z_x = beta0*TA_z_x;
			Tbetaf_z_x = Tbetaf_x*(n*deltaf*TDamIndex_z) + Tbetaf*(n*deltaf*TDamIndex_z_x);
			Phi_z_x = TA_z_x - (TPow_z*Tbetak*Tbetaf_x + TPow_z*Tbetak_x*Tbetaf + TPow*Tbetaf_z_x*Tbetak + TPow*Tbetaf_z*Tbetak_x + TPow*Tbetaf_x*Tbetak_z + TPow*Tbetaf*Tbetak_z_x)*Psi;
			numerF_z_x = Phi_z_x - (a_x*f_z*geps*pow(Phi,2) + a*f_z*geps_x*pow(Phi,2) + 2*a*f_z*geps*Phi*Phi_x) - (a_x*f*geps_z*pow(Phi,2) + a*f*geps_z_x*pow(Phi,2) + 2*a*f*geps_z*Phi*Phi_x);
			denomF_z_x = 2*(1 + a*f*geps*Phi)*(a_x*f*geps*Phi + a*f*geps_x*Phi + a*f*geps*Phi_x);
			F_z_x = - numerF_z/denomF_z - dStrain*(numerF_z_x*denomF_z - numerF_z*denomF_z_x)/pow(denomF_z,2);
		
		// Compute tangent	
			double DzDx = -(F_x*F_z - F*F_z_x)/(pow(F_z,2));
			Ttangent = alpha*k + (1-alpha)*k*DzDx;
			
		}
		else {
			Ttangent = alpha*k + (1-alpha)*k;
		}
		
	

	}
	
	
    return 0;
}

double 
BoucWenInfill::getStress(void)
{
    return Tstress;
}

double 
BoucWenInfill::getInitialTangent(void)
{
    return (alpha*k + (1-alpha)*k);
}


double 
BoucWenInfill::getTangent(void)
{
    return Ttangent;
}

double 
BoucWenInfill::getStrain(void)
{
    return Tstrain;
}

int 
BoucWenInfill::commitState(void)
{
    // Commit trial history variables
	xmaxp = xmax;
	Cstrain = Tstrain;
	Cz = Tz;
	Ce = Te;

    return 0;
}

int 
BoucWenInfill::revertToLastCommit(void)
{
	xmax = xmaxp;
    return 0;
}

int 
BoucWenInfill::revertToStart(void)
{
    xmaxp = 0.0;
	Tstrain = 0.0;
	Cstrain = 0.0;
	Tz = 0.0;
	Cz = 0.0;
	Te = 0.0;
	Ce = 0.0;
	Tstress = 0.0;
	Ttangent = alpha*k + (1-alpha)*k;

    return 0;
}

UniaxialMaterial *
BoucWenInfill::getCopy(void)
{
	BoucWenInfill*theCopy =
	new BoucWenInfill(this->getTag(), mass, alpha, beta0, eta0, n, k, xy,
	deltak, deltaf, psi, Zs, As, epsp, tolerance, maxNumIter);
    theCopy->xmax = xmax;
	theCopy->xmaxp = xmaxp;
    theCopy->Tstrain = Tstrain;
    theCopy->Cstrain = Cstrain;
    theCopy->Tz = Tz;
    theCopy->Cz = Cz;
    theCopy->Te = Te;
    theCopy->Ce = Ce;
    theCopy->Tstress = Tstress;
	theCopy->Ttangent = Ttangent;

    return theCopy;
}

int 
BoucWenInfill::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(1+15+4);
  data(0) = this->getTag();
  
  data(1) = mass;
  data(2) = alpha;
  data(3) = beta0;
  data(4) = eta0;
  data(5) = n;
  data(6) = k;
  data(7) = xy;
  data(8) = deltak;
  data(9) = deltaf;
  data(10) = psi;
  data(11) = Zs;
  data(12) = As;
  data(13) = epsp;
  data(14) = tolerance;
  data(15) = maxNumIter;

  data(16) = xmaxp;
  data(17) = Cstrain;
  data(18) = Cz;
  data(19) = Ce;
  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "BoucWenInfill::sendSelf() - failed to send data" << endln;

  return res;
}

int 
BoucWenInfill::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(20);
  
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) {
    opserr << "BoucWenInfill::recvSelf() - failed to receive Vector" << endln;
    return res;
  }
  
  this->setTag(int(data(0)));

  mass = data(1);
  alpha = data(2);
  beta0 = data(3);
  eta0 = data(4);
  n = data(5);
  k = data(6);
  xy = data(7);
  deltak = data(8);
  deltaf = data(9);
  psi = data(10);
  Zs = data(11);
  As = data(12);
  epsp = data(13);
  tolerance = data(14);
  maxNumIter = (int)data(15);

  xmaxp = data(16);
  Cstrain = data(17);
  Cz = data(18);
  Ce = data(19);
  
  xmax = xmaxp;
  Tstrain = Cstrain;
  Tz = Cz;
  Te = Ce;
  
  this->revertToLastCommit();
    
  return res;
}

void 
BoucWenInfill::Print(OPS_Stream &s, int flag)
{
    s << "BoucWenInfill, tag: " << this->getTag() << endln;
	s << "  mass: " << mass << endln;
	s << "  alpha: " << alpha << endln;
	s << "  beta0: " << beta0 << endln;
    s << "  eta0: " << eta0 << endln;
    s << "  n: " << n << endln;
    s << "  k: " << k << endln;
    s << "  xy: " << xy << endln;
   	s << "  deltak: " << deltak << endln;
	s << "  deltaf: " << deltaf << endln;
	s << "  psi: " << psi << endln;
	s << "  Zs: " << Zs << endln;
	s << "  As: " << As << endln;
	s << "  epsp: " << epsp << endln;
}
