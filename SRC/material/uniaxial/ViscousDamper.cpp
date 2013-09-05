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
                                                                        
// $Revision: 1.7 $
// $Date: 2009/03/05 00:52:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ViscousDamper.cpp,v $
                                                                        
// Written: Sarven Akcelyan and Dimitrios G. Lignos, PhD, McGill University
// Created: January 2013
// Updated: September 2013
// Revision: A
//
// Description: This file contains the class interface for 
// Viscous Damper Model Relationship of the form F = K*u + C*pow(V,a)
// Reference: Kasai K, Oohara K. (2001). “Algorithm and Computer Code To Simulate Response of Nonlinear Viscous Damper”. 
// Proceedings Passively Controlled Structure Symposium 2001, Yokohama, Japan.
// Variables:
// K: axial stiffness of a damper
// C: Velocity constant of a damper
// Alpha: Exponent of velocity of a damper

#include <math.h>

#include <elementAPI.h>
#include <ViscousDamper.h>
#include <Vector.h>
#include <Channel.h>

#include <OPS_Globals.h>

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
  double dData[3];
  int numData = 1;
	// Check tag
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial  ViscousDamper tag" << endln;
    return 0;
  }
  // Check if we have 3 input variables for K, C, Alpha
  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid Args want: uniaxialMaterial ViscousDamper tag? K? C? Alpha? "<< endln;
    
    return 0;	
  }
  
  // Parsing was successful, allocate the material with zero index
  theMaterial = new ViscousDamper(iData[0], 
				  dData[0], dData[1], dData[2]);
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type ViscousDamper Material\n";
    return 0;
  }
  
  return theMaterial;
}


ViscousDamper::ViscousDamper(int tag, double k, double c, double a)
:UniaxialMaterial(tag,MAT_TAG_ViscousDamper), K(k), C(c), Alpha(a)
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
 K(0.0), C(0.0), Alpha(0.0)
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
	
	// Determine Incremental strain
	double dStrain = (strain - Cstrain);

	// Determine the incremental strain rate from incremental strain
	
	double dVel = dStrain/ops_Dt;
    
	// Determine the average velocity from past and current step
	
	double dVelAv = 0.5*(dVel + CdVel);
	
	double a0 = Tfd/C;
	// 4th Order Runge-Kutta with Classical Coefficients
	
	double k0 = (TdVel - sgn(a0) * pow(fabs(a0),1.0/Alpha)) * K * ops_Dt;
	
	double a1 = (Tfd + k0 * 0.5)/C;
	
	double k1 = (dVelAv - sgn(a1) * pow(fabs(a1),1.0/Alpha)) * K * ops_Dt;
	
	double a2 = (Tfd + k1 * 0.5)/C;
	
	double k2 = (dVelAv - sgn(a2) * pow(fabs(a2),1.0/Alpha)) * K * ops_Dt;
	
	double a3 = (Tfd + k2)/C;
	
	double k3 = (dVel - sgn(a3) * pow(fabs(a3),1.0/Alpha)) * K * ops_Dt;

	
	// Total Stress computation based on 4-th Order R-K
	
	Tfd = Tfd + (1.0/6.0) * (k0 + 2.0 * k1 + 2.0 * k2 + k3);
	
	// Total strain in the elastic part of the damper 
	double Tstrains = Tfd/K;
	
	// Total strain in the viscous part of the damper
	//double Tstraind = Tstrain - Tstrains;
	double Tstraind =  strain - Tstrains;
	// Total Stress of the damper
	
	Tstress = Tfd;
	Tstrain = strain;
	TdVel = dVel;
	Ttangent = K;
	
    return 0;
}

double ViscousDamper::getStress(void)
{

	return  Tstress;
}

double ViscousDamper::getTangent(void)
{
	return 0.0;
}

double ViscousDamper::getInitialTangent(void)
{
    return K;
}

double ViscousDamper::getDampTangent(void)
{
	//double DTangent = Tstress/Tstrain;
    return 0.0;	
}


double 
ViscousDamper::getStrain(void)
{
    return Tstrain;
}

double 
ViscousDamper::getStrainRate(void)
{
    return 0;
}

int 
ViscousDamper::commitState(void)
{
	//commit trial  variables
    Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;
	CdVel = TdVel;
	Cfd = Tfd;
	
	return 0;
}

int 
ViscousDamper::revertToLastCommit(void)
{
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
	TdVel = CdVel;
	Tfd = Cfd;
	
    return 0;
}

int 
ViscousDamper::revertToStart(void)
{
    // Initialize state variables
	Tstrain=0.0;
	Tstress=0.0;
	Ttangent = K;
	TdVel = 0.0;
	Tfd = 0.0;
	
	Cstrain=0.0;
	Cstress = 0.0;
	Ctangent = K;
	CdVel = 0.0;
	Cfd = 0.0;
	
    return 0;
}

int
ViscousDamper::sgn(double dVariable){ 
    if (dVariable<-1e-5){
		return -1;
	}else{
		return 1;
	}
}

UniaxialMaterial *
ViscousDamper::getCopy(void)
{
    ViscousDamper *theCopy = new ViscousDamper(this->getTag(), K, C, Alpha);
    // Converged state variables
	theCopy->Cstrain = Cstrain;
	theCopy->Cstress = Cstress;
	theCopy->Ctangent = Ctangent;
	theCopy->CdVel = CdVel;
	theCopy->Cfd = Cfd;
	
	// Trial state variables
    theCopy->Tstrain = Tstrain;
    theCopy->Tstress = Tstress;	
	theCopy->Ttangent = Ttangent;
	theCopy->TdVel = TdVel;
	theCopy->Tfd = Tfd;
	
    return theCopy;
}

int 
ViscousDamper::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  static Vector data(9);
  data(0) = this->getTag();

  // Material properties
  data(1) = K;
  data(2) = C;
  data(3) = Alpha;
  
  // State variables from last converged state
  data(4) = Cstrain;
  data(5) = Cstress;
  data(6) = Ctangent;
  data(7) = CdVel;
  data(8) = Cfd;
	
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
  static Vector data(9);
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
	
	// State variables from last converged state 
	Cstrain = data(4);
	Cstress = data(5);
	Ctangent = data(6);
	CdVel = data(7);
	Cfd = data(8);
	  
	//Copy converged state values into trial values
	//Tstrain = Cstrain;
	//Tstress = Cstress;
	//Ttangent = Ctangent;
	//TdVel = CdVel;
	//Tfd = Cfd;
	  
  }
    
  return res;
}

void 
ViscousDamper::Print(OPS_Stream &s, int flag)
{
    s << "ViscousDamper tag: " << this->getTag() << endln;
    s << "  K: " << K << endln;	
    s << "  C: " << C << endln;
    s << "  Alpha: " << Alpha << endln;
	
}


