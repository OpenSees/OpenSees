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
// $Date: 2006-01-17 20:44:32 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/Isolator2spring.cpp,v $

// Written: K. Ryan
// Created: September 2003
// Updates: November 2005
//
// Description: This file contains the class implementation for a "two-spring isolator" 
// material.  This material is based on the two-spring model originally developed by 
// Koh and Kelly to represent the buckling behavior of an elastomeric bearing.  The 
// material model has been modified to include material nonlinearity and optional 
// strength degradation.

#include <Isolator2spring.h>           
#include <math.h>

#include <Channel.h>
#include <Matrix.h>
#include <Vector.h>
#include <elementAPI.h>

void* OPS_Isolator2spring()
{
    if (OPS_GetNumRemainingInputArgs() < 8) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section Iso2spring tag? tol? k1? Fy? k2? kv? hb? Pe? <Po?>" << endln;
	return 0;
    }    
	  
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid Iso2spring tag" << endln;
	return 0;
    }

    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 8) numdata = 8;
    double data[8] = {0,0,0,0,0,0,0,0};
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double inputs\n";
	opserr << "section Iso2spring: " << tag << endln;
	return 0;
    }
    double tol = data[0];
    double k1 = data[1];
    double Fy = data[2];
    double kb = data[3];
    double kvo = data[4];
    double hb = data[5];
    double Pe = data[6];
    double Po = data[7];

    return new Isolator2spring(tag, tol, k1, Fy, kb, kvo, hb, Pe, Po);
}

Vector Isolator2spring::s(2);
Vector Isolator2spring::s3(3);
Vector Isolator2spring::f0(5);
Matrix Isolator2spring::df(5,5);
ID Isolator2spring::code(3);

Isolator2spring::Isolator2spring
(int tag, double tol_in, double k1_in, double Fy_in, double kb_in, double kvo_in, double hb_in, double Pe_in, 
 double po_in) :
 SectionForceDeformation(tag, SEC_TAG_Isolator2spring),
 tol(tol_in), k1(k1_in), Fyo(Fy_in), kbo(kb_in), kvo(kvo_in), h(hb_in), Pe(Pe_in), po(po_in), x0(5), ks(3,3)
{
  this->revertToStart();
  
  pcr = sqrt(Pe*kbo*h);
  H = k1*kbo/(k1 - kbo);
  
  code(0) = SECTION_RESPONSE_P;
  code(1) = SECTION_RESPONSE_VY;
  code(2) = SECTION_RESPONSE_MZ;
  
}

Isolator2spring::Isolator2spring():
 SectionForceDeformation(0, SEC_TAG_Isolator2spring),
 tol(1.0e-12), k1(0.0), Fyo(0.0), kbo(0.0), kvo(0.0), h(0.0), Pe(0.0), po(0.0),
 x0(5), ks(3,3)
{

        this->revertToStart();

	//pcr = sqrt(Pe*kbo*h);
	//H = k1*kbo/(k1 - kbo);

	code(0) = SECTION_RESPONSE_P;
	code(1) = SECTION_RESPONSE_VY;
	code(2) = SECTION_RESPONSE_MZ;
}

Isolator2spring::~Isolator2spring()
{
	// Nothing to do here
}

int
Isolator2spring::setTrialSectionDeformation(const Vector &e)
{
  utpt[0] = e(0);
  utpt[1] = e(1);
  return 0;
}

const Matrix&
Isolator2spring::getSectionTangent(void)
{
  
  // ks is computed in getStressResultant 
  return ks;
}

const Matrix&
Isolator2spring::getInitialTangent(void)
{
  // Initial tangent uses nominal properties of the isolator.  
  ks(0,0) = k1;
  ks(1,1) = kvo;
  ks(0,1) = ks(1,0) = 0.0;
  ks(0,2) = ks(1,2) = ks(2,2) = ks(2,1) = ks(2,0) = 0.0;
  
	return ks;
}

const Vector&
Isolator2spring::getStressResultant(void)
{

  double Fy;
  if (po < 1.0e-10) {
    // No strength degradation
    Fy = Fyo;
  } else {
    // Strength degradation based on bearing axial load
    double p2 = x0(1)/po;
    if (p2<0) {
      p2 = 0.0;
    }
    Fy = Fyo*(1-exp(-p2));
  }
  
  
  // Material stresses using rate independent plasticity, return mapping algorithm
  
  // Compute trial stress using elastic tangent
  double fb_try = k1*(x0(2)-sP_n);
  double xi_try = fb_try - q_n;
  
  // Yield function
  double Phi_try = fabs(xi_try) - Fy;
  
  double fspr;
  double dfsds;
  double del_gam;
  int sign;
  
  // Elastic step
  if (Phi_try <= 0.0) {
    // Stress and tangent, update plastic deformation and back stress
    fspr = fb_try;
    dfsds = k1;
    sP_n1 = sP_n;
    q_n1 = q_n;
  }
  
  // Plastic step
  else {
    // Consistency parameter
    del_gam = Phi_try/(k1+H);
    
    sign = (xi_try < 0) ? -1 : 1;
    // Return stress to yield surface
    fspr = fb_try - del_gam*k1*sign;
    dfsds = kbo;
    // Update plastic deformation and back stress
    sP_n1 = sP_n + del_gam*sign;
    q_n1 = q_n + del_gam*H*sign;
  }
  
  // Nonlinear equilibrium and kinematic equations; want to find the 
  // zeros of these equations.
  f0(0) = x0(0) - fspr + x0(1)*x0(3);
  f0(1) = x0(0)*h - Pe*h*x0(3) + x0(1)*(x0(2)+h*x0(3));
  f0(2) = x0(1) - kvo*x0(4);
  f0(3) = utpt[0] - x0(2) - h*x0(3);
  f0(4) = -utpt[1] - x0(2)*x0(3) - h/2.0*x0(3)*x0(3) - x0(4);
  
  int iter = 0;
  double normf0 = f0.Norm();
  static Matrix dfinverse(5,5);
  
  // Solve nonlinear equations using Newton's method
  while (normf0 > tol) {
    
    iter += 1;
    
    // Formulate Jacobian of nonlinear equations
    df(0,0) = 1.0;
    df(0,1) = x0(3);
    df(0,2) = -dfsds;
    df(0,3) = x0(1);
    df(0,4) = 0.0;
    
    df(1,0) = h;
    df(1,1) = x0(2) + h*x0(3);
    df(1,2) = x0(1);
    df(1,3) = (x0(1) - Pe)*h;
    df(1,4) = 0.0;
    
    df(2,0) = 0.0;
    df(2,1) = 1.0;
    df(2,2) = 0.0;
    df(2,3) = 0.0;
    df(2,4) = -kvo;
    
    df(3,0) = 0.0;
    df(3,1) = 0.0;
    df(3,2) = -1.0;
    df(3,3) = -h;
    df(3,4) = 0.0;
    
    df(4,0) = 0.0;
    df(4,1) = 0.0;
    df(4,2) = -x0(3);
    df(4,3) = -(x0(2) + h*x0(3));
    df(4,4) = -1.0;
    
    df.Invert(dfinverse);
    // Compute improved estimate of solution x0
    x0 -= dfinverse*f0;
    
    if (po > 1.0e-10) { // Update strength according to axial load
      double p2 = x0(1)/po;
      if (p2<0) {
	p2 = 0.0;
      }
      Fy = Fyo*(1-exp(-p2));
    }
    
    // Apply plasticity theory again, return mapping algorithm 
    fb_try = k1*(x0(2) - sP_n);
    xi_try = fb_try - q_n;
    
    Phi_try = fabs(xi_try) - Fy;
    // Elastic step
    if (Phi_try <= 0.0) {
      fspr = fb_try;
      dfsds = k1;
      sP_n1 = sP_n;
      q_n1 = q_n;
    }
    
    // Plastic step
    else {
      del_gam = Phi_try/(k1+H);
      sign = (xi_try < 0) ? -1 : 1;
      fspr = fb_try - del_gam*k1*sign;
      dfsds = kbo;
      sP_n1 = sP_n + del_gam*sign;
      q_n1 = q_n + del_gam*H*sign;
    }
    
    // Estimate the residual
    f0(0) = x0(0) - fspr + x0(1)*x0(3);
    f0(1) = x0(0)*h - Pe*h*x0(3) + x0(1)*(x0(2)+h*x0(3));
    f0(2) = x0(1) - kvo*x0(4);
    f0(3) = utpt[0] - x0(2) - h*x0(3);
    f0(4) = -utpt[1] - x0(2)*x0(3) - h/2.0*x0(3)*x0(3) - x0(4);
    
    normf0 = f0.Norm();
    
    if (iter > 19) {
      opserr << "WARNING! Iso2spring: Newton iteration failed. Norm Resid: " << normf0  << endln;
      break;
    }
  }
  
  // Compute stiffness matrix by three step process
  double denom = h*dfsds*(Pe - x0(1)) - x0(1)*x0(1);
  static Matrix fkin(3,2);
  fkin(0,0) = 1.0;
  fkin(1,0) = h;
  fkin(2,0) = 0.0;
  fkin(0,1) = -x0(3);
  fkin(1,1) = -(x0(2) + h*x0(3));
  fkin(2,1) = -1.0;
  
  static Matrix feq(3,3);
  feq(0,0) = (Pe-x0(1))*h/denom;
  feq(0,1) = feq(1,0) = x0(1)/denom;
  feq(1,1) = dfsds/denom;
  feq(0,2) = feq(1,2) = feq(2,0) = feq(2,1) = 0.0;
  feq(2,2) = 1.0/kvo;
  
  static Matrix ftot(2,2);
  static Matrix ktot(2,2);
  ftot.Zero();
  ftot.addMatrixTripleProduct(0.0,fkin,feq,1.0);
  ftot.Invert(ktot);
  
  ks(0,0) = ktot(0,0);
  ks(1,0) = ktot(1,0);
  ks(0,1) = ktot(0,1);
  ks(1,1) = ktot(1,1);
  ks(0,2) = ks(1,2) = ks(2,2) = ks(2,1) = ks(2,0) = 0.0;
  
  
  // Compute force vector
  s3(0) = x0(0);
  s3(1) = -x0(1);
  s3(2) = (x0(1)*utpt[0] + x0(0)*h)/2.0;
  return s3;
}

const Vector&
Isolator2spring::getSectionDeformation(void)
{
  // Write to static variable for return
  s(0) = utpt[0];
  s(1) = utpt[1];
  
  return s;
}

int
Isolator2spring::commitState(void)
{
        sP_n = sP_n1;
      	q_n = q_n1;

	return 0;
}

int
Isolator2spring::revertToLastCommit(void)
{
	return 0;
}

int
Isolator2spring::revertToStart(void)
{
	for (int i = 0; i < 2; i++) {
		utpt[i]  = 0.0;
	}

	sP_n = 0.0;
	sP_n1 = 0.0;
	q_n  = 0.0;
	q_n1 = 0.0;

	x0.Zero(); 
	ks.Zero();

	return 0;
}

SectionForceDeformation*
Isolator2spring::getCopy(void)
{
  	Isolator2spring *theCopy =
		new Isolator2spring (this->getTag(), tol, k1, Fyo, kbo, kvo, h, Pe, po);
 
        for (int i = 0; i < 2; i++) {
	 	theCopy->utpt[i]  = utpt[i];
	}
	
	theCopy->sP_n = sP_n; 
	theCopy->sP_n1  = sP_n1; 
	theCopy->q_n = q_n;
	theCopy->q_n1 = q_n1;
	theCopy->pcr = pcr;
	theCopy->H = H;

	theCopy->x0 = x0;
	theCopy->ks = ks;
        
	return theCopy;  
}

const ID&
Isolator2spring::getType(void)
{
	return code;
}

int
Isolator2spring::getOrder(void) const
{
	return 3;
}

int 
Isolator2spring::sendSelf(int cTag, Channel &theChannel)
{
        int res = 0;
	
	static Vector data(13);
	    
	data(0) = this->getTag();
	data(1) = tol;
	data(2) = k1;
	data(3) = Fyo;
	data(4) = kbo;
	data(5) = kvo;
	data(6) = h;
	data(7) = Pe;
	data(8) = po;
	data(9) = sP_n;
	data(10) = q_n;
	data(11) = H;
	data(12) = pcr;
	
	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) 
	  opserr << "Isolator2spring::sendSelf() - failed to send data\n";
	return res;
}

int 
Isolator2spring::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
        int res = 0;
 
	static Vector data(13);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);
	
	if (res < 0) {
	  opserr << "Isolator2spring::recvSelf() - failed to receive data\n";
	  this->setTag(0);      
	}
	else {
	  this->setTag((int)data(0));
	  tol = data(1);
	  k1 = data(2);
	  Fyo = data(3);
	  kbo = data(4);
	  kvo = data(5);
	  h = data(6);
	  Pe = data(7);
	  po  = data(8);
	  sP_n  = data(9);
	  q_n = data(10);
	  H = data(11);
	  pcr = data(12);
	  
	  // Set the trial state variables
	  revertToLastCommit();
	}
	
	return res;
}

void
Isolator2spring::Print(OPS_Stream &s, int flag)
{

	s << "Isolator2spring, tag: " << this->getTag() << endln;
	s << "\tol:    " << tol << endln;
	s << "\tk1:    " << k1 << endln; 
	s << "\tFy:    " << Fyo << endln; 
	s << "\tk2:    " << kbo << endln; 
	s << "\tkv:    " << kvo << endln;
	s << "\th:     " << h << endln; 
	s << "\tPe:    " << Pe << endln; 
	s << "\tPo:    " << po << endln; 
}
