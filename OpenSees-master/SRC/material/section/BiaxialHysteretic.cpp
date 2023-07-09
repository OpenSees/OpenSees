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
                                                                        
// $Revision: 1.23 $
// $Date: 2009-10-02 20:48:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/BiaxialHysteretic.cpp,v $
                                                                        
                                                                        
// File: ~/section/BiaxialHysteretic.C
//
// Written: MHS
// Created: June 2000
// Revision: A
//
// Purpose: This file contains the implementation for the BiaxialHysteretic class. 

#include <stdlib.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Vector.h>
#include <Matrix.h>
#include <MatrixUtil.h>
#include <classTags.h>
#include <BiaxialHysteretic.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <ID.h>

#include <string.h>
#include <math.h>

#include <classTags.h>
#include <elementAPI.h>
#include <vector>

double BiaxialHysteretic::sqrtpi = sqrt(3.1415926535897932384626);
double BiaxialHysteretic::sqrttwo = sqrt(2.0);

void* OPS_BiaxialHysteretic()
{
    if (OPS_GetNumRemainingInputArgs() < 6) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section BiaxialHysteretic tag? k? fc? fn? alp? als? <eta? r0? rp? rs? rc? rn? Rs? sig? lmbda? code1? code2?>" << endln;
	return 0;
    }
  
    int tag;
  
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid BiaxialHysteretic tag" << endln;
	return 0;
    }
  
    double data[] = {0.0,0.0,0.0,0.0,0.0, 0.6,0.0,0.0,0.0,0.0,0.0,0.0,0.1,0.0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 14) numdata = 14;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid BiaxialHysteretic input" << endln;
	return 0;
    }

    int code[2];
    code[0] = SECTION_RESPONSE_MZ;
    code[1] = SECTION_RESPONSE_MY;
    if (numdata == 14) {
      for (int i = 0; i < 2; i++) {
	const char *type = OPS_GetString();
	code[i] = 0;
	if (strcmp(type,"Mz") == 0) 
	  code[i] = SECTION_RESPONSE_MZ;
	else if (strcmp(type,"P") == 0)
	  code[i] = SECTION_RESPONSE_P;
	else if (strcmp(type,"Vy") == 0)
	  code[i] = SECTION_RESPONSE_VY;
	else if (strcmp(type,"My") == 0)
	  code[i] = SECTION_RESPONSE_MY;
	else if (strcmp(type,"Vz") == 0)
	  code[i] = SECTION_RESPONSE_VZ;
	else if (strcmp(type,"T") == 0)
	  code[i] = SECTION_RESPONSE_T;
	else {
	  opserr << "WARNING invalid code" << endln;
	  opserr << "\nsection BiaxialHysteretic: " << tag << endln;
	  return 0;
	}
      }
    }
    
    SectionForceDeformation *sec = new BiaxialHysteretic(tag, data[0], data[1], data[2],
							 data[3], data[4], data[5],
							 data[6], data[7], data[8],
							 data[9], data[10], data[11],
							 data[12], data[13], code[0],
							 code[1]);
    return sec;
}

// #define maxOrder 10

// Vector BiaxialHysteretic::s(2);
// Matrix BiaxialHysteretic::ks(2,2);
// Matrix BiaxialHysteretic::fs(2,2);
// ID BiaxialHysteretic::code(2);

// constructors:
BiaxialHysteretic::BiaxialHysteretic(int tag, double _k, double _fc, double _fn,
				     double _alp, double _als, double _eta, double _r0,
				     double _rp, double _rs, double _rc, double _rn,
				     double _Rs, double _sig, double _lmbda,
				     int code1, int code2):
    SectionForceDeformation(tag, SEC_TAG_BiaxialHysteretic),
    Fh(0.0), kh(0.0), Fp(0.0), kp(0.0), ku(0.0),
    Et(_fn*_fn/_k), Eh(0.0), r0(_r0), rp(_rp), rs(_rs), rc(_rc), rn(_rn),
    fn(_fn), fc(_fc), k(_k), alp(_alp), als(_als), eta(_eta),
    Rs(_Rs), sig(_sig), lmbda(_lmbda), ufx(2), ufy(2),
    ui(2), u(2), Li(2), Fi(2), L(2), F(2), sF(2),
    du(2), loading(2), loadingprev(2),
    uxmax(2), uymax(2), Kt(2,2), code(2),
    otherDbTag(0), parameterID(0), dedh(2)
{
    if (k <= 0) {
	opserr << "WARNING: k <= 0 ::BiaxialHysteretic\n";
	return;
    }
    if (fn <= 0) {
	opserr << "WARNING: fn <= 0 ::BiaxialHysteretic\n";
	return;
    }
    du[0] = Vector(1);
    du[1] = Vector(1);

    loading(0) = 1.0;
    loadingprev(0) = -1.0;

    updateSprings();
    code(0) = code1;
    code(1) = code2;

}

// constructor for blank object that recvSelf needs to be invoked upon
BiaxialHysteretic::BiaxialHysteretic():
    SectionForceDeformation(0, SEC_TAG_BiaxialHysteretic),
    Fh(0.0), kh(0.0), Fp(0.0), kp(0.0), ku(0.0),
    Et(0.0), Eh(0.0), r0(0.0), rp(0.0), rs(0.0), rc(0.0), rn(0.0),
    fn(0.0), fc(0.0), k(0.0), alp(0.0), als(0.0), eta(0.0),
    Rs(0.0), sig(0.0), lmbda(0.0), ufx(2), ufy(2),
    ui(2), u(2), Li(2), Fi(2), L(2), F(2), sF(2),
    du(2), loading(2), loadingprev(2),
    uxmax(2), uymax(2), Kt(2,2), code(2),
    otherDbTag(0), parameterID(0), dedh(2)
{
    du[0] = Vector(1);
    du[1] = Vector(1);
    loading(0) = 1.0;
    loadingprev(0) = -1.0;

    updateSprings();
}

// destructor:
BiaxialHysteretic::~BiaxialHysteretic()
{

}

int BiaxialHysteretic::setTrialSectionDeformation (const Vector &def)
{
    // set new displacement
    u = def;

    // update internal states
    Kt.Zero();
    for (int a=0; a<2; ++a) {

	// update loading state
	if (updateLoadingState(a) < 0) {
	    opserr << "WARNING: failed to update loading state\n";
	    return -1;
	}

	// update zero force point
	if (loading(a) == -1) {
	    if (updateZeroForcePoint(a) < 0) {
		opserr << "WARNING: failed to update zero force point\n";
		return -1;
	    }
	    if (updateLoadingState(a) < 0) {
		opserr << "WARNING: failed to update loading state\n";
		return -1;
	    }
	}

	// update forces
	if (updateForce(a) < 0) {
	    opserr << "WARNING: failed to update force\n";
	    return -1;
	}

	// update tangent
	if (updateTangent(a) < 0) {
	    opserr << "WARNING: failed to update tangent\n";
	    return -1;
	}
    }

    // update tangent for spring 3
    if (updateTangent(2) < 0) {
	opserr << "WARNING: failed to update tangent\n";
	return -1;
    }
  
    return 0;
}

const Vector &
BiaxialHysteretic::getSectionDeformation(void)
{
    return u;
}

const Matrix &
BiaxialHysteretic::getSectionTangent(void)
{
    return Kt;
}

const Matrix &
BiaxialHysteretic::getInitialTangent(void)
{
    Kt(0,0) = k;
    Kt(1,1) = k;
    Kt(0,1) = 0.0;
    Kt(1,0) = 0.0;

    return Kt;
}

const Matrix &
BiaxialHysteretic::getInitialFlexibility(void)
{
    Kt(0,0) = 1.0/k;
    Kt(1,1) = 1.0/k;
    Kt(0,1) = 0.0;
    Kt(1,0) = 0.0;
    
    return Kt;
}

const Vector &
BiaxialHysteretic::getStressResultant(void)
{
    sF(0) = 0.0;
    sF(1) = 0.0;

    for (int a=0; a<2; ++a) {
	if (L(a) > 0) {
	    sF(0) += F(a)*(u(0)-ufx(a))/L(a);
	    sF(1) += F(a)*(u(1)-ufy(a))/L(a);
	}
    }

    sF(0) += ku*u(0);
    sF(1) += ku*u(1);

    return sF;
}

SectionForceDeformation *
BiaxialHysteretic::getCopy(void)
{
    BiaxialHysteretic *theCopy = 0;
    
    theCopy = new BiaxialHysteretic(this->getTag(),
				    k,fc,fn,alp,als,eta,
				    r0,rp,rs,rc,rn,
				    Rs,sig,lmbda,code(0),code(1));
    return theCopy;
}

const ID&
BiaxialHysteretic::getType ()
{
    return code;
}

int
BiaxialHysteretic::getOrder () const
{
    return 2;
}

int
BiaxialHysteretic::commitState(void)
{
    // maximum inelastic displacements
    getStressResultant();
    double uxdd = u(0) - sF(0)/k;
    double uydd = u(1) - sF(1)/k;
    if (u(0)>0 && uxdd>uxmax(0)) {
	uxmax(0) = uxdd;
    } else if (u(0)<0 && uxdd<uxmax(1)) {
	uxmax(1) = uxdd;
    }
    if (u(1)>0 && uydd>uymax(0)) {
	uymax(0) = uydd;
    } else if (u(1)<0 && uydd<uymax(1)) {
	uymax(1) = uydd;
    }

    // update states
    updateEnergy();
    updateSprings();

    ui = u;
    Fi = F;
    Li = L;
    du[0] = Vector(1);
    du[1] = Vector(1);
    loadingprev = loading;
    loading = ID(2);
    
    return 0;
}

int
BiaxialHysteretic::revertToLastCommit(void)
{
    return 0;
}	

int
BiaxialHysteretic::revertToStart(void)
{
    return 0;
}


int
BiaxialHysteretic::sendSelf(int cTag, Channel &theChannel)
{
    int res = 0;

    return res;
}


int
BiaxialHysteretic::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    int res = 0;

    return res;
}

void
BiaxialHysteretic::Print(OPS_Stream &s, int flag)
{
    s << "BiaxialHysteretic, tag: " << this->getTag() << endln;
}


// AddingSensitivity:BEGIN ////////////////////////////////////
int
BiaxialHysteretic::setParameter(const char **argv, int argc, Parameter &param)
{
    if (argc < 1)
	return -1;

    // if (strstr(argv[0],"k") != 0 || strstr(argv[0],"E") != 0) {
    // 	param.setValue(k);
    // 	return param.addObject(1, this);
    // }
    // if (strstr(argv[0],"My") != 0 || strstr(argv[0],"Fy") != 0) {
    // 	param.setValue(My);
    // 	return param.addObject(2, this);
    // }
    // if (strstr(argv[0],"a") != 0 || strstr(argv[0],"b") != 0) {
    // 	param.setValue(a);
    // 	return param.addObject(3, this);
    // }

    return -1;
}

const Vector &
BiaxialHysteretic::getSectionDeformationSensitivity(int gradIndex)
{
    return u;
}

const Vector &
BiaxialHysteretic::getStressResultantSensitivity(int gradIndex,
						 bool conditional)
{
    return sF;
}

const Matrix &
BiaxialHysteretic::getSectionTangentSensitivity(int gradIndex)
{
    return Kt;
}

const Matrix&
BiaxialHysteretic::getInitialTangentSensitivity(int gradIndex)
{
    return Kt;
}

int
BiaxialHysteretic::commitSensitivity(const Vector& defSens,
				     int gradIndex, int numGrads)
{
    return 0;
}

// AddingSensitivity:END ///////////////////////////////////

const Vector&
BiaxialHysteretic::getdedh(void)
{
    return dedh;
}


void
BiaxialHysteretic::updateSprings()
{
    if (loadingprev(0)*loading(0) == -1 ||
	loadingprev(1)*loading(1) == -1) {
	double energyRatio = -Eh/Et;

	double phi0 = exp(r0*energyRatio);
	double phip = exp(rp*energyRatio);
	double phis = exp(rs*energyRatio);
	double phic = exp(rc*energyRatio);
	double phin = exp(rn*energyRatio);

	double k0 = phi0*k;
	double kp0 = phip*alp*k;
    
	// ultimate spring
	ku = phis*als*k;
    
	double fc0 = phic*fc;
	double fn0 = phin*fn;
    
    
	double uc = fc0/k0;
	double un = uc + (fn0-fc0)/kp0;
        
	double k1 = k0 - ku;
	double k2 = kp0 - ku;

	double fn1 = fn0 - ku*un;
	double fc1 = fc0 - ku*uc;

	double uc1 = fc1/k1;

	// hysteretic spring
	Fh = fc1 - k2*uc1;
	kh = k1 - k2;

	// post-yielding spring
	Fp = fn1-Fh;
	kp = k2;

    }
}

int
BiaxialHysteretic::updateLoadingState(int a)
{
    if (a < 0 || a > 1) return -1;
    
    // previous and current distance to zero force point
    Li(a) = sqrt((ui(0)-ufx(a))*(ui(0)-ufx(a))+(ui(1)-ufy(a))*(ui(1)-ufy(a)));
    L(a) = sqrt((u(0)-ufx(a))*(u(0)-ufx(a))+(u(1)-ufy(a))*(u(1)-ufy(a)));

    // direction of n
    double nx=0, ny=0;
    if (fabs(u(0)-ui(0)) < 1e-12) {
	nx = ui(0);
	ny = ufy(a);
    } else {
	double m = (u(1)-ui(1))/(u(0)-ui(0));
	nx = (ufx(a)+ui(0)*m*m-ui(1)*m+ufy(a)*m)/(m*m+1);
	ny = (ufy(a)*m*m+ui(1)-ui(0)*m+ufx(a)*m)/(m*m+1);
    }

    double n = (ui(0)-nx)*(u(0)-nx)+(ui(1)-ny)*(u(1)-ny);

    // get du
    if (n >= 0) {
	du[a] = Vector(1);
	du[a](0) = L(a) - Li(a);
	loading(a) = sign(du[a](0));
    } else {
	double L3 = sqrt((nx-ufx(a))*(nx-ufx(a))+(ny-ufy(a))*(ny-ufy(a)));
	du[a] = Vector(2);
	du[a](0) = L3 - Li(a);
	du[a](1) = L(a) - L3;
	loading(a) = 1;
	if (du[a](0)>0 || du[a](1)<0) {
	    opserr << "WARNING: material is unloading then loading but calculated du shows the opposite\n";
	    return -1;
	}
    }

    
    return 0;
}


int
BiaxialHysteretic::updateZeroForcePoint(int a)
{
    if (a < 0 || a > 1) return -1;
    
    // last is zero force point
    if (Li(a) == 0) return 0;

    // spring 1
    double Lhat = 0.0;
    if (a == 0) {
	double b = Fh/(1.0-2.0*eta);
	double s = Rs*sqrt((uxmax(0)-uxmax(1))*(uxmax(0)-uxmax(1))+
			   (uymax(0)-uymax(1))*(uymax(0)-uymax(1)));
	Lhat = Li(a)-b/kh*log(b/(b-Fi(a)))
	    +s*(erf(lmbda*Fh/(sqrttwo*sig*Fh))-erf((Fi(a)+lmbda*Fh)/(sqrttwo*sig*Fh)));
    } else if (a == 1) {
	Lhat = Li(a) - Fi(a)/kp;
    }

    ufx(a) += (ui(0)-ufx(a))*Lhat/Li(a);
    ufy(a) += (ui(1)-ufy(a))*Lhat/Li(a);
    
    return 0;
}

int
BiaxialHysteretic::updateForce(int a)
{
    if (a == 0) {
	if (du[a].Size() == 1) {
	    double b = Fh/(eta*sign(du[a](0))+1-eta);

	    if (Rs == 0.0) {
		F(a) = b-(b-Fi(a))*exp(-du[a](0)*kh/b);
	    } else {
		double x = Fi(a);
		if (newton(x,du[a](0),Fi(a)) < 0) {
		    opserr << "WARNING: failed to converge to get force\n";
		    return -1;
		}
		F(a) = x;
	    }
	} else {

	    if (Rs == 0.0) { 
		double b = Fh/(eta*sign(du[a](0))+1-eta);
		double Fhat = b - (b-Fi(a))*exp(-du[a](0)*kh/b);
		b = Fh/(eta*sign(du[a](1))+1-eta);
		F(a) = b - (b-Fhat)*exp(-du[a](1)*kh/b);
	    } else {
		double x = Fi(a);
		if (newton(x,du[a](0),Fi(a)) < 0) {
		    opserr << "WARNING: failed to converge to get force\n";
		    return -1;
		}
		double Fhat = x;
		x = Fhat;
		if (newton(x,du[a](1),Fhat) < 0) {
		    opserr << "WARNING: failed to converge to get force\n";
		    return -1;
		}
		F(a) = x;
	    }
	}
    } else if (a == 1) {

	if (du[a].Size() == 1) {
	    F(a) = Fi(a) + kp*du[a](0);
	} else {
	    double Fhat = Fi(a) + kp*du[a](0);
	    if (Fhat > Fp) Fhat = Fp;
	    if (Fhat < -Fp) Fhat = -Fp;
	    F(a) = Fhat + kp*du[a](1);
	}

	if (F(a) > Fp) F(a) = Fp;
	if (F(a) < -Fp) F(a) = -Fp;
    }

    return 0;
}

int
BiaxialHysteretic::updateTangent(int a)
{
    if (a == 0) {
	double qx = u(0) - ufx(a);
	double qy = u(1) - ufy(a);
	if (L(a) == 0) return 0;
	double La2 = L(a)*L(a);
	double La3 = La2*L(a);

	double b = Fh/(eta*sign(du[a](0))+1-eta);
	double s =Rs*sqrt((uxmax(0)-uxmax(1))*(uxmax(0)-uxmax(1))+
			  (uymax(0)-uymax(1))*(uymax(0)-uymax(1)));
	double Q = (F(a)-lmbda*Fh*sign(du[a](0)))/(sqrttwo*sig*Fh);
	double dFdu = 1.0/(s*sqrttwo*exp(-Q*Q)/(sqrtpi*sig*Fh)+b/(kh*(b-F(a))));

	Kt(0,0) += dFdu*qx*qx/La2 + F(a)*qy*qy/La3;
	Kt(0,1) += dFdu*qx*qy/La2 - F(a)*qx*qy/La3;
	Kt(1,0) += dFdu*qx*qy/La2 - F(a)*qx*qy/La3;
	Kt(1,1) += dFdu*qy*qy/La2 + F(a)*qx*qx/La3;
	
    } else if (a == 1) {

	double qx = u(0) - ufx(a);
	double qy = u(1) - ufy(a);
	if (F(a) < Fp) {
	    Kt(0,0) += kp;
	    Kt(1,1) += kp;
	} else {
	    Kt(0,0) += qy*qy;
	    Kt(0,1) += -qx*qy;
	    Kt(1,0) += -qx*qy;
	    Kt(1,1) += qx*qx;
	}
    } else if (a == 2) {
	Kt(0,0) += ku;
	Kt(1,1) += ku;
    }
    
    return 0;
}

void
BiaxialHysteretic::updateEnergy()
{
    bool adding = true;
    if (loadingprev(0)*loading(0) == -1) {
	// first spring
	double b = Fh/(eta*sign(du[0](0))+1-eta);
	double s =Rs*sqrt((uxmax(0)-uxmax(1))*(uxmax(0)-uxmax(1))+
			  (uymax(0)-uymax(1))*(uymax(0)-uymax(1)));
	double Q = (F(0)-lmbda*Fh*sign(du[0](0)))/(sqrttwo*sig*Fh);
	double dFdu = 1.0/(s*sqrttwo*exp(-Q*Q)/(sqrtpi*sig*Fh)+b/(kh*(b-F(0))));
	Eh -= 0.5*F(0)*F(0)/dFdu;
	if (Eh < 0) Eh = 0;
	adding = false;
    }

    if (loadingprev(1)*loading(1) == -1) {
	
	// second spring
	Eh -= 0.5*F(1)*F(1)/kp;
	if (Eh < 0) Eh = 0;
	adding = false;
    }

    if (adding) {
	for (int a=0; a<2; ++a) {
	    double f = 0.5*(F(a)+Fi(a));
	    Eh += fabs(f*du[a](0));
	}
    }
}

int
BiaxialHysteretic::sign(double v)
{
    if (v > 0) return 1;
    if (v < 0) return -1;
    return 0;
}

double
BiaxialHysteretic::spring1(double du1, double F1, double F)
{
    double b = Fh/(eta*sign(du1)+1-eta);
    double s = Rs*sqrt((uxmax(0)-uxmax(1))*(uxmax(0)-uxmax(1))+
		       (uymax(0)-uymax(1))*(uymax(0)-uymax(1)));
    double Qi = (F1-lmbda*Fh*sign(du1))/(sqrttwo*sig*Fh);
    double Q = (F-lmbda*Fh*sign(du1))/(sqrttwo*sig*Fh);

    return exp(-kh/b*(du1-s*(erf(Q)-erf(Qi))))-(b-F)/(b-F1);
}

double
BiaxialHysteretic::dspring1(double du1, double F1, double F)
{
    double b = Fh/(eta*sign(du1)+1-eta);
    double s = Rs*sqrt((uxmax(0)-uxmax(1))*(uxmax(0)-uxmax(1))+
		       (uymax(0)-uymax(1))*(uymax(0)-uymax(1)));
    double Qi = (F1-lmbda*Fh*sign(du1))/(sqrttwo*sig*Fh);
    double Q = (F-lmbda*Fh*sign(du1))/(sqrttwo*sig*Fh);

    return exp(-kh/b*(du1-s*(erf(Q)-erf(Qi))))*(kh*s*sqrttwo*exp(-Q*Q))/(b*sqrtpi*sig*Fh)+1.0/(b-F1);
}

int
BiaxialHysteretic::newton(double& x, double du1, double F1, 
			  double tol, int maxiter)
{
    double R = spring1(du1,F1,x);
    int niter = 0;

    while (niter <= maxiter) {
	if (fabs(R) < tol) {
	    return 0;
	}
	double Kt = dspring1(du1,F1,x);
	double dx = R/(-Kt);

	x += dx;

	R = spring1(du1,F1,x);
	++niter;
    }
    
    return -1;
}
