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
                                                                        
// Written: Kathryn Petek	
//          December 2004
// Modified: Chris McGann
//           January 2011
                                                                      
// Description: This file contains the implementation for the BoundingCamClay class.

#include <BoundingCamClay.h>
#include <BoundingCamClay3D.h>
#include <BoundingCamClayPlaneStrain.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

const double BoundingCamClay::one3   = 1.0/3.0 ;
const double BoundingCamClay::two3   = 2.0/3.0;
const double BoundingCamClay::root23 = sqrt(2.0/3.0);

double BoundingCamClay::mElastFlag = 1;

#include <elementAPI.h>

static int numBoundingCamClayMaterials = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_BoundingCamClayMaterial)
{
    if (numBoundingCamClayMaterials == 0) {
      numBoundingCamClayMaterials++;
      opserr << "BoundingCamClay nDmaterial - Written: C.McGann, K.Petek, P.Arduino, U.Washington\n";
    }
  
    NDMaterial *theMaterial = 0;
  
    int numArgs = OPS_GetNumRemainingInputArgs();
  
    if (numArgs < 10) {
      opserr << "Want: nDMaterial BoundingCamClay tag? massDensity? C? bulk? OCR? mu_o? alpha? lambda? h? m?" << endln;
      return 0;	
    }
    
    int tag;
    double dData[9];
  
    int numData = 1;
    if (OPS_GetInt(&numData, &tag) != 0) {
      opserr << "WARNING invalid nDMaterial BoundingCamClay material tag" << endln;
      return 0;
    }
    numData = 9;
    if (OPS_GetDouble(&numData, dData) != 0) {
      opserr << "WARNING invalid material data for nDMaterial BoundingCamClay material with tag: " << tag << endln;
      return 0;
    }
  
    theMaterial = new BoundingCamClay(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], 
                                              dData[6], dData[7], dData[8]);
    
    if (theMaterial == 0) {
      opserr << "WARNING ran out of memory for nDMaterial BoundingCamClay material with tag: " << tag << endln;
    }
  
    return theMaterial;
}

// full constructor
BoundingCamClay::BoundingCamClay(int tag, int classTag, double mDen, double C, double bulk, double OCR, double mu_o,
				                                        double Alpha, double lambda, double h, double m)
  : NDMaterial(tag,ND_TAG_BoundingCamClay),
    mEpsilon(6), 
    mEpsilon_P(6),
    mEpsilon_n_P(6),
    mSigma(6),
    mSigma_n(6),
    mSIGMAo(6),
    mSIGMAo_n(6),
    mM(6,6),
    mCe(6,6),
    mCep(6,6),
    mI1(6),
    mIIco(6,6),
	mIIcon(6,6),
	mIImix(6,6),
    mIIvol(6,6),
	mIIdevCon(6,6),
	mIIdevMix(6,6),
	mState(7)
{
	massDen = mDen;
	iC = C;
	mBulk = bulk;
	iOCR = OCR;
	imu_o = mu_o;
	ialpha = Alpha;
	ilambda = lambda;
	ih = h;
	im = m;

    this->initialize();
}

   
// null constructor
BoundingCamClay ::BoundingCamClay() 
    : NDMaterial(),
	mEpsilon(6), 
    mEpsilon_P(6),
    mEpsilon_n_P(6),
    mSigma(6),
    mSigma_n(6),
    mSIGMAo(6),
    mSIGMAo_n(6),
    mM(6,6),
    mCe(6,6),
    mCep(6,6),
    mI1(6),
    mIIco(6,6),
	mIIcon(6,6),
	mIImix(6,6),
    mIIvol(6,6),
	mIIdevCon(6,6),
	mIIdevMix(6,6),
	mState(7)
{
	massDen = 0.0;
	iC = 1.0;
	mBulk = 1.0;
	iOCR = 1.0;
	imu_o = 0.0;
	ialpha = 0.0;		
	ilambda = 1.0;
	ih = 0.0;
	im = 1.0;

    this->initialize();
}

// destructor
BoundingCamClay::~BoundingCamClay()
{
}

// zero internal variables
void BoundingCamClay::initialize()
{
	// strain and stress terms
	mEpsilon.Zero();
	mEpsilon_P.Zero();
	mEpsilon_n_P.Zero();
	mSigma.Zero();
	mSigma_n.Zero();
	mSIGMAo.Zero();
	mSIGMAo_n.Zero();

	// loading/bounding function relation scalar
	mKappa_n = iOCR - 1.0;
	mKappa = mKappa_n;

    // initialize surface radius variables
	mR_n = 1.0;
	mr_n = mR_n/iOCR;

    // initial stress ratio
	mStressRatio = 1.0;

	// initialize remaining variables
	ikappa = 0.0001;
	iepsE_vo = 0.0;

	// critical state hardening law constant
	mTHETA = 1.0/(ilambda - ikappa);

	// boolean variable
	flagReversal = false;
	    
	// 2nd-order Identity Tensor
	mI1.Zero();
	mI1(0) = 1.0;
	mI1(1) = 1.0;
	mI1(2) = 1.0;

	// 4th-order mixed variant identity
	mIImix.Zero();
	for (int i = 0; i<6; i++) {
		mIImix(i,i) = 1.0;
	}

	// 4th-order covariant identity
	mIIco = mIImix;
	mIIco(3,3) = 2.0;
	mIIco(4,4) = 2.0;
	mIIco(5,5) = 2.0;

	// 4th-order contravariant identity
	mIIcon = mIImix;
	mIIcon(3,3) = 0.5;
	mIIcon(4,4) = 0.5;
	mIIcon(5,5) = 0.5;

	// 4th-order Volumetric Tensor, IIvol = I1 tensor I1
	mIIvol.Zero();
	for (int i = 0; i<3; i++) {
		mIIvol(i,0) = 1.0;
		mIIvol(i,1) = 1.0;
		mIIvol(i,2) = 1.0;
	}

	// 4th-order contravariant deviatoric tensor
	mIIdevCon = mIIcon - one3*mIIvol;

	// 4th-order mixed variant deviatoric tensor
	mIIdevMix = mIImix - one3*mIIvol;

	// 4th-order tensor M --> needs covariant formulation
	mM = mIIco - (one3 - (pow((iC/3.0),2.0)))*mIIvol;

	// state parameter vector for recorders
	mState.Zero();

	initializeState = true;
}

NDMaterial*
BoundingCamClay::getCopy(const char *type)
{
	if (strcmp(type,"PlanStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		BoundingCamClayPlaneStrain *clone;
		clone = new BoundingCamClayPlaneStrain(this->getTag(), massDen, iC, mBulk, iOCR, imu_o, ialpha, ilambda, ih, im);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		BoundingCamClay3D *clone;
     	clone = new BoundingCamClay3D(this->getTag(), massDen, iC, mBulk, iOCR, imu_o, ialpha, ilambda, ih, im);
	 	return clone;
  	} else {
	  	opserr << "BoundingCamClay::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int 
BoundingCamClay::commitState(void)
{
	// update state variables for next step
	mEpsilon_n_P  = mEpsilon_P;
	mSigma_n      = mSigma;
	mSIGMAo_n     = mSIGMAo;

	mr_n = mr;
	mR_n = mR;
	mKappa_n = mKappa;

	return 0;
}
 
int BoundingCamClay::revertToLastCommit (void)
{
    return 0;
}

int BoundingCamClay::revertToStart(void)
{
	// added for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->initialize();
	}

    return 0;
}

NDMaterial*
BoundingCamClay::getCopy (void)
{
	opserr << "BoundingCamClay::getCopy -- subclass responsibility\n"; 
  	exit(-1);
  	return 0;
}

const char*
BoundingCamClay::getType (void) const
{
    opserr << "BoundingCamClay::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
BoundingCamClay::getOrder (void) const
{
    opserr << "BoundingCamClay::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}

//--------------------Plasticity-------------------------------------
// plasticity integration routine
void BoundingCamClay::plastic_integrator() 
{
	double f, ev, es, p, q;
	double kappa, r, R;
	double norm_e = 0;
	Vector SIGMAo(6);
	Vector epsilonET(6);
	Vector epsilonE(6);
	Vector sigma(6);
	Vector alpha(6);
	Vector xi(6);
	Vector df_dSigma(6);
	Vector e(6);
	Vector n(6);

	// initialize working variables
	kappa = mKappa_n;
	R = mR_n;
	r = mr_n;
    p = 0.0;
    q = 0.0;

	mEpsilon_P = mEpsilon_n_P;
    SIGMAo = mSIGMAo_n;

	// trial elastic strain
	epsilonET = mEpsilon - mEpsilon_P;
		
	// trial elastic volumetric strain
	ev = GetTrace(epsilonET);
	// trial elastic deviatoric strain tensor (contravariant)
	e = mIIdevCon*epsilonET;
	// norm of trial elastic deviatoric strain (contravariant-type norm)
	norm_e = GetContraNorm(e);
	// trial second invariant of deviatoric strain tensor
	es = root23*norm_e;

	// force elastic response if initialization analysis is designated===================================
	if (mElastFlag == 0) {
		//opserr << "force elastic response" << endln;
		double qCalc;

		// elastic strain = trial elastic strain
		epsilonE = epsilonET;

		// set elastic tangent tensor Ce
		mCe = GetElasticOperator(p, ev, es, n);
		// consistent elastoplastic tangent = elastic tangent
		mCep = mCe;

		// compute stress 
		sigma = mCe*epsilonE;

		// update stress invariants for recorders
		p = one3*GetTrace(sigma);
		q = sqrt(3.0/2.0)*GetContraNorm(sigma - p*mI1);
		qCalc = root23*q;

		if (p != 0) {

    		// compute reference pressure, p_o
    		if (q > 1.0e-10) {
    			double c = sqrt(pow(p,2.0) + pow(qCalc,2.0));
    			double phi = atan(qCalc/p);
    			double ang = 3.141592654 + 2.0*phi;
    			r = c/sqrt(2 - 2*cos(ang));
    		} else {
    			r = -p/2.0;
    		}
    		mp_o = -2.0*r/iC;
    
    		// designate compressibility index kappa
    		ikappa = -p/mBulk;
			// critical state hardening law constant
			mTHETA = 1.0/(ilambda - ikappa);
			mStressRatio = p/mp_o;
    
    		// compute R based on overconsolidation ratio
    		R = iOCR*r;
    
    		// surface relation scalar
    		kappa = R/r - 1.0;
    		
    		initializeState = false;
		}

    // proceed with full algorithm if initialization analysis is not designated==========================
	} else if (mElastFlag == 1) {

    	// set initial volumetric strain
    	if (!initializeState) {
    		iepsE_vo = ikappa*log(mStressRatio) + GetTrace(mEpsilon);
    		initializeState = true;
			flagReversal = true;
    	}
    
    	// trial mean normal stress
    	p = (1.0 + (3.0/2.0)*ialpha*es*es/ikappa)*mp_o*exp((iepsE_vo - ev)/ikappa);
    			
    	// trial second invariant of deviatoric stress tensor
    	q = 3.0*(imu_o - ialpha*mp_o*exp((iepsE_vo - ev)/ikappa))*es;
    
    	// compute trial deviatoric normal (contravariant)
        if (norm_e < 1.0e-13) { 
    		// check for numerical accuracy
        	n.Zero();
        } else {
    		n = e/norm_e;
    	}
    
    	// trial stress
    	sigma = p*mI1 + root23*q*n;
    	// trial backstress
    	alpha = (kappa*SIGMAo - (1.0/iC)*mI1)/(1.0 + kappa)*R;
    	// computational variable, stress difference
    	xi = sigma - alpha;
    
    	// trial yield function value
    	f = DoubleDot2_2(DoubleDot2_4(xi, mM), xi) - r*r;
    
    	// tolerance for checking trial yield function
    	double tolerance = -1.0e-7;
    	
    	// check trial state
    	if (f <= tolerance) {
			//opserr << "elastic/unloading step " << f << endln;
    		// elastic/unloading step, trial state = true state
    
        	// recenter ellipse if stress reversal
        	if (!flagReversal) {
                //opserr << "re-centering -------------- " << endln;
           	   	mSIGMAo = mSigma_n/mR_n;
        		flagReversal = true;
        	}
        		
        	// variables for kappa update
        	Vector sigma_o(6); 
        	Vector beta(6);
        	Vector temp(6);
        	double a, b, c;
        
        	// projection center
        	sigma_o = R*mSIGMAo;
        	// backstress on bounding function
        	beta = -(R/iC)*mI1;
        
        	// terms for quadratic eqn in kappa
        	temp = DoubleDot2_4((sigma - sigma_o), mM);
        	a = DoubleDot2_2(temp, (sigma - sigma_o));
        	b = 2.0*(DoubleDot2_2(temp, (sigma - beta)));
        	temp = DoubleDot2_4((sigma - beta), mM);
        	c = (DoubleDot2_2(temp, (sigma - beta))) - R*R;
        
    		// update kappa
        	kappa = (-b + sqrt((b*b - 4.0*a*c)))/(2.0*a);
        			
        	// update hardening response variable r
        	r = R/(1.0 + kappa);
        	// check axis ratio -> can create stiff system if r/R < 0.1
        	double rOverR = r/R;
        	if (rOverR < 0.1) {
        		r = 0.1*R;
    			kappa = R/r - 1.0;
        	}
    
    		// elastic strain = trial elastic strain
    		epsilonE = epsilonET;
    
    		// set elastic tangent tensor Ce
    		mCe = GetElasticOperator(p, ev, es, n);
    		// consistent elastoplastic tangent = elastic tangent
    		mCep = mCe;

    	} else if (f > tolerance) {
			//opserr << "update of loading function required " << f << endln;
    		// update of loading function required
    		double Deps_vp, rho, eta, nu;
    		double denom;
    		double dgamma;
    		Vector x(8);
    		Vector dx(8);
    		Vector Resid(8);
    			
    		// elastic strain = trial elastic strain
    		epsilonE = epsilonET;
    			
    		// set initial consistency parameter
    		dgamma = 0.0;
    		// compute df_dSigma (covariant)
    		df_dSigma = 2.0*DoubleDot4_2(mM, xi);
    		// initial change in plastic volumetric strain 
    		Deps_vp = 0.0;
    	
    		// initialize computational terms
            denom = 1.0 + mTHETA*Deps_vp;
    		rho = mTHETA*R/denom;
    		eta = mTHETA*(ih*pow(kappa,im) + r)/denom;
    		nu  = mTHETA*im*ih*pow(kappa,(im - 1.0))*Deps_vp/denom;
    			
    		// initialize residual vector {b}
    		Resid.Zero();
    		Resid(7) = f;
    	
    		// initialize vector of unknowns {x}
    		for (int i = 0; i < 6; i++) {
    			x(i) = epsilonE(i);
    		}
    		x(6) = kappa;
    		x(7) = dgamma;
    			
    		// iteration terms
    		double resNorm = 1.0;
    		double TOL = 1e-10;
    		int k = 0;
    	
    		// computational terms for Jacobian of residual
    		Matrix A(8,8);
    		Matrix Ainv(8,8);
    		Matrix Temp1(6,6);
    		Vector Temp2(6);
    		Matrix A11(6,6);
    		Vector A12(6);
    		Vector A21(6);
    		Vector A31(6);
    		double A32;
    
    		// iterative newton loop
    		while ((resNorm > TOL) && (k < 25)) {
    	
    			k += 1;
    
    			// set elastic tangent
    			mCe = GetElasticOperator(p, ev, es, n);
    
    			// Compute consistent Jacobian of the residuals [A]
    			Temp1 = Dyadic2_2(alpha, mI1);
    			Temp2 = SIGMAo + (1.0/iC)*mI1;
    
    			A11 = mIImix + 2.0*dgamma*(DoubleDot4_4(mM, mCe)) - 2.0*dgamma*rho/R*(DoubleDot4_4(mM, Temp1));
    			A12 = -2.0*dgamma*R/((1.0 + kappa)*(1 + kappa))*(DoubleDot4_2(mM, Temp2));
    			A21 = (rho - (1.0 + kappa)*eta)*mI1;
    			A31 = DoubleDot2_4(df_dSigma, mCe) - rho/R*(DoubleDot2_4(df_dSigma, Temp1)) - 2.0*eta*r*mI1;
    			A32 = 2.0*r*nu - r/(1.0 + kappa)*(DoubleDot2_2(df_dSigma, Temp2));
    
    			for (int i = 0; i < 6; i++) {
    				for (int j = 0; j < 6; j++) {
    					A(i,j) = A11(i,j);      // A11
    				}
    				A(i,6) = A12(i);            // A12
    				A(i,7) = df_dSigma(i);      // A13
    				A(6,i) = A21(i);            // A21
    				A(7,i) = A31(i);            // A31
    			}
    			A(6,6) = (1.0 + kappa)*nu - r;  // A22
    			A(6,7) = 0.0;                   // A23
    			A(7,6) = A32;                   // A32
    			A(7,7) = 0.0;                   // A33
    
    			// compute inverse of Jacobian of residuals
    			A.Invert(Ainv);
    
    			// compute incremental change in vector of unknowns
    			dx = Ainv*Resid;
    			x -= dx;
    
    			// update unknown terms using new vector {x}
    			for (int i = 0; i < 6; i++) {
    				epsilonE(i) = x(i);
    			}
    			kappa = x(6);
    			dgamma = x(7);
    
    			// update volumetric strain invariant
    			ev = GetTrace(epsilonE);
    			// update elastic deviatoric strain tensor (contravariant)
    			e = mIIdevCon*epsilonE;
    			// update norm of elastic deviatoric strain (contravariant-type norm)
    			norm_e = GetContraNorm(e);
    			// update second invariant of deviatoric strain tensor
    			es = root23*norm_e;
    
    			// update deviatoric normal (contravariant)
        		if (norm_e < 1.0e-13) { 
    				// check for numerical accuracy
        			n.Zero();
        		} else {
    				n = e/norm_e;
    			}
    
    			// update stress invariants
    			p = (1.0 + (3.0/2.0)*ialpha/ikappa*es*es)*mp_o*exp((iepsE_vo - ev)/ikappa);
    			q = 3.0*(imu_o - ialpha*mp_o*exp((iepsE_vo - ev)/ikappa))*es;
    
    			// update stress tensor
    			sigma = p*mI1 + root23*q*n; 
    			// update backstress
    			alpha = (kappa*SIGMAo - (1.0/iC)*mI1)/(1.0 + kappa)*R;
    			// update stress difference term
    			xi = sigma - alpha;
    
    			// update df_dSigma
    			df_dSigma = 2.0*DoubleDot4_2(mM, xi); 
    			// compute change in plastic volumetric strain
    			Deps_vp = GetTrace((epsilonET - epsilonE));
    			
    			// update hardening response variables (r and R)
    			denom = 1.0 + mTHETA*Deps_vp;
    			R = mR_n/denom;
    			r = (mr_n - mTHETA*ih*pow(kappa,im)*Deps_vp)/denom;
    
    			// update computational terms
    			rho = mTHETA*R/denom;
    			eta = mTHETA*(ih*pow(kappa, im) + r)/denom;
    			nu  = mTHETA*im*ih*pow(kappa, (im-1.0))*Deps_vp/denom;
    
    			// compute new yield function value
    			f = DoubleDot2_2((DoubleDot2_4(xi, mM)), xi) - r*r;
    
    			// update residual vector {b}
    			for (int i = 0; i<6; i++) {
    				Resid(i) = epsilonE(i) - epsilonET(i) + dgamma*df_dSigma(i);
    			}
    			Resid(6) = R - (1.0 + kappa)*r;
    			Resid(7) = f;
    
    			// compute norm of residual vector
    			resNorm = Resid.Norm();
    			//opserr << "iteration " << k << " residual norm " << resNorm << endln;
    		}
    
    		// update plastic strain
    		mEpsilon_P = mEpsilon - epsilonE;
			//opserr << "updating plastic strain " << endln;
			//opserr << "mEpsilon " << mEpsilon << endln;
			//opserr << "epsilonE " << epsilonE << endln;
    
    		// compute elastic compliance tensor
    		Matrix De(6,6);
    		De = GetComplianceOperator(p, ev, es, epsilonE);
    		// compute elastoplastic tangent operator
    		mCep = GetCep(kappa, r, R, dgamma, rho, eta, nu, SIGMAo, xi, df_dSigma, De);
    
    		// update remaining variables
    		flagReversal = false;
    		mSIGMAo = SIGMAo;
    	}
	}
	
	// update state variables
	mr = r;
	mR = R;
	mKappa = kappa;
    mSigma = sigma;

	// update recorder variables
	double evp = GetTrace(mEpsilon_P);
	Vector esp(6);
	esp = mIIdevMix*mEpsilon_P;
	double norm_ep = GetCovariantNorm(mEpsilon_P);
	double norm_dev_ep = GetCovariantNorm(esp);

	mState(0) = p;
	mState(1) = q;
	mState(2) = evp;
	mState(3) = norm_dev_ep;
	mState(4) = norm_ep;
	mState(5) = mr;
	mState(6) = mR;

	return;
}

Matrix
BoundingCamClay::GetElasticOperator(double p, double ev, double es, Vector n) 
// returns the elastic tangent operator
{
	Matrix Ce(6,6);
	Matrix temp(6,6);
	double Omega;
	double De11, De22, De12;

	Omega = (iepsE_vo - ev)/ikappa;

	if (mElastFlag == 0) {
		De11 = mBulk;
	} else {
		De11 = -p/ikappa;
	}
	De22 = 3.0*(imu_o - ialpha*mp_o*exp(Omega));
	De12 = 3.0*mp_o*ialpha*es*exp(Omega)/ikappa;

	temp = Dyadic2_2(mI1,n) + Dyadic2_2(n,mI1);

	// elastic tangent operator (contravariant)
	Ce = two3*De22*mIIcon + (De11 - (2.0/9.0)*De22)*mIIvol + root23*De12*temp;

	//opserr << "elastic tangent " << Ce << endln;

	return Ce;
}

Matrix
BoundingCamClay::GetComplianceOperator(double p, double ev, double es, Vector strain)
// returns the tangential elastic compliance tensor, inv(Ce)
{
	Vector e(6);
	Vector n(6);
	Matrix D(6,6);
	Matrix temp(6,6);
	double norm_e;
	double Omega;
	double De11, De12, De22, det;
	double D22inv;
	double Ee11, Ee12, Ee22;

	Omega = (iepsE_vo - ev)/ikappa;

	// terms from Ce
	if (mElastFlag == 0) {
		De11 = mBulk;
	} else {
		De11 = -p/ikappa;
	}
	De22 = 3.0*(imu_o - ialpha*mp_o*exp(Omega));
	De12 = 3.0*mp_o*ialpha*es*exp(Omega)/ikappa;

	det = De11*De22 - De12*De12;
	D22inv = 1.0/De22;

	// 2x2 matrix E is inv(D)
	Ee11 = De22/det;
	Ee22 = De11/det;
	Ee12 = -De12/det;

	// compute covariant deviatoric normal
	e = mIIdevMix*strain;
	norm_e = GetCovariantNorm(e);
	if (norm_e < 1.0e-13) { 
		// check for numerical accuracy
    	n.Zero();
    } else {
		n = e/norm_e;
	}

	temp = Dyadic2_2(mI1,n) + Dyadic2_2(n,mI1);

	// elastic compliance operator (covariant)
	D = 1.5*D22inv*mIIco + (Ee11/9.0 - 0.5*D22inv)*mIIvol + Ee12/(sqrt(6.0))*temp 
	    + 1.5*(Ee22 - D22inv)*Dyadic2_2(n,n);

	//opserr << "elastic compliance " << D << endln;

	return D;
}

Matrix 
BoundingCamClay::GetCep(double kappa, double r, double R, double dgamma, double rho, double eta,
                        double nu, Vector SIGMAo, Vector xi, Vector df_dSigma, Matrix De)
// returns the consistent elastoplastic tangent operator (contravariant)
{
	Vector alpha_R(6);
	Vector alpha_k(6);
	Vector phiAlpha_R(6);
	Vector phiAlpha_k(6);
	Matrix Q(4,4);
	Matrix Qinv(4,4);
	Matrix D(6,6);
	Matrix Cinv(6,6);
	Matrix Cep(6,6);
	double temp;
	double I1phiAlpha_R;
	double I1phiAlpha_k;

	temp = 1.0/(1.0 + kappa);
	
	// derivatives of backstress wrt R and kappa
	alpha_R = temp*(kappa*SIGMAo - (1.0/iC)*mI1);
	alpha_k = temp*temp*R*(SIGMAo + (1.0/iC)*mI1);

	// intermediate computational variables for tensor U
	phiAlpha_R = 2.0*DoubleDot4_2(mM, alpha_R);
	phiAlpha_k = 2.0*DoubleDot4_2(mM, alpha_k);

	// intermediate computational variables for matrix [Q]
	I1phiAlpha_R = DoubleDot2_2(mI1, phiAlpha_R);
	I1phiAlpha_k = DoubleDot2_2(mI1, phiAlpha_k);

	// set non-singular array [Q]
	Q.Zero();
	// first row
	Q(0,0) = 1.0 - dgamma*rho*I1phiAlpha_R;
	Q(0,2) = -(dgamma*rho*I1phiAlpha_k);
	Q(0,3) = rho*GetTrace(df_dSigma);
	// second row
	Q(1,0) = -(dgamma*eta*I1phiAlpha_R);
	Q(1,1) = 1.0;
	Q(1,2) = nu - dgamma*eta*I1phiAlpha_k;
	Q(1,3) = eta*GetTrace(df_dSigma);
	// third row
	Q(2,0) = 1.0;
	Q(2,1) = -1.0 - kappa;
	Q(2,2) = -r;
	// fourth row
	Q(3,0) = -1.0*DoubleDot2_2(df_dSigma, alpha_R);
	Q(3,1) = -2.0*r;
	Q(3,2) = -1.0*DoubleDot2_2(df_dSigma, alpha_k);

	// invert array Q
	Q.Invert(Qinv);
	
	// computational tensors L and U
    Vector L1(6);
	Vector L2(6);
	Vector L3(6);
	Vector L4(6);
	Vector U1(6);
	Vector U2(6);
	Vector U3(6);
	Vector U4(6);

	L1 = 2.0*dgamma*rho*DoubleDot2_4(mI1, mM);
	L2 = 2.0*dgamma*eta*DoubleDot2_4(mI1, mM);
	L3.Zero();
	L4 = df_dSigma;

	U1 = -dgamma*phiAlpha_R;
	U2.Zero();
	U3 = -dgamma*phiAlpha_k;
	U4 = df_dSigma;

	// compute inverse of Cep
	Cinv = De + 2.0*dgamma*mM 
		    - Qinv(0,0)*Dyadic2_2(U1,L1)
			- Qinv(0,1)*Dyadic2_2(U1,L2)
			- Qinv(0,2)*Dyadic2_2(U1,L3)
			- Qinv(0,3)*Dyadic2_2(U1,L4)
			- Qinv(1,0)*Dyadic2_2(U2,L1)
			- Qinv(1,1)*Dyadic2_2(U2,L2)
			- Qinv(1,2)*Dyadic2_2(U2,L3)
			- Qinv(1,3)*Dyadic2_2(U2,L4)
			- Qinv(2,0)*Dyadic2_2(U3,L1)
			- Qinv(2,1)*Dyadic2_2(U3,L2)
			- Qinv(2,2)*Dyadic2_2(U3,L3)
			- Qinv(2,3)*Dyadic2_2(U3,L4)
			- Qinv(3,0)*Dyadic2_2(U4,L1)
			- Qinv(3,1)*Dyadic2_2(U4,L2)
			- Qinv(3,2)*Dyadic2_2(U4,L3)
			- Qinv(3,3)*Dyadic2_2(U4,L4);

	// invert to get Cep
	Cinv.Invert(Cep);

	//opserr << "Cinv " << Cinv << endln;
	//opserr << "Cep " << Cep << endln;

	return Cep;
}

double
BoundingCamClay::GetContraNorm(Vector v)
// computes contravariant (stress-type) norm of input 6x1 tensor
{
	double result = 0.0;
	
	for (int i = 0; i < 3; i++) {
		result += v(i)*v(i);
	}
	for (int i = 3; i < 6; i++) {
		result += 2.0*v(i)*v(i);
	}

    return sqrt(result);
}
		
double
BoundingCamClay::GetCovariantNorm(Vector v) 
// computes the norm of the input argument (for strain-type storage)
{
	if (v.Size() != 6) {
		opserr << "ERROR! BoundingCamClay::NormEngStrain requires vector of size(6)!" << endln;
	}
    double result = 0.0;

	for (int i = 0; i < 3; i++) {
		result += v(i)*v(i);
	}
	for (int i = 3; i < 6; i++) {
		result += 0.5*v(i)*v(i);
	}

    return sqrt(result);
}

double
BoundingCamClay::GetTrace(Vector v) 
// computes the trace of the input argument
{
	if (v.Size() != 6)
		opserr << "ERROR! BoundingCamClay::GetTrace requires vector of size(6)!" << endln;

	return (v(0) + v(1) + v(2));
}

double
BoundingCamClay::DoubleDot2_2(Vector v1, Vector v2)
// computes doubledot product for vector-vector arguments
{
	double result = 0.0;
	
	if (v1.Size() != v2.Size()) {
		opserr << "ERROR! BoundingCamClay::DoubleDot2_2 function requires vectors of equal size!" << endln;
	}

	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}

Vector
BoundingCamClay::DoubleDot2_4(Vector v1, Matrix m1)
// computes doubledot product for vector-matrix arguments
{
	Vector result(6);
	result.Zero();

	if (v1.Size() != m1.noRows()) {
		opserr << "ERROR! BoundingCamClay::DoubleDot2_4 function requires Size(v1) = noRows(m1) " << endln;
	}

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(j) += v1(i) * m1(i,j);
	}

	return result;
}

Vector
BoundingCamClay::DoubleDot4_2(Matrix m1, Vector v1)
// computes doubledot product for matrix-vector arguments
{
	Vector result(6);
	result.Zero();

	if (m1.noCols() != v1.Size()) {
		opserr << "ERROR! BoundingCamClay::DoubleDot4_2 function requires noCols(m1) = Size(v1) " << endln;
	}

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(i) += m1(i,j) * v1(j);
	}

	return result;
}

Matrix
BoundingCamClay::DoubleDot4_4(Matrix m1, Matrix m2)
// computes doubledot product for matrix-matrix arguments
{
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < m1.noRows(); i ++) {
		for (int j = 0; j < m2.noCols(); j++) {
			for (int k = 0; k<m1.noRows(); k++) 
				result(i,j) += m1(i,k) * m2(k,j);
		}
	}
	return result;
}

Matrix 
BoundingCamClay::Dyadic2_2(Vector v1, Vector v2)
// computes dyadic product for two vector-storage arguments
{
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}

Vector
BoundingCamClay::GetState()
// returns vector of state parameters for recorders
{
	return mState;
}

Vector
BoundingCamClay::GetCenter()
// returns the projection for recorders
{
	return mSIGMAo;
}

Response*
BoundingCamClay::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->GetState());
	else if (strcmp(argv[0], "center") == 0)
		return new MaterialResponse(this, 4, this->GetCenter());
	else
		return 0;
}

int
BoundingCamClay::getResponse(int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getStrain();
			return 0;
		case 3:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = GetState();
			return 0;
		case 4:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = GetCenter();
			return 0;
		default:
			return -1;
	}
}

int
BoundingCamClay::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(8);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = iC;
  data(cnt++) = iOCR;
  data(cnt++) = ikappa;
  data(cnt++) = imu_o;
  data(cnt++) = ialpha;
  data(cnt++) = ilambda;
  data(cnt++) = ih;
  data(cnt++) = im;
  data(cnt++) = iepsE_vo;

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "BoundingCamClay::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  

  return 0;
 
}

int 
BoundingCamClay::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(7);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "BoundingCamClay::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));
  iC = data(cnt++);
  iOCR = data(cnt++);
  ikappa = data(cnt++);
  imu_o = data(cnt++);
  ialpha = data(cnt++);
  ilambda = data(cnt++);
  ih = data(cnt++);
  im = data(cnt++);
  iepsE_vo = data(cnt++);

  return 0;

}

void BoundingCamClay::Print(OPS_Stream &s, int flag )
{
  s << "BoundingCamClay" << endln;
}

int
BoundingCamClay::setParameter(const char **argv, int argc, Parameter &param)
{
    if (strcmp(argv[0],"materialState") == 0) {
        // switch elastic/plastic state
        return param.addObject(5,this);
    } else {
        // invalid parameter type
        opserr << "WARNING: invalid parameter command for BoundingCamClay nDMaterial with tag: " << this->getTag() << endln;
        return -1;
    }

    return -1;
}

int
BoundingCamClay::updateParameter(int responseID, Information &info)
{
	if (responseID == 5) {
        // materialState called - update mElasticFlag
		mElastFlag = (int)info.theDouble;
	}
	
	return 0;
}
