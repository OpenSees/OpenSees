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
                                                                        
// $Revision: 1.3 $
// $Date: 2010-03-04 19:10:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoundingCamClay.cpp,v $

// Written: kap	
// Created: 12/04
                                                                      
// Description: This file contains the implementation for the BoundingCamClay class.
//				
//		This implementation of Drucker-Prager allows for non-associative flow 
//		through the use of the rho_bar parameter.  It also includes linear 
//		kinematic hardening and linear and nonlinear isotropic hardening.
//
//

#include <BoundingCamClay.h>
#include <BoundingCamClay3D.h>
#include <Information.h>
#include <MaterialResponse.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

//#include <myDebug.h>

const double BoundingCamClay :: one3   = 1.0 / 3.0 ;
const double BoundingCamClay :: two3   = 2.0 / 3.0 ;
const double BoundingCamClay :: root23 = sqrt( 2.0 / 3.0 ) ;

#include <elementAPI.h>
static int numBoundingCamClayMaterials = 0;
#define OPS_Export extern "C"

OPS_Export void *
OPS_NewBoundingCamClayMaterial(void)
{
  if (numBoundingCamClayMaterials == 0) {
    numBoundingCamClayMaterials++;
    opserr << "BoundingCamClay nDmaterial - Written by Kathryn Petek and Pedro Arduino - Copyright@2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 12) {
    opserr << "Want: nDMaterial BoundingCamClay tag? C? r? R? p_o? kappa? mu_o? alpha? lambda? h? m? epsE_vo?" << endln;
    return 0;	
  }
  
  int tag;
  double dData[11];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid tag for  BoundingCamClay material" << endln;
    return 0;
  }
  numData = 11;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial BoundingCamClay material  with tag: " << tag << endln;
    return 0;
  }

  theMaterial = new BoundingCamClay(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], 
				    dData[5], dData[6], dData[7], dData[8], dData[9], dData[10]);

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial BoundingCamClay material  with tag: " << tag << endln;
  }

  return theMaterial;
}




//full constructor
BoundingCamClay::BoundingCamClay(int tag, int classTag, double C, double r, double R, double p_o, double kappa,
				 double mu_o, double alpha, double lambda, double h, 
				 double m, double epsE_vo)
    : NDMaterial(tag,ND_TAG_BoundingCamClay),
      mEpsilon(6), 
      //mEpsilon_n(6),
      mEpsilon_P(6),
      mEpsilon_n_P(6),
      mSigma(6),
      mSigma_n(6),
      mSIGMAo(6),
      mSIGMAo_n(6),
      //mBeta(6),
      //mBeta_n(6),
      mM(6,6),
      mphi(6,6),
      mCe(6,6),
      mCep(6,6),
      mI1(6),
      mII(6,6),
      mIIvol(6,6)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::BoundingCamClay(...)" << endln;
#endif

		iC = C;
		mr_n = r;
		mR_n = R;
		ip_o = p_o;
		ikappa = kappa;
		imu_o = mu_o;
		ialpha = alpha;
		ilambda = lambda;
		ih = h;
		im = m;
		iepsE_vo = epsE_vo;

        this->initialize();
}

   
//null constructor
BoundingCamClay ::BoundingCamClay  () 
    : NDMaterial(),
    mEpsilon(6), 
	//mEpsilon_n(6),
	mEpsilon_P(6),
	mEpsilon_n_P(6),
	mSigma(6),
    mSigma_n(6),
	mSIGMAo(6),
	mSIGMAo_n(6),
	//mBeta(6),
	//mBeta_n(6),
    mM(6,6),
	mphi(6,6),
	mCe(6,6),
	mCep(6,6),
	mI1(6),
	mII(6,6),
    mIIvol(6,6)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::BoundingCamClay()" << endln;
#endif

		iC = 1.0;			// watch out for 1/C
		ip_o = 0.0;
		ikappa = 0.0;       
		imu_o = 0.0;
		ialpha = 0.0;		
		ilambda = 1.0;      // watch out for 1/(lambda - kappa)
		ih = 0.0;
		im = 1.0;
		iepsE_vo = 0.0;

        this->initialize();
}

//destructor
BoundingCamClay::~BoundingCamClay  ()
{
}

//zero internal variables
void BoundingCamClay::initialize( )
{
#ifdef DEBUG
        opserr << "BoundingCamClay::initialize()" << endln;
#endif

		mKappa_n = mR_n / mr_n - 1.0;
		mKappa = mKappa_n;
		mTHETA = 1.0/(ilambda - ikappa);

		mEpsilon.Zero();
		//mEpsilon_n.Zero();
		mEpsilon_P.Zero();
		mEpsilon_n_P.Zero();
		mSigma.Zero();
		mSigma_n.Zero();
		mSIGMAo.Zero();
		mSIGMAo_n.Zero();

		flagReversal = false;
	    
		// 2nd order Identity Tensor
		mI1.Zero();
		mI1(0) = 1.0;
		mI1(1) = 1.0;
		mI1(2) = 1.0;

		// 4th order Identity Tensor
		mII.Zero();
		for (int i = 0;  i<6; i++)
			mII(i,i) = 1.0;

		// 4th order Volumetric Tensor
		// IIvol = I1 tensor I1
		mIIvol.Zero();
		for (int i = 0; i<3; i++){
			mIIvol(i,0) = 1.0;
			mIIvol(i,1) = 1.0;
			mIIvol(i,2) = 1.0;
		}

		mM.Zero();
		mM = mII + (-one3 + pow((iC*one3),2.0)) * mIIvol;
		mphi = 2*mM;

}


NDMaterial * BoundingCamClay::getCopy (const char *type)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::getCopy(..)" << endln;
#endif

  if (strcmp(type,"ThreeDimensional")==0 || strcmp(type, "3D") ==0) 
  {  BoundingCamClay3D *clone;
     clone = new BoundingCamClay3D(this->getTag(), iC, mr_n, mR_n, ip_o, ikappa, imu_o,
		 ialpha, ilambda, ih, im, iepsE_vo);
	 return clone;
  }
  else {
	  opserr << "BoundingCamClay::getCopy failed to get copy: " << type << endln;
	  return 0;
  }
}

int BoundingCamClay::commitState (void)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::commitState()" << endln;
#endif

		//mEpsilon_n    = mEpsilon;
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
#ifdef DEBUG
        opserr << "BoundingCamClay::revertToLastCommit()" << endln;
#endif

    return 0;
}

int BoundingCamClay::revertToStart(void)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::revertToStart()" << endln;
#endif

    this->initialize();

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

//plasticity integration routine
void BoundingCamClay:: plastic_integrator() 
{
#ifdef DEBUG
        opserr << "BoundingCamClay::plastic_integrator()" << endln;
#endif

		double f, ev, es, p, q;
		double kappa, r, R;
		Vector epsilonET(6);      // Trial Elastic strain
		Vector sigma(6);
		Vector alpha(6);
		Vector sigmaMinusAlpha(6);
		Vector dfdSigma(6);
		Vector SIGMAo(6);
		Vector s(6);
		Vector n(6);

		// initialize working variables
		kappa = mKappa_n;
		r = mr_n;
		R = mR_n;

		mEpsilon_P = mEpsilon_n_P;
		SIGMAo = mSIGMAo_n;

		epsilonET = mEpsilon - mEpsilon_P;  // trial elastic strain

		ev = GetVolInv(epsilonET);
		es = GetDevInv(epsilonET);

		p = (1.0 + 3.0/2.0*ialpha/ikappa*es*es)*ip_o*exp((iepsE_vo-ev)/ikappa);
		q = 3.0*(imu_o - ialpha*ip_o*exp((iepsE_vo-ev)/ikappa))*es;

		s = epsilonET - one3*ev*mI1;
		n = s;
        if (Norm_EngStrain(s) == 0) 
            n = n/1e-15;
        else 
			n = n / (Norm_EngStrain(s));

		sigma = p*mI1 + root23*q*n; 
		alpha = (kappa*SIGMAo - 1.0/iC*mI1)/(1.0+kappa)*R;
		sigmaMinusAlpha = sigma - alpha;
		dfdSigma = DoubleDot4_2(mphi, sigmaMinusAlpha);

		// Calculate f
		f = DoubleDot2_2(DoubleDot2_4(sigmaMinusAlpha, mM), sigmaMinusAlpha) - r*r;

		// elastic step
		if (f <= 0.0) {
opserr << "f <= 0.0 : elastic step or unloading " << endln;
			// recenter ellipse if stress reversal
			if (!flagReversal) {
                SIGMAo = mSigma_n / R;		// = sigma_n / R_n
				flagReversal = true;
			}

			// update kappa
			Vector sigma_o(6);
			Vector beta(6);
			Vector temp(6);
			double a, b, c;

			sigma_o = R * SIGMAo;
			beta = - R/iC*mI1;

			temp = DoubleDot2_4((sigma - sigma_o), mM);
			a = DoubleDot2_2(temp, (sigma - sigma_o));
			b = 2* DoubleDot2_2(temp, (sigma - beta));
			temp = DoubleDot2_4((sigma - beta), mM);
			c = DoubleDot2_2(temp, (sigma - beta))-R*R;

			kappa = (-b + sqrt((b*b - 4*a*c)))/ 2*a;

			// update r
			r = R / (1.0 + kappa);

			double rOverR = 0.1;
			if ( (r/R) < rOverR) {
				r = rOverR * R;
			}

			// mCep = elastic modulus
			mCe = GetElasticModulus(p, q, ev, es, n);
			mCep = mCe;
		
		} else {
			double dgamma = 0.0;
			double Deps_vp, rho, eta, nu;
			double denom;
			Vector epsilonE(6);
			Vector x(8);
			Vector dx(8);
			Vector Resid(8);

			epsilonE = epsilonET;

			// initialize all terms for Newton Iteration
			Deps_vp = 0.0;  // = tr(dgamma * dfdsigma) = tr(0)
            denom = (1.0 + mTHETA * Deps_vp);
			
			rho = mTHETA*R / denom;
			eta = mTHETA * (ih * pow(kappa, im) + r) / denom;
			nu  = mTHETA*im*ih*pow(kappa, (im-1.0)) * Deps_vp / denom;

			Resid.Zero();
			Resid(7) = f;

			for (int i = 0; i < 6; i++)
				x(i) = epsilonE(i);
			x(6) = kappa;
			x(7) = dgamma;
			
			double RESK, RES0, RESPrev;
			RESK = 1.0;
			RES0 = 1.0;
			RESPrev = 1.0;

			int k = 0;

			// Newton loop
			double TOL = 1e-6;
			//while (Resid.Norm() > TOL) {
			//while ((RESK > TOL) && (RESK <= RESPrev)) {
			while (RESK > TOL) {
				Matrix A(8,8);
				Matrix Ainv(8,8);
				Matrix Temp1(6,6);
				Vector Temp2(6);
				Matrix A11(6,6);
				Vector A12(6);
				Vector A21(6);
				Vector A31(6);
				//Matrix Ce(6,6);

				k += 1;

				Temp1 = Dyadic2_2(alpha, mI1);
				Temp2 = SIGMAo + 1.0/iC*mI1;
				mCe = GetElasticModulus(p, q, ev, es, n);

				// Compute A and Ainv
				A11 = mII + dgamma*DoubleDotStrain4_4(mphi, mCe) - dgamma*rho/R*DoubleDot4_4(mphi, Temp1);
				A12 = -dgamma*R/(1.0+kappa)/(1+kappa) * DoubleDot4_2(mphi, Temp2);
				A21 = (rho - (1.0+kappa)*eta)*mI1;
				A31 = DoubleDotStrain2_4(dfdSigma, mCe) - rho/R*DoubleDot2_4(dfdSigma, Temp1) - 2.0*eta*r*mI1;

				for (int i = 0; i < 6; i++) {
					for (int j = 0; j < 6; j++) {
						A(i,j) = A11(i,j);   // A11
					}
					A(i,6) = A12(i);         // A12
					A(i,7) = dfdSigma(i);    // A13
					A(6,i) = A21(i);         // A21
					A(7,i) = A31(i);         // A31
				}
				A(6,6) = (1.0+kappa)*nu - r;  // A22
				A(6,7) = 0.0;                 // A23
				A(7,6) = 2.0*r*nu - r/(1+kappa) * DoubleDot2_2(dfdSigma, Temp2);  // A32
				A(7,7) = 0.0;                 // A33

				A.Invert(Ainv);

				dx = -1*Ainv*Resid; // Note: should be dx = -Ainv*Resid
				x += dx;  // Note: should be x += dx  .  Corrected for above

				// update all terms
				for (int i = 0; i < 6; i++) 
					epsilonE(i) = x(i);

				kappa  = x(6);
				dgamma = x(7);

				// update all terms 
				ev = GetVolInv(epsilonE);
				es = GetDevInv(epsilonE);

				p = (1.0 + 3.0/2.0*ialpha/ikappa*es*es)*ip_o*exp((iepsE_vo-ev)/ikappa);
				q = 3.0*(imu_o - ialpha*ip_o*exp((iepsE_vo-ev)/ikappa))*es;

				s = epsilonE - one3*ev*mI1;
				n = s;
				if (Norm_EngStrain(s) == 0) 
                    n = n/1e-15;
                else 
					n = n / (Norm_EngStrain(s));

				sigma = p*mI1 + root23*q*n; 
				alpha = (kappa*SIGMAo - 1.0/iC*mI1)/(1.0+kappa)*R;
				sigmaMinusAlpha = sigma - alpha;
				dfdSigma = DoubleDot4_2(mphi, sigmaMinusAlpha); 

				Deps_vp = GetVolInv((epsilonET - epsilonE));
				denom = (1 + mTHETA * Deps_vp);

				R = mR_n/denom;  // update R
				r = (mr_n - mTHETA*ih*pow(kappa,im)*Deps_vp)/denom;  // update r

				rho = mTHETA*R / denom;
				eta = mTHETA * (ih * pow(kappa, im) + r) / denom;
				nu  = mTHETA*im*ih*pow(kappa, (im-1)) * Deps_vp / denom;

//opserr << "A   = " << A;
//opserr << "Ainv = " << Ainv << endln; 
opserr << "r (k)   = " << r << endln;
opserr << "R (k)   = " << R << endln;
opserr << "mTHETA = " << mTHETA << endln;
opserr << "ih = " << ih << endln;
opserr << "kappa = " << kappa << endln;
opserr << "im = " << im << endln;
opserr << "Deps_vp = " << Deps_vp << endln;
opserr << "denom = " << denom << endln;
opserr << " pow(kappa,im)  = " << pow(kappa,im)  << endln;
opserr << "eta = " << eta << endln;
opserr << endln;
opserr << "dfdSigma = " << dfdSigma;
opserr << "mCe = " << mCe;
opserr << "DoubleDotStrain2_4(dfdSigma, mCe) = " << DoubleDotStrain2_4(dfdSigma, mCe);
opserr << endln;

				// Calculate f
				f = DoubleDot2_2(DoubleDot2_4(sigmaMinusAlpha, mM), sigmaMinusAlpha) - r*r;

				// Compute residual
				for (int i = 0; i<6; i++) 
					Resid(i) = epsilonE(i) - epsilonET(i) + dgamma*dfdSigma(i);
				Resid(6) = R - (1.0+kappa)*r;
				Resid(7) = f;

				RESPrev = RESK;
				RESK = Resid.Norm();
				if ( k == 1)
					RES0 = RESK;
				//RESK = RESK/RES0;

opserr << "Resid = " << Resid;
opserr << "RESK = " << RESK << endln;
opserr << "END OF " << k << " ITERATION on R -------------------------- " << endln;
opserr << endln;
			}
opserr << "CONVERGED in " << k << " steps____________________________________________" << endln;


			// update terms
			mEpsilon_P = mEpsilon - epsilonE;
			//mbeta = alpha*(1+kappa) - kappa*SIGMAo*R;


			// Compute Cep
			mCep = GetConsistentModulus(kappa, dgamma, rho, eta, nu, r, R, p, ev, es, SIGMAo, sigmaMinusAlpha, dfdSigma, n);
opserr << "mCep = " << mCep;
			flagReversal = false;
		}



		// update remaining state variables
		mr = r;
		mR = R;
		mKappa = kappa;
		mSIGMAo = SIGMAo;
		mSigma = sigma;

opserr << "END OF SetTrialStrain() *************************************************************" << endln;
opserr << endln;
		

	return;
}

////plasticity integration routine
//void BoundingCamClay:: plastic_integrator() 
//{
//#ifdef DEBUG
//        opserr << "BoundingCamClay::plastic_integrator()" << endln;
//#endif
//
//		double f, ev, es, p, q;
//		double kappa, r, R;
//		Vector epsilonET(6);      // Trial Elastic strain
//		Vector sigma(6);
//		Vector alpha(6);
//		Vector sigmaMinusAlpha(6);
//		Vector dfdSigma(6);
//		Vector SIGMAo(6);
//		Vector s(6);
//		Vector n(6);
//
//		// initialize working variables
//		kappa = mKappa_n;
//		r = mr_n;
//		R = mR_n;
//
//		mEpsilon_P = mEpsilon_n_P;
//		SIGMAo = mSIGMAo_n;
//
//		epsilonET = mEpsilon - mEpsilon_P;  // trial elastic strain
//
//		ev = GetVolInv(epsilonET);
//		es = GetDevInv(epsilonET);
//
//		p = (1.0 + 3.0/2.0*ialpha/ikappa*es*es)*ip_o*exp((iepsE_vo-ev)/ikappa);
//		q = 3.0*(imu_o - ialpha*ip_o*exp((iepsE_vo-ev)/ikappa))*es;
//
//		s = epsilonET - one3*ev*mI1;
//		n = s;
//        if (Norm_EngStrain(s) == 0) 
//            n = n/1e-15;
//        else 
//			n = n / (Norm_EngStrain(s));
//
//		sigma = p*mI1 + root23*q*n; 
//		alpha = (kappa*SIGMAo - 1.0/iC*mI1)/(1.0+kappa)*R;
//		sigmaMinusAlpha = sigma - alpha;
//		
//		// Calculate f
//		f = DoubleDot2_2(DoubleDot2_4(sigmaMinusAlpha, mM), sigmaMinusAlpha) - r*r;
//
//		// elastic step
//		if (f <= 0.0) {
//opserr << "f <= 0.0 : elastic step or unloading " << endln;
//			// recenter ellipse if stress reversal
//			if (!flagReversal) {
//                SIGMAo = mSigma_n / R;		// = sigma_n / R_n
//				flagReversal = true;
//			}
//
//			// update kappa
//			Vector sigma_o(6);
//			Vector beta(6);
//			Vector temp(6);
//			double a, b, c;
//
//			sigma_o = R * SIGMAo;
//			beta = - R/iC*mI1;
//
//			temp = DoubleDot2_4((sigma - sigma_o), mM);
//			a = DoubleDot2_2(temp, (sigma - sigma_o));
//			b = 2* DoubleDot2_2(temp, (sigma - beta));
//			temp = DoubleDot2_4((sigma - beta), mM);
//			c = DoubleDot2_2(temp, (sigma - beta))-R*R;
//
//			kappa = (-b + sqrt((b*b - 4*a*c)))/ 2*a;
//
//			// update r
//			r = R / (1.0 + kappa);
//
//			double rOverR = 0.1;
//			if ( (r/R) < rOverR) {
//				r = rOverR * R;
//			}
//
//			// mCep = elastic modulus
//			mCe = GetElasticModulus(p, q, ev, es, n);
//			mCep = mCe;
//		
//		} else {
//			double dgamma = 0.0;
//			double Deps_vp, rho, eta, nu;
//			double denom;
//			Vector epsilonE(6);
//			Vector x(8);
//			Vector dx(8);
//			Vector Resid(8);
//
//			epsilonE = epsilonET;
//
//			mCe = GetElasticModulus(p, q, ev, es, n);
//
//			dfdSigma = DoubleDot4_2(mphi, sigmaMinusAlpha);
//
//			Deps_vp = GetVolInv((dgamma*dfdSigma));
//
//			denom = (1.0 + mTHETA * Deps_vp);
//			rho = mTHETA*R / denom;
//			eta = mTHETA * (ih * pow(kappa, im) + r) / denom;
//			nu  = mTHETA*im*ih*pow(kappa, (im-1.0)) * Deps_vp / denom;
//
//			f = DoubleDot2_2((DoubleDot2_4(sigmaMinusAlpha, mM)), sigmaMinusAlpha) - r*r;
//
//			Resid.Zero();
//			Resid(7) = f;
//
//			for (int i = 0; i < 6; i++)
//				x(i) = epsilonE(i);
//			x(6) = kappa;
//			x(7) = dgamma;
//			
//			double RESK, RES0, RESPrev;
//			RESK = 1.0;
//			RES0 = 1.0;
//			RESPrev = 1.0;
//
//			int k = 0;
//
//			// Newton loop
//			double TOL = 1e-6;
//			//while (Resid.Norm() > TOL) {
//			//while ((RESK > TOL) && (RESK <= RESPrev)) {
//			while (RESK > TOL) {
//				Matrix A(8,8);
//				Matrix Ainv(8,8);
//				Matrix Temp1(6,6);
//				Vector Temp2(6);
//				Matrix A11(6,6);
//				Vector A12(6);
//				Vector A21(6);
//				Vector A31(6);
//				//Matrix Ce(6,6);
//
//				k += 1;
//
//				Temp1 = Dyadic2_2(alpha, mI1);
//				Temp2 = SIGMAo + 1.0/iC*mI1;
//				
//				// Compute A and Ainv
//				A11 = mII + dgamma*DoubleDotStrain4_4(mphi, mCe) - dgamma*rho/R*DoubleDot4_4(mphi, Temp1);
//				A12 = -dgamma*R/(1.0+kappa)/(1+kappa) * DoubleDot4_2(mphi, Temp2);
//				A21 = (rho - (1.0+kappa)*eta)*mI1;
//				A31 = DoubleDotStrain2_4(dfdSigma, mCe) - rho/R*DoubleDot2_4(dfdSigma, Temp1) - 2.0*eta*r*mI1;
//
//				for (int i = 0; i < 6; i++) {
//					for (int j = 0; j < 6; j++) {
//						A(i,j) = A11(i,j);   // A11
//					}
//					A(i,6) = A12(i);         // A12
//					A(i,7) = dfdSigma(i);    // A13
//					A(6,i) = A21(i);         // A21
//					A(7,i) = A31(i);         // A31
//				}
//				A(6,6) = (1.0+kappa)*nu - r;  // A22
//				A(6,7) = 0.0;                 // A23
//				A(7,6) = 2.0*r*nu - r/(1+kappa) * DoubleDot2_2(dfdSigma, Temp2);  // A32
//				A(7,7) = 0.0;                 // A33
//
//				A.Invert(Ainv);
//
//				dx = -1*Ainv*Resid; // Note: should be dx = -Ainv*Resid
//				x += dx;  // Note: should be x += dx  .  Corrected for above
//
//				// update all terms
//				for (int i = 0; i < 6; i++) 
//					epsilonE(i) = x(i);
//
//				kappa  = x(6);
//				dgamma = x(7);
//
//				// update all terms 
//				ev = GetVolInv(epsilonE);
//				es = GetDevInv(epsilonE);
//
//				p = (1.0 + 3.0/2.0*ialpha/ikappa*es*es)*ip_o*exp((iepsE_vo-ev)/ikappa);
//				q = 3.0*(imu_o - ialpha*ip_o*exp((iepsE_vo-ev)/ikappa))*es;
//
//				s = epsilonE - one3*ev*mI1;
//				n = s;
//				if (Norm_EngStrain(s) == 0) 
//                    n = n/1e-15;
//                else 
//					n = n / (Norm_EngStrain(s));
//
//				sigma = p*mI1 + root23*q*n; 
//				alpha = (kappa*SIGMAo - 1.0/iC*mI1)/(1.0+kappa)*R;
//				sigmaMinusAlpha = sigma - alpha;
//
//				mCe = GetElasticModulus(p, q, ev, es, n);
//				dfdSigma = DoubleDot4_2(mphi, sigmaMinusAlpha); 
//
//				//Deps_vp = GetVolInv((dgamma*dfdSigma));
//				Deps_vp = GetVolInv(epsilonET - epsilonE);
//				denom = (1 + mTHETA * Deps_vp);
//				R = mR_n/denom;  // update R
//				r = (mr_n - mTHETA*ih*pow(kappa,im)*Deps_vp)/denom;  // update r
//
//				rho = mTHETA*R / denom;
//				eta = mTHETA * (ih * pow(kappa, im) + r) / denom;
//				nu  = mTHETA*im*ih*pow(kappa, (im-1)) * Deps_vp / denom;
//
////opserr << "A   = " << A;
////opserr << "Ainv = " << Ainv << endln; 
//opserr << "r (k)   = " << r << endln;
//opserr << "R (k)   = " << R << endln;
//opserr << "mTHETA = " << mTHETA << endln;
//opserr << "ih = " << ih << endln;
//opserr << "kappa = " << kappa << endln;
//opserr << "im = " << im << endln;
//opserr << "Deps_vp = " << Deps_vp << endln;
//opserr << "denom = " << denom << endln;
//opserr << " pow(kappa,im)  = " << pow(kappa,im)  << endln;
//opserr << "eta = " << eta << endln;
//opserr << endln;
//opserr << "dfdSigma = " << dfdSigma;
//opserr << "mCe = " << mCe;
//opserr << "DoubleDotStrain2_4(dfdSigma, mCe) = " << DoubleDotStrain2_4(dfdSigma, mCe);
//opserr << endln;
//
//				// Calculate f
//				f = DoubleDot2_2(DoubleDot2_4(sigmaMinusAlpha, mM), sigmaMinusAlpha) - r*r;
//
//				// Compute residual
//				for (int i = 0; i<6; i++) 
//					Resid(i) = epsilonE(i) - epsilonET(i) + dgamma*dfdSigma(i);
//				Resid(6) = R - (1.0+kappa)*r;
//				Resid(7) = f;
//
//				RESPrev = RESK;
//				RESK = Resid.Norm();
//				if ( k == 1)
//					RES0 = RESK;
//				//RESK = RESK/RES0;
//
//opserr << "Resid = " << Resid;
//opserr << "RESK = " << RESK << endln;
//opserr << "END OF " << k << " ITERATION on R -------------------------- " << endln;
//opserr << endln;
//			}
//opserr << "CONVERGED in " << k << " steps____________________________________________" << endln;
//
//
//			// update terms
//			mEpsilon_P = mEpsilon - epsilonE;
//			//mbeta = alpha*(1+kappa) - kappa*SIGMAo*R;
//
//
//			// Compute Cep
//			mCep = GetConsistentModulus(kappa, dgamma, rho, eta, nu, r, R, p, ev, es, SIGMAo, sigmaMinusAlpha, dfdSigma, n);
//opserr << "mCep = " << mCep;
//			flagReversal = false;
//		}
//
//
//
//		// update remaining state variables
//		mr = r;
//		mR = R;
//		mKappa = kappa;
//		mSIGMAo = SIGMAo;
//		mSigma = sigma;
//
//opserr << "END OF SetTrialStrain() *************************************************************" << endln;
//opserr << endln;
//		
//
//	return;
//}

Matrix BoundingCamClay:: GetElasticModulus(double p, double q, double ev, double es, Vector n) {
#ifdef DEBUG
        opserr << "BoundingCamClay::GetElasticModulus(...)" << endln;
#endif

	Matrix Ce(6,6);
	Matrix temp(6,6);
	double Omega;
	double C1, C2, C3;

	Omega = -(ev - iepsE_vo)/ikappa;

	C1 = 2.0*(imu_o - ialpha*ip_o*exp(Omega));

	C2 = -p/ikappa - C1 / 3.0;

	C3 = root23 * 3.0*ip_o*ialpha*es/ikappa*exp(Omega);

	temp = Dyadic2_2(mI1, n) + Dyadic2_2(n,mI1);

	Ce = C1*mII  + C2*mIIvol  +  C3*temp;

	return Ce;
}

Matrix BoundingCamClay:: 
GetConsistentModulus(double kappa, double dgamma, double rho, double eta, double nu, double r, double R, double p, double ev, double es, Vector SIGMAo, Vector sigmaMinusAlpha, Vector dfdSigma, Vector n){
#ifdef DEBUG
        opserr << "BoundingCamClay::GetConsistentModulus(...)" << endln;
#endif
			Vector alpha_R(6);
			Vector alpha_k(6);
			Vector phiAlpha_R(6);
			Vector phiAlpha_k(6);
			Matrix Q(4,4);
			Matrix Qinv(4,4);
			Matrix D(6,6);
			Matrix Cinv(6,6);
			Matrix Cep(6,6);

			double temp = 1.0/(1.0 + kappa);
			alpha_R = temp * (kappa*SIGMAo - 1.0/iC*mI1);
			alpha_k = temp*temp*R * (SIGMAo + 1.0/iC*mI1);

			dfdSigma = DoubleDot4_2(mphi, sigmaMinusAlpha);
			
			phiAlpha_R = DoubleDot4_2(mphi, alpha_R);
			phiAlpha_k = DoubleDot4_2(mphi, alpha_k);

			double I1phiAlpha_R = DoubleDot2_2(mI1, phiAlpha_R);
			double I1phiAlpha_k = DoubleDot2_2(mI1, phiAlpha_k);

			Q.Zero();
			Q(0,0) = 1.0 - dgamma*rho*I1phiAlpha_R;
			//Q(0,1) = 0;
			Q(0,2) =   -(dgamma*rho*I1phiAlpha_k);
			Q(0,3) = rho* GetVolInv(dfdSigma);

			Q(1,0) =   - dgamma*eta*I1phiAlpha_R;
			Q(1,1) = 1.0;
			Q(1,2) = nu - dgamma*eta*I1phiAlpha_k;
			Q(1,3) = eta* GetVolInv(dfdSigma);

			Q(2,0) = 1.0;
			Q(2,1) = -1.0 - kappa;
			Q(2,2) = -r;
			//Q(2,3) = 0;

			Q(3,0) = -DoubleDot2_2(dfdSigma, alpha_R);
			Q(3,1) = -2.0*r;
			Q(3,2) = -DoubleDot2_2(dfdSigma, alpha_k);
			//Q(3,3) = 0;

			Q.Invert(Qinv);

            Vector L1(6);
			Vector L2(6);
			//Vector L4(6);
			Vector U1(6);
			Vector U3(6);
			//Vector U4(6);

			L1 = dgamma*rho*mI1;
			L2 = dgamma*eta*DoubleDot2_4(mI1, mphi);
			//L3.Zero();
			//L4 = dfdSigma;

			U1 = -dgamma*phiAlpha_R;
			// U2.Zero();
			U3 = -dgamma*phiAlpha_k;
			// U4 = dfdSigma;

			D = GetCompliantModulus(p, ev, es, n);

			Cinv = D + dgamma*mphi 
				- Qinv(0,0)*Dyadic2_2(U1,L1)
				- Qinv(0,3)*Dyadic2_2(U1,dfdSigma)
				- Qinv(2,0)*Dyadic2_2(U3,L1)
				- Qinv(2,1)*Dyadic2_2(U3,L2)
				- Qinv(3,0)*Dyadic2_2(dfdSigma,L1)
				- Qinv(3,1)*Dyadic2_2(dfdSigma,L2);

			Cinv.Invert(Cep);

	return Cep;
}

Matrix BoundingCamClay:: GetCompliantModulus(double p, double ev, double es, Vector n){
#ifdef DEBUG
        opserr << "BoundingCamClay::GetCompliantModulus(...)" << endln;
#endif
	Matrix D(2,2);
	double Omega;
	double D11, D12, D22, det, C1, C2, C3, C4;

	Omega = -(ev - iepsE_vo)/ikappa;

	D11 = -p / ikappa;
	D12 = (3.0*ip_o*ialpha*es)/ikappa*exp(Omega);
	D22 = 3.0*(imu_o - ialpha*ip_o*exp(Omega));
	det = D11*D22-D12*D12;

	// E11 = D22/det;
	// E12 = E21 = -D12/det;
	// E22 = D11/det;

	C1 = 3.0/2.0/D22;                      // = 3/2 * D22^-1
	C2 = 1.0/9.0*D22/det - 0.5/D22;        // = 1/9*E11 - 0.5*D22^-1
	C3 = -D12/det/sqrt(6.0);               // = 1/root6*E12
	C4 = 3.0/2.0*(D11/det - 1/D22);        // = 3/2* (D22 - D22^-1)

	D = C1 * mII + C2*mIIvol + C3*(Dyadic2_2(mI1, n) + Dyadic2_2(n, mI1)) + C4*(Dyadic2_2(n, n));

	return D;
}
		
double
BoundingCamClay::Norm_EngStrain(Vector v) 
{
	if (v.Size() != 6)
		opserr << "\n ERROR! BoundingCamClay::NormEngStrain requires vector of size(6)!" << endln;
    double result;

	for (int i = 0; i < 3; i++) 
		result += v(i)*v(i);
	for (int i = 3; i < 6; i++)
		result += 0.5*v(i)*v(i);

  return sqrt(result);
}
										  

double BoundingCamClay:: GetVolInv(Vector v) {
	// Note: This function is only for vectors of size 6!
	if (v.Size() != 6)
		opserr << "\n ERROR! BoundingCamClay::GetVolInv requires vector of size(6)!" << endln;

	return (v(0) + v(1) + v(2));
}

double BoundingCamClay:: GetDevInv(Vector v) {
	// Note: This function is only for vectors of size 6!
	if (v.Size() != 6)
		opserr << "\n ERROR! BoundingCamClay::GetDevInv requires vector of size(6)!" << endln;

	double m, DevInv;
	Vector d(6);

	// Calculate Deviatoric Part
	d = v;
	m = one3*(v(0) + v(1) + v(2));
	for (int i = 0; i < 3; i++) 
		d(i) = d(i) - m;

	// Calculate Deviatoric Invariant
	//DevInv = root23*(d.Norm());
	DevInv = root23*(Norm_EngStrain(d));

	return DevInv;
}



double BoundingCamClay:: DoubleDot2_2(Vector v1, Vector v2){
	double result = 0.0;
	
	if (v1.Size() != v2.Size())
		opserr << "\n ERROR! BoundingCamClay::DoubleDot2_2 function requires vectors of equal size!" << endln;

	for (int i = 0; i < v1.Size(); i++) 
		result += v1(i) * v2(i);

	return result;
}

Vector BoundingCamClay:: DoubleDot2_4(Vector v1, Matrix m1){
	Vector result(6);
	result.Zero();

	if (v1.Size() != m1.noRows())
		opserr << "\n ERROR! BoundingCamClay::DoubleDot2_4 function requires Size(v1) = noRows(m1) " << endln;

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(j) += v1(i) * m1(i,j);
	}

	return result;
}

Vector BoundingCamClay:: DoubleDotStrain2_4(Vector v1, Matrix m1){
	Vector result(6);
	result.Zero();

	if (v1.Size() != m1.noRows())
		opserr << "\n ERROR! BoundingCamClay::DoubleDot2_4 function requires Size(v1) = noRows(m1) " << endln;

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) {
			if (i <=2) 
                result(j) += v1(i) * m1(i,j);
			else
                result(j) += 0.5 * v1(i) * m1(i,j);
		}
	}
	return result;
}

Vector BoundingCamClay:: DoubleDot4_2(Matrix m1, Vector v1){
	Vector result(6);
	result.Zero();

	if (m1.noCols() != v1.Size())
		opserr << "\n ERROR! BoundingCamClay::DoubleDot4_2 function requires noCols(m1) = Size(v1) " << endln;

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(i) += m1(i,j) * v1(j);
	}

	return result;
}

Matrix BoundingCamClay:: DoubleDot4_4(Matrix m1, Matrix m2){
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

Matrix BoundingCamClay:: DoubleDotStrain4_4(Matrix m1, Matrix m2){
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < m1.noRows(); i ++) {
		for (int j = 0; j < m2.noCols(); j++) {
			for (int k = 0; k<m1.noRows(); k++) {
				if ( k <=2)
                    result(i,j) += m1(i,k) * m2(k,j);
				else 
					result(i,j) += 0.5 * m1(i,k) * m2(k,j);
			}
		}
	}
    return result;
}

Matrix BoundingCamClay:: Dyadic2_2(Vector v1, Vector v2){
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}



Response*
BoundingCamClay::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getStress());
	else
		return 0;
}

int BoundingCamClay::getResponse (int responseID, Information &matInfo)
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
				*(matInfo.theVector) = getStress();
			return 0;
		default:
			return -1;
	}
}

int BoundingCamClay::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(8);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = iC;
  data(cnt++) = ip_o;
  data(cnt++) = ikappa;
  data(cnt++) = imu_o;
  data(cnt++) = ialpha;
  data(cnt++) = ilambda;
  data(cnt++) = ih;
  data(cnt++) = im;
  data(cnt++) = iepsE_vo;
// NEED TO ADD STATE INFO!


  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "BoundingCamClay::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  

  return 0;
 
}


int BoundingCamClay::recvSelf(int commitTag, Channel &theChannel, 
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
  ip_o = data(cnt++);
  ikappa = data(cnt++);
  imu_o = data(cnt++);
  ialpha = data(cnt++);
  ilambda = data(cnt++);
  ih = data(cnt++);
  im = data(cnt++);
  iepsE_vo = data(cnt++);
  // NEED TO ADD STATE INFO!!!

  return 0;

}

void BoundingCamClay::Print(OPS_Stream &s, int flag )
{
  s << "BoundingCamClay" << endln;
}

