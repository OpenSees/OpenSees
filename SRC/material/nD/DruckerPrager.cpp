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

// $Revision: 1.2 $
// $Date: 2010-02-04 20:50:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DruckerPrager.cpp,v $
                                                                      
// Written: Kathryn Petek, Peter Mackenzie-Helnwein, and Pedro Arduino
// Created: 12/04
// Modified: 06/09
//
// Description: This file contains the implementation for the DruckerPrager class.
// In its original version this file was called DruckerPragerTensionCutoff.
//				
//		This implementation of Drucker-Prager allows for non-associative flow 
//		through the use of the rho_bar parameter.  It also includes linear 
//		kinematic hardening and linear and nonlinear isotropic hardening.
//

#include <DruckerPrager.h>
#include <DruckerPrager3D.h>
#include <DruckerPragerPlaneStrain.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>

//#include <myDebug.h>
#ifdef DEBUG_LEVEL
#undef DEBUG_LEVEL
#endif

const double DruckerPrager :: one3   = 1.0 / 3.0 ;
const double DruckerPrager :: two3   = 2.0 / 3.0 ;
const double DruckerPrager :: root23 = sqrt( 2.0 / 3.0 ) ;

//  0 = elastic+no param update; 1 = elastic+param update; 2 = elastoplastic+no param update (default)
double		DruckerPrager::mElastFlag = 2;


#include <elementAPI.h>
#define OPS_Export 

static int numDruckerPragerMaterials = 0;

OPS_Export void *
OPS_NewDruckerPragerMaterial(void)
{
  if (numDruckerPragerMaterials == 0) {
    numDruckerPragerMaterials++;
    OPS_Error("DruckerPrager nDmaterial - Written: K.Petek, P.Mackenzie-Helnwein, P.Arduino, U.Washington\n", 1);
  }

  // Pointer to a uniaxial material that will be returned
  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 13) {
    opserr << "Want: nDMaterial DruckerPrager tag? massDensity? K? G? sigma_y? rho? rho_bar? Kinf? Ko? delta1? delta2? H? theta? <atm?>" << endln;
    return 0;	
  }
  
  int tag;
  double dData[13];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial DruckerPrager material  tag" << endln;
    return 0;
  }
  if (numArgs == 13)
    numData = 12;
  else
    numData = 13;

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial DruckerPrager material  with tag: " << tag << endln;
    return 0;
  }

  if (numArgs  == 13)
    theMaterial = new DruckerPrager(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
				    dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);
  else
    theMaterial = new DruckerPrager(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
				    dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);

  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial DruckerPrager material  with tag: " << tag << endln;
  }

  return theMaterial;
}



//full constructor
DruckerPrager::DruckerPrager(int tag, int classTag, double bulk, double shear, double s_y, double r,
			                                        double r_bar, double Kinfinity, double Kinit, double d1,
			                                        double d2, double H, double t, double mDen, double atm)
  : NDMaterial(tag,ND_TAG_DruckerPrager),
    mEpsilon(6), 
    mEpsilon_n_p(6),
    mEpsilon_n1_p(6),
    mSigma(6),
    mBeta_n(6),
    mBeta_n1(6),
    mCe(6,6),
    mCep(6,6),
    mI1(6),
    mIIvol(6,6),
    mIIdev(6,6),
    mState(5)
{
	massDen =  mDen;
    mKref    =  bulk;
    mGref    =  shear;
    mPatm	 =  atm;
    mK       =  bulk;
    mG       =  shear;
    msigma_y =  s_y;
    mrho     =  r;
    mrho_bar =  r_bar;
	mKinf    =  Kinfinity;
	mKo      =  Kinit;
	mdelta1  =  d1; 
	mdelta2  =  d2;
    mHard    =  H;
    mtheta   =  t;
	if (mrho == 0.0) { 
		mTo = 1e10;
	}
	else { 
		mTo      =  root23 * msigma_y / mrho;  // add something incase rho = 0
	}
	// Use these values to deactivate yield surface 1 - Create Pure Tension Cutoff
	//msigma_y = 1e10;
	//mTo      = 100;

    this->initialize();
}
   
//null constructor
DruckerPrager ::DruckerPrager  () 
    : NDMaterial(),
    mEpsilon(6), 
    mEpsilon_n_p(6),
    mEpsilon_n1_p(6),
    mSigma(6),
    mBeta_n(6),
    mBeta_n1(6),
	mCe(6,6),
	mCep(6,6),
	mI1(6),
    mIIvol(6,6),
    mIIdev(6,6),
	mState(5)
{
	massDen =  0.0;
    mKref    =  0.0;
    mGref    =  0.0;
    mPatm	 =  101.0;
    mK       =  0.0;
    mG       =  0.0;
    msigma_y =  1e+10;
    mrho     =  0.0;
    mrho_bar =  0.0;
	mKinf    =  0.0;
	mKo      =  0.0;
	mdelta1  =  0.0;
	mdelta2  =  0.0;
    mHard    =  0.0;
    mtheta   =  0.0;
	mTo      =  0.0;

    this->initialize();
}

//destructor
DruckerPrager::~DruckerPrager  ()
{
}

//zero internal variables
void DruckerPrager::initialize( )
{
    mEpsilon.Zero();
    mEpsilon_n_p.Zero();
    mEpsilon_n1_p.Zero();

    mSigma.Zero();
    mBeta_n.Zero();
    mBeta_n1.Zero();

	mAlpha1_n = 0.0;
	mAlpha1_n1 = 0.0;
	mAlpha2_n = 0.0;
	mAlpha2_n1 = 0.0;
    mFlag = 1;

 	mHprime = (1-mtheta)*mHard;

    // 2nd order Identity Tensor
    mI1.Zero();
    mI1(0) = 1;
    mI1(1) = 1;
    mI1(2) = 1;

    // 4th order Volumetric Tensor
    // IIvol = I1 tensor I1
    mIIvol.Zero();
    mIIvol(0,0) = 1;
    mIIvol(0,1) = 1;
    mIIvol(0,2) = 1;
    mIIvol(1,0) = 1;
    mIIvol(1,1) = 1;
    mIIvol(1,2) = 1;
    mIIvol(2,0) = 1;
    mIIvol(2,1) = 1;
    mIIvol(2,2) = 1;

    // 4th order Deviatoric Tensor
    // Note:  this is the contravariant form!
    //        useable for s^a = 2G * IIdev^ab * epsilon_b
    // (Need a different form for s^a = IIdev ^a_b * sigma^a)
    mIIdev.Zero();
    mIIdev(0,0) =   two3;
    mIIdev(0,1) = - one3;
    mIIdev(0,2) = - one3;
    mIIdev(1,0) = - one3;
    mIIdev(1,1) =   two3;
    mIIdev(1,2) = - one3;
    mIIdev(2,0) = - one3;
    mIIdev(2,1) = - one3;
    mIIdev(2,2) =   two3;
    mIIdev(3,3) = 0.5;
    mIIdev(4,4) = 0.5;
    mIIdev(5,5) = 0.5;

    mCe  = mK * mIIvol + 2*mG*mIIdev;
    mCep = mCe;
	mState.Zero();
}


NDMaterial * DruckerPrager::getCopy (const char *type)
{
  	if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		DruckerPragerPlaneStrain *clone;
		clone = new DruckerPragerPlaneStrain(this->getTag(), mK, mG, msigma_y, mrho, mrho_bar, mKinf, mKo,
		                                                     mdelta1, mdelta2, mHard, mtheta, massDen, mPatm);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type, "3D") ==0) {  
		DruckerPrager3D *clone;
     	clone = new DruckerPrager3D(this->getTag(),  mK, mG, msigma_y, mrho, mrho_bar, mKinf, mKo,
		                                             mdelta1, mdelta2, mHard, mtheta, massDen, mPatm);
	 	return clone;
  	} else {
	  	opserr << "DruckerPrager::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int DruckerPrager::commitState (void)
{
    mEpsilon_n_p = mEpsilon_n1_p;
    mAlpha1_n    = mAlpha1_n1; 
	mAlpha2_n    = mAlpha2_n1; 
    mBeta_n      = mBeta_n1;

    return 0;
}
 
int DruckerPrager::revertToLastCommit (void)
{
    return 0;
}

int DruckerPrager::revertToStart(void)
{
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->initialize();
	}

    return 0;
}

NDMaterial*
DruckerPrager::getCopy (void)
{
  opserr << "DruckerPrager::getCopy -- subclass responsibility\n"; 
  exit(-1);
  return 0;
}

const char*
DruckerPrager::getType (void) const
{
    opserr << "DruckerPrager::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
DruckerPrager::getOrder (void) const
{
    opserr << "DruckerPrager::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}

//--------------------Plasticity-------------------------------------
//plasticity integration routine
void DruckerPrager:: plastic_integrator( ) 
{
		bool okay;		// boolean variable to ensure satisfaction of multisurface kuhn tucker conditions
		double f1;
		double f2;
		double norm_eta;
		double Invariant_1;
		double Invariant_ep;
		double norm_ep;
		double norm_dev_ep;
		Vector epsilon_e(6);
		Vector s(6);
		Vector eta(6);
		Vector dev_ep(6);
		Vector Jact(2);

		double fTOL;
		double gTOL;
		fTOL = 0.0;
		gTOL = -1.0e-10;
		
        double NormCep;

		double alpha1;			// hardening parameter for DP surface
		double alpha2;			// hardening parameter for tension cut-off
		Vector n(6);			// normal to the yield surface in strain space
		Vector R(2);			// residual vector
		Vector gamma(2);		// vector of consistency parameters
		Vector dgamma(2);		// incremental vector of consistency parameters
		Matrix g(2,2);			// jacobian of the corner region (return map)
		Matrix g_contra(2,2);	// inverse of jacobian of the corner region

        // set trial state:

        // epsilon_n1_p_trial = ..._n1_p  = ..._n_p
        mEpsilon_n1_p = mEpsilon_n_p;

		// alpha1_n+1_trial
		mAlpha1_n1 = mAlpha1_n;
		// alpha2_n+1_trial
		mAlpha2_n1 = mAlpha2_n;

        // beta_n+1_trial
        mBeta_n1 = mBeta_n;

        // epsilon_elastic = epsilon_n+1 - epsilon_n_p
		epsilon_e = mEpsilon - mEpsilon_n1_p;

        // trial stress
		mSigma = mCe*epsilon_e;

        // deviator stress tensor: s = 2G * IIdev * epsilon_e
        //I1_trial
		Invariant_1 = ( mSigma(0) + mSigma(1) + mSigma(2) );

        // s_n+1_trial
		s = mSigma - (Invariant_1/3.0)*mI1;

        //eta_trial = s_n+1_trial - beta_n;
		eta = s - mBeta_n;
		
		// compute yield function value (contravariant norm)
        norm_eta = sqrt(eta(0)*eta(0) + eta(1)*eta(1) + eta(2)*eta(2) + 2*(eta(3)*eta(3) + eta(4)*eta(4) + eta(5)*eta(5)));

        // f1_n+1_trial
		f1 = norm_eta + mrho*Invariant_1 - root23*Kiso(mAlpha1_n1);

		// f2_n+1_trial
		f2 = Invariant_1 - T(mAlpha2_n1);
		
		// update elastic bulk and shear moduli 
 		this->updateElasticParam();

		// check trial state
		int count = 1;
		if ((f1<=fTOL) && (f2<=fTOL) || mElastFlag < 2) {

			okay = true;
			// trial state = elastic state - don't need to do any updates.
			mCep = mCe;
			count = 0;

			// set state variables for recorders
            Invariant_ep = 	mEpsilon_n1_p(0)+mEpsilon_n1_p(1)+mEpsilon_n1_p(2);

			norm_ep  = sqrt(mEpsilon_n1_p(0)*mEpsilon_n1_p(0) + mEpsilon_n1_p(1)*mEpsilon_n1_p(1) + mEpsilon_n1_p(2)*mEpsilon_n1_p(2)
                           + 0.5*(mEpsilon_n1_p(3)*mEpsilon_n1_p(3) + mEpsilon_n1_p(4)*mEpsilon_n1_p(4) + mEpsilon_n1_p(5)*mEpsilon_n1_p(5)));
			
			dev_ep = mEpsilon_n1_p - one3*Invariant_ep*mI1;

            norm_dev_ep  = sqrt(dev_ep(0)*dev_ep(0) + dev_ep(1)*dev_ep(1) + dev_ep(2)*dev_ep(2)
                           + 0.5*(dev_ep(3)*dev_ep(3) + dev_ep(4)*dev_ep(4) + dev_ep(5)*dev_ep(5)));

			mState(0) = Invariant_1;
			mState(1) = norm_eta;
       		mState(2) = Invariant_ep;
        	mState(3) = norm_dev_ep;
			mState(4) = norm_ep;
			return;
		}
		else {
			// plastic correction required
			okay = false;

			// determine number of active surfaces.  size & fill Jact
			if ( (f1 > fTOL ) && (f2 <= fTOL) ) {
				// f1 surface only
				Jact(0) = 1;
				Jact(1) = 0;
			}
			else if ( (f1 <= fTOL ) && (f2 > fTOL) ) {
				// f2 surface only
				Jact(0) = 0;
				Jact(1) = 1;
			}
			else if ( (f1 > fTOL ) && (f2 > fTOL) ) {
				// both surfaces active
				Jact(0) = 1;
				Jact(1) = 1;
			}
		} 

		//-----------------MultiSurface Placity Return Map--------------------------------------
		while (!okay) {

			alpha1 = mAlpha1_n;
			alpha2 = mAlpha2_n;
	
			//  n = eta / norm_eta;  (contravaraint)
			if (norm_eta < 1.0e-13) {
				n.Zero();
			} else {
				n = eta/norm_eta;
			}
			
			// initialize R, gamma1, gamma2, dgamma1, dgamma2 = 0
			R.Zero();  
			gamma.Zero(); 
			dgamma.Zero();
			// initialize g such that det(g) = 1
			g(0,0) = 1; 
			g(1,1) = 1;
			g(1,0) = 0;
			g(0,1) = 0;

			// Newton procedure to compute nonlinear gamma1 and gamma2
			//initialize terms
			for (int i = 0; i < 2; i++) {
				if (Jact(i) == 1) {
					R(0) = norm_eta - (2*mG + two3*mHprime)*gamma(0) + mrho*Invariant_1 
						   - 9*mK*mrho*mrho_bar*gamma(0) - 9*mK*mrho*gamma(1) - root23*Kiso(alpha1);
					g(0,0) = -2*mG - two3*(mHprime + Kisoprime(alpha1)) - 9*mK*mrho*mrho_bar;
				} else if (Jact(i) == 2) {
					R(1) = Invariant_1 - 9*mK*mrho_bar*gamma(0) - 9*mK*gamma(1) - T(alpha2);
					g(1,1) = -9*mK + mdelta2*T(alpha2);
				}
			}
			if (Jact(0) == 1 && Jact(1) == 1) {
				g(0,1) = -9*mK*mrho;
				g(1,0) = mrho_bar*(-9*mK + mdelta2*T(alpha2));
			} 
			g.Invert(g_contra);

			// iteration counter
			int m = 0;

			//iterate
			while ((fabs(R.Norm()) > 1e-10) && (m < 10)) {

				dgamma = -1*g_contra * R;
				gamma += dgamma;

				//update alpha1 and alpha2
				alpha1 = mAlpha1_n + root23*gamma(0);
				alpha2 = mAlpha2_n + mrho_bar*gamma(0) + gamma(1);

				// reset g & R matrices
				g(0,0) = 1;
				g(1,1) = 1;
				g(1,0) = 0;
				g(0,1) = 0;
				R.Zero();
				for (int i = 0; i < 2; i++) {
				if (Jact(i) == 1) {
					R(0) = norm_eta - (2*mG + two3*mHprime)*gamma(0) + mrho*Invariant_1 
						   - 9*mK*mrho*mrho_bar*gamma(0) - 9*mK*mrho*gamma(1) - root23*Kiso(alpha1);
					g(0,0) = -2*mG - two3*(mHprime + Kisoprime(alpha1)) - 9*mK*mrho*mrho_bar;
				} else if (Jact(i) == 2) {
					R(1) = Invariant_1 - 9*mK*mrho_bar*gamma(0) - 9*mK*gamma(1) - T(alpha2);
					g(1,1) = -9*mK + mdelta2*T(alpha2);
				}
				}
				if (Jact(0) == 1 && Jact(1) == 1) {
					g(0,1) = -9*mK*mrho;
					g(1,0) = mrho_bar*(-9*mK + mdelta2*T(alpha2));
				} 
				g.Invert(g_contra);

				m++;
			}

			// check maintain Kuhn-Tucker conditions
			f1 = norm_eta - (2*mG + two3*mHprime)*gamma(0) + mrho*Invariant_1 
							-9*mK*mrho*mrho_bar*gamma(0) - 9*mK*mrho*gamma(1) - root23*Kiso(alpha1);
			f2 = Invariant_1 - 9*mK*mrho_bar*gamma(0) - 9*mK*gamma(1) - T(alpha2);

			if ( count > 100 ) {
				okay = true;
                break;
			}

			// check active surfaces
			if ((Jact(0) == 1) && (Jact(1) == 0)) {
				// f2 may be > or < f2_tr because of softening of f2 related to alpha1
				if (f2 >= fTOL) {
					// okay = false;
					Jact(0) = 1;
					Jact(1) = 1;
					count += 1;

				} else {
					okay = true;

				}
			} else if ((Jact(0) == 0) && (Jact(1) == 1)) {
				// f1 will always be less than f1_tr
				okay = true;
			} else if ((Jact(0) == 1) && (Jact(1) == 1)) {
				if ((gamma(0) <= gTOL) && (gamma(1) > gTOL)){
					// okay = false;
					Jact(0) = 0;
					Jact(1) = 1;
					count += 1;
				} else if ((gamma(0) > gTOL) && (gamma(1) <= gTOL)){
					// okay = false;
					Jact(0) = 1;
					Jact(1) = 0;
					count += 1;
				} else if ((gamma(0) > gTOL) && (gamma(1) > gTOL)) {
					okay = true;
				}
			}
							
			if ( (count > 3) && (!okay) ) {
				Jact(0) = 1;
				Jact(1) = 1;
				count += 100;
			}

			if ( count > 3 ) {
				opserr << "Jact = " << Jact;
				opserr << "count = " << count << endln;
			}

		} // end of while(!okay) loop


		//update everything and exit!

		Vector b1(6);
		Vector b2(6);
		Vector n_covar(6);
		Vector temp1(6);
		Vector temp2(6);

		// update alpha1 and alpha2
		mAlpha1_n1 = alpha1;
		mAlpha2_n1 = alpha2;

		//update epsilon_n1_p
		//first calculate n_covar
		// n_a = G_ab * n^b = covariant
		n_covar(0) = n(0);
		n_covar(1) = n(1);
		n_covar(2) = n(2);
		n_covar(3) = 2*n(3);
		n_covar(4) = 2*n(4);
		n_covar(5) = 2*n(5);
		mEpsilon_n1_p = mEpsilon_n_p + (mrho_bar*gamma(0) + gamma(1))*mI1 + gamma(0)*n_covar;

           
        Invariant_ep = 	mEpsilon_n1_p(0)+mEpsilon_n1_p(1)+mEpsilon_n1_p(2);

		norm_ep  = sqrt(mEpsilon_n1_p(0)*mEpsilon_n1_p(0) + mEpsilon_n1_p(1)*mEpsilon_n1_p(1) + mEpsilon_n1_p(2)*mEpsilon_n1_p(2)
                           + 0.5*(mEpsilon_n1_p(3)*mEpsilon_n1_p(3) + mEpsilon_n1_p(4)*mEpsilon_n1_p(4) + mEpsilon_n1_p(5)*mEpsilon_n1_p(5)));
			
		dev_ep = mEpsilon_n1_p - one3*Invariant_ep*mI1;

        norm_ep  = sqrt(dev_ep(0)*dev_ep(0) + dev_ep(1)*dev_ep(1) + dev_ep(2)*dev_ep(2)
                     + 0.5*(dev_ep(3)*dev_ep(3) + dev_ep(4)*dev_ep(4) + dev_ep(5)*dev_ep(5)));

		// update sigma
		mSigma -= (3*mK*mrho_bar*gamma(0) + 3*mK*gamma(1))*mI1 + 2*mG*gamma(0)*n;

		s            -= 2*mG*gamma(0) * n;
		Invariant_1  -= 9*mK*mrho_bar*gamma(0) + 9*mK*gamma(1);
		//mSigma        = s + Invariant_1/3.0 * mI1;

		//update beta_n1
		mBeta_n1 = mBeta_n - (two3*mHprime*gamma(0))*n;
			
		// update Cep
		// note: Cep is contravariant
		if ((Jact(0) == 1) && (Jact(1) == 0)) {
			b1 = 2*mG*n + 3*mK*mrho*mI1;
			b2.Zero();
		} else if ((Jact(0) == 0) && (Jact(1) == 1)){
			b1.Zero();
			b2 = 3*mK*mI1;
		} else if ((Jact(0) == 1) && (Jact(1) == 1)){
			b1 = 2*mG*n + 3*mK*mrho*mI1;
			b2 = 3*mK*mI1;
		}

		temp1 = g_contra(0,0)*b1 + g_contra(0,1)*b2;  
		temp2 = mrho_bar*temp1 + g_contra(1,0)*b1 + g_contra(1,1)*b2;

		NormCep = 0.0;
		for (int i = 0; i < 6; i++){
			for (int j = 0; j < 6; j++) {
				mCep(i,j) = mCe(i,j)
						  + 3*mK * mI1(i)*temp2(j)  
						  + 2*mG * n(i)*temp1(j)
						  - 4*mG*mG/norm_eta*gamma(0) * (mIIdev(i,j) - n(i)*n(j));
				NormCep += mCep(i,j)*mCep(i,j);
			}
		}

		if ( NormCep < 1e-10){
			mCep = 1.0e-3 * mCe;
			opserr << "NormCep = " << NormCep << endln;
		}

		mState(0) = Invariant_1;
		mState(1) = norm_eta;
        mState(2) = Invariant_ep;
        mState(3) = norm_dev_ep;
		mState(4) = norm_ep;

	return;
}

int DruckerPrager::updateElasticParam( )
{
    double Sigma_mean = 0.0;
	if ( mElastFlag == 1 && mFlag == 1){
		Sigma_mean = -one3*(mSigma(0)+mSigma(1)+mSigma(2));
        if (Sigma_mean < 0.0) Sigma_mean = 0.0;  // prevents modulus update for cases where tension exists 
        mK = mKref * pow(1+(Sigma_mean/mPatm), 0.5);
        mG = mGref * pow(1+(Sigma_mean/mPatm), 0.5);
   		mCe  = mK * mIIvol + 2*mG*mIIdev;
        mFlag = 0;
		//opserr << "Plastic Integrator -->" << "K = " << mK  << "  G =" << mG << endln;
	} else if ( mElastFlag != 1 ){
        mFlag = 1;
	}

	return 0;
}

double DruckerPrager::Kiso(double alpha1)
{
	return msigma_y + mtheta * mHard * alpha1 + (mKinf - mKo) * (1 - exp(-mdelta1 * alpha1));
}


double DruckerPrager::Kisoprime(double alpha1)
{
	return mtheta * mHard + (mKinf - mKo) * mdelta1*  exp(-mdelta1 * alpha1);
}

double DruckerPrager:: T(double alpha2) 
{
}


double DruckerPrager::deltaH(double dGamma)
{
	return mHprime * root23 * dGamma;
}


Vector DruckerPrager::getState()
{
	return mState;
}


Response*
DruckerPrager::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse =0;
	const char *matType = this->getType();

	output.tag("NdMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());

	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int DruckerPrager::getResponse (int responseID, Information &matInfo)
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
				*(matInfo.theVector) = getState();
			return 0;
		default:
			return -1;
	}
}

int DruckerPrager::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 2)
    	return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0],"updateMaterialStage") == 0) {
			return param.addObject(1, this);
		} else if (strcmp(argv[0],"shearModulus") == 0) {
    		return param.addObject(10, this);
		} else if (strcmp(argv[0],"bulkModulus") == 0) {
    		return param.addObject(11, this);
		}
	}

    return -1;
}

int
DruckerPrager::updateParameter(int responseID, Information &info)
{
	// updateMaterialStage called
	if (responseID == 1) {
		mElastFlag = info.theInt;
	}
    // materialState called
	if (responseID == 5) {
		mElastFlag = info.theDouble;
	}
	if (responseID == 10) {
    	mElastFlag = info.theDouble;
	}
	if (responseID==11) {
    	mElastFlag=info.theDouble;
  	}

	return 0;
}

int DruckerPrager::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state
  // int a vector object
  static Vector data(8);
  int cnt = 0;
  data(cnt++) = this->getTag();
  data(cnt++) = mK;
  data(cnt++) = mG;
  data(cnt++) = mrho;
  data(cnt++) = mrho_bar;
  data(cnt++) = mHard;
  data(cnt++) = mtheta;
// NEED TO ADD STATE INFO!


  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "DruckerPrager::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  
  return 0;
 
}

int DruckerPrager::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(7);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "DruckerPrager::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag( int(data(cnt++)) );
  mK = data(cnt++);
  mG = data(cnt++);
  mrho = data(cnt++);
  mrho_bar = data(cnt++);
  mHard = data(cnt++);
  mtheta = data(cnt++);

  return 0;

}

void DruckerPrager::Print(OPS_Stream &s, int flag )
{
  s << "DruckerPrager" << endln;
}

