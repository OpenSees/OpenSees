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
                                                                        
//----------------------------------------------------------------------------------------------------------------------------
// Developed by:
// Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
// Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
// Adam M. Knaack (adam.knaack@schaefer-inc.com) 
// Schaefer-Inc, Cincinnati, Ohio, USA
// Yahya C. Kurama (ykurama@nd.edu)
// Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Created: 2019
// Last updated: 2019
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Description: This file contains the source code of TDConcreteMC10NL. 
// TDConcreteMC10NL is a time-dependent concrete material model that calculates
// creep and shrinkage strains.
/*-------------------------------
! Concrete Compression - Non-linear with linear softening
! Concrete Tension - Tamai, S., Shima, H., Izumo, J., Okamura, H. 1988. Average Stress-Strain Relationship in Post Yield Range of Steel Bar in Concrete, Concrete Library of JSCE, No. 11, 117-129.
! Concrete Creep - Linear superposition of creep coefficient, Model Code 2010 time function
! Concrete Shrinkage - Model Code 2010 time function
-------------------------------*/
// Detailed descriptions of the model and its implementation can be found in the following:
// (1) Knaack, A.M., Kurama, Y.C. 2018. Modeling Time-Dependent Deformations: Application for Reinforced Concrete Beams with 
//     Recycled Concrete Aggregates. ACI Structural J. 115, 175�190. doi:10.14359/51701153
// (2) Knaack, A.M., 2013. Sustainable concrete structures using recycled concrete aggregate: short-term and long-term behavior
//     considering material variability. PhD Dissertation, Civil and Environmental Engineering and Earth Sciences, University of Notre Dame, Notre Dame, Indiana, USA, 680 pp.
// A manual describing the use of the model and sample files can be found at:
// <https://data.mendeley.com/datasets/z4gxnhchky/3>
//----------------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------------------------------------------------------
// Disclaimer: This software is provided �as is�, without any warranties, expressed or implied. In no event shall the developers be liable for any claim, damages, or liability arising from or in connection with this software.
//----------------------------------------------------------------------------------------------------------------------------



#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>


#include "TDConcreteMC10NL.h"
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <iostream> //Added by AMK for debugging
#include <elementAPI.h> //Added by AMK to use methods for parsing data line;
#include <Domain.h> //Added by AMK to get current Domain time;
using namespace std; //Added by AMK for debugging
#include <MaterialResponse.h>
#include <Vector.h>


//Added by AMK to use dylib:
//-----------------------------------------------------------------------
	#ifdef _USRDLL
	#define OPS_Export extern "C" _declspec(dllexport)
	#elif _MACOSX
	#define OPS_Export extern "C" __attribute__((visibility("default")))
	#else
	#define OPS_Export extern "C"
	#endif

	static int numTDConcreteMC10NL = 0;

//	OPS_Export void * //ntosic: eliminated AMK code
	void * //ntosic: new code over AMK
	OPS_TDConcreteMC10NL() {
		// Print description of material model:
		if (numTDConcreteMC10NL == 0) {
			opserr << "Time-Dependent Concrete Material Model - Written by Nikola Tosic, 2019 \n";
			numTDConcreteMC10NL = 1;
		}

		// Pointer to a uniaxial material that will be returned:
			UniaxialMaterial *theMaterial = 0;
		
		// Parse the input line for the material parameters:
			int iData;
			int numData;
			int numArgs;
		
			numArgs = OPS_GetNumRemainingInputArgs();
			//ntosic
			if (numArgs == 19) {
				//TDConcreteMC10NL(int tag, double _fc, double _epsc0, double _fcu,
				//double _epscu, double _tcr, double _ft, double _Ets, double _Ec, double _age, double _epsshu)
				double dData[18];
			
				//Collect material tag:
				numData = 1;
				if (OPS_GetIntInput(&numData, &iData) != 0) {
					opserr << "WARNING: invalid uniaxialMaterial TDConcreteMC10NL tag\n";
					return 0;
				}
			
				//Collect input parameters: 
				numData = 18; //ntosic
				if (OPS_GetDoubleInput(&numData, dData) != 0) {
					opserr << "WARNING: invalid material property definition\n";
					return 0;
				}
			
				//Create a new materiadouble 
				//ntosic
				theMaterial = new TDConcreteMC10NL(iData,dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],dData[7],dData[8],dData[9],dData[10],dData[11], dData[12], dData[13], dData[14], dData[15], dData[16], dData[17]);
                if (theMaterial == 0) {
					opserr << "WARNING: could not create uniaxialMaterial of type TDConcreteMC10NL \n";
					return 0;
				}
			
				//Return new material:
				return theMaterial;
			}
	}

//-----------------------------------------------------------------------


TDConcreteMC10NL::TDConcreteMC10NL(int tag, double _fc, double _fcu, double _epscu, double _ft, double _Ec, double _Ecm, double _beta, double _age, double _epsba, double _epsbb, double _epsda, double _epsdb, double _phiba, double _phibb, double _phida, double _phidb, double _tcast, double _cem): 
  UniaxialMaterial(tag, 0), //Changed by ntosic (from TDConcrete by AMK)
  fc(_fc), fcu(_fcu), epscu(_epscu), ft(_ft), Ec(_Ec), Ecm(_Ecm), beta(_beta), age(_age), epsba(_epsba), epsbb(_epsbb), epsda(_epsda), epsdb(_epsdb), phiba(_phiba), phibb(_phibb), phida(_phida), phidb(_phidb), tcast(_tcast), cem(_cem)
{
  ecminP = 0.0;
  ecmaxP = 0.0; //ntosic
  deptP = 0.0;

	sigCr = fabs(sigCr); //ntosic: CHANGE?
  eP = Ec; //Added by AMK
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
	e = Ec; //Added by AMK
	Et = Ec;
	count = 0; //Added by AMK
	epsInit = 0.0; //Added by AMK
	sigInit = 0.0; //Added by AMK
	eps_total = 0.0; //Added by AMK
	epsP_total = 0.0; //Added by AMK
	//ntosic: SPLIT INTO BASIC AND DRYING!
	eps_m = 0.0; //Added by AMK
	eps_crb = 0.0; //Added by ntosic
	eps_crd = 0.0; //Added by ntosic
	eps_shb = 0.0; //Added by ntosic
	eps_shd = 0.0; //Added by ntosic
	epsP_crb = 0.0; //Added by ntosic
	epsP_crd = 0.0; //Added by ntosic
	epsP_shb = 0.0; //Added by ntosic
	epsP_shd = 0.0; //Added by ntosic
	epsP_m = 0.0; //Added by AMK
	
	t_load = -1.0; //Added by AMK
	crack_flag = 0;
    iter = 0;
	
	
	
	//Change inputs into the proper sign convention: ntosic: changed
		fc = -1.0*fabs(fc);
        fcu = -1.0*fabs(fcu);
        epscu = -1.0*fabs(epscu);		
		epsba = -1.0*fabs(epsba);
		epsda = -1.0*fabs(epsda);
		phiba = 1.0*fabs(phiba);
		phida = 1.0*fabs(phida);
}

TDConcreteMC10NL::TDConcreteMC10NL(void):
  UniaxialMaterial(0, 0) //Changed by AMK (from MAT_TAG_Concrete02 to 0)
{
 
}

TDConcreteMC10NL::~TDConcreteMC10NL(void)
{
  // Does nothing
}

UniaxialMaterial*
TDConcreteMC10NL::getCopy(void)
{
  TDConcreteMC10NL *theCopy = new TDConcreteMC10NL(this->getTag(), fc, fcu, epscu, ft, Ec, Ecm, beta, age, epsba, epsbb, epsda, epsdb, phiba, phibb, phida, phidb, tcast, cem); //ntosic
  
  return theCopy;
}

double
TDConcreteMC10NL::getInitialTangent(void)
{
	return Ec; //Added by AMK
}

double
TDConcreteMC10NL::getCurrentTime(void)
{
	double currentTime;
	Domain * theDomain = ops_TheActiveDomain;

	if (theDomain != 0) {
		currentTime = theDomain->getCurrentTime();
	}
    
    return currentTime;
}	
//ntosic
double
TDConcreteMC10NL::setCreepBasicStrain(double time, double stress)
{
    double creepBasic;
    double runSum = 0.0;
    
    DTIME_i[count] = ops_Dt;
 
	for (int i = 1; i<=count; i++) {
                PHIB_i[i] = setPhiBasic(time,TIME_i[i]); //Determine PHI //ntosic: PHIB
                runSum += PHIB_i[i]*DSIG_i[i]/Ecm; //CONSTANT STRESS within Time interval //ntosic: changed to Ecm from Ec (according to Model Code formulation of phi basic)
    }
    
    phib_i = PHIB_i[count];
    creepBasic = runSum;
    return creepBasic;
    
}
//ntosic
double
TDConcreteMC10NL::setCreepDryingStrain(double time, double stress)
{
	double creepDrying;
	double runSum = 0.0;

	DTIME_i[count] = ops_Dt;

	for (int i = 1; i <= count; i++) {
		PHID_i[i] = setPhiDrying(time, TIME_i[i]); //Determine PHI //ntosic: PHID
		runSum += PHID_i[i] * DSIG_i[i] / Ecm; //CONSTANT STRESS within Time interval //ntosic: changed to Ecm from Ec (according to Model Code formulation of phi drying)
	}

	phid_i = PHID_i[count];
	creepDrying = runSum;
	return creepDrying;

}
//ntosic
double
TDConcreteMC10NL::setPhiBasic(double time, double tp)
{
	// ntosic: Model Code 2010 Equations
	double tmtp = time - tp;
	double tpa = tp * pow(9.0 / (2.0 + pow(tp, 1.2)) + 1.0, cem);
	double phiBasic = phiba * log(pow(30.0 / tpa + 0.035, 2.0) * (tmtp / phibb) + 1.0);
	return phiBasic;
}
//ntosic
double
TDConcreteMC10NL::setPhiDrying(double time, double tp)
{
	// ntosic: Model Code 2010 Equations
	double tmtp = time - tp;
	double tpa = tp * pow(9.0 / (2.0 + pow(tp, 1.2)) + 1.0, cem);
	double phiDrying = phida / (0.1 + pow(tpa, 0.2)) * pow(tmtp, 1.0 / (2.3 + 3.5 / pow(tpa, 0.5))) / pow(phidb + tmtp, 1.0 / (2.3 + 3.5 / pow(tpa, 0.5)));
	return phiDrying;
}
//ntosic
double 
TDConcreteMC10NL::setShrinkBasic(double time)
{
    double shrinkBasic = epsba * (1 - exp(-0.2* epsbb * pow(time, 0.5))); //ntosic: Model Code 2010 Equations
	return shrinkBasic;
}

//ntosic
double
TDConcreteMC10NL::setShrinkDrying(double time)
{
	double tD = age; //Age at initiation of drying
	double shrinkDrying = 0.0; //ntosic
	if (time - (tD) < 0) {
		shrinkDrying = 0.0;
	}
	else {
		shrinkDrying = epsda * pow(time - (tD), 0.5) / pow(epsdb + time - (tD), 0.5); //ntosic: Model Code 2010 Equations
	}
	return shrinkDrying;
}

int
TDConcreteMC10NL::setTrialStrain(double trialStrain, double strainRate)
{
	double t = getCurrentTime();
    double tol = 1.e-4; // 9/13
    double test = 10.0; // 9/13
    double sigI = 0.0;  // 9/13
    int niter = 500;  // 9/13
	
    //opserr<<"\n trialStrain = "<<trialStrain;
    
	// Need to initialize count and initial stress/strain:
    /*
    if (ops_Creep == 0) {
		count = 0;
	} else if (ops_Creep==1 && count == 0){
		count = 1;
    }
    */
	
	// Check casting age:
	if (t-tcast<(2.0-0.0001)) { //Assumed that concrete can only carry load once hardened at 2 days following casting
		eps_crb = 0.0; //ntosic
		eps_crd = 0.0; //ntosic
		eps_shb = 0.0; //ntosic
		eps_shd = 0.0; //ntosic
		eps_m = 0.0;
		eps_total = trialStrain;
		sig = 0.0;
	} else { // Concrete has hardened and is ready to accept load
		// Initialize total strain:
        	eps_total = trialStrain;
	
		// Calculate shrinkage Strain:
            if (iter < 1) {
                eps_shb = setShrinkBasic(t); //ntosic
				eps_shd = setShrinkDrying(t); //ntosic
            }

    	// Calculate creep and mechanical strain, assuming stress remains constant in a time step:
    	if (ops_Creep == 1) {
        	if (fabs(t-TIME_i[count]) <= 0.0001) { //If t = t(i-1), use creep/shrinkage from last calculated time step
            	eps_crb = epsP_crb; //ntosic
				eps_crd = epsP_crd; //ntosic
            	eps_shb = epsP_shb; //ntosic
				eps_shd = epsP_shd; //ntosic
            	eps_m = eps_total - eps_crb - eps_crd - eps_shb - eps_shd; //ntosic
            	sig = setStress(eps_m, e);
            
        	} else { // if the current calculation is a new time step
        		//if (crackP_flag == 1 && sigP >= 0.0) { //if cracking occurs and previous stress is positive, 
        		// creep strain should not increase
        		//	eps_cr = epsP_cr;
        		//} else {
        		//	eps_cr = setCreepStrain(t,sig);
        		//}
        		//if (t < tcast) {
        		//cout << "\nWARNING: TDConcrete loaded before tcast, creep and shrinkage not calculated";
        		//	eps_sh = epsP_sh;
        		//	eps_cr = epsP_cr;
        		//	eps_m = eps_total - eps_cr - eps_sh;
        		//	sig = setStress(eps_m, e);
        		//} else {
				if (iter < 1) {
                    eps_crb = setCreepBasicStrain(t,sig); 
					eps_crd = setCreepDryingStrain(t, sig);
				}
        		eps_m = eps_total - eps_crb - eps_crd - eps_shb - eps_shd; //ntosic
        		sig = setStress(eps_m, e);
				//}
        	}
    	} else { //Static Analysis using previously converged time-dependent strains
        	    eps_crb = epsP_crb; //ntosic
				eps_crd = epsP_crd; //ntosic
            	eps_shb = epsP_shb; //ntosic
				eps_shd = epsP_shd; //ntosic
            	eps_m = eps_total - eps_crb - eps_crd - eps_shb - eps_shd; //ntosic
    	        sig = setStress(eps_m, e);
    	}
		//
		//opserr<<"\n   eps_cr = "<<eps_cr;
		//opserr<<"\n   eps_sh = "<<eps_sh;
		//opserr<<"\n   eps_m = "<<eps_m;
		//opserr<<"\n   sig = "<<sig;
	}
    iter ++;
	return 0;
}

double
TDConcreteMC10NL::setStress(double strain, double &stiff)
{
// Determine proper load path (comp load, comp unload, tens load, tens unload):
    double stress=0.0;
    crack_flag = crackP_flag;
    ecmin = ecminP; //Initialized as ecmin = 0; ecmin should never be positive
    ecmax = ecmaxP; //Initialized as ecmax = 0; ecmax should never be negative
    
    if (strain <= ecmin) { // Concrete in compression loading
        this->Compr_Envlp(strain,stress,stiff);
        ecmin = strain;			// reset ecmin
        crack_flag = 0;			// concrete in compression, no cracking
    } else { // Concrete in either: Comp Unload, Tens Load, or Tens Unload/reload
    	if (strain < 0.0) { // Compression Unloading
    		//stiff = Ec;
    		//stress = strain * stiff;
    		this->Compr_Envlp(strain,stress,stiff);
    	} else { // either Tens Load, Tens Unload, or Tens reload
    		double et0 = ft/Ec;
    		if (strain >= ecmax) { //Tens Load or reload if strain is larger than before
    		//Need to check whether cracking has occured or not
    		//If cracked, then reloading occurs along Et
    		//If not cracked, then loading occurs according to Tens_Envlp
    			ecmax = strain; // reset ecmax
    			this->Tens_Envlp(strain, stress, stiff);
    			if (strain >= et0) {//cracking has occured, set cracking flag
    				crack_flag = 1;
    			}
    		} else { //Tens Unload or Tens Reload
    			if (strain<=et0 && ecmax<=et0) { //Linear unloading/reloading, i.e, not cracked
    					this->Tens_Envlp(strain,stress,stiff);
    			} else { // Nonlinear unloading/reloading, i.e., cracked
    				stress = Et*strain;
    				stiff = Et;
    			}
    		}
    	}
    }
    return stress;
}

double
TDConcreteMC10NL::getStrain(void)
{
	return eps_total; //Added by AMK
  //return eps;
}
//ntosic
double
TDConcreteMC10NL::getPHIB_i(void)
{
	return phib_i;
}
//ntosic
double
TDConcreteMC10NL::getPHID_i(void)
{
	return phid_i;
}

double 
TDConcreteMC10NL::getStress(void)
{
	return sig;
}

double 
TDConcreteMC10NL::getTangent(void)
{
	return e;
}
//ntosic
double
TDConcreteMC10NL::getCreepBasic(void)
{
	return eps_crb;
}
//ntosic
double
TDConcreteMC10NL::getCreepDrying(void)
{
	return eps_crd;
}
//ntosic
double
TDConcreteMC10NL::getShrinkBasic(void)
{
	return eps_shb;
}
//ntosic
double
TDConcreteMC10NL::getShrinkDrying(void)
{
	return eps_shd;
}

double
TDConcreteMC10NL::getMech(void)
{
	return eps_m;
}

int 
TDConcreteMC10NL::commitState(void)
{
  iter = 0;
  ecminP = ecmin;
  ecmaxP = ecmax;
  deptP = dept;
  
  dsig_i[count]=sig-sigP;
  /* 5/8/2013: commented the following lines so that the DSIG_i[count+1]=sig-sigP;*/
  //if (crack_flag == 1) {// DSIG_i will be different depending on how the fiber is cracked
  //	if (sig < 0 && sigP > 0) { //if current step puts concrete from tension to compression, DSIG_i will be only the comp. stress
  //		DSIG_i[count+1] = sig;
  //	}
  //	if (sig > 0) {// Concrete should not creep when crack is opened
  //		DSIG_i[count+1] = 0.0;
  //	}
  //	if (sig > 0 && sigP < 0) {//if current step goes from compression to tension, DSIG_i will be the stress difference
  //		DSIG_i[count+1] = sig-sigP;
  //	}
  //} else { //concrete is uncracked, DSIG = sig - sigP
  //	DSIG_i[count+1] = sig-sigP;
  //}
  DSIG_i[count+1] = sig-sigP;
  
  //Secant Stiffness for determination of creep strain:
      if (fabs(eps_m/sig)>Ec) {  //ntosic: originaly was eps_m/sig
          E_i[count+1] = Ec;
      } else {
          E_i[count+1] = fabs(sig/eps_m); //ADDED 7/22
      }
  
      if (isnan(E_i[count+1])) {
          E_i[count+1] = Ec;
      }
  
  
  TIME_i[count+1] = getCurrentTime();
    
  eP = e;
  sigP = sig;
  epsP = eps;
	
 //Added by AMK:
	epsP_total = eps_total; //Added by AMK;
    epsP_shb = eps_shb; //ntosic
	epsP_shd = eps_shd; //ntosic
	epsP_crb = eps_crb; //ntosic
	epsP_crd = eps_crd; //ntosic
	epsP_m = eps_m;
	//ntosic: strain compression limit changed to 0.4fc/Ec; Include nonlinear creep coefficient?
    if (eps_m < 0 && fabs(eps_m)>0.40*fabs(fc/Ec)) {
        double s = fabs(eps_m/fc)*Ec;
		cout<<"\n          Strain Compression Limit Exceeded: "<<s<<"fc'";
    }
		//Cracking flags:
		crackP_flag = crack_flag;
		
	//cracked reloading/unloading stiffness:
	if (crackP_flag==1) {
		if (sig/eps_m<Et) {
			Et = sig/eps_m;
		}
	}
	
	if (count==0) {
		epsInit = epsP_total;
		sigInit = sigP;
	}
	
	if (sigInit<0.0 && t_load<0.0) {
		t_load = getCurrentTime();
		sigInit = sigP;
		epsInit = epsP_m;
	} else if (sigInit>0.0 && sigP<0.0 && t_load<0.0) {
		t_load = getCurrentTime();
		sigInit = sigP;
		epsInit = epsP_m;
	}
	
	//if (ops_Creep==1) {
	//	count++;
	//}
	count++;
	
  return 0;
}

int 
TDConcreteMC10NL::revertToLastCommit(void)
{
  	eps_total = epsP_total; //Added by AMK;
    eps_shb = epsP_shb; //ntosic
	eps_shd = epsP_shd; //ntosic
	eps_crb = epsP_crb; //ntosic
	eps_crd = epsP_crd; //ntosic
	eps_m = epsP_m;  
    
  ecmin = ecminP;
  ecmax = ecmaxP; //ntosic
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int 
TDConcreteMC10NL::revertToStart(void)
{
  ecminP = 0.0;
  ecmaxP = 0.0; //ntosic
  deptP = 0.0;

	eP = Ec;
  epsP = 0.0;
  sigP = 0.0;
  eps = 0.0;
  sig = 0.0;
	e = Ec;
  
	if (ops_Creep==0) {
		count = 0;
	} else {
		count = 1;
	}
	
  return 0;
}

int 
TDConcreteMC10NL::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(24); //ntosic
  data(0) =fc;
  data(1) =fcu;
  data(2) = epscu;
  data(3) =ft;    
  data(4) =Ec;
  data(5) =Ecm;  //ntosic 
  data(6) =beta;   
  data(7) =age; 
  data(8) =epsba; //ntosic   
  data(9) =epsbb; //ntosic   
  data(10) =epsda; //ntosic
  data(11) =epsdb; //ntosic  
  data(12) =phiba; //ntosic
  data(13) =phibb; //ntosic
  data(14) =phida; //ntosic
  data(15) =phidb; //ntosic
  data(16) =cem; //ntosic
  data(17) = ecminP; //ntosic
  data(18) = ecmaxP; //ntosic
  data(19) = deptP; //ntosic
  data(20) = epsP; //ntosic
  data(21) = sigP; //ntosic
  data(22) = eP; //ntosic
  data(23) = this->getTag();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "TDConcreteMC10NL::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
TDConcreteMC10NL::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(24); //ntosic

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "TDConcreteMC10NL::recvSelf() - failed to recvSelf\n";
    return -1;
  }
  
  fc = data(0); //ntosic
  fcu = data(1); //ntosic
  epscu = data(2); //ntosic
  ft = data(3);   
  Ec = data(4);
  Ecm = data(5); //ntosic
  beta = data(6);   
  age = data(7); 
  epsba = data(8);  //ntosic 
  epsbb = data(9);  //ntosic 
  epsda = data(10); //ntosic
  epsdb = data(11); //ntosic     
  phiba = data(12);  //ntosic
  phibb = data(13);  //ntosic 
  phida = data(14); //ntosic
  phidb = data(15); //ntosic
  cem = data(16); //ntosic
  ecminP = data(17); //ntosic
  ecmaxP = data(18); //ntosic
  deptP = data(19); //ntosic
  epsP = data(20); //ntosic
  sigP = data(21); //ntosic
  eP = data(22); //ntosic
  this->setTag(data(23));

  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
TDConcreteMC10NL::Print(OPS_Stream &s, int flag)
{
  s << "TDConcreteMC10NL:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}




void
TDConcreteMC10NL::Tens_Envlp (double epsc, double &sigc, double &Ect)
{
/*-----------------------------------------------------------------------
! monotonic envelope of concrete in tension (positive envelope)
!
!   ft    = concrete tensile strength
!   Ec0   = initial tangent modulus of concrete 
!   Ets   = tension softening modulus
!   eps   = strain
!
!   returned variables
!    sigc  = stress corresponding to eps
!    Ect  = tangent concrete modulus
!-----------------------------------------------------------------------*/
  
	double Ec0 = Ec;
	double eps0 = ft / Ec0;
	double epsu = ft * (1.0 / Ets + 1.0 / Ec0);
	double b = beta;
	// USE THIS ONE
	if (epsc <= eps0) {
		sigc = epsc * Ec0;
		Ect = Ec0;
	}
	else {
		Ect = -b * eps0*ft / pow(epsc, 2)*pow(eps0 / epsc, b - 1.0); //original
		sigc = ft * pow(eps0 / epsc, b);
	}
   
  //THiS IS FOR TESTING LINEAR
  //sigc = epsc*Ec0;
  //Ect = Ec0;
   /*
    if (epsc<=epsu) {
      Ect  = -Ets;
      sigc = ft-Ets*(epsc-eps0);
    } else {
      Ect  = 1.0e-10;
      sigc = 1.0e-10;
    }
    */
    
  return;
}

//ntosic  
void
TDConcreteMC10NL::Compr_Envlp(double epsc, double &sigc, double &Ect)
{
	//Linear
	//Ect = Ec;
	//sigc = Ect*epsc;
	//Non-linear with linear softening
	/*-----------------------------------------------------------------------
	! monotonic envelope of concrete in compression (negative envelope)
	!
	!   fc    = concrete compressive strength
	!   fcu   = stress at ultimate (crushing) strain
	!   epscu = ultimate (crushing) strain
	!   Ec0   = initial concrete tangent modulus
	!   epsc  = strain
	!
	!   returned variables
	!   sigc  = current stress
	!   Ect   = tangent concrete modulus
	-----------------------------------------------------------------------*/

	double Ec0 = Ec; //ntosic
	double epsc0 = 2.0*fc / Ec0; //ntosic

	double ratLocal = epsc / epsc0;
	if (epsc >= epsc0) {
		sigc = fc * ratLocal*(2.0 - ratLocal);
		Ect = Ec0 * (1.0 - ratLocal);
	}
	else {

		//   linear descending branch between epsc0 and epscu
		if (epsc > epscu) {
			sigc = (fcu - fc)*(epsc - epsc0) / (epscu - epsc0) + fc;
			Ect = (fcu - fc) / (epscu - epsc0);
		}
		else {

			// flat friction branch for strains larger than epscu

			sigc = fcu;
			Ect = 1.0e-10;
			//       Ect  = 0.0  
			return;
		}
	}
}

int
TDConcreteMC10NL::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

/* Methods added by AMK: */

Response* 
TDConcreteMC10NL::setResponse(const char **argv, int argc,
							  OPS_Stream &theOutput)
{	
	Response *theResponse = 0;
	
	theOutput.tag("UniaxialMaterialOutput");
	theOutput.attr("matType", this->getClassType());
	theOutput.attr("matTag", this->getTag());
	
	// stress
	if (strcmp(argv[0],"stress") == 0) {
		theOutput.tag("ResponseType", "sigma11");
		theResponse =  new MaterialResponse(this, 1, this->getStress());
	}  
	// tangent
	else if (strcmp(argv[0],"tangent") == 0) {
		theOutput.tag("ResponseType", "C11");
		theResponse =  new MaterialResponse(this, 2, this->getTangent());
	}
	
	// strain
	else if (strcmp(argv[0],"strain") == 0) {
		theOutput.tag("ResponseType", "eps11");
		theResponse =  new MaterialResponse(this, 3, this->getStrain());
	}
	
	// strain
	else if ((strcmp(argv[0],"stressStrain") == 0) || 
			 (strcmp(argv[0],"stressANDstrain") == 0) ||
			 (strcmp(argv[0],"stressAndStrain") == 0)) {
		theOutput.tag("ResponseType", "sig11");
		theOutput.tag("ResponseType", "eps11");
		theResponse =  new MaterialResponse(this, 4, Vector(2));
	}
	
	else if (strcmp(argv[0],"CreepStressStrainTangent")==0) {
		theOutput.tag("ResponseType", "sig11");
		theOutput.tag("ResponseType", "eps11");
		theOutput.tag("ResponseType", "C11");
		theOutput.tag("ResponseType", "CreepBasicStrain"); //ntosic
		theOutput.tag("ResponseType", "CreepDryingStrain"); //ntosic
		theOutput.tag("ResponseType", "MechStrain");
		theOutput.tag("ResponseType", "ShrinkBasicStrain"); //ntosic
		theOutput.tag("ResponseType", "ShrinkDryingStrain"); //ntosic
		theOutput.tag("ResponseType", "t_load");
		theResponse = new MaterialResponse(this, 6, Vector(8));
	}
	
	else if ((strcmp(argv[0],"stressStrainTangent") == 0) || 
			 (strcmp(argv[0],"stressANDstrainANDtangent") == 0)) {
		theOutput.tag("ResponseType", "sig11");
		theOutput.tag("ResponseType", "eps11");
		theOutput.tag("ResponseType", "C11");
		theResponse =  new MaterialResponse(this, 5, Vector(3));
	}
	
	// stress sensitivity for local sensitivity recorder purpose.  Quan 2009
	// limit:  no more than 10000 random variables/sensitivity parameters
	else if (strstr(argv[0],"stressSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradient = atoi(token);
		theOutput.tag("ResponseType", "sigsens11");
		theResponse =  new MaterialResponse(this, gradient+10000, this->getStress());
	}
	// strain sensivitiy
	else if (strstr(argv[0],"strainSensitivity") != 0) {
		char *token = strtok((char *) argv[0], " ");
		if (token != NULL) token = strtok(NULL, " ");
		int gradient = atoi(token);
		theOutput.tag("ResponseType", "epssens11");
		theResponse =  new MaterialResponse(this, gradient+20000, this->getStrain());
	}
	
	
	theOutput.endTag();
	return theResponse;
	
}

int 
TDConcreteMC10NL::getResponse(int responseID, Information &matInfo)
{
	static Vector stressStrain(2);
	static Vector stressStrainTangent(3);
	static Vector CreepStressStrainTangent(8); //Added by AMK //ntosic: modified number
	// each subclass must implement its own stuff   
	
	// added for sensitivity recorder. Quan 2009
	if ((responseID>10000)&&(responseID<20000)){
		matInfo.setDouble(this->getStressSensitivity(responseID-10000,false));
		return 0;
	}
	else if (responseID>20000){
		matInfo.setDouble(this->getStrainSensitivity(responseID-20000));
		return 0;
	}
	
	switch (responseID) {
		case 1:
			matInfo.setDouble(this->getStress());
			return 0;
			
		case 2:
			matInfo.setDouble(this->getTangent());
			return 0;      
			
		case 3:
			matInfo.setDouble(this->getStrain());
			return 0;      
			
		case 4:
			stressStrain(0) = this->getStress();
			stressStrain(1) = this->getStrain();
			matInfo.setVector(stressStrain);
			return 0;
			
		case 5:
			stressStrainTangent(0) = this->getStress();
			stressStrainTangent(1) = this->getStrain();
			stressStrainTangent(2) = this->getTangent();
			matInfo.setVector(stressStrainTangent);
			return 0;
		
		case 6:
			CreepStressStrainTangent(0) = this->getStress();
			CreepStressStrainTangent(1) = this->getStrain();
			CreepStressStrainTangent(2) = this->getTangent();
			CreepStressStrainTangent(3) = this->getCreepBasic(); //ntosic
			CreepStressStrainTangent(4) = this->getCreepDrying(); //ntosic
			CreepStressStrainTangent(5) = this->getMech();
			CreepStressStrainTangent(6) = this->getShrinkBasic(); //ntosic
			CreepStressStrainTangent(7) = this->getShrinkDrying(); //ntosic
			matInfo.setVector(CreepStressStrainTangent);
			return 0;
			
		default:      
			return -1;
	}
}


