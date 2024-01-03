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
 // Adam M. Knaack (adam.knaack@schaefer-inc.com) 
 // Schaefer-Inc, Cincinnati, Ohio, USA
 // Nikola D. Tosic (ntosic@imk.grf.bg.ac.rs)
 // Department for Materials and Structure, Faculty of Civil Engineering, University of Belgrade, Serbia
 // Yahya C. Kurama (ykurama@nd.edu)
 // Department of Civil and Environmental Engineering and Earth Sciences, College of Engineering, University of Notre Dame, Notre Dame, Indiana, USA
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Created: 2012
 // Last updated: 2019
 //----------------------------------------------------------------------------------------------------------------------------

 //----------------------------------------------------------------------------------------------------------------------------
 // Description: This file contains the source code of TDConcrete. 
 // TDConcrete is a time-dependent concrete material model that calculates
 // creep and shrinkage strains.
 /*-------------------------------
 ! Concrete Compression - Linear
 ! Concrete Tension - Tamai, S., Shima, H., Izumo, J., Okamura, H. 1988. Average Stress-Strain Relationship in Post Yield Range of Steel Bar in Concrete, Concrete Library of JSCE, No. 11, 117-129.
 ! Concrete Creep - Linear superposition of creep coefficient, ACI 209 time function
 ! Concrete Shrinkage - ACI 209 time function
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


#include "TDConcrete.h" //Changed by AMK
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>
#include <Information.h>
#include <elementAPI.h> //Added by AMK to use methods for parsing data line;
#include <Domain.h> //Added by AMK to get current Domain time;
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

	static int numTDConcrete = 0;

//	OPS_Export void *
	void *
	OPS_TDConcrete() {
		// Print description of material model:
		if (numTDConcrete == 0) {
			opserr << "Time-Dependent Concrete Material Model - Written by Adam Knaack, University of Notre Dame, 2012 \n";
			numTDConcrete = 1;
		}

		// Pointer to a uniaxial material that will be returned:
			UniaxialMaterial *theMaterial = 0;
		
		// Parse the input line for the material parameters:
			int iData;
			int numData;
			int numArgs;
		
			numArgs = OPS_GetNumRemainingInputArgs();
		
			if (numArgs == 13) {
				//TDConcrete(int tag, double _fc, double _epsc0, double _fcu,
				//double _epscu, double _tcr, double _ft, double _Ets, double _Ec, double _age, double _epsshu)
				double dData[12];
			
				//Collect material tag:
				numData = 1;
				if (OPS_GetIntInput(&numData, &iData) != 0) {
					opserr << "WARNING: invalid uniaxialMaterial TDConcrete tag\n";
					return 0;
				}
			
				//Collect input parameters:
				numData = 12;
				if (OPS_GetDoubleInput(&numData, dData) != 0) {
					opserr << "WARNING: invalid material property definition\n";
					return 0;
				}
			
				//Create a new materiadouble
				theMaterial = new TDConcrete(iData,dData[0],dData[1],dData[2],dData[3],dData[4],dData[5],dData[6],dData[7],dData[8],dData[9],dData[10],dData[11]);
                if (theMaterial == 0) {
					opserr << "WARNING: could not create uniaxialMaterial of type TDConcrete \n";
					return 0;
				}
			
				//Return new material:
				return theMaterial;
			}

			return 0;
	}

//-----------------------------------------------------------------------


TDConcrete::TDConcrete(int tag, double _fc, double _ft, double _Ec, double _beta, double _age, double _epsshu, double _epssha, double _tcr, double _epscru, double _epscra, double _epscrd, double _tcast): 
  UniaxialMaterial(tag, MAT_TAG_TDConcrete),
  fc(_fc), ft(_ft), Ec(_Ec), beta(_beta), age(_age), epsshu(_epsshu), epssha(_epssha), tcr(_tcr), epscru(_epscru), epscra(_epscra), epscrd(_epscrd), tcast(_tcast)
{
  ecminP = 0.0;
  deptP = 0.0;

  //sigCr = fabs(sigCr);
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
	
	eps_m = 0.0; //Added by AMK
	eps_cr = 0.0; //Added by AMK
	eps_sh = 0.0;
	epsP_cr = 0.0; //Added by AMK
	epsP_sh = 0.0; 
	epsP_m = 0.0; //Added by AMK
	
	t_load = -1.0; //Added by AMK
	crack_flag = 0;
    iter = 0;
	
	
	
	//Change inputs into the proper sign convention:
		fc = -fabs(fc); 
		epsshu = -fabs(epsshu);
		epscru = fabs(epscru); 
}

TDConcrete::TDConcrete(void):
  UniaxialMaterial(0, MAT_TAG_TDConcrete)
{
 
}

TDConcrete::~TDConcrete(void)
{
  // Does nothing
}

UniaxialMaterial*
TDConcrete::getCopy(void)
{
  TDConcrete *theCopy = new TDConcrete(this->getTag(), fc, ft, Ec, beta, age, epsshu, epssha, tcr, epscru, epscra, epscrd, tcast); 
  
  return theCopy;
}

double
TDConcrete::getInitialTangent(void)
{
	return Ec; //Added by AMK
}

double
TDConcrete::getCurrentTime(void)
{
	double currentTime;
	Domain * theDomain = ops_TheActiveDomain;

	if (theDomain != 0) {
		currentTime = theDomain->getCurrentTime();
	}
    
    return currentTime;
}	

double
TDConcrete::setCreepStrain(double time, double stress)
{
    double creep;
    double runSum = 0.0;
    
    DTIME_i[count] = ops_Dt;
    
    for (int i = 1; i<=count; i++) {
                PHI_i[i] = setPhi(time,TIME_i[i]); //Determine PHI
                runSum += PHI_i[i]*DSIG_i[i]/Ec; //CONSTANT STRESS within Time interval
    }
    
    phi_i = PHI_i[count];
    creep = runSum;
    return creep;
    
}

double 
TDConcrete::setPhi(double time, double tp)
{	
	// ACI Equation:
	double tmtp = time-tp;
	double f1 = pow((4+0.85*tp)/tp,0.5);
	double f2 = pow(tmtp,epscra)/(epscrd+pow(tmtp,epscra))*epscru;
	double f3 = (1.25*pow((tp-tcast),-0.118))/(1.25*pow(tcr,-0.118));
	double phi = f2*f3;
	return phi;
}

double 
TDConcrete::setShrink(double time)
{
	double tD = age; //Age at initiation of drying
    double shrink = 0.0;
    if (time-(tD) < 0) {
        shrink = 0.0;
    } else {
        shrink = (time-(tD)) / (epssha + (time - (tD))) * epsshu;
    }
	return shrink;
}

int
TDConcrete::setTrialStrain(double trialStrain, double strainRate)
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
		eps_cr = 0.0;
		eps_sh = 0.0;
		eps_m = 0.0;
		eps_total = trialStrain;
		sig = 0.0;
	} else { // Concrete has hardened and is ready to accept load
		// Initialize total strain:
        	eps_total = trialStrain;
	
		// Calculate shrinkage Strain:
            if (iter < 1) {
                eps_sh = setShrink(t);
            }

    	// Calculate creep and mechanical strain, assuming stress remains constant in a time step:
    	if (ops_Creep == 1) {
        	if (fabs(t-TIME_i[count]) <= 0.0001) { //If t = t(i-1), use creep/shrinkage from last calculated time step
            	eps_cr = epsP_cr;
            	eps_sh = epsP_sh;
            	eps_m = eps_total - eps_cr - eps_sh;
            	sig = setStress(eps_m, e);
            
        	} else { // if the current calculation is a new time step
        		//if (crackP_flag == 1 && sigP >= 0.0) { //if cracking occurs and previous stress is positive, 
        		// creep strain should not increase
        		//	eps_cr = epsP_cr;
        		//} else {
        		//	eps_cr = setCreepStrain(t,sig);
        		//}
        		//if (t < tcast) {
        		//opserr << "\nWARNING: TDConcrete loaded before tcast, creep and shrinkage not calculated" << endln;
        		//	eps_sh = epsP_sh;
        		//	eps_cr = epsP_cr;
        		//	eps_m = eps_total - eps_cr - eps_sh;
        		//	sig = setStress(eps_m, e);
        		//} else {	
                if (iter < 1) {
                    eps_cr = setCreepStrain(t,sig); 
                }
        		eps_m = eps_total - eps_cr - eps_sh;
        		sig = setStress(eps_m, e);
        		//}
        	}
    	} else { //Static Analysis using previously converged time-dependent strains
        	    eps_cr = epsP_cr;
            	eps_sh = epsP_sh;
            	eps_m = eps_total-eps_cr-eps_sh;
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
TDConcrete::setStress(double strain, double &stiff)
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
    		//Need to check whether cracking has occurred or not
    		//If cracked, then reloading occurs along Et
    		//If not cracked, then loading occurs according to Tens_Envlp
    			ecmax = strain; // reset ecmax
    			this->Tens_Envlp(strain, stress, stiff);
    			if (strain >= et0) {//cracking has occurred, set cracking flag
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
TDConcrete::getStrain(void)
{
	return eps_total; //Added by AMK
  //return eps;
}

double
TDConcrete::getPHI_i(void)
{
	return phi_i;
}

double 
TDConcrete::getStress(void)
{
	return sig;
}

double 
TDConcrete::getTangent(void)
{
	return e;
}

double
TDConcrete::getCreep(void)
{
	return eps_cr;
}

double
TDConcrete::getShrink(void)
{
	return eps_sh;
}

double
TDConcrete::getMech(void)
{
	return eps_m;
}

int 
TDConcrete::commitState(void)
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
      if (fabs(eps_m/sig)>Ec) {
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
    epsP_sh = eps_sh;
	epsP_cr = eps_cr;
	epsP_m = eps_m;
    if (eps_m < 0 && fabs(eps_m)>0.50*fabs(fc/Ec)) {
        double s = fabs(eps_m/fc)*Ec;
	s = 0.5*fabs(fc/Ec);
	opserr << "Strain Compression Limit Exceeded: " << eps_m << ' ' << -s << endln;
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
TDConcrete::revertToLastCommit(void)
{
  	eps_total = epsP_total; //Added by AMK;
    eps_sh = epsP_sh;
	eps_cr = epsP_cr;
	eps_m = epsP_m;  
    
  ecmin = ecminP;;
  dept = deptP;
  
  e = eP;
  sig = sigP;
  eps = epsP;
  return 0;
}

int 
TDConcrete::revertToStart(void)
{
  ecminP = 0.0;
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
TDConcrete::sendSelf(int commitTag, Channel &theChannel)
{
  static Vector data(14);
  data(0) =ft;    
  data(1) =Ec; 
  data(2) =beta;   
  data(3) =age; 
  data(4) =epsshu;   
  data(5) =epssha;    
  data(6) =tcr;   
  data(7) =epscru;
  data(8) =epscra; 
  data(9) =epscrd;     
  data(10) = this->getTag();
  data(11) = fc;
  data(12) = tcast;
  data(13) = count;
  
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "TDConcrete::sendSelf() - failed to sendSelf\n";
    return -1;
  }
  return 0;
}

int 
TDConcrete::recvSelf(int commitTag, Channel &theChannel, 
	     FEM_ObjectBroker &theBroker)
{

  static Vector data(14);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "TDConcrete::recvSelf() - failed to recvSelf\n";
    return -1;
  }
  
  ft = data(0);   
  Ec = data(1);
  beta = data(2);   
  age = data(3); 
  epsshu = data(4);   
  epssha = data(5);    
  tcr = data(6);   
  epscru = data(7);
  epscra = data(8); 
  epscrd = data(9);   
  this->setTag(data(10));
  fc = data(11);
  tcast = data(12);
  count = (int)data(13);
  
  e = eP;
  sig = sigP;
  eps = epsP;
  
  return 0;
}

void 
TDConcrete::Print(OPS_Stream &s, int flag)
{
  s << "TDConcrete:(strain, stress, tangent) " << eps << " " << sig << " " << e << endln;
}




void
TDConcrete::Tens_Envlp (double epsc, double &sigc, double &Ect)
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
  double eps0 = ft/Ec0;
  double epsu = ft*(1.0/Ets+1.0/Ec0);
  double b = beta;
  // USE THIS ONE
  if (epsc<=eps0) {
    sigc = epsc*Ec0;
    Ect  = Ec0;
  } else {
    Ect = -b*eps0*ft/pow(epsc,2)*pow(eps0/epsc,b-1.0);
    sigc = ft*pow(eps0/epsc,b);
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

  
void
TDConcrete::Compr_Envlp (double epsc, double &sigc, double &Ect) 
{
//Linear
Ect = Ec;
sigc = Ect*epsc;
  return;
}

int
TDConcrete::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}

/* Methods added by AMK: */

Response* 
TDConcrete::setResponse(const char **argv, int argc,
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
		theOutput.tag("ResponseType", "CreepStrain");
		theOutput.tag("ResponseType", "MechStrain");
		theOutput.tag("ResponseType", "ShrinkStrain");
		theOutput.tag("ResponseType", "t_load");
		theResponse = new MaterialResponse(this, 6, Vector(6));
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
TDConcrete::getResponse(int responseID, Information &matInfo)
{
	static Vector stressStrain(2);
	static Vector stressStrainTangent(3);
	static Vector CreepStressStrainTangent(6); //Added by AMK
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
			CreepStressStrainTangent(3) = this->getCreep();
			CreepStressStrainTangent(4) = this->getMech();
			CreepStressStrainTangent(5) = this->getShrink();
			matInfo.setVector(CreepStressStrainTangent);
			return 0;
			
		default:      
			return -1;
	}
}


