
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete06.cpp,v $
// Created: 07/06
// Modified by: LMS 
// Description: This file contains the class implementation for Concrete06. Based on Concrete01.
// Envelope:	Compression:	Thorenfeldt base curve as presented by Collins and Porasz (1989)
//				Tension:		linear pre-cracking and exponent (stiffening) post-cracking  by Belarbi and Hsu (1994)
// Cyclic:		Compression:	linear unloading and reloading path with stiffness degradation (similar to Yassin (1994)) 
//								Unloading path is defined with half stiffness value of loading path.
//								Both rules are connected to each other through paths with initial elastic stiffness (in compression).
//								plastic strain eq. to closely match Palermo (2003) alphaC=0.32
//				Tension:		same linear loading and unloading path with stiffness degradation.
//								plastic strain eq. to closely match Palermo (2003) alphaT=0.08
  
#include <Concrete06.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <math.h>
#include <float.h>
#include <elementAPI.h>
#include <string.h>

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))
#define MAX(a, b)  (((a) > (b)) ? (a) : (b))

void * OPS_ADD_RUNTIME_VPV(OPS_Concrete06)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 10) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: uniaxialMaterial Concrete06 ";
	opserr << "tag? fc? eo? r? k? alphaC? fcr? ecr? b? alphaT?\n";
	return 0;
    }

    int tag;
    numdata = 1;
    if (OPS_GetIntInput(&numdata,&tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[9];
    numdata = 9;
    if (OPS_GetDoubleInput(&numdata,data)) {
	opserr << "WARNING invalid double data\n";
	return 0;
    }

    UniaxialMaterial* mat = new Concrete06(tag,data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],data[8]);
    if (mat == 0) {
	opserr << "WARNING: failed to create Concrete06 material\n";
	return 0;
    }

    return mat;
}

Concrete06::Concrete06
(int tag, double FC, double EO, double R, double K, double ALPHAC, double FCR, double ECR, double B, double ALPHAT)
  :UniaxialMaterial(tag, MAT_TAG_Concrete06),
   fc(FC), eo(EO), r(R), k(K), alphaC(ALPHAC), fcr(FCR), ecr(ECR), b(B), alphaT(ALPHAT),
   Cecmax(0.0), Cscmax(0.0), Cet(0.0), CetAccum(0.0), Cet1(0.0), Cet2(0.0), 
	Cstrain(0.0), Cstress(0.0)
{
	// Make all concrete parameters negative
	if (fc > 0.0)
		fc = -fc;
		
	if (eo > 0.0)
		eo = -eo;
 
	Cecmax = -0.00000001; // to avoid Ec=0
    envelopeC(Cecmax);
	Cscmax= Tstress;
	
	Cstmax = fcr;
	Cetmax= ecr;
	CEt = fcr/ecr;
   
	// Initial tangent
	double Ec0 = fc/eo*r/(r-1.0); // initial stiffness compression
	Ctangent = Ec0;
	Ttangent = Ec0;
	CEr1 = Ec0;
	CEr2 = Ec0;
	
	Eci = Ec0;
	Eti = CEt;

	/*
	envelopeC(5.0*eo);			// reference point for degrading stiffness (reloading) at 5.0*eo
	double s5eo=Tstress;
	envelopeT(5.0*ecr);			// reference point for degrading stiffness (reloading or unloading) at 5.0*ecr
	double s5ecr=Tstress;

	ecref = (s5eo - alphaC*Eci*5.0*eo)/(Eci - alphaC*Eci);
	scref = Eci*ecref;
	etref = (s5ecr - alphaT*Eti*5.0*ecr)/(Eti - alphaT*Eti);
	stref = Eti*etref;
*/
 	// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Concrete06::Concrete06():UniaxialMaterial(0, MAT_TAG_Concrete06),
    fc(0.0), eo(0.0), r(0.0), k(0.0), alphaC(0.0), fcr(0.0), ecr(0.0), b(0.0), alphaT(0.0),
     Cecmax(0.0), Cscmax(0.0), Cet(0.0), CetAccum(0.0), Cet1(0.0), Cet2(0.0), 
	Cstrain(0.0), Cstress(0.0)
{
	// Set trial values
	this->revertToLastCommit();

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

Concrete06::~Concrete06 ()
{
   // Does nothing
}

int Concrete06::setTrialStrain (double strain, double strainRate)
{
  // Set trial strain
  this->revertToLastCommit();
  Tstrain = strain;

  if ((Tstrain) - DBL_EPSILON <= ecmax) {  //compressive envelope
	  ecmax=Tstrain;
	  envelopeC(Tstrain);
	  scmax=Tstress;
	  et1=MAX((1.0-exp(-ecmax/eo*alphaC))*ecmax, ecmax-scmax/Eci);
	  Er1=MIN(scmax/(ecmax-et1), Eci);
	  //Er1=(scmax-scref)/(ecmax-ecref);
	  //et1=ecmax-scmax/Er1;
	  et=et1;
  }
  else if (fabs(et - et2) <= DBL_EPSILON) { 
	  if ((Tstrain - DBL_EPSILON >= ecmax) && (Tstrain + DBL_EPSILON <= et)) { 
		  DefLoop(Er2);

	  }
	  else if ((Tstrain-et) + DBL_EPSILON <= stmax/Et) {
		  Tstress=Et*(Tstrain-et);
		  Ttangent=Et;
	  }
	  else {								// tension envelope
		  etmax=Tstrain;
		  envelopeT(Tstrain-et2+etAccum);
		  stmax=Tstress;
		  
		  et2=MIN((1.0-exp(-etmax/ecr*alphaT))*etmax, etmax-stmax/Eti);
	      Et=MIN(stmax/(etmax-et2), Eti);
		  
		 // Et=(stmax-stref)/(etmax-etref+etAccum-et2);
		 // et2=etmax-stmax/Et;
		  etAccum=et2-et+etAccum;
		  Er2=scmax/(ecmax-et2);
		  et=et2;
	  }
	  
  }
  else { //et=et1
	  if ((Tstrain - DBL_EPSILON >= ecmax) && (Tstrain + DBL_EPSILON <= et)) { 
		  DefLoop(Er1);

	  }
	  else if ((Tstrain-et) + DBL_EPSILON <= stmax/Et) {
		  Tstress=Et*(Tstrain-et);
		  Ttangent=Et;

	  }
	  else {								// tension envelope
		  etmax=Tstrain;
		  envelopeT(Tstrain-et1+etAccum);
 		  stmax=Tstress;
 	
		  et2=MIN((1.0-exp(-etmax/ecr*alphaT))*etmax, etmax-stmax/Eti);
	      Et=MIN(stmax/(etmax-et2), Eti);

		  //Et=(stmax-stref)/(etmax-etref+etAccum-et1);
		  //et2=etmax-stmax/Et;
		  etAccum=et2-et1+etAccum;
		  Er2=scmax/(ecmax-et2);
		  et=et2;
	  }
  }
  return 0;
}


void Concrete06::envelopeC (double e)
{
	double x = e/eo;
	
	if (e > eo) {
		double xr = pow(x,r);
		Tstress = fc*(r*x/(r-1.0+xr));
		Ttangent = fc/eo*(r/(r-1.0+xr)-r*r*xr/((r-1.0+xr)*(r-1.0+xr)));
	}
	else {
		double xkr = pow(x,k*r);
		Tstress = fc*(r*x/(r-1.0+xkr));
		Ttangent = fc/eo*(r/(r-1.0+xkr)-k*r*r*xkr/((r-1.0+xkr)*(r-1.0+xkr)));
	}
}


void Concrete06::envelopeT (double e)
{
	if (e + DBL_EPSILON >= ecr) { 
		double x = ecr/e;
		double xb = pow(x,b);
		Tstress = fcr*xb;
		Ttangent = -fcr*xb*b/e;
	}
	else {
		double Eti=fcr/ecr;
		Tstress = Eti*e;
		Ttangent = Eti;
	}
}


void Concrete06::DefLoop (double Erj)
{
	double sSup = scmax+Erj*(Tstrain-ecmax);
	//double sInf = 0.5*Erj*(Tstrain-et);
	double sInf = MIN(Erj, 0.071*Eci)*(Tstrain-et); //Seckin(1981) from Palermo(2003)

	double sTrial = Cstress + Eci*(Tstrain-Cstrain);

	if ((sSup - DBL_EPSILON <= sTrial) && (sInf + DBL_EPSILON >= sTrial)) {	
		Tstress = sTrial;
		Ttangent = Eci;	
	}
	else if (sInf - DBL_EPSILON <= sTrial) { 
		Tstress = sInf;
		Ttangent = MIN(Erj, 0.071*Eci);	
		//Ttangent = 0.5*Erj;
	}
	else {
		Tstress = sSup;
		Ttangent = Erj;	
	}
}


double Concrete06::getStress ()
{
   return Tstress;
}

double Concrete06::getStrain ()
{
   return Tstrain;
}

double Concrete06::getTangent ()
{
   return Ttangent;
}

int Concrete06::commitState ()
{
   // History variables
	Cecmax=ecmax;
	Cet=et;
	CetAccum=etAccum;
	Cscmax=scmax;
	Cet1=et1;
	Cet2=et2;
	Cstmax=stmax;
	Cetmax=etmax;
	CEt=Et;
	CEr1=Er1;
	CEr2=Er2;

   // State variables
	Cstrain = Tstrain;
	Cstress = Tstress;
	Ctangent = Ttangent;

   return 0;
}

int Concrete06::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
	ecmax=Cecmax;
	et=Cet;
	etAccum=CetAccum;
	scmax=Cscmax;
	et1=Cet1;
	et2=Cet2;
	stmax=Cstmax;
	etmax=Cetmax;
	Et=CEt;
	Er1=CEr1;
	Er2=CEr2;

   // Recompute trial stress and tangent
	Tstrain = Cstrain;
	Tstress = Cstress;
	Ttangent = Ctangent;
 
   return 0;
}

int Concrete06::revertToStart ()
{

// Initial tangent
	double Ec0 = fc/eo*r/(r-1.0); // initial stiffness compression
	
	Cecmax = -0.00000001; 
    envelopeC(Cecmax);
	Cscmax= Tstress;

	Cet=0.0;
	CetAccum=0.0;
	Cet1=0.0;
	Cet2=0.0;
	Cstmax=fcr;
	Cetmax=ecr;
	CEt=fcr/ecr;
	CEr1=Ec0;
	CEr2=Ec0;

	Eci = Ec0;
    Eti = CEt;
/*
	envelopeC(5.0*eo);			// reference point for degrading stiffness (reloading) at 5.0*eo
	double s5eo=Tstress;
	envelopeT(5.0*ecr);			// reference point for degrading stiffness (reloading or unloading) at 5.0*ecr
	double s5ecr=Tstress;

	ecref = (s5eo - alphaC*Eci*5.0*eo)/(Eci - alphaC*Eci);
	scref = Eci*ecref;
	etref = (s5ecr - alphaT*Eti*5.0*ecr)/(Eti - alphaT*Eti);
	stref = Eti*etref;
*/
   // History variables

   // State variables
	Cstrain = 0.0;
	Cstress = 0.0;
	Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   return 0;
}

UniaxialMaterial* Concrete06::getCopy ()
{
   Concrete06* theCopy = new Concrete06(this->getTag(),
                                     fc, eo, r, k, alphaC, fcr, ecr, b, alphaT);

   // Converged history variables
	
   theCopy-> Cecmax = Cecmax;
   theCopy-> Cet = Cet ;
   theCopy-> CetAccum = CetAccum ;
   theCopy-> Cscmax = Cscmax ;
   theCopy-> Cet1 = Cet1 ;
   theCopy-> Cet2 = Cet2 ;
   theCopy-> Cstmax = Cstmax ;
   theCopy-> Cetmax = Cetmax ;
   theCopy-> CEt = CEt ;
   theCopy-> CEr1 = CEr1 ;
   theCopy-> CEr2 = CEr2 ;

   // Converged state variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   return theCopy;
}

int Concrete06::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(24);
   data(0) = this->getTag();

   // Material properties
   data(1) = fc;
   data(2) = eo;
   data(3) = r;
   data(4) = k;
   data(5) = ecr;
   data(6) = fcr;
   data(7) = b;
   data(8) = alphaC;
   data(9) = alphaT;


   // History variables from last converged state

   data(10) = Cecmax;
   data(11) = Cet;
   data(12) = CetAccum;
   data(13) = Cscmax;
   data(14) = Cet1;
   data(15) = Cet2;
   data(16) = Cstmax;
   data(17) = Cetmax;
   data(18) = CEt;
   data(19) = CEr1;
   data(20) = CEr2;


   // State variables from last converged state
   data(21) = Cstrain;
   data(22) = Cstress;
   data(23) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "Concrete06::sendSelf() - failed to send data\n";

   return res;
}

int Concrete06::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(24);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "Concrete06::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
    fc = data(1);
    eo = data(2);
    r = data(3) ;
    k = data(4) ;
    ecr = data(5);
    fcr = data(6);
    b = data(7); 
    alphaC = data(8);
    alphaT = data(9);

      // History variables from last converged state
    Cecmax = data(10);
    Cet = data(11) ;
    CetAccum = data(12) ;
    Cscmax = data(13) ;
    Cet1 = data(14) ;
    Cet2 = data(15);
    Cstmax = data(16);
    Cetmax = data(17) ;
    CEt = data(18);
    CEr1 = data(19);
    CEr2 = data(20);

      // State variables from last converged state
    Cstrain = data(21);
    Cstress = data(22);
    Ctangent = data(23);

      // Set trial state variables
    Tstrain = Cstrain;
    Tstress = Cstress;
    Ttangent = Ctangent;
   }

   return res;
}

void Concrete06::Print (OPS_Stream& s, int flag)
{
   s << "Concrete06, tag: " << this->getTag() << endln;
   s << "  fc: " << fc << endln;
   s << "  eo: " << eo << endln;
   s << "  r: " << r << endln;
   s << "  k: " << k  << endln;
   s << "  ecr : " << ecr   << endln;
   s << "  fcr : " << fcr   << endln;
   s << "  b: " << b  << endln;
   s << "  alphaC: " << alphaC  << endln;
   s << "  alphaT: " << alphaT  << endln;

}

int
Concrete06::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = eo;
    return 0;
  } else
    return -1;
}
