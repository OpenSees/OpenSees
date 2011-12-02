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

// $Revison
// $Date
// $Source

// written by: Davide Fugazza

// Description: This file contains the class implementation for SMAMaterial.

#include <SMAMaterial.h> 
#include <Vector.h>
#include <Channel.h>
#include <math.h>

SMAMaterial::SMAMaterial(int tag, double e, double el, double s_as_s, double s_as_f, 
			 double s_sa_s, double s_sa_f)
  :UniaxialMaterial(tag,MAT_TAG_SMA), 
   E(e), eps_L(el), sig_AS_s(s_as_s), sig_AS_f(s_as_f), sig_SA_s(s_sa_s), sig_SA_f(s_sa_f), Tcsi(0.0),
   Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{
   // Initialize varibles

   this->revertToStart();
}

SMAMaterial::SMAMaterial()
  :UniaxialMaterial(0,MAT_TAG_SMA), 
   E(0.0), eps_L(0.0), sig_AS_s(0.0), sig_AS_f(0.0), sig_SA_s(0.0), sig_SA_f(0.0), Tcsi(0.0),
   Tstrain(0.0), Tstress(0.0), Ttangent(0.0)
{
   // Initialize variables

   this->revertToStart();
}

SMAMaterial::~SMAMaterial()
{
   // does nothing
}

int 
SMAMaterial::setTrialStrain(double strain, double strainRate)
{
  
///////////////////////////////////////////// 
 
Tstrain = strain;

// Trick to avoid not-converged solutions

//if (Tstrain*Cstrain <= 0) {

//	 if (Tstrain > 0) {

//	 Cstrain = 1.0e-9;

//	 }

//	 else {

//   Cstrain = -1.0e-9;

//	 }

//Ccsi = 0;
//Tstress = 0;
//Cstress = E*Cstrain; 
//Tcsi = 0;
//Ttangent = 0;

//}

// If absolute value of current strain is <= 1.0e-9 (tolerance) then perform elastic step

//if (fabs(Tstrain) <= 1.0e-9) {

//	 Tcsi = 0.0;
//	 Tstress = E*Tstrain;
//	 Ttangent = E;

//}
	
 // Statements written to avoid defining for each quadrant " >= " or " <= " when 
 // evaluating the sign of (Tstrain - Cstrain)

//if (fabs(Tstrain - Cstrain) <= 1.0e-9) {

//	 if (Tstrain < 0) {

//	 Tcsi = Ccsi;
//	 Tstress = E*(Tstrain + Tcsi*eps_L); 
//	 Ttangent = E;
  
//	 }

//	 else {

//    Tcsi = Ccsi;
//	 Tstress = E*(Tstrain - Tcsi*eps_L); 
//	 Ttangent = E;
 
//	 }
   
//}

//-------------------------------------------------------------------------------------
// FIRST QUADRANT (sigma > 0 and epsilon > 0): transformation A ---> S (loading)
//-------------------------------------------------------------------------------------

if ((Tstrain - Cstrain) > 0 && Tstrain > 0) {  

         double eps_AS_s = sig_AS_s/E + Ccsi*eps_L;
	 double eps_AS_f = sig_AS_f/E +      eps_L;

	 if (Tstrain <= eps_AS_s) {

	 Tcsi = Ccsi;
	 Tstress = E*(Tstrain - Tcsi*eps_L);
	 Ttangent = E;

         }

	 else  if (Tstrain > eps_AS_s && Tstrain < eps_AS_f) { 

		   if (Cstress <= sig_AS_s) {

		   Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_AS_f - E*Tstrain + sig_AS_s)/(-sig_AS_f + Ccsi*E*eps_L - E*eps_L + sig_AS_s);
		   Tstress = E*(Tstrain - Tcsi*eps_L);
		   double H = -(1 - Ccsi)*E/((1 - Ccsi)*(-E*eps_L) + sig_AS_s - sig_AS_f);
		   Ttangent = E*(1 - eps_L*H);

		   }

		   else {

	           Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_AS_f - E*Tstrain + Cstress)/(-sig_AS_f + Ccsi*E*eps_L - E*eps_L + Cstress);
                   Tstress = E*(Tstrain - Tcsi*eps_L);
                   double H = -(1 - Ccsi)*E/((1 - Ccsi)*(-E*eps_L) + Cstress - sig_AS_f);
                   Ttangent = E*(1 - eps_L*H);
         
		   }

         }

	 else {

         Tcsi = 1;
	 Tstress = E*(Tstrain - Tcsi*eps_L);
	 Ttangent = E;

	 }

}

//if ((Tstrain - Cstrain) > 0 && Tstrain < 0 ) { 
//	Tcsi = 0.0 ;
//	Tstress = 0.0 ;
//	Ttangent = 0.0 ;
//}
//-------------------------------------------------------------------------------------
// FIRST QUADRANT (sigma > 0 and epsilon > 0): transformation S ---> A (unloading)
//-------------------------------------------------------------------------------------

if ((Tstrain - Cstrain) < 0 && Tstrain > 0) { 

//	Tcsi = 0.0 ;
//	Tstress = 0.0 ;
//	Ttangent = 0.0 ;
         double eps_SA_s = sig_SA_s/E + Ccsi*eps_L;
	 double eps_SA_f = sig_SA_f/E;

	 if (Tstrain >= eps_SA_s) {

	 Tcsi = Ccsi;
	 Tstress = E*(Tstrain - Tcsi*eps_L);
	 Ttangent = E;

         }

	 else  if (Tstrain < eps_SA_s && Tstrain >= eps_SA_f) { 

		   if (Cstress > sig_SA_s) {

	  	   Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_SA_f)/(-sig_SA_f + Ccsi*E*eps_L + sig_SA_s);
		   Tstress = E*(Tstrain - Tcsi*eps_L);
		   double H = (Ccsi*E)/(- Ccsi*(-E*eps_L) + sig_SA_s - sig_SA_f);
		   Ttangent = E*(1 - eps_L*H);

		   }

		   else {

		  Tcsi = (Ccsi*E*Tstrain - Ccsi*sig_SA_f)/(-sig_SA_f + Ccsi*E*eps_L + Cstress);
          Tstress = E*(Tstrain - Tcsi*eps_L);
          double H = (Ccsi*E)/(-Ccsi*(-E*eps_L) + Cstress - sig_SA_f);
          Ttangent = E*(1 - eps_L*H);
         
		   }

         }

	 else {

         Tcsi = 0;
	 Tstress = E*(Tstrain);
	 Ttangent = E;

	 }

}

//if ((Tstrain - Cstrain) < 0 && Tstrain > 0) { 

//	Tcsi = 0;
//	Tstress = E*(Tstrain);
//	Ttangent = E;
//}
//-------------------------------------------------------------------------------------
// THIRD QUADRANT (sigma < 0 and epsilon < 0): transformation S ---> A (loading)
//-------------------------------------------------------------------------------------

if ((Tstrain - Cstrain) < 0 && Tstrain < 0.0) { 
	

//       double eps_AS_s = -sig_AS_s/E - Ccsi*eps_L ;
//	 double eps_AS_f = -sig_AS_f/E - eps_L ;

//	 if (Tstrain >= eps_AS_s) {

//	 Tcsi = Ccsi;
//	 Tstress = E*(Tstrain + Tcsi*eps_L);
//	 Ttangent = E;

//       }

//	 else  if (Tstrain < eps_AS_s && Tstrain > eps_AS_f) { 

//		   if (Cstress >= -sig_AS_s) {

//		   Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_AS_f - E*Tstrain - sig_AS_s)/(sig_AS_f - Ccsi*E*eps_L + E*eps_L - sig_AS_s);
//		   Tstress = E*(Tstrain + Tcsi*eps_L);
//		   double H = (1 - Ccsi)*E/((1 - Ccsi)*(E*eps_L) - sig_AS_s + sig_AS_f);
//		   Ttangent = E*(1 - eps_L*H);

//		   }

//		   else {

//                 Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_AS_f - E*Tstrain + Cstress)/(sig_AS_f - Ccsi*E*eps_L + E*eps_L + Cstress);
//                 Tstress = E*(Tstrain + Tcsi*eps_L);
//                 double H = (1 - Ccsi)*E/((1 - Ccsi)*(E*eps_L) + Cstress + sig_AS_f);
//                 Ttangent = E*(1 - eps_L*H);
         
//		   }

//       }

//	 else {

//       Tcsi = 1;
//	 Tstress = E*(Tstrain + Tcsi*eps_L);
//	 Ttangent = E;

//	 }

//}

	 Tcsi = 0.0;
	 Tstress = 0.0;
	 Ttangent = E;

}
//-------------------------------------------------------------------------------------
// THIRD QUADRANT (sigma < 0 and epsilon < 0): transformation S ---> A (unloading)
//-------------------------------------------------------------------------------------

if  ((Tstrain - Cstrain) > 0 && Tstrain < 0.0) { 

//if (Tstrain < 0) {
//       double eps_SA_s = -sig_SA_s/E - Ccsi*eps_L ;
//	 double eps_SA_f = -sig_SA_f/E;

//	 if (Tstrain <= eps_SA_s) {

//	 Tcsi = Ccsi;
//	 Tstress = E*(Tstrain + Tcsi*eps_L);
//	 Ttangent = E;

//       }

//	 else  if (Tstrain > eps_SA_s && Tstrain <= eps_SA_f) { 

//		   if (Cstress < -sig_SA_s) {

//		   Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_SA_f)/(sig_SA_f - Ccsi*E*eps_L - sig_SA_s);
//		   Tstress = E*(Tstrain + Tcsi*eps_L);
//		   double H = -(Ccsi*E)/(-Ccsi*(E*eps_L) - sig_SA_s + sig_SA_f);
//		   Ttangent = E*(1 - eps_L*H);

//		   }

//		   else {

//		   Tcsi = (Ccsi*E*Tstrain + Ccsi*sig_SA_f)/(sig_SA_f - Ccsi*E*eps_L + Cstress);
//                 Tstress = E*(Tstrain + Tcsi*eps_L);
//                 double H = -(Ccsi*E)/(-Ccsi*(E*eps_L) + Cstress + sig_SA_f);
//                 Ttangent = E*(1 - eps_L*H);
         
//		   }

//       }

//	 else {

//       Tcsi = 0;
//	 Tstress = E*(Tstrain); 
//	 Ttangent = E;

//	 }

//}

	Tcsi = 0.0;
	Tstress = 0;
	Ttangent = E;

}



// Trick to avoid not-converged solutions. This line enforce the tangent to assume
// the elastic value (i.e. Young's modulus)

//if (Tstrain*Cstrain < 1.0e-9) {

//	 Ttangent = E;

//}

//if (Tstrain*Cstrain <= 0) {


//	 else {

//   Cstrain = -1.0e-9;

//	 }

//Ccsi = 0;
//Tstress = 0;
//Cstress = E*Cstrain; 
//Tcsi = 0;
//Ttangent = 0;

//}

//////////////////////////////////////////// 

  return 0;
}

double 
SMAMaterial::getStress(void)
{
  return Tstress;
}

double 
SMAMaterial::getTangent(void)
{
  return Ttangent;
}

double 
SMAMaterial::getInitialTangent(void)
{
  return 0.0;
}

double 
SMAMaterial::getStrain(void)
{
  return Tstrain;
}

int 
SMAMaterial::commitState(void)
{

  // Commit trial state variables

  Cstrain = Tstrain ;
  Cstress = Tstress ;
  Ccsi = Tcsi ;

  return 0;
}

int 
SMAMaterial::revertToLastCommit(void)
{
  return 0;
}

int 
SMAMaterial::revertToStart(void)
{
  // Reset committed history variables

  Cstrain = 0.0 ;
  Cstress = 0.0 ;
  Ccsi = 0.0;

  // Reset trial history variables

  // Initialize state variables

  Tcsi = 0.0;
  Tstrain = 0.0 ;
  Tstress = 0.0 ;
  Ttangent = E ;

  return 0;
}

UniaxialMaterial *
SMAMaterial::getCopy(void)
{
  SMAMaterial *theCopy = new SMAMaterial(this->getTag(), E, eps_L, sig_AS_s, sig_AS_f, sig_SA_s, sig_SA_f);

  // Copy committed history variables
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ccsi = Ccsi;

  // Copy trial state variables
  theCopy->Tcsi = Tcsi;
  theCopy->Tstrain = Tstrain;
  theCopy->Tstress = Tstress;
  theCopy->Ttangent = Ttangent;

  return theCopy;
}

int 
SMAMaterial::sendSelf(int cTag, Channel &theChannel)
{
  return -1;
}

int 
SMAMaterial::recvSelf(int cTag, Channel &theChannel, 
			      FEM_ObjectBroker &theBroker)
{
  return -1;
}

void 
SMAMaterial::Print(OPS_Stream &s, int flag)
{
  s << "SMAMaterial, tag: " << this->getTag() << endln;
  s << " E:         " << E << endln;
  s << " eps_L:     " << eps_L << endln;
  s << " sig_AS_s:  " << sig_AS_s << endln;
  s << " sig_AS_f:  " << sig_AS_f << endln;
  s << " sig_SA_s:  " << sig_SA_s << endln;
  s << " sig_SA_f:  " << sig_SA_f << endln;
  return;
}


