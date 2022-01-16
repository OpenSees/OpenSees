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
                                                                        
// $Revision: 1.5 $
// $Date: 2010-05-06 00:10:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/ConfinedConcrete01.cpp,v $
                                                                        
// Description: This file contains the class definition for ConfinedConcrete01.
// Compressive envelope curve:     concrete model for confined concrete proposed by F. Braga, R. Gigliotti, M. Laterza in  
//                                 "Analytical Stress-Strain Relationship for Concrete Confined by Steel Stirrups
//                                 and/or FRP Jackets", ASCE, Journal of Structural Engineering, Vol. 132, No. 9, 
//                                 September 1, 2006 
// Unloading and reloanding curve: stress-strain model proposed by D. Karsan, and J. Jirsa in 
//                                 "Behavior of concrete under compressive loadings"
//                                 Journal of the Structural Division, ASCE, Vol. 95, No. ST12, December, 1969
// Tensile envelope strength:      model has no tensile strength
//                                 


// Written: Michele D'Amato, University of Basilicata, Potenza, Italy -email: damato.mic@gmail.com
//          Newton Le, University of California, Davis, CA, USA -email: newle@ucdavis.edu
// Created: 07/2008

#include <vector>
#include <ConfinedConcrete01.h>

#include <Matrix.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>
#include <stdio.h>
#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#define OPS_Export 
#define PI 3.14159265358979323846

/*
USER INPUT PARAMETERS

----Unconfined Concrete Properties----
fpc: f'c
stRatio: fc/f'c
epscuOption: 0-direct input of epscu
    1-Scott
    2-gamma
epscu: input dependent upon epscuOption
nuOption: 0-constant
		  1-variable, upper bound 0.5
          2-variable, no upper bound 
nuc: poisson's ratio (if nuOption is 0)

---Section Geometric Properties---
secType: see Fig. 8 in the paper (Journal of Structural Engineering, page 1408)
dim: the number of square sections (determined from secType)
semiLength[]: semilengths of the square section 

----Transverse Reinforcement Properties----

phis[]: transverse bar diameter or wrapping thickness
S[] : spacing of transverse reinforcement
fyh[]: yield strength of transverse reinforcement
mueps[]: ductility factor of transverse reinforcement
Es0[]: initial elastic modulus of transverse reinforcement
haRatio[]: hardening ratio of transverse reinforcement
wrappingArea: full area of wrap (unrolled)

----Longitudinal Reinforcement Properties----
longDiameter: longitudinal bar diameter

----Attard and Setunge Properties----
aggregateType: 0-crushed
               1-gravel
concreteType: 0-without silica fume
              1-with silica fume 

*/

ConfinedConcrete01::ConfinedConcrete01()
        :UniaxialMaterial(0, MAT_TAG_ConfinedConcrete01),
         CminStrain(0.0), CendStrain(0.0),
         Cstrain(0.0), Cstress(0.0)
{
}


ConfinedConcrete01::ConfinedConcrete01(int tag, std::vector<double> *eps, std::vector<double> *sigmac) 
  :UniaxialMaterial(tag, MAT_TAG_ConfinedConcrete01),
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0) 
{
  this->sigmac = sigmac;
  this->eps = eps;
}

ConfinedConcrete01::ConfinedConcrete01(int tag, int secType, int dim, std::vector<double> semiLength, 
				       std::vector<double> phis, std::vector<double> S, 
				       std::vector<double> fyh, std::vector<double> Es0, 
				       std::vector<double> haRatio, std::vector<double> mueps, 
				       std::vector<double> As, std::vector<double> Is, 
				       double rhos, double fpc, double stRatio, double Ec, 
				       int epscuOption, double epscu, double epscuLimit, int nuOption, double nuc, 
				       double phiLon, int concrType, int aggrType, double tol, int maxNumIter)
  :UniaxialMaterial(tag, MAT_TAG_ConfinedConcrete01),
   CminStrain(0.0), CendStrain(0.0),
   Cstrain(0.0), Cstress(0.0) 
{
  int ii;
  /*
  opserr << tag << " " <<  secType << " " << dim << " " << rhos << " " << fpc << " " << stRatio << " " << Ec << endln;
  opserr << epscuOption << " " <<  epscu << " " << epscuLimit << " " << nuOption << " " << nuc << " " << phiLon << " " ;
  opserr << concrType << " " << aggrType << " " << tol << " " << maxNumIter << endln;

  opserr << "1: "; for(ii=0; ii < semiLength.size(); ii++) opserr << semiLength[ii] << " "; opserr << endln;
  opserr << "2: "; for(ii=0; ii < phis.size(); ii++) opserr << phis[ii] << " ";  opserr << endln;
  opserr << "3: "; for(ii=0; ii < S.size(); ii++) opserr << S[ii] << " ";  opserr << endln;
  opserr << "4: "; for(ii=0; ii < fyh.size(); ii++) opserr << fyh[ii] << " ";  opserr << endln;
  opserr << "5: "; for(ii=0; ii < Es0.size(); ii++) opserr << Es0[ii] << " ";  opserr << endln;
  opserr << "6: "; for(ii=0; ii < haRatio.size(); ii++) opserr << haRatio[ii] << " ";  opserr << endln;
  opserr << "7: "; for(ii=0; ii < mueps.size(); ii++) opserr << mueps[ii] << " ";  opserr << endln;
  opserr << "8: "; for(ii=0; ii < As.size(); ii++) opserr << As[ii] << " ";  opserr << endln;
  opserr << "9: "; for(ii=0; ii < Is.size(); ii++) opserr << Is[ii] << " ";  opserr << endln;
  */
  
  double epsc, epsic, fic, ft, fpl, alpha, Eti;
  double fc = fpc * stRatio;
  
  setupAttardSetunge(fpc, stRatio, Ec, aggrType, concrType, 
		     epsc, fc, epsic, fic, ft, fpl,alpha, Eti);


  opserr << fpc << " " << stRatio << " " << Ec << " " << aggrType << " " << concrType<< " " << epsc << " " << fc << " " << epsic << " " << fic << " " << ft << " " << fpl << " " << alpha << Eti << endln;

  bglModel (semiLength, epscu, epscuOption, epscuLimit, nuOption, nuc, epsc, fc, epsic,
	    fic, ft, fpl, alpha, Eti, phis, As, Is, S, fyh, mueps, Es0, haRatio, 
	    phiLon, secType, dim, tol, maxNumIter);

  opserr  << epscu  << " " << epscuOption << " " << epscuLimit;
  opserr << " " << nuOption << " " << nuc << " " << epsc << " ";
  opserr << fc << " " << epsic << " " << fic << " " << ft << " ";
  //  opserr << fpl << " " << alpha << " " << Eti << " ";
  // opserr << phis << " " << As << " " << Is << " " << fyh << " " << mueps << " " << Es0 << " " << haRatio << " " 
  opserr << phiLon << " " << secType << " " << dim << " " << tol << " " << maxNumIter << endln;

  
  // Make all concrete parameters negative
  this->epscu = -epscu;
  this->fpcu = -(*sigmac)[((int) sigmac->size()) - 1];
  
  for (int i = 0; i < (int) eps->size(); i++) {
    (*eps)[i] = -(*eps)[i];
    if ((*sigmac)[i] != (*sigmac)[i]) {
      (*sigmac)[i] = 0.0;
    } else {
      (*sigmac)[i] = -(*sigmac)[i];
    }
  }
  // Find maximum strength
  double max = 0;
  for (int i = 0; i < (int) eps->size(); i++) {
    //printf("%f\t%f\n", (*eps)[i], (*sigmac)[i]);
    //opserr << (*eps)[i] << "\t" << sigmac->at(i) << "\n";
    if (sigmac->at(i) > max) {
      max = sigmac->at(i);
      epsc0 = -eps->at(i);
    }
  }
  
  // Initial tangent
  Ctangent = Eti;
  CunloadSlope = Eti;
  Ttangent = Eti;
  
  // Set trial values
  this->revertToLastCommit();
  
  // AddingSensitivity:BEGIN /////////////////////////////////////
  parameterID = 0;
  SHVs = 0;
  // AddingSensitivity:END //////////////////////////////////////
}

ConfinedConcrete01::~ConfinedConcrete01 ()
{
  // Does nothing
}

void ConfinedConcrete01::setupAttardSetunge(double fpc, double stRatio, double Ec, double aggrType, 
						double concrType, double &epsc, double &fc, 
						double &epsic, double &fic, double &ft, double &fpl,
						double &alpha, double &Eti) {

	//---Unconfined concrete----------------------
	//Stress peak
	fc = fpc*stRatio;

	//Stress where ec is calculated
	fpl = 0.45 * fc;

	//Concrete elastic modulus
//	rho = 2306;   //surface dry unite weight in kg/m3
//	if (fpc <= 50) {
//		ec = 0.043*pow(rho,1.5)*pow(fpc,0.5); // MPa (ACI 318 for NSC)
//	} else {
  //      ec = (3320*pow(fpc,0.5)+6900)*pow(rho/2320,1.5); // MPa (ACI 363 for HSC)
	//}

    //Initial tangent modulus
	if (fpc >= 100) {
		alpha = 1.0;
	} else if (fpc <= 20) {
		alpha = 1.17;
	} else {
		alpha = 1.17-0.17*(fpc-20)/(100-20);
	}
	
	Eti = alpha*Ec;
	this->Ec0 = Eti;

	//Strain at stress peak
	if (aggrType == 0) {
		epsc = (fpc/Ec)*4.26/pow(fpc,0.25);
	} else {
		epsc = (fpc/Ec)*3.78/pow(fpc,0.25);
	}

    //Stress and strain at point of inflection
	epsic = (2.5-0.3*log(fc))*epsc;
    fic = (1.41-0.17*log(fc))*fc;   

    //Tensile strength of concrete
	if (concrType == 0) {
		ft = 0.9*(0.32*pow(fc,0.67)); // without silica fume
	} else {
		ft = 0.9*(0.62*pow(fc,0.5));  // with silica fume
	}
}

void ConfinedConcrete01::bglModel(std::vector<double> semiLength, double & epscu, int epscuOption, double epscuLimit,
				  int nuOption, double nuc, double epsc, double fc, 
				  double epsic, double fic, double ft, double fpl, 
				  double alpha, double Eti, std::vector<double> phis, 
				  std::vector<double> As, std::vector<double> Is,
				  std::vector<double> S, 
				  std::vector<double> fyh, std::vector<double> mueps, 
				  std::vector<double> Es0, std::vector<double> haRatio, double phiLon, 
				  int secType, int dim, double tol, double maxNumIter) {
  
  std::vector<double> Ec, nu;
  /*
    eps: axial strain of confined concrete
    sigmac: corresponding stress in confined concrete
    Ec: secant modulus of " "
    nu: poisson ratio of unconfined concrete
  */
  
  std::vector<std::vector<double> > sigmaTr(1, std::vector<double>(dim,0.0));
  std::vector<std::vector<double> > epsTr(1, std::vector<double>(dim,0.0));
  std::vector<std::vector<double> > A(1, std::vector<double>(dim,0.0));
  std::vector<std::vector<double> > B(1, std::vector<double>(dim,0.0));
  std::vector<std::vector<double> > Es;
  std::vector<std::vector<double> > frm(1, std::vector<double>(dim+1,0.0));
  
  /* 
     sigmaTr: stress in transverse reinforcement
     epsTr: strain in transverse reinf
     Es: secant modulus of " "
     A: Airy's constant
     B: Airy's constant
     frm: confining pressure (dim+1 columns)
  */
  double gamma, DEc;
  double fcu;
  
  //	int np; //number of points of each curve (the size of vectors defined above)
  
  int i,j;
  
  double fcc, epscc, dEpsZ, k;
  std::vector<double> ksl;
  
  eps = new std::vector<double>();
  eps->push_back(0.0);
  sigmac = new std::vector<double>();
  sigmac->push_back(0.0);
  
  /*
    for (j = 0; j < dim; j++) {
    printf("semilength %f\n", semiLength[j]);
    printf("phi %f\n", phis[j]);
    printf("As %f\n", As[j]);
    printf("Is %f\n", Is[j]);
    printf("S %f\n", S[j]);
    }
    printf("epscu %f\n", epscu);
  */
  //---------------------Confinement along the column----------------------	
  for (j = 0; j < dim; j++) {
    
    k = confAlongCol(semiLength[j], phis[j], S[j], dim, phiLon);
    ksl.push_back(k);
    //printf("ksl %d %f\n", j, ksl[j]);
  }
  
  //----------Calculation of parameters depending on axial strain---------  
		 fcc = 0.0;
		 fcu = 0.0;
	 if (epscuOption == 2) { //If you assign gamma ratio. 
		 gamma = epscu;
		 epscu = epscuLimit; 
		 epscc = 0.0;
	 }
      
//-------Axial strain Increment for each step---------
      //(suggested 0.0005)
      dEpsZ = 0.0005;
      
// ----Setting of initial value of some parameters
      
     Ec.push_back(Eti); // Initial concrete modulus
	 Es.push_back(Es0); // Trasnversal reinforcement modulus
	 int totalIter = 0;

	  // --------------Beginning of increment of axial strain procedure--------
	  for (i = 0; eps->at(i) <= epscu + dEpsZ/10.0 && (epscuOption != 2 || eps->at(i) <= epscuLimit + dEpsZ/10.0); i++) {
		  int iter = 0;
		  nu.push_back(0.0);
		  double difference = 0.0;
		  bool c4 = 0;
		  //--------------Beginning of iteration method for each step-------------
		  do { // Beginning of the iteration loop

			  //---- Modulus of Poisson of concrete (by using the results obtained by Kupfer)
			  if (nuOption == 0)
				  nu[i] = nuc;
			  else {
				  double nu0 = 0.2;
				  double numax = 0.5;
				  double rn = eps->at(i)/epsc;
				  double nuTemp = nu0*(1+0.2*rn-pow(rn,2)+1.55*pow(rn,3));
				  if (nuOption == 1 && nuTemp >= numax)
					  nuTemp = numax;
				  nu[i] = nuTemp;
			  }

			  //------------- (Secant) Modulus of transverse reinforcement-----------
			  // ---- External transverse reinforcement
			  if (secType == 5) {  // S5 without wrapping  
				  trReinfModSqSec(fyh[0], mueps[0], Es0[0], As[0], S[0], semiLength[0], 
					  haRatio[0],  A[i][0], B[i][0], epsTr[i][0], 
					  sigmaTr[i][0], Es[i][0]);            
				  trReinfModCircSec(fyh[1],mueps[1], Es0[1], As[1], semiLength[1], 
					  haRatio[1], B[i][1], epsTr[i][1], sigmaTr[i][1], 
					  Es[i][1]);

				  if (dim == 3) { //S5 with Wrapping
					  trReinfModSqSec(fyh[2], mueps[2], Es0[2], As[2], S[2], semiLength[2], 
						  haRatio[2],  A[i][2], B[i][2], epsTr[i][2], 
						  sigmaTr[i][2], Es[i][2]); 
				  }
			  } else if (secType == 6) { //S6
				  for (j = 0; j < dim; j++) {
					  trReinfModCircSec(fyh[j], mueps[j], Es0[j], As[j], semiLength[j], 
						  haRatio[j], B[i][j], epsTr[i][j], sigmaTr[i][j], 
						  Es[i][j]);
				  }
			  } else { //Other Sections
				  for (j = 0; j < dim; j++) {
					  trReinfModSqSec(fyh[j],mueps[j],Es0[j],As[j],S[j],semiLength[j], haRatio[j],
						  A[i][j], B[i][j], epsTr[i][j], sigmaTr[i][j], Es[i][j]);
				  }
			  }
			  // ----------------------------- Airy's Constants -----------------------
			  if (secType == 5) { //S5 without wrapping
				  airyConSqSec (semiLength[0], Ec[i], nu[i], As[0], Is[0], Es[i][0],
					  S[0], eps->at(i), A[i][0], B[i][0]);
				  B[i][1] = airyConCircSec (semiLength[1], Ec[i], nu[i], As[1], Es[i][1],S[1],
					  eps->at(i));

				  if (dim == 3) //S5 with wrapping
					  airyConSqSec (semiLength[2], Ec[i], nu[i], As[2], Is[2], Es[i][2],
						S[2], eps->at(i), A[i][2], B[i][2]);
			  } else if (secType == 6) { //S6
				  for (j = 0; j < dim; j++) {
					  B[i][j] = airyConCircSec (semiLength[j], Ec[i], nu[i], As[j], Es[i][j],
						  S[j], eps->at(i));
				  }
			  } else { //Other sections
				  for (j = 0; j < dim; j++) {
					  airyConSqSec(semiLength[j], Ec[i], nu[i], As[j], Is[j], Es[i][j],
						  S[j], eps->at(i), A[i][j], B[i][j]);
				  }
			  }

			  // ---------------------------- Confining pressure ----------------------   
			  //---- External transverse reinforcement
			  if (secType == 5) {   // S5 without wrapping     
				  frm[i][0] = confPressSqSec (semiLength[0], B[i][0], ksl[0]);           
				  frm[i][1] = confPressCircSec (S[1], B[i][1], ksl[1]);

				  if (dim == 3) //S5 with wrapping
					  frm[i][2] = confPressSqSec (semiLength[2], B[i][2], ksl[2]);   
			  } else if (secType == 6) { //S6
				  for (j = 0; j < dim; j++) {
					  frm[i][j] = confPressCircSec (S[j], B[i][j], ksl[j]);  
				  }
			  } else { //Other sections
				  for (j = 0; j < dim; j++) {
					  frm[i][j] = confPressSqSec (semiLength[j], B[i][j], ksl[j]);
				  }
			  }

			  // ---- Superposition of confining pressures
			  // ---- Confining pressures are already multiplied by ksl (for wrapping ksl=1)
			  superPosConfPress (B, semiLength, frm, secType, dim, i);

			  // ---- (Secant) Modulus and stress of confined concrete by using Attard e Setunge      
			  double EcInp = Ec[i]; // Saving of concrete modulus used. Ec^k

			  double EcOut;

			  attSet (epsc, fc, epsic, fic, ft, fpl, alpha, Eti, eps->at(i), frm[i][dim], 
				  (*sigmac)[i], Ec[i], fcu, epscu, epscuOption, epscuLimit);

			 // if (epscuOption == 2 && i == 0) {
			//	  epscu = 0.2;
			//  }

			  // ------------Checking of tolerance for iteration method----------------     
			  EcOut = Ec[i]; // Saving of concrete modulus obtained Ec*^k
			  DEc = (EcOut-EcInp)/EcInp; // DEc(i)^k=(Ec^*k-Ec^k)/Ec^k
			  //if (sigmac->at(i) < fcc)
			//		Ec[i] = (EcInp + EcOut)/2;
			
			  double prevDiff = difference;
				difference = EcOut - EcInp;
				if (iter > 0 && ( (prevDiff < 0 && difference > 0) || (prevDiff > 0 && difference < 0)))
					c4 = 1;

				if (c4)
					Ec[i] = (EcInp + EcOut)/2;
			  iter++;
		  } while (fabs(DEc) > tol && iter < maxNumIter);

		  if (iter == maxNumIter) {
			  printf("WARNING: Does not converge with specified maximum number of iterations.\n");
		  }
		  totalIter += iter;
// -----In the next iteration Ec^(k+1)=E*c^k
//-------------------End of iteration method for each step--------------

// --- Calculation of Confined Concrete Peak Stress 
		  if (sigmac->at(i) > fcc) {
			  fcc = sigmac->at(i);  // Strength peak of confined concrete
			  this->fpc = -fcc;
			  epscc = eps->at(i);    // Corresponding strain
			  this->epsc0 = -epscc;
			  if (epscuOption == 2)
				fcu = gamma * fcc;   // Ultimate stress of confined concrete
		  }                
          //printf("fcc: %f\n", fcc);
		  // ------------Increment of axial strain for the next step---------------    
		  if (eps->at(i) + dEpsZ <= epscu + dEpsZ/10.0 && (epscuOption != 2 || eps->at(i) + dEpsZ <= epscuLimit + dEpsZ/10.0)) {
			  eps->push_back(0.0);
			  sigmac->push_back(0.0);
			  Ec.push_back(0.0);
			  sigmaTr.resize(i+2);
			  epsTr.resize(i+2);
			  Es.resize(i+2);
			  A.resize(i+2);
			  B.resize(i+2);
			  frm.resize(i+2);

			  sigmaTr[i+1].resize((int) sigmaTr[0].size());
			  epsTr[i+1].resize((int) epsTr[0].size());
			  Es[i+1].resize((int) Es[0].size());
			  A[i+1].resize((int) A[0].size());
			  B[i+1].resize((int) B[0].size());
			  frm[i+1].resize((int) frm[0].size());

			  (*eps)[i+1] = eps->at(i)+dEpsZ;  // Axial strain for the next step

			  for (j=0; j < dim; j++) {        
				  Es[i+1][j] = Es[i][j];   // Steel modulus for the next step
				  Ec[i+1] = Ec[i];       // Concrete modulus for the next step
			  }

			  if (eps->at(i+1) > epscu)
				  (*eps)[i+1] = epscu;

			 // printf("Point added: %f %f.Iterations: %d\n", eps->at(i), sigmac->at(i), iter);
		  } else {
			 // printf("Total iterations: %d\n", totalIter);
		    //---------------------Writing of output file----------------------------  
			 // printf("Eci\n");
 
			  /*
			  for(int n = 0; n < Ec.size(); n++) {
					printf("Ec %d %f\n", n, Ec[n]);
			  }
			  
			  for(int g = 0; g < dim; g++) {
                 printf("i %d\n", g);
                 printf("epsTr\tsigmaTr\n");
			   for(int n = 0; n < epsTr.size(); n++) {
			     printf("%f\t%f\n", epsTr[n][g], sigmaTr[n][g]);
                 }
			   }
			   */
		           
              /*
			  printf("eps\tfrm\n");
			  for (int n=0; n < frm.size(); n++) {
				  printf("%f\t%f\t%f\t%f\t%f\t%f\%f\n", (*eps)[n], frm[n][0], frm[n][1], frm[n][2], frm[n][3], frm[n][4], frm[n][5]);
			  }
              */
			  /*
			  printf("--------------------------------------------\n");
			  printf("eps\tA\tB\n");
			  for (int n=0; n < eps->size(); n++) {
				  printf("%f\t%e\t%e\n", (*eps)[n], A[n][0], B[n][0]);
			  }
			  */
			  return;
		  }
		  // -----------------End of increment of axial strain procedure-----------
	}
}


// *******************CONFINING ALONG THE COLUMN*************************
double ConfinedConcrete01::confAlongCol (double semiLength, double phis, double S, int dim, double phiLon) {

      double k, kc, ksil, beta, ksist;
					
       kc = pow((1-S/(4*semiLength)),2); //Sheink and Uzumeri, 1980

	   if (kc > 1 )
       kc = 1;

       ksil = phiLon/S;
       beta = phis/phiLon;
       ksist = phis/semiLength;
      
       k = 45.0*pow(ksil,3)/(45.0*pow(ksil,3)+beta*ksist);
	   
       if (k < kc)
         k = kc;
	
       //FMK
       if (phis == 0.0)
         k =1.0;

      return k;
}

void ConfinedConcrete01::trReinfModSqSec (double syh, double duc, double Eel, double Atr,
					  double Str, double DimSec, double hard, double CA,
					  double CB, double& etr, double& sgtr, double& Essec) {
						  
	double etry, etru;

       //  Essec is the secant modulus calculated in previous step!!!!

    etry = syh/Eel; //Yield strain of transversal steel
	etr = Str*(pow(DimSec,3))*(CA+3*CB)/(3*Essec*Atr);
    etru = duc*etry;
	
	if (etr > etru)
		return;
         
//     Transversal strain is calculated by using Airy's constants and 
//     modulus of concrete and trasversal reinf. obtained in the previous
//     iteration

	if (etr >= etry) {
		Essec = (syh+hard*Eel*(etr-etry))/etr; //Steel modulus calculated in this step.
//It will be used for the next iteration            
        sgtr = syh+hard*Eel*(etr-etry);
	} else {
		sgtr = etr*Essec;
	}
}

// **(SECANT) MODULUS OF TRANSVERSAL REINFORCEMENT FOR CIRCULAR SECTION Es(epsz)
void ConfinedConcrete01::trReinfModCircSec (double syh, double duc, double Eel, double Atr,
						double DimSec, double hard, double CB, double & etr,
						double & sgtr, double & Essec) {
	double etry, etru;

    //Essec is the secant modulus calculated in previous step!!!!
	etry = syh/Eel; //Yield strain of transverse steel
    etr = CB*DimSec/(Essec*Atr); 
    etru = duc*etry;

    if (etr > etru) 
		return;

//     Transversal strain is calculated by using Airy's constants and 
//     modulus of concrete and trasversal reinf. obtained in the previous iteration

	if (etr >= etry) {
		Essec = (syh+hard*Eel*(etr-etry))/etr; //Steel modulus calculated in this step.
											   //It will be used for the next iteration
        sgtr = syh+hard*Eel*(etr-etry);
	} else {
		sgtr = etr*Essec;
        //EsOut = Essec;
	}    
}

// **************AIRY'S CONSTANTS FOR SQUARE SECTION*********************
void ConfinedConcrete01::airyConSqSec (double semiLength, double Ecz, double nuz, double As, 
				   double Is, double Esz, double S, double epsz, double & Az,
				   double & Bz) {
	
	double    NCA, DCA, NCB, DCB, DCA1, DCA2, DCA3;    
	
	NCA = 21*S*pow(Ecz,2)*Esz*As*nuz*semiLength;
    DCA1 = 25*pow(S,2)*pow(Ecz,2)*pow(semiLength,4);
	DCA2 = 6*S*Ecz*Esz*semiLength*(315*Is*(nuz+1)+2*pow(semiLength,2)*As*(2*nuz+5));
	DCA3 = -1890*pow(Ecz,2)*Is*As*(pow(nuz,2)-1);
    DCA = DCA1+DCA2+DCA3;
    Az = NCA/DCA*epsz;

    NCB = 18*Ecz*Esz*As*nuz*(S*Ecz*pow(semiLength,3)+ 105*Esz*Is*(nuz+1));
    DCB = pow(semiLength,2)*DCA;
    Bz = NCB/DCB*epsz;  
}

// *************AIRY'S CONSTANTS FOR CIRCULAR SECTION *******************
double ConfinedConcrete01::airyConCircSec (double R, double Ecz, double nuz, double As, 
					 double Esz, double S, double epsz) {
	
	double Nq, Dq;

	Nq = Ecz*Esz*As*nuz*S;
	Dq = R*Ecz*S+Esz*As*(1-nuz)*(nuz*epsz+1);

	return Nq/Dq*epsz;        
}

// **************CONFINING PRESSURE FOR SQUARE SECTION*******************
double ConfinedConcrete01::confPressSqSec (double semiLength, double Bz, double k) {  
      return k*(Bz*pow(semiLength,2.0)); // Absolute value
}

// **********CONFINING PRESSURE FOR CIRCULAR SECTION*********************
double ConfinedConcrete01::confPressCircSec (double S, double q, double k) {  
      return k*(q/S); // Absolute value     
}

// **********SUPERPOSITION OF CONFINING PRESSURE*************************

void ConfinedConcrete01::superPosConfPress (std::vector<std::vector<double> > B, std::vector<double> semiLength, 
						  std::vector<std::vector<double> > & frm, int secType, int dim, int i) {

    double  fr1, fr2;

	if ((secType == 1) || (secType == 6)) {  //S1 and C
		if (dim == 1)
			frm[i][dim] = frm[i][0];
		else
            frm[i][dim] = pow(semiLength[0]/semiLength[1],2)*frm[i][0]+frm[i][1];
	} else if ((secType == 2) || (secType == 3)) { //S2 or S3 
        if (dim == 2) // Without wrapping
            frm[i][dim] = frm[i][0]+frm[i][1]/2;
        else // With wrapping dim=3
            frm[i][dim] = (frm[i][0]+0.5*frm[i][1])*pow(semiLength[0]/semiLength[2],2)+frm[i][2];
	} else if (secType == 41) { // 2l(1)=2l(3)=L, 2l(2)=a, 2l(4)=b
		fr1 = frm[i][1]*pow(semiLength[1],2)+frm[i][0]*(semiLength[0]-semiLength[1])*semiLength[1]; 
		fr1 = fr1/(semiLength[0]*semiLength[1]);
		 
        fr2 = frm[i][3]*pow(semiLength[3],2)+frm[i][2]*(semiLength[2]-semiLength[3])*semiLength[3];
        fr2 = fr2/(semiLength[2]*semiLength[3]);

		if (dim == 4) //Without wrapping 
			frm[i][dim] = (1+semiLength[3]/semiLength[0])*fr1+semiLength[3]/semiLength[0]*fr2;   
		else // With wrapping dim=5
			frm[i][dim] = fr1*semiLength[0]*(semiLength[0]+semiLength[3])/pow(semiLength[4],2) +
				fr2*semiLength[3]*semiLength[0]/pow(semiLength[4],2)+frm[i][4];  
	} else if (secType == 42) { // 2l(1)=2l(2)=2l(3)=L, 2l(4)=b
		fr1 = frm[i][0];
		fr2 = (frm[i][2]*pow(semiLength[2],2)+frm[i][1]*(semiLength[1]-semiLength[2])*semiLength[2]);
		fr2 = fr2/(semiLength[1]*semiLength[2]);

		if (dim == 3)
			frm[i][dim] = fr1+(2*semiLength[2]/semiLength[0])*fr2;
		else
			frm[i][dim] = fr1*pow(semiLength[0]/semiLength[3],2) + 2*(semiLength[2]*semiLength[0]/pow(semiLength[3],2))*fr2+frm[i][3];
	} else if (secType == 5) {// S5
		if (dim == 2) // Without wrapping
			frm[i][dim] = frm[i][1]*PI/4.0+frm[i][0];
		else // With wrapping
			frm[i][dim] = frm[i][0]*pow(semiLength[0]/semiLength[2],2)+frm[i][1]*(PI/4.0)*pow(semiLength[0]/semiLength[2],2)+frm[i][2];
	} else if (secType == 7) {   // l(1)=d, l(2)=c, l(3)=l(1)+cover, l(4)=l(2)+cover 
		if (dim == 2) { // Without wrapping
			frm[i][dim] = (frm[i][1]*pow(semiLength[1],2)+frm[i][0]*(semiLength[0]-semiLength[1])*
				semiLength[1])/(semiLength[0]*semiLength[1]);
		} else { // With wrapping
                fr1 = frm[i][1]*pow(semiLength[1],2)+frm[i][0]*(semiLength[0]-semiLength[1])*semiLength[1];
                fr1 = fr1 / (semiLength[0]*semiLength[1]);
                fr2 = frm[i][3]*pow(semiLength[3],2)+frm[i][2]*
					(semiLength[2]-semiLength[3])*semiLength[3];
                fr2 = fr2 / (semiLength[2]*semiLength[3]);
                frm[i][dim] = (semiLength[0]*semiLength[1])/
					(semiLength[2]*semiLength[3])*fr1+fr2;
		}
	}
}
// ************************STRESS-STRAIN RELATIONSHIP******************** 
// ********FOR CONFINED CONCRETE PROPOSED BY ATTARD AND SETUNGE**********
		void ConfinedConcrete01::attSet (double epsc, double fc, double epsic, double fic,
			double ft, double fpl, double alpha, double Eti, double epsz, double fr,
		  double & sigmaci, double &Eczi, double fcu, double &epscu, double epscuOption, double epscuLimit) {  
//     Attard&Setunge's variables 
      double    Aa, Ba1, Ba2, Ba, Ca, Da, Ad, Bd, Cd, Dd; // Sargin's costants    
      double    f2ic, k, eps0, f0, epsi, fi, f2i, eps2i, Ei, E2i; 
      double    X, Y; // Adimensional point of stress-strain curve
      double    DX;   // If you set gamma ratio these variables are used for finding epscu corresponding to fcu

	  if (fr == 0) { // If the confining pressure is equal zero
// ------------------------------Sargin's constant------------------------         
// -----Unconfined Concerete
//      Ascending branch
         Aa = Eti*epsc/fc;
         Ba1 = pow((Aa-1.0),2.0)/(alpha*(1.0-fpl/fc));
         Ba2 = pow(Aa,2.0)*(1-alpha)/(pow(alpha,2)*fpl/fc*(1-fpl/fc));         
         Ba = Ba1+Ba2-1.0;
         Ca = Aa-2.0;
         Da = Ba+1.0;

//      Descending branch      
         Ad = fic/(epsc*epsic)*pow(epsic-epsc,2)/(fc-fic);
         Bd = 0;
         Cd = Ad-2.0;
         Dd = 1.0;
        
         f0 = fc;     // Peak Strength of Unconfined Concrete
         eps0 = epsc; // Corresponding Strain
/*
      X = epsz/eps0;
	  if (epsz <= eps0)
         Y = (Aa*X+Ba*pow(X,2))/(1+Ca*X+Da*pow(X,2));
      else
         Y = (Ad*X+Bd*pow(X,2))/(1+Cd*X+Dd*pow(X,2));

	  sigmaci = Y*f0;
*/

	  } else {
// -----Confined Concerete
     
// -----Confined Peak Stress and Corresponding Strain
         k = 1.25*(1.0+0.062*fr/fc)*pow(fc,(-0.21)); // MPa        
         f0 = pow(fr/ft+1.0,k)*fc;               // MPa         
         eps0 = (1.0+(17.0-0.06*fc)*fr/fc)*epsc;   // MPa

// -----Stress and Strain at Point of Inflection
         fi = ((fic/fc-1.0)/(5.06*pow((fr/fc),0.57)+1.0)+1.0)*f0;
         epsi = ((epsic/epsc-2.0)/(1.12*pow(fr/fc,0.26)+1.0)+2.0)*eps0;       

// -----Stress and Strain at Second Point         
         f2ic = (1.45-0.25*log(fc))*fc;
         f2i = ((f2ic/fc-1.0)/(6.35*pow(fr/fc,0.62)+1.0)+1.0)*f0;        
         eps2i = 2.0*epsi-eps0;

//------Sargin's constant     
//      Ascending branch
         Aa = Eti*eps0/f0;
         Ba1 = pow(Aa-1.0,2.0)/(alpha*(1.0-fpl/f0));
         Ba2 = pow(Aa,2.0)*(1.0-alpha)/(pow(alpha,2.0)*fpl/f0*(1.0-fpl/f0));    
         Ba = Ba1+Ba2-1.0;

		 if (Ba <= -1.0) 
            Ba = -1.0;

         Ca = Aa-2.0;
         Da = Ba+1.0;
//      Descending branch      
         Ei = fi/epsi;
         E2i = f2i/eps2i;
         
         Ad = (eps2i-epsi)/eps0*(eps2i*Ei/(f0-fi)-4*epsi*E2i/(f0-f2i));
         Bd = (epsi-eps2i)*(Ei/(f0-fi)-4*E2i/(f0-f2i));
         Cd = Ad-2.0;
         Dd = Bd+1.0; 
	  }

//-----------Stress corresponding to axial strain for each step----------
      X = epsz/eps0;
	  if (epsz <= eps0)
         Y = (Aa*X+Ba*pow(X,2))/(1+Ca*X+Da*pow(X,2));
      else
         Y = (Ad*X+Bd*pow(X,2))/(1+Cd*X+Dd*pow(X,2));

//	  if (Y*f0 > sigmaci)
		sigmaci = Y*f0;
      
//-----Ultimate strain corresponding to ultimate stress of confined concrete
//-----It works only if you set gamma
	  if (epscuOption == 2) {
         Y = fcu/f0;
         DX = pow(Ad*(1-Y)+2*Y,2)+4*(Bd*(1-Y)-Y)*Y;
         X = -(Ad*(1-Y)+2*Y)-pow(DX,0.5);
         X = X/(2*(Bd*(1-Y)-Y));
         epscu = X*eps0;
         if (X <= 0)
            epscu = epscuLimit;
	  }

//-----------------Secant modulus of concrete for each  step-------------      
      if (sigmaci == 0)
         Eczi = Eti;
      else
         Eczi = sigmaci/epsz;
}

int ConfinedConcrete01::setTrialStrain (double strain, double strainRate)
{

   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tstrain = Cstrain;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON)
    return 0;

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    return 0;
  }
  
  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain < Cstrain) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
  return 0;
}



int 
ConfinedConcrete01::setTrial (double strain, double &stress, double &tangent, double strainRate)
{
	 // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;
   Tstress = Cstress;
   Ttangent = Ctangent;
   Tstrain = Cstrain;

  // Determine change in strain from last converged state
  double dStrain = strain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON) {
    stress = Tstress;
    tangent = Ttangent;
    return 0;
  }

  // Set trial strain
  Tstrain = strain;
  
  // check for a quick return
  if (Tstrain > 0.0) {
    Tstress = 0;
    Ttangent = 0;
    stress = 0;
    tangent = 0;
    return 0;
  }
  
  
  // Calculate the trial state given the change in strain
  // determineTrialState (dStrain);
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*Tstrain - TunloadSlope*Cstrain;
  
  // Material goes further into compression
  if (strain <= Cstrain) {
    TminStrain = CminStrain;
    TendStrain = CendStrain;
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
  //opserr << "Concrete01::setTrial() " << strain << " " << tangent << " " << strain << endln;
  
  stress = Tstress;
  tangent =  Ttangent;
  
  return 0;
}

void ConfinedConcrete01::determineTrialState (double dStrain)
{  
  TminStrain = CminStrain;
  TendStrain = CendStrain;
  TunloadSlope = CunloadSlope;
  
  double tempStress = Cstress + TunloadSlope*dStrain;
  
  // Material goes further into compression
  if (Tstrain <= Cstrain) {
    
    reload ();
    
    if (tempStress > Tstress) {
      Tstress = tempStress;
      Ttangent = TunloadSlope;
    }
  }
  
  // Material goes TOWARD tension
  else if (tempStress <= 0.0) {
    Tstress = tempStress;
    Ttangent = TunloadSlope;
  }
  
  // Made it into tension
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
  
}

void ConfinedConcrete01::reload ()
{
  if (Tstrain <= TminStrain) {
    
    TminStrain = Tstrain;
    
    // Determine point on envelope
    envelope ();
    
    unload ();
  }
  else if (Tstrain <= TendStrain) {
    Ttangent = TunloadSlope;
    Tstress = Ttangent*(Tstrain-TendStrain);
  }
  else {
    Tstress = 0.0;
    Ttangent = 0.0;
  }
}

void ConfinedConcrete01::envelope ()
{
	//find Tstrain between two points on curve and do linear interpolation

	if (Tstrain > 0.0 || Tstrain < epscu) {
		Ttangent = 0.0;
		Tstress = 0.0;
	} else
		for (int i = 0; i < (int) eps->size(); i++) {
				if (Tstrain > eps->at(i)) {
				Ttangent = (sigmac->at(i) - sigmac->at(i-1))/(eps->at(i) - eps->at(i-1));
				Tstress = (Tstrain - eps->at(i-1)) * Ttangent + sigmac->at(i-1);
				break;
			}
		}
/*
  if (Tstrain > epsc0) {
    double eta = Tstrain/epsc0;
    Tstress = fpc*(2*eta-eta*eta);
    double Ec0 = 2.0*fpc/epsc0;
    Ttangent = Ec0*(1.0-eta);
  }
  else if (Tstrain > epscu) {
    Ttangent = (fpc-fpcu)/(epsc0-epscu);
    Tstress = fpc + Ttangent*(Tstrain-epsc0);
  }
  else {
    Tstress = fpcu;
    Ttangent = 0.0;
  }
  */
}

void ConfinedConcrete01::unload ()
{
  double tempStrain = TminStrain;
  
  if (tempStrain < epscu)
    tempStrain = epscu;
  
  double eta = tempStrain/epsc0;
  
  double ratio = 0.707*(eta-2.0) + 0.834;
  
  if (eta < 2.0)
    ratio = 0.145*eta*eta + 0.13*eta;
  
  TendStrain = ratio*epsc0;
  
  double temp1 = TminStrain - TendStrain;
  
  //double Ec0 = 2.0*fpc/epsc0;


  double temp2 = Tstress/Ec0;
  
  if (temp1 > -DBL_EPSILON) {	// temp1 should always be negative
    TunloadSlope = Ec0;
  }
  else if (temp1 <= temp2) {
    TendStrain = TminStrain - temp1;
    TunloadSlope = Tstress/temp1;
  }
  else {
    TendStrain = TminStrain - temp2;
    TunloadSlope = Ec0;
  }
}

double ConfinedConcrete01::getStress ()
{
   return Tstress;
}

double ConfinedConcrete01::getStrain ()
{
   return Tstrain;
}

double ConfinedConcrete01::getTangent ()
{
   return Ttangent;
}

int ConfinedConcrete01::commitState ()
{
   // History variables
   CminStrain = TminStrain;
   CunloadSlope = TunloadSlope;
   CendStrain = TendStrain;

   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int ConfinedConcrete01::revertToLastCommit ()
{
   // Reset trial history variables to last committed state
   TminStrain = CminStrain;
   TendStrain = CendStrain;
   TunloadSlope = CunloadSlope;

   // Recompute trial stress and tangent
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int ConfinedConcrete01::revertToStart ()
{
//	double Ec0 = 2.0*fpc/epsc0;

   // History variables
   CminStrain = 0.0;
   CunloadSlope = Ec0;
   CendStrain = 0.0;

   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = Ec0;

   // Reset trial variables and state
   this->revertToLastCommit();

   // Quan April 2006---
   if (SHVs !=0) {SHVs->Zero();}
   parameterID=0;

   return 0;
}

UniaxialMaterial* ConfinedConcrete01::getCopy ()
{
   ConfinedConcrete01* theCopy = new ConfinedConcrete01(this->getTag(), eps, sigmac);

  theCopy->fpc = fpc;    // Compressive strength
  theCopy->epsc0 = epsc0;  // Strain at compressive strength
  theCopy->fpcu = fpcu;   // Crushing strength
  theCopy->epscu = epscu;  // Strain at crushing strength
  theCopy->Ec0 = Ec0;


  // Converged history variables
  theCopy->CminStrain = CminStrain;
  theCopy->CunloadSlope = CunloadSlope;
  theCopy->CendStrain = CendStrain;
  
  // Converged state variables
  theCopy->Cstrain = Cstrain;
  theCopy->Cstress = Cstress;
  theCopy->Ctangent = Ctangent;

  return theCopy;
}

int ConfinedConcrete01::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(11);
   data(0) = this->getTag();

   // Material properties
   data(1) = fpc;
   data(2) = epsc0;
   data(3) = fpcu;
   data(4) = epscu;

   // History variables from last converged state
   data(5) = CminStrain;
   data(6) = CunloadSlope;
   data(7) = CendStrain;

   // State variables from last converged state
   data(8) = Cstrain;
   data(9) = Cstress;
   data(10) = Ctangent;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "ConfinedConcrete01::sendSelf() - failed to send data\n";

   return res;
}

int ConfinedConcrete01::recvSelf (int commitTag, Channel& theChannel,
                                 FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(11);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);

   if (res < 0) {
      opserr << "ConfinedConcrete01::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

      // Material properties 
      fpc = data(1);
      epsc0 = data(2);
      fpcu = data(3);
      epscu = data(4);

      // History variables from last converged state
      CminStrain = data(5);
      CunloadSlope = data(6);
      CendStrain = data(7);

      // State variables from last converged state
      Cstrain = data(8);
      Cstress = data(9);
      Ctangent = data(10);

      // Set trial state variables
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }

   return res;
}

void ConfinedConcrete01::Print (OPS_Stream& s, int flag)
{
   s << "ConfinedConcrete01, tag: " << this->getTag() << endln;
   s << "  fpc: " << fpc << endln;
   s << "  epsc0: " << epsc0 << endln;
   s << "  fpcu: " << fpcu << endln;
   s << "  epscu: " << epscu << endln;
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
ConfinedConcrete01::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"fc") == 0) {// Compressive strength
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"epsco") == 0) {// Strain at compressive strength
    return param.addObject(2, this);
  }
  else if (strcmp(argv[0],"fcu") == 0) {// Crushing strength
    return param.addObject(3, this);
  }
  else if (strcmp(argv[0],"epscu") == 0) {// Strain at crushing strength
    return param.addObject(4, this);
  }
  
  return -1;
}
    
                            

int
ConfinedConcrete01::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
		this->fpc = info.theDouble;
		break;
	case 2:
		this->epsc0 = info.theDouble;
		break;
	case 3:
		this->fpcu = info.theDouble;
		break;
	case 4:
		this->epscu = info.theDouble;
		break;
	default:
		break;
	}
        
	// Make all concrete parameters negative
	if (fpc > 0.0)
		fpc = -fpc;

	if (epsc0 > 0.0)
		epsc0 = -epsc0;

	if (fpcu > 0.0)
		fpcu = -fpcu;

	if (epscu > 0.0)
		epscu = -epscu;

	// Initial tangent
//	double Ec0 = 2*fpc/epsc0;
	Ctangent = Ec0;
	CunloadSlope = Ec0;
	Ttangent = Ec0;
   	TunloadSlope = CunloadSlope;

	return 0;
}




int
ConfinedConcrete01::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double
ConfinedConcrete01::getStressSensitivity(int gradNumber, bool conditional)
{
	// Initialize return value
	double TstressSensitivity = 0.0;
	double dktdh = 0.0;
	double TstrainSensitivity = 0.0;


	// Pick up sensitivity history variables
	double CminStrainSensitivity = 0.0;
	double CunloadSlopeSensitivity = 0.0;
	double CendStrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	double CstrainSensitivity = 0.0;
	if (SHVs != 0) {
		CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
		CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
		CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
		CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
		CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
	}


	// Assign values to parameter derivatives (depending on what's random)
	double fpcSensitivity = 0.0;
	double epsc0Sensitivity = 0.0;
	double fpcuSensitivity = 0.0;
	double epscuSensitivity = 0.0;

	if (parameterID == 1) {
		fpcSensitivity = 1.0;
	}
	else if (parameterID == 2) {
		epsc0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		fpcuSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		epscuSensitivity = 1.0;
	}


	// Strain increment 
	double dStrain = Tstrain - Cstrain;

	// Evaluate stress sensitivity 
	if (dStrain < 0.0) {					// applying more compression to the material

		if (Tstrain < CminStrain) {			// loading along the backbone curve

			if (Tstrain > epsc0) {			//on the parabola
				
				TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
					      + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
						  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
				
				dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
					  * (1.0-Tstrain/epsc0)
					  - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
					  / (epsc0*epsc0);
			}
			else if (Tstrain > epscu) {		// on the straight inclined line
//cerr << "ON THE STRAIGHT INCLINED LINE" << endl;

				dktdh = ( (fpcSensitivity-fpcuSensitivity)
					  * (epsc0-epscu) 
					  - (fpc-fpcu)
					  * (epsc0Sensitivity-epscuSensitivity) )
					  / ((epsc0-epscu)*(epsc0-epscu));

				double kt = (fpc-fpcu)/(epsc0-epscu);

				TstressSensitivity = fpcSensitivity 
					      + dktdh*(Tstrain-epsc0)
						  + kt*(TstrainSensitivity-epsc0Sensitivity);
			}
			else {							// on the horizontal line
//cerr << "ON THE HORIZONTAL LINES" << endl;
				TstressSensitivity = fpcuSensitivity;
				dktdh = 0.0;
			
			}
		}
		else if (Tstrain < CendStrain) {	// reloading after an unloading that didn't go all the way to zero stress
//cerr << "RELOADING AFTER AN UNLOADING THAT DIDN'T GO ALL THE WAY DOWN" << endl;
			TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
				      + CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);

			dktdh = CunloadSlopeSensitivity;
		}
		else {

			TstressSensitivity = 0.0;
			dktdh = 0.0;

		}
	}
	else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
//cerr << "UNLOADING, BUT NOT ALL THE WAY DOWN" << endl;
		TstressSensitivity = CstressSensitivity 
			               + CunloadSlopeSensitivity*dStrain
				           + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);

		dktdh = CunloadSlopeSensitivity;
	}
	else {									// unloading all the way down to zero stress
//cerr << "UNLOADING ALL THE WAY DOWN" << endl;

		TstressSensitivity = 0.0;
		dktdh = 0.0;

	}

	return TstressSensitivity;
}



int
ConfinedConcrete01::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{

	// Initialize unconditaional stress sensitivity
	double TstressSensitivity = 0.0;
	double dktdh = 0.0;


	// Assign values to parameter derivatives (depending on what's random)
	double fpcSensitivity = 0.0;
	double epsc0Sensitivity = 0.0;
	double fpcuSensitivity = 0.0;
	double epscuSensitivity = 0.0;

	if (parameterID == 1) {
		fpcSensitivity = 1.0;
	}
	else if (parameterID == 2) {
		epsc0Sensitivity = 1.0;
	}
	else if (parameterID == 3) {
		fpcuSensitivity = 1.0;
	}
	else if (parameterID == 4) {
		epscuSensitivity = 1.0;
	}


	// Pick up sensitivity history variables
	double CminStrainSensitivity = 0.0;
	double CunloadSlopeSensitivity = 0.0;
	double CendStrainSensitivity = 0.0;
	double CstressSensitivity = 0.0;
	double CstrainSensitivity = 0.0;
	
	if (SHVs == 0) {
		SHVs = new Matrix(5,numGrads);
		CunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
	}
	else {
		CminStrainSensitivity   = (*SHVs)(0,(gradNumber-1));
		CunloadSlopeSensitivity = (*SHVs)(1,(gradNumber-1));
		CendStrainSensitivity   = (*SHVs)(2,(gradNumber-1));
		CstressSensitivity      = (*SHVs)(3,(gradNumber-1));
		CstrainSensitivity      = (*SHVs)(4,(gradNumber-1));
	}


	// Strain increment 
	double dStrain = Tstrain - Cstrain;

	// Evaluate stress sensitivity 
	if (dStrain < 0.0) {					// applying more compression to the material

		if (Tstrain < CminStrain) {			// loading along the backbone curve

			if (Tstrain > epsc0) {			//on the parabola
				
				TstressSensitivity = fpcSensitivity*(2.0*Tstrain/epsc0-(Tstrain/epsc0)*(Tstrain/epsc0))
					      + fpc*( (2.0*TstrainSensitivity*epsc0-2.0*Tstrain*epsc0Sensitivity)/(epsc0*epsc0) 
						  - 2.0*(Tstrain/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)/(epsc0*epsc0));
				
				dktdh = 2.0*((fpcSensitivity*epsc0-fpc*epsc0Sensitivity)/(epsc0*epsc0))
					  * (1.0-Tstrain/epsc0)
					  - 2.0*(fpc/epsc0)*(TstrainSensitivity*epsc0-Tstrain*epsc0Sensitivity)
					  / (epsc0*epsc0);
			}
			else if (Tstrain > epscu) {		// on the straight inclined line

				dktdh = ( (fpcSensitivity-fpcuSensitivity)
					  * (epsc0-epscu) 
					  - (fpc-fpcu)
					  * (epsc0Sensitivity-epscuSensitivity) )
					  / ((epsc0-epscu)*(epsc0-epscu));

				double kt = (fpc-fpcu)/(epsc0-epscu);

				TstressSensitivity = fpcSensitivity 
					      + dktdh*(Tstrain-epsc0)
						  + kt*(TstrainSensitivity-epsc0Sensitivity);
			}
			else {							// on the horizontal line

				TstressSensitivity = fpcuSensitivity;
				dktdh = 0.0;
			
			}
		}
		else if (Tstrain < CendStrain) {	// reloading after an unloading that didn't go all the way to zero stress

			TstressSensitivity = CunloadSlopeSensitivity * (Tstrain-CendStrain)
				      + CunloadSlope * (TstrainSensitivity-CendStrainSensitivity);

			dktdh = CunloadSlopeSensitivity;
		}
		else {

			TstressSensitivity = 0.0;
			dktdh = 0.0;

		}
	}
	else if (Cstress+CunloadSlope*dStrain<0.0) {// unloading, but not all the way down to zero stress
	
		TstressSensitivity = CstressSensitivity 
			               + CunloadSlopeSensitivity*dStrain
				           + CunloadSlope*(TstrainSensitivity-CstrainSensitivity);

		dktdh = CunloadSlopeSensitivity;
	}
	else {									// unloading all the way down to zero stress

		TstressSensitivity = 0.0;
		dktdh = 0.0;

	}

	// Commit some history variables
	(*SHVs)(3,(gradNumber-1)) = TstressSensitivity;
	(*SHVs)(4,(gradNumber-1)) = TstrainSensitivity;





	// Possibly update history variables for the three ordinary history variable derivatives
	double epsTemp, epsTempSensitivity;
	double eta, etaSensitivity;
	double ratio, ratioSensitivity;
	double temp1, temp1Sensitivity;
	double temp2, temp2Sensitivity;
	double TminStrainSensitivity;
	double TunloadSlopeSensitivity;
	double TendStrainSensitivity;

	if (dStrain<0.0 && Tstrain<CminStrain) {

		TminStrainSensitivity = TstrainSensitivity;

		if (Tstrain < epscu) {

			epsTemp = epscu; 

			epsTempSensitivity = epscuSensitivity;

		}
		else {

			epsTemp = Tstrain;

			epsTempSensitivity = TstrainSensitivity;
		}

		eta = epsTemp/epsc0;

		etaSensitivity = (epsTempSensitivity*epsc0-epsTemp*epsc0Sensitivity) / (epsc0*epsc0);

		if (eta < 2.0) {

			ratio = 0.145 * eta*eta + 0.13*eta;

			ratioSensitivity = 0.29 * eta * etaSensitivity + 0.13 * etaSensitivity;

		}
		else {

			ratio = 0.707*(eta-2.0) + 0.834;

			ratioSensitivity = 0.707 * etaSensitivity;
		}

		temp1 = Tstrain - ratio * epsc0;

		temp1Sensitivity = TstrainSensitivity - ratioSensitivity * epsc0
			                                  - ratio * epsc0Sensitivity;

		temp2 = Tstress * epsc0 / (2.0*fpc); 
		
		temp2Sensitivity = (2.0*fpc*(TstressSensitivity*epsc0+Tstress*epsc0Sensitivity)
			-2.0*Tstress*epsc0*fpcSensitivity) / (4.0*fpc*fpc);

		if (temp1 == 0.0) {

			TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
		}
		else if (temp1 < temp2) {

			TendStrainSensitivity = TstrainSensitivity - temp1Sensitivity;

			TunloadSlopeSensitivity = (TstressSensitivity*temp1-Tstress*temp1Sensitivity) / (temp1*temp1);

		}
		else {

			TendStrainSensitivity = TstrainSensitivity - temp2Sensitivity;

			TunloadSlopeSensitivity = (2.0*fpcSensitivity*epsc0-2.0*fpc*epsc0Sensitivity) / (epsc0*epsc0);
		}
	}
	else {
		TminStrainSensitivity = CminStrainSensitivity;
		TunloadSlopeSensitivity = CunloadSlopeSensitivity;
		TendStrainSensitivity = CendStrainSensitivity;
	}



	(*SHVs)(0,(gradNumber-1)) = TminStrainSensitivity;
	(*SHVs)(1,(gradNumber-1)) = TunloadSlopeSensitivity;
	(*SHVs)(2,(gradNumber-1)) = TendStrainSensitivity;

	return 0;
}

int
ConfinedConcrete01::getVariable(const char *varName, Information &theInfo)
{
  if (strcmp(varName,"ec") == 0) {
    theInfo.theDouble = epsc0;
    return 0;
  } else
    return -1;
}


static int numConfinedConcrete01Materials = 0;

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_ConfinedConcrete01Material)
{
  if (numConfinedConcrete01Materials == 0) {
    numConfinedConcrete01Materials++;
    opserr << "ConfinedConceret01 unaxial material - Written by M.D'Amato, University of Basilicata, Italy 2009\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  /*
    USER INPUT PARAMETERS

    ----Unconfined Concrete Properties----
    1.  fc
    3.  epscuOption:	
       0-direct input of epscu
       1-Scott
       2-gamma
    4.  epscu: input dependent upon epscuOption
    5.  nuOption: 0-constant
			  1-variable, upper bound 0.5
			  2-variable, no upper bound 
    6.  nuc: poisson's ratio (if nuOption is 0)

    ---Section Geometric Properties---
    7.  secType: see Fig. 8 in the paper (Journal of Structural Engineering, page 1408)
    8.  dim: the number of square sections (determined from secType)


    ----Longitudinal Reinforcement Properties----
    17. phiLon: longitudinal bar diameter

    ----Attard and Setunge Properties----
    18. aggrType: 0-crushed
		  1-gravel
    19. concrType: 0-without silica fume
		   1-with silica fume

    ----Transverse Reinforcement Properties----

    9.  semiLength[]: semilengths of the square section 
    10. phis[]: transverse bar diameter or wrapping thickness
    11. S[] : spacing of transverse reinforcement
    12. fyh[]: yield strength of transverse reinforcement
    13. mueps[]: ductility factor of transverse reinforcement
    14. Es0[]: initial elastic modulus of transverse reinforcement
    15. haRatio[]: hardening ratio of transverse reinforcement

    ---TCL Command Reference----------------------------

    OLD:
    uniaxialMaterial ConfinedConcrete01 $tag $secType $fpc <-stRatio $stRatio> $Ec
    <$epscu> <-scott> <-gamma $gamma <-epscu $epscuLimit>>   ;one option required
    <$nu> <-varub> <-varnoub>	;one option required
    $L1 <$L2> <$L3>	  ;Geometric properties
    $phis $S $fyh $Es0 $haRatio <-mu $mu>   ;External transverse
    <-phi $phisi> <-S $Si> <-fyh $fyhi> <-Es0 $Es0i> <-haRatio $haRatioi> <-mu $mui>	;Internal transverse
    <-wrap $cover $Am $Sw $ful $Es0w>  ;External wrapping
    $phiLon
    <-gravel> <-silica> <-tol $tol> <-maxNumIter $maxNumIter> ;Optional
  
    PROPOSED:
    uniaxialMaterial ConfinedConcrete01 $tag $secType $fpc $Ec 
    <-epscu $epscu> <-scott> <-gamma $gamma>  ;one option required
    <-nu $nu> <-varub> <-varnoub> ;one option required
    $L1 ($L2) ($L3) ; provision of L2 and L3 depend on previous args passed
    $phis $S $fyh $Es0 $haRatio $mu
    $phiLon
    FOLLOWING ARE OPTIONAL 
    <-internal $phisi $Si $fyhi $Es0i $haRatioi $mui 
    <-wrap $cover $Am $Sw $ful $Es0w>			   
    <-gravel> <-silica> <-tol $tol> <-maxNumIter $maxNumIter>  <-epscuLimit $epscuLimit> 
    <-stRatio $stRatio>

*/
  int argLoc = 4;
  int numReq = 12;
  double temp;
  int tag;
  int secType;
  int dim = 1;
  std::vector<double> semiLength, phis, S, fyh, Es0, haRatio, mueps, As, Is;
  double rhos;
  double fpc;
  double stRatio = 1.0;
  double Ec;
  int epscuOption; //0-direct input of epscu
  // 1-Scott
  // 2-gamma
  double epscu; // input dependent upon epscuOption
  int nuOption = 0;	// 0-constant
  // 1-variable, upper bound 0.5
  // 2-variable, no upper bound 
  double nuc = 0;			// poisson's ratio (if nuOption is 0)
  double phiLon;
  int concrType = 0;
  int aggrType = 0;

  double tol = 0.000001;
  int maxNumIter = 500;
  
  int argc = OPS_GetNumRemainingInputArgs();

  if (argc < numReq) {
    opserr << "WARNING insufficient arguments\n";
    return 0;
  }

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" << endln;
    return 0;
  }

  // ==Read required ConfinedConcrete01 material parameters===========================
  
  // --Parse section type-----------------------------------------------------
  //  char argvS[5];
  const char *argvS;
  argvS = OPS_GetString();
  /*
  if (OPS_GetString(argvS, 5) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" << endln;
    return 0;
  }
  */

  if (strcmp(argvS, "S1") == 0 || strcmp(argvS, "s1") == 0) {
    secType = 1;
  } else if (strcmp(argvS, "S2") == 0 || strcmp(argvS, "s2") == 0) {
    secType = 2;
  } else if (strcmp(argvS, "S3") == 0 || strcmp(argvS, "s3") == 0) {
    secType = 3;
  } else if (strcmp(argvS, "S4a") == 0 || strcmp(argvS, "s4a") == 0) {
    secType = 41;
  } else if (strcmp(argvS, "S4b") == 0 || strcmp(argvS, "s4b") == 0) {
    secType = 42;
  } else if (strcmp(argvS, "S5") == 0 || strcmp(argvS, "s5") == 0) {
    secType = 5;
  } else if (strcmp(argvS, "C") == 0 || strcmp(argvS, "c") == 0) {
    secType = 6;
  } else if (strcmp(argvS, "R") == 0 || strcmp(argvS, "r") == 0) {
    secType = 7;
  } else {
    opserr << "WARNING invalid section type, should be:(S1-S3, S4a, S4b, S5, R, C)\n";
    opserr << "ConfinedConcrete01 material: " << tag << endln;
    return 0;
  }
  
  // --Parse concrete properties-------------------------------------------------
  if (OPS_GetDouble(&numData, &fpc) != 0) {
    opserr << "WARNING invalid fpc\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  
  if (fpc < 0.0)
    fpc = -fpc;
  
  //printf("fpc: %f\n", fpc);
  
  if (OPS_GetDouble(&numData, &Ec) != 0) {
    opserr << "WARNING invalid Ec\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  
  //printf("Ec: %f\n", Ec);
  
  // --Parse epsilon_cu: 3 options-----------------------------------------------
  double epscuLimit = 0.05;

  //char argvEPSCU[10];
  const char *argvEPSCU = OPS_GetString();
  argLoc++;

  /*
  if (OPS_GetString(argvEPSCU, 10) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" << endln;
    return 0;
  } else

  */

  if (strcmp(argvEPSCU, "-scott") == 0 || strcmp(argvEPSCU, "-Scott") == 0) {
    epscuOption = 1;
  } else if (strcmp(argvEPSCU, "-gamma") == 0) {
    // Update number of arguments required
    numReq++;
    if (argc < numReq) {
      opserr << "WARNING insufficient arguments\n";
      return 0;
    }
    epscuOption = 2;
    if (OPS_GetDouble(&numData, &epscu) != 0) {
      opserr << "WARNING invalid gamma\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;
    } else 
      argLoc++;

  } else { //Specify epscu directly
    epscuOption = 0;
    if (OPS_GetDouble(&numData, &epscu) != 0) {
      opserr << "WARNING invalid epscu\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;	
    } else 
      argLoc++;    
  }
  

  if (epscuOption != 1 && epscu < 0.0)
    epscu = -epscu;
  
  // --Parse nu----------------------------------------------------------------

  //  char argvNu[10];
  const char *argvNu = OPS_GetString();
  /*
  if (OPS_GetString(argvNu, 10) != 0) {
    opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" << endln;
    return 0;
  } else

  */
  argLoc++;    

  if (strcmp(argvNu, "-varUB") == 0 || strcmp(argvNu, "-varub") == 0) {
    nuOption = 1;
  } else if (strcmp(argvNu, "-varNoUB") == 0 || strcmp(argvNu, "-varnoub") == 0) {
    nuOption = 2;
  } else {
    nuOption = 0;
    if (OPS_GetDouble(&numData, &nuc) != 0) {
      opserr << "WARNING invalid nu\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;	
    } else
      argLoc++;
  }

  // --Parse section geometry--------------------------------------------------
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid L1\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  } else
    argLoc++;

  semiLength.push_back(temp/2);
  
  if (secType == 2 || secType == 3) {
    dim++;
    semiLength.push_back(semiLength[0]*sqrt(2.0)/2);
  } else if (secType == 5) {
    dim++;
    semiLength.push_back(semiLength[0]);
  }
  
  if (secType == 42) {
    dim++;
    semiLength.push_back(semiLength[0]);
  }
  
  if (secType == 41 || secType == 42 || secType == 7) {
    dim++;
    if (OPS_GetDouble(&numData, &temp) != 0) {
      opserr << "WARNING invalid L2\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;	
    }
    semiLength.push_back(temp/2);
  }
  
  if (secType == 41) {
    dim += 2;
    semiLength.push_back(semiLength[0]); //semiLength[2] = semiLength[0]
    if (OPS_GetDouble(&numData, &temp) != 0) {
      opserr << "WARNING invalid L3\n";
      opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;	
    }
    semiLength.push_back(temp/2);
  }


// --Parse transverse reinforcement properties---------------------------------
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid phis\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;
  }
  phis.push_back(temp);
  As.push_back(PI*pow(phis[0],2)/4);
  Is.push_back(PI*pow(phis[0],4)/64);
  
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid S\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  S.push_back(temp);
  
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid fyh\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  fyh.push_back(temp);
  
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid Es0\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  Es0.push_back(temp);
  
  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid haRatio\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  haRatio.push_back(temp);

  if (OPS_GetDouble(&numData, &temp) != 0) {
    opserr << "WARNING invalid mu\n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }
  mueps.push_back(temp);
  
  if (OPS_GetDouble(&numData, &phiLon) != 0) {
    opserr << "WARNING invalid phiLon \n";
    opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
    return 0;	
  }

  // Preload transverse properties in other dimensions
  for (int i = 1; i < dim; i++) {
    if (secType == 2) {
      As.push_back(As[0]*sqrt(2.0)/2.0/2.0);
      Is.push_back(0.0);
    } else {
      As.push_back(As[0]);
      Is.push_back(Is[0]);
    }
    phis.push_back(phis[0]);
    S.push_back(S[0]); 
    fyh.push_back(fyh[0]);
    mueps.push_back(mueps[0]);
    Es0.push_back(Es0[0]);
    haRatio.push_back(haRatio[0]);
  }


  // LATER 
  argc = OPS_GetNumRemainingInputArgs();
  while (argc > 0) {
    const char * argvLoc = OPS_GetString();;
    /*
    if (OPS_GetString(argvLoc, 10) != 0) {
      opserr << "WARNING invalid uniaxialMaterial ConfinedConcrete01 tag" << endln;
      return 0;
    }
    */

    if (strcmp(argvLoc, "-stRatio") == 0) {
      if (OPS_GetDouble(&numData, &stRatio) != 0 || stRatio > 1.0 || stRatio < 0.0) {
	opserr << "WARNING invalid stRatio\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
    }

    else if (strcmp(argvLoc, "-epscuLimit") == 0) {
      if (OPS_GetDouble(&numData, &epscuLimit) != 0) {
	opserr << "WARNING invalid epscuLimit\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;
      }
    }

    // internal transverse
    else if ((strcmp(argvLoc, "-internalT") == 0) || (strcmp(argvLoc, "-internal") == 0)) {
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid phi (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      phis[dim-1] = temp;
      if (secType == 2) {
	As[dim-1] = temp*sqrt(2.0)/2.0/2.0;
	Is[dim-1] = 0.0;
      } else {
	As[dim-1] = PI*pow(temp,2)/4;
	Is[dim-1] = PI*pow(temp,4)/64;
      }

      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid S (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      S[dim-1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid fyh (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      fyh[dim-1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid Es0 (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      Es0[dim-1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid haRatio (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      haRatio[dim-1] = temp;

      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid mu (stirrups)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      mueps[dim-1] = temp;
    
      // One internal transverse property parsed
      // Check wrapping for allowable sections
      if (!((2 <= secType && secType <= 5) || secType == 41 || secType == 42)) {
	opserr << "WARNING Stirrups only valid for section types S2, S3, S4a, S4b, S5\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      
      if (secType == 7) {
	phis[dim-2] = phis[dim-1];
	As[dim-2] = As[dim-1];
	Is[dim-2] = Is[dim-1];
	S[dim-2] = S[dim-1];
	fyh[dim-2] = fyh[dim-1];
	Es0[dim-2] = Es0[dim-1];
	haRatio[dim-2] = haRatio[dim-1];
	mueps[dim-2] = mueps[dim-1];
      }
    }

    else if (strcmp(argvLoc, "-wrap") == 0) {
      if (argc < 6) {
	opserr << "WARNING insufficient arguments with -wrap option\n";
	return 0;
      }
      
      dim++;
      phis.push_back(0.0);
      Is.push_back(0.0);
      mueps.push_back(1.0);
      haRatio.push_back(1.0);
      argLoc++;
      
      // Get the cover, update semiLength vector
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid cover\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      semiLength.push_back(semiLength[0] + phis[0]/2 + temp);
      
      if (secType == 7) {
	dim++;
	semiLength.push_back(semiLength[1] + phis[1]/2 + temp);
      }
      
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid As (wrapping)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      As.push_back(temp);
      
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid S (wrapping) \n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
      return 0;	
      }
      S.push_back(temp);
      
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid fyh (wrapping)\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      fyh.push_back(temp);
      
      if (OPS_GetDouble(&numData, &temp) != 0) {
	opserr << "WARNING invalid Es0 (wrapping) \n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      Es0.push_back(temp);
      
      if (secType == 7) {
	phis.push_back(phis[dim-2]);
	As.push_back(As[dim-2]);
	Is.push_back(Is[dim-2]);
	S.push_back(S[dim-2]); 
	fyh.push_back(fyh[dim-2]);
	mueps.push_back(mueps[dim-2]);
	Es0.push_back(Es0[dim-2]);
	haRatio.push_back(haRatio[dim-2]);
      }
    }

    else if (strcmp(argvLoc, "-gravel") == 0) {
      aggrType = 1;
      
    } else if (strcmp(argvLoc, "-silica") == 0) {
      concrType = 1;
      
    } else if (strcmp(argvLoc, "-tol") == 0) {
      
      if (OPS_GetDouble(&numData, &tol) != 0) {
	opserr << "WARNING invalid tol\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      
    } else if (strcmp(argvLoc, "-maxNumIter") == 0) {
      if (OPS_GetInt(&numData, &maxNumIter) != 0) {
	opserr << "WARNING invalid maxNumIter\n";
	opserr << "uniaxialMaterial ConfinedConcrete01: " << tag << endln;
	return 0;	
      }
      
    } else {
      opserr << "WARNING invalid argument(s) :" << argvLoc << endln;
      return 0;
    }

    argc = OPS_GetNumRemainingInputArgs();
  }

  
  //--Transverse ratio--------------------------------------------------------
  if (secType == 1) { //S1
    rhos = 2*As[0]/(S[0]*semiLength[0]); // Add. External wrapping is neglected
  } else if (secType == 2) {   //S2
    rhos =(2*As[0]+As[1])/(S[0]*semiLength[0]); 
  } else if (secType == 3) { //S3
    rhos = (2*As[0]+sqrt(2.0)*As[1])/(S[0]*semiLength[0]);
  } else if (secType == 41) { //S41
    rhos = (3*As[0]+As[1])/(S[0]*semiLength[0]);
  } else if (secType == 42) { //S42
    rhos = 2*(As[0]+As[1])/(S[0]*semiLength[0]);
  } else if (secType == 5) { //S5
    rhos = 2*(As[0]+(1+sqrt(2.0))/3.0*As[1])/(S[0]*semiLength[0]);
  } else if (secType == 6) { //S6
    rhos = 2*As[0]/(S[0]*semiLength[0]); 
    // Add. External wrapping is neglected
  } else { //S7
    rhos = As[0]*(semiLength[0]+semiLength[1])/(S[0]*semiLength[0]*semiLength[1]);  
    // Add. External wrapping is neglected
  }
  
  //printf("rhos %f\n", rhos);
  

  
  //Parse additional optional parameters




  
  if (epscuOption == 1) {
    epscu = 0.004+0.9*rhos*(fyh[0]/300.0); // Scott et al. 1982
  }
  
  /*
    printf("SemiLengths: \n");
    
    for (int n = 0; n < semiLength.size(); n++) {
    printf("%f\n", semiLength[n]);
    }
    
    printf("Dimension: %d\n", dim);
    printf("Nu: %f\n", nuc);
  */
  
  // Parsing was successful, allocate the material
  theMaterial = new ConfinedConcrete01(tag, secType, dim, semiLength, phis, S, fyh, 
				       Es0, haRatio, mueps, As, Is, rhos, fpc, stRatio, Ec, epscuOption, epscu, epscuLimit,
				       nuOption, nuc, phiLon, concrType, aggrType, tol, maxNumIter);
  
  return theMaterial;
}

		
