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
                                                                        
// $Revision: 1.19 $
// $Date: 2007/11/30 23:34:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/Concrete01.cpp,v $
                                                                        
// Written: Murat 
// Created: 03/08
// 
// Description: This file contains the class implementation for 
// Concrete with . 
//


#include <NDMaterial.h>
#include <ConcreteMcftNonLinear5.h>
#include <Matrix.h>
#include <Vector.h>
#include <math.h> 
#include <Parameter.h>
#include <Information.h>

#include <elementAPI.h>

void *OPS_ConcreteMcftNonlinear5()
{
  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 9) {
    opserr << "ERROR not enough input args: nDMaterial ConcreteMcftNonlinear5 tag? fcu? ecu? Ec? fcr? Esv? fyv? alphaV? RoV?" << endln;
    return 0;
  }

  int tag;
  double dData[8];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "ERROR nDMaterial ConcreteMcftNonlinear5 - unable to read matTag" << endln;
    return 0;
  }

  numData = 8;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "ERROR nDMaterial ConcreteMcftNonlinear5 - unable to read inputs" << endln;    
    return 0;
  }
  
  NDMaterial *theMaterial = new ConcreteMcftNonLinear5(tag,
						       dData[0], dData[1], dData[2], dData[3],
						       dData[4], dData[5], dData[6], dData[7]);
  if (theMaterial == 0) {
    opserr << "ERROR - could not create nDMaterial ConcreteMcftNonlinear5 with tag " << tag << endln;
    return 0;
  }
  
  return theMaterial;
}
	
// constructor: 
ConcreteMcftNonLinear5::ConcreteMcftNonLinear5 
(int tag, double fcui, double ecui, double Eci, double fcri, double Esvi, double fyvi, double alphaVi, double RoVi)
   :NDMaterial(tag, ND_TAG_ConcreteMcftNonLinear5), 
   fcu(fcui), ecu(ecui), Ec(Eci), fcr(fcri), Esv(Esvi), fyv(fyvi), alphaV(alphaVi), RoV(RoVi), epsf(2), dDri(2,2),
	TT(3,3), dDdfcu(2,2), dDdRoV(2,2), Dr(2,2), sigf(2), fx(0.0), fxy(0.0), fy(0.0), dsigfdfcu(2), dsigfdRoV(2)
	

{
exP = 0.0;
exyP = 0.0;
fxP = 0.0; 
fxyP =0.0;

Dr(0,0) = Ec;
Dr(0,1) = 0.0;
Dr(1,0) = 0.0;
Dr(1,1) = Ec/2;

dfxdfcuP  = 0.0;
dfxydfcuP = 0.0;
dfxdRoVP  = 0.0;
dfxydRoVP = 0.0;

//Sensitivity
parameterID = 0;

//opserr << " check1 " << endln;
}


ConcreteMcftNonLinear5::ConcreteMcftNonLinear5(void)
  :NDMaterial(0, ND_TAG_ConcreteMcftNonLinear5), epsf(2), dDri(2,2),
	TT(3,3), dDdfcu(2,2), dDdRoV(2,2), Dr(2,2), sigf(2),  fx(0.0), fxy(0.0), fy(0.0), 
	dsigfdfcu(2), dsigfdRoV(2)
 {
exP = 0.0;
exyP = 0.0;
fxP = 0.0; 
fxyP =0.0;

Dr(0,0) = Ec;
Dr(0,1) = 0.0;
Dr(1,0) = 0.0;
Dr(1,1) = Ec/2;

dfxdfcuP  = 0.0;
dfxydfcuP = 0.0;
dfxdRoVP  = 0.0;
dfxydRoVP = 0.0;

//Sensitivity
parameterID = 0;

//opserr << " check2 " << endln;
}

ConcreteMcftNonLinear5::~ConcreteMcftNonLinear5(void)
{
	
}


NDMaterial*
ConcreteMcftNonLinear5::getCopy (void)
{
	   
	ConcreteMcftNonLinear5 *theCopy;

theCopy = new ConcreteMcftNonLinear5 (this->getTag(), fcu, ecu, Ec, fcr, Esv, fyv, alphaV, RoV);
//opserr << " check3 " << endln;
	return theCopy;


}

NDMaterial*
ConcreteMcftNonLinear5::getCopy (const char *type)
{
//opserr << " check4 " << endln;
	return this->getCopy();


}




int
ConcreteMcftNonLinear5::setTrialStrain (const Vector &strain)
{
  epsf = strain;
  
//opserr << " check5 " << endln;
 // //opserr<< "Strain" << epsf <<  endln;
  return 0;
}

int
ConcreteMcftNonLinear5::setTrialStrain (const Vector &strain, const Vector &rate)
{
   epsf = strain;
//opserr << " check6 " << endln;
  return 0;
  

}

int
ConcreteMcftNonLinear5::setTrialStrainIncr (const Vector &strain)
{
  epsf += strain;
  //opserr << " check7 " << endln;
  return 0;
  

}

int
ConcreteMcftNonLinear5::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  epsf += strain;
  //opserr << " check8 " << endln;
  return 0;
  

}


const Vector&
ConcreteMcftNonLinear5::getStrain (void)
{
	//opserr << " check9 " << endln;
	return epsf;


}


const Matrix&
ConcreteMcftNonLinear5::getTangent (void)
{
//opserr << " check10 " << endln;
return Dr;


}

const Matrix&
ConcreteMcftNonLinear5::getInitialTangent (void)
{

    Dr(0,0) = Ec;
    Dr(0,1) = 0.0;
	Dr(1,0) = 0.0;
	Dr(1,1) = Ec/2;

//opserr << " check11 " << endln;
 return Dr;


}

const Vector&
ConcreteMcftNonLinear5::getStress (void)
{
///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

	double ex  = epsf(0);
	double 	exy = epsf(1);

	double exdelta = ex - exP;
	double exydelta = exy - exyP;
//opserr << exdelta <<endln;
	
	Dr.Zero();
	dDdfcu.Zero();
	dDdRoV.Zero();
	dDri.Zero();
	TT.Zero();

		InitCrackAngle = 0.000001;
	int counter1 = 0;
    double pi = 3.141592654;
	double degreetorad = 3.141592654*2/360.0;
	double nE  = Ec/(Ec-fcu/ecu);


	//angle search algorithm boundaries
	double bound1 = InitCrackAngle*degreetorad;
	double bound2 = (90-InitCrackAngle)*degreetorad;
	double ResiStress = 1;
	double ResiStressCheck = 1000.0 ;
	double toleranceEpsy = 0.0001;
	double stepsize;
	double nstep = 90;


	// temporary parameters
	double theta ;
	double e1, e2, ey, sig1, sig2;

	//final values
	double e1f, e2f, eyf, fxf, fyf, fxyf, thetaf;


	// stepsize and increments in angle
	stepsize = (bound2-bound1)/nstep;
	theta = bound1+ counter1*stepsize;
//----------------------------------------------------------------------------------------------------------





     //compatibility equations depend on the sign of exy strain because e1 is always defined as Tension
	 //where e2 is always defined as Compression

	double countStep = 1;
	if (exy != 0 && abs(exy/ex) > 0.01 ) {

			//Update the principle axis Angle LOOP
			while((abs(ResiStress) > toleranceEpsy) ) {
			//opserr << "theta" << theta << endln;
			if ( exy > 0) {
				//CASE 1 formulation 
				//opserr << "exy + " << endln;
				e2 = ex - exy * tan(theta) / 2;
		
				
			} else if (exy < 0 ) {
				//CASE 2 formulation
				//opserr << "exy - " << endln;
				e2 = ex + exy * tan(theta) / 2;

			}
			
			if (e2<0) {

			double tan2 = tan(theta)*tan(theta);
			e1 = (ex - e2 + ex*tan2)/ tan2;
			ey = e1 + e2 - ex;

			
			//Stresses
			if (e1 <= fcr/Ec) {
				sig1 = Ec*e1;
			} else {
			sig1 = fcr / (1+sqrt(500*e1));
			}
			sig2 = (e2/ecu)*fcu*nE/(nE-1+pow(e2/ecu,nE));

		



			
			

			fxy = (sig1 - sig2)/2 * sin(2*theta);
			fx = sig2 + fxy*tan(theta);
			fy = sig1 - fxy*tan(theta);

			//Force equilibrium in vertical direction

			ResiStress = Esv * RoV * ey + fy; 

				
			
			
//CheckPoint 1

			if (countStep  > 2 && ResiStress * ResiStressCheck < 0 ) {
					bound1 = theta - stepsize;
					bound2 = theta + stepsize;
					nstep = 10;
					stepsize = (bound2-bound1)/nstep;
					counter1 = 0;
			}

//CheckPoint 2

			if ( ResiStress < toleranceEpsy ) {
					FinalAnglex = theta;
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;
			}

//CheckPoint 3



			if( (countStep > 2) && (ResiStressCheck < 0 && ResiStress < 0) && (ResiStress < ResiStressCheck)) {
				
					e1 = e1f;
					e2 = e2f;
					ey = eyf;
					fx = fxf;
					fy = fyf;
					fxy= fxyf;
					theta = thetaf;
					
					FinalAnglex = theta;
					if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
					Strain1 = e1;
					Strain2 = e2;
					Sigma1 = sig1;
					Sigma2 = sig2;
					epsy = ey;

					
			sigf(0) = fx;
			sigf(1) = fxy;

					break;
 
			}


//CheckPoint 4

		
				

			if (countStep == 90) {
				break;
			}


			counter1++;		//count theta stepsize
			
			if(abs(ResiStress) > toleranceEpsy)
			theta = bound1 + counter1*stepsize;
		

			ResiStressCheck = ResiStress;	
			countStep++;	//count computation steps
			

			// record previous step parameters
			e1f = e1;
			e2f = e2;
			eyf = ey;
			fxf = fx;
			fyf = fy;
			fxyf= fxy;
			thetaf = theta;
	

} else if (e2>0) {

			counter1++;		//count theta stepsize
			theta = bound1 + counter1*stepsize;


			ResiStressCheck = ResiStress;	
			countStep++;	
			
		}

}

//CHECK tangent matrix Dr11 values
		if(exy>0) {
			
			Dr(0,0) = this->c1tmd00(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(0,1) = this->c1tmd01(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,0) = this->c1tmd10(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,1) = this->c1tmd11(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
		
			dDdfcu(0,0) = this->c1dd00dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(0,1) = this->c1dd01dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,0) = this->c1dd10dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,1) = this->c1dd11dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdRoV(0,0) = this->c1dd00dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(0,1) = this->c1dd01dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,0) = this->c1dd10dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,1) = this->c1dd11dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);


		} else if (exy<0) {

			Dr(0,0) = this->c2tmd00(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(0,1) = this->c2tmd01(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,0) = this->c2tmd10(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,1) = this->c2tmd11(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdfcu(0,0) = this->c2dd00dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(0,1) = this->c2dd01dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,0) = this->c2dd10dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,1) = this->c2dd11dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdRoV(0,0) = this->c2dd00dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(0,1) = this->c2dd01dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,0) = this->c2dd10dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,1) = this->c2dd11dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);


		}


			sigf(0) = fx;
			sigf(1) = fxy;

			
	} else {
						
			if (ex<0) {

				e2 = ex;
				e1 = 0.0;
				ey = 0.0;
				fx = (e2/ecu)*fcu*nE/(nE-1+pow(e2/ecu,nE));
				
				if (e1 <= fcr/Ec) {
					fy = Ec*e1;
				} else {
					fy = fcr / (1+sqrt(500*e1));
				}
				
				
				fy = e1* Ec;
				sig2 = fx;
				sig1 = fy; 
				fxy = 0.0;
				FinalAnglex = 0.001;
				if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
				Strain1 = e1;
				Strain2 = e2;
				Sigma1 = sig1;
				Sigma2 = sig2;
				epsy = ey;
			
			
			sigf(0) = fx;
			sigf(1) = fxy;

			Dr(0,0) = this->c2tmd00(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(0,1) = this->c2tmd01(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,0) = this->c2tmd10(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,1) = this->c2tmd11(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdfcu(0,0) = this->c2dd00dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(0,1) = this->c2dd01dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,0) = this->c2dd10dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,1) = this->c2dd11dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdRoV(0,0) = this->c2dd00dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(0,1) = this->c2dd01dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,0) = this->c2dd10dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,1) = this->c2dd11dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			} else if (ex>0) {

				e1 = ex;
				e2 = 0.0;
				ey = e2;
				if (e1 <= fcr/Ec) {
					fx = Ec*e1;
				} else {
					fx = fcr / (1+sqrt(500*e1));
				}
				
				fx = e1* Ec;
				fy = (e2/ecu)*fcu*nE/(nE-1+pow(e2/ecu,nE));
				fxy = 0.0;
				sig1= fx;
				sig2= fy;
				FinalAnglex = 89.999;
				if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
				Strain1 = e1;
				Strain2 = e2;
				Sigma1 = sig1;
				Sigma2 = sig2;
				epsy = ey;


			sigf(0) = fx;
			sigf(1) = fxy;


			Dr(0,0) = this->c1tmd00(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(0,1) = this->c1tmd01(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,0) = this->c1tmd10(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			Dr(1,1) = this->c1tmd11(ex,exy,FinalAnglex,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);


			dDdfcu(0,0) = this->c1dd00dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(0,1) = this->c1dd01dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,0) = this->c1dd10dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdfcu(1,1) = this->c1dd11dfcu(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			dDdRoV(0,0) = this->c1dd00dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(0,1) = this->c1dd01dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,0) = this->c1dd10dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);
			dDdRoV(1,1) = this->c1dd11dRoV(ex,exy,theta,Ec,nE,fcu,ecu,e1,fcr,Esv,RoV);

			} else { 
				e1 = 0.0;
				e2 = 0.0;
				ey = 0.0;
				fx = 0.000;
				fy = 0.000;
				fxy =0.000;
				sig1=fx;
				sig2=fy;
				FinalAnglex = 0.001;
				if(e1 >= fcr/Ec) { 
						crackLabel = 1;
					} else {
						crackLabel = 0;
					}
				Strain1 = e1;
				Strain2 = e2;
				Sigma1 = sig1;
				Sigma2 = sig2;
				epsy = ey;

			
			sigf(0) = fx;
			sigf(1) = fxy;

			Dr(0,0) = Ec;
			Dr(0,1) = 0.0;
			Dr(1,0) = 0.0;
			Dr(1,1) = Ec/2;

			dDdfcu(0,0) = 0.0;
			dDdfcu(0,1) = 0.0;
			dDdfcu(1,0) = 0.0;
			dDdfcu(1,1) = 0.0;

			dDdRoV(0,0) = 0.0;
			dDdRoV(0,1) = 0.0;
			dDdRoV(1,0) = 0.0;
			dDdRoV(1,1) = 0.0;
			}
			
	}


//---------------------------------------------------------------------------------------------------------


	t00 = Dr(0,0);
	t01 = Dr(0,1);
	t10 = Dr(1,0);
	t11 = Dr(1,1);
  

	sigf(0) =  fxP  + Dr(0,0) * exdelta + Dr(0,1) * exydelta;
	sigf(1) =  fxyP + Dr(1,0) * exdelta + Dr(1,1) * exydelta;
//opserr << "check 12a" << endln;
	
	dsigfdfcu(0) =	 dfxdfcuP   + 	dDdfcu(0,0) * exdelta + 	dDdfcu(0,1) * exydelta;
	dsigfdfcu(1) =   dfxydfcuP  + 	dDdfcu(1,0) * exdelta + 	dDdfcu(1,1) * exydelta;

	dsigfdRoV(0) =	 dfxdRoVP   + 	dDdRoV(0,0) * exdelta + 	dDdRoV(0,1) * exydelta;
	dsigfdRoV(1) =   dfxydRoVP  + 	dDdRoV(1,0) * exdelta + 	dDdRoV(1,1) * exydelta;


//opserr << " check12 " << endln;
//opserr << " Dr = "  << Dr << endln;

//opserr << " sensfcu = "  << dsigfdfcu << endln;
//getchar();
	return sigf;
	

}



/////////////////////////////////////////////////////////////////////////////////////////////
int
ConcreteMcftNonLinear5::commitState (void)
{
	exP = epsf(0);
	exyP = epsf(1);
	fxP = sigf(0);
	fxyP = sigf(1);

	dfxdfcuP = dsigfdfcu(0);
	dfxydfcuP = dsigfdfcu(1);
	
	dfxdRoVP = dsigfdRoV(0);
	dfxydRoVP = dsigfdRoV(1);
	
	//opserr << " check13 " << endln;
	return 0;
	

}

int
ConcreteMcftNonLinear5::revertToLastCommit (void)
{
//opserr << " check14 " << endln;
	return 0;
	

}

int
ConcreteMcftNonLinear5::revertToStart (void)
{
epsf.Zero();

	//opserr << " check15 " << endln;
	return 0;
	

}


const char*
ConcreteMcftNonLinear5::getType (void) const
{
	
	//opserr << " check16 " << endln;
    return "BeamFiber2d";
	

}

int 
ConcreteMcftNonLinear5::getOrder(void) const
{
	
  //opserr << " check17 " << endln;
  return 2;
  

}


void
ConcreteMcftNonLinear5::Print (OPS_Stream &s, int flag)
{
	
	//opserr << " check18 " << endln;
	return;
	
}

int
ConcreteMcftNonLinear5::sendSelf(int commitTag, Channel &theChannel)
{
	
	//opserr << " check19 " << endln;
	return 0;
	

}

int
ConcreteMcftNonLinear5::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	
	//opserr << " check20 " << endln;
	return 0;
	

}




Response *
ConcreteMcftNonLinear5::setResponse (const char **argv, int argc, 
				   OPS_Stream &s)
{

//opserr << " check21 " << endln;
	Response *theRes = NDMaterial::setResponse(argv, argc, s);
	if (theRes != 0)
		return theRes;
	else {

		if (strcmp(argv[0], "crackAngle") == 0) {
			return new MaterialResponse(this, 10, Vector(5));
		}
		else if (strcmp(argv[0], "fiberStress") == 0) {
			return new MaterialResponse(this, 11, Vector(8));
		}
		else
			return 0;
	}
	

}


// AddingSensitivity:BEGIN ///////////////////////////////////
int
ConcreteMcftNonLinear5::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"fcu") == 0) {// Compressive strength
	//  opserr << " check22 " << endln;
    return param.addObject(1, this);
  }
  else if (strcmp(argv[0],"RoV") == 0) {// Vertical Reinforcement Ratio
	//  opserr << " check22 " << endln;
    return param.addObject(2, this);
  }  
  return -1;
   

}
    
                            

int
ConcreteMcftNonLinear5::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case 1:
	//	opserr << " check23 " << endln;
		this->fcu = info.theDouble;
		break;

	case 2:
	//	opserr << " check23 " << endln;
		this->RoV = info.theDouble;
		break;
	default:
		break;
	}
        


	return 0;
	

}

int
ConcreteMcftNonLinear5::activateParameter(int passedParameterID)
{  
	
	//opserr << " check24 " << endln;
	parameterID = passedParameterID;

	return 0;
	

}

const Vector&
ConcreteMcftNonLinear5::getStressSensitivity(int gradNumber, bool conditional)
{
  static Vector zerodsigdh(2);

  if (  parameterID == 1 ) {
    //opserr << " check25 " << endln;
    return dsigfdfcu;
  }
  if ( parameterID == 2 ) {
    //  opserr << " check25 " << endln;
    return dsigfdRoV;
  }

  return zerodsigdh;
}



int
ConcreteMcftNonLinear5::commitSensitivity(const Vector &strainGradient, int gradNumber, int numGrads)
{
//opserr << " check27 " << endln;
	return 0;
	

}

// AddingSensitivity:END /////////////////////////////////////////////








int 
ConcreteMcftNonLinear5::getResponse (int responseID, Information &matInformation)
{
//opserr << " check28 " << endln;
	static Vector crackInfo(5);
	if (responseID == 10) {
		
		crackInfo(0) = epsf(0);
		crackInfo(1) = epsf(1);
		crackInfo(2) = epsy;
		crackInfo(3) = FinalAnglex;
		crackInfo(4) = crackLabel;



		matInformation.setVector(crackInfo);
	} 
	static Vector prinStress(8);
	if (responseID == 11) {
		
		prinStress(0) = Sigma1;
		prinStress(1) = Sigma2;
		prinStress(2) = fx;
		prinStress(3) = fxy;
		prinStress(4) = t00;
		prinStress(5) = t01;
		//prinStress(6) = t10;
		//prinStress(7) = t11;
		prinStress(6) = dsigfdfcu(0);
		prinStress(7) = dsigfdfcu(1);
		matInformation.setVector(prinStress);
	}

	return 0;
	

}


double ConcreteMcftNonLinear5::c1tmd00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{									  
	//checkpoint
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d00;

	if(e1 > fcr/Ec) {

	 d00 = -((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2))) + 
   (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
      )/2. - ((Esv*RoV - (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))\
         - (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
                (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/
         2. + cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 
              2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
               pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),
               2))))/2. - (pow(sect,2)*sin(2*theta)*
         (-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2.
        - cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	d00 =-((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2))) + 
   (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*
           (Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/
         (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
              ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*
           (Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
              pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/
            (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/
          (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
	}

		return d00;
}

double ConcreteMcftNonLinear5::c1tmd01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{


	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d01;	

   if(e1 > fcr/Ec) {

	 d01 = (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
   (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
            pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
      )/2. - ((Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))\
         - (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
                (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/
         2. + cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 
              2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
               pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),
               2))))/2. - (pow(sect,2)*sin(2*theta)*
         (-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2.
        - cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	d01 = (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
      pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
   (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((Ec*cott)/2. - 
        (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - 
             (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)
       *((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/
         (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
              ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*
           (Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
              pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/
            (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/
          (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
   }


		return d01;
}

double ConcreteMcftNonLinear5::c1tmd10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d10;

	if(e1 > fcr/Ec) {

	 d10 = (sin(2*theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
      )/2. - ((Esv*RoV - (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))\
         - (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 
              2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
               pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),
               2))))/2. - (pow(sect,2)*sin(2*theta)*
         (-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2.
        - cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	 d10 = (sin(2*theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*
           (Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*
              ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
              pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/
            (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/
          (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
	
	}
	
	return d10;
}

double ConcreteMcftNonLinear5::c1tmd11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d11;
	

if(e1 > fcr/Ec) {

		d11 = (sin(2*theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
            pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
      )/2. - ((Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))\
         - (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*
                     ((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 
              2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
               pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),
               2))))/2. - (pow(sect,2)*sin(2*theta)*
         (-((fcu*nE*(ex - (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2.
        - cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

		d11 = (sin(2*theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
           pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   (((Ec*cott)/2. + Esv*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - 
             (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)
       *((sin(2*theta)*(Ec*pow(cott,2)*
              ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*(-(exy*pow(sect,2))/2. + 
         pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            ((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
              pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/
            (2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
            ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          ((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/
          (ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));

}

if(exy<0) 
	d11 = -1*d11;
	

		return d11;
}

double ConcreteMcftNonLinear5::tm3(double e2, double ecu, double fcu, double nE)
{
	double df2de2 = -((e2*pow(e2/ecu,-1 + nE)*fcu*pow(nE,2))/
      (pow(ecu,2)*pow(-1 + pow(e2/ecu,nE) + nE,2))) + 
   (fcu*nE)/(ecu*(-1 + pow(e2/ecu,nE) + nE));

	return df2de2;
}

double ConcreteMcftNonLinear5::c2tmd00(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d00;

if(e1 > fcr/Ec) {

	d00 =-((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2))) + 
   (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.\
    - ((Esv*RoV - (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2.)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
            pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2. + (pow(sect,2)*sin(2*theta)*
           (-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*
          (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/
       2. - (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	d00= -((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2))) + 
   (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*
           (Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
            pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
              (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
              (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
        ))/(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 
         2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + 
         pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - (pow(sect,2)*sin(2*theta)*
         (Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));
}

	//if(d00 < 0.0)
	//	d00 =0.0;

	return d00;
}

double ConcreteMcftNonLinear5::c2tmd01(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d01;

if(e1 > fcr/Ec) {

	 d01 = -(fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
   (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.\
    - ((Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2.)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
            pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2. + (pow(sect,2)*sin(2*theta)*
           (-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*
          (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/
       2. - (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	d01= -(fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
   (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + 
        (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + 
             (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
            pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
              (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*
              (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
        ))/(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 
         2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + 
         pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - (pow(sect,2)*sin(2*theta)*
         (Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));


}
	return d01;
}

double ConcreteMcftNonLinear5::c2tmd10(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d10;
	if(e1 > fcr/Ec) {

	 d10 = (sin(2*theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.\
    - ((Esv*RoV - (5*sqrt(5.0)*fcr)/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2.)*((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2. + cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*
          (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/
       2. - (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	 d10= (sin(2*theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Esv*RoV - (sin(2*theta)*tan(theta)*
           (Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*
              (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
        ))/(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 
         2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + 
         pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - (pow(sect,2)*sin(2*theta)*
         (Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));

	}
	//if(d10 < 0.0)
	//	d10 =0.0;
	return d10;

}

double ConcreteMcftNonLinear5::c2tmd11(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Esv, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double d11;

if(e1 > fcr/Ec) {

	 d11 = (sin(2*theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.\
    - ((Esv*RoV*(-cott/2. + tan(theta)/2.) + 
        (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2.)*((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*
                   (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)))
           )/2. + cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Esv*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*
          (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*
                 (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/
       2. - (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
              (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/
            (ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	 d11 = (sin(2*theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
           pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((-(Ec*cott)/2. + Esv*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + 
             (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*
              (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/
              (2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
        ))/(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 
         2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Esv*RoV*((exy*pow(sect,2))/2. + 
         pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*
            (-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
              pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - (pow(sect,2)*sin(2*theta)*
         (Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))
         )/2. - cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*
          (-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))) ;

}
if(exy<0) 
	d11 = -1*d11;

	return d11;

}







double ConcreteMcftNonLinear5::c1dd00dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	//checkpoint
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd00dfcu;

	if(e1 > fcr/Ec) {
	
	dd00dfcu = -((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2))) + nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   (((exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	dd00dfcu = -((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2))) + nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))) + 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))));
	
	}

		return dd00dfcu;
}

double ConcreteMcftNonLinear5::c1dd01dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{


	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd01dfcu;	

   if(e1 > fcr/Ec) {

	 dd01dfcu=(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. - (((exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        )*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	dd01dfcu = (pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. - (((exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        )*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))));
   }


		return dd01dfcu;
}

double ConcreteMcftNonLinear5::c1dd10dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd10dfcu;

	if(e1 > fcr/Ec) {

	 dd10dfcu = (sin(2*theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   ((-((nE*cos(2*theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	 dd10dfcu=(sin(2*theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-((nE*cos(2*theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))));
	
	}
	
	return dd10dfcu;
}

double ConcreteMcftNonLinear5::c1dd11dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd11dfcu;
	

if(e1 > fcr/Ec) {

		dd11dfcu = (sin(2*theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. + (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   ((-((nE*cos(2*theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        )*((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

		dd11dfcu=(sin(2*theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. + (((nE*pow(sect,2)*sin(2*theta)*(ex - (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   ((-((nE*cos(2*theta)*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*(-(exy*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))) + 
   (sin(2*theta)*tan(theta)*(-(pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))
        )*((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))));

}

if(exy<0) 
	dd11dfcu = -1*dd11dfcu;
	

		return dd11dfcu;
}


double ConcreteMcftNonLinear5::c2dd00dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd00dfcu;

if(e1 > fcr/Ec) {

	dd00dfcu = -((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2))) + nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((-(exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	dd00dfcu = -((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
      (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2))) + nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))) + 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))));
}

	
	return dd00dfcu;
}

double ConcreteMcftNonLinear5::c2dd01dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd01dfcu;

if(e1 > fcr/Ec) {

	dd01dfcu = -(pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. - ((-(exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))
        )*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

		dd01dfcu = -(pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
    (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
   ((-(exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))) + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))));


}
	return dd01dfcu;
}

double ConcreteMcftNonLinear5::c2dd10dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd10dfcu;
	if(e1 > fcr/Ec) {

	 dd10dfcu = (sin(2*theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   ((-((nE*cos(2*theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      ((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	dd10dfcu=(sin(2*theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   ((Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-((nE*cos(2*theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - nE/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))));

	}
	
	return dd10dfcu;

}

double ConcreteMcftNonLinear5::c2dd11dfcu(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd11dfcu;

if(e1 > fcr/Ec) {

	dd11dfcu = (sin(2*theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))
        ))/2. + (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   ((-((nE*cos(2*theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))
        )*((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (2.*(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
         (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))));

	} else {

	dd11dfcu = (sin(2*theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
   (((nE*pow(sect,2)*sin(2*theta)*(ex + (exy*tan(theta))/2.))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (nE*cos(2*theta)*tan(theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
        (sin(2*theta)*tan(theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   ((-((nE*cos(2*theta)*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
        (sin(2*theta)*((exy*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))) + 
   (sin(2*theta)*tan(theta)*((pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (2.*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
        Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))));

}
if(exy<0) 
	dd11dfcu = -1*dd11dfcu;

	return dd11dfcu;

}


// RoV Parameter
double ConcreteMcftNonLinear5::c1dd00dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	//checkpoint
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd00dRoV;

	if(e1 > fcr/Ec) {
	
	dd00dRoV = (Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	dd00dRoV =(Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
	
	}

		return dd00dRoV;
}

double ConcreteMcftNonLinear5::c1dd01dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{


	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd01dRoV;	

   if(e1 > fcr/Ec) {

	 dd01dRoV= (Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*(cott/2. - tan(theta)/2.)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	dd01dRoV = (Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*(cott/2. - tan(theta)/2.)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
   }


		return dd01dRoV;
}

double ConcreteMcftNonLinear5::c1dd10dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd10dRoV;

	if(e1 > fcr/Ec) {

	 dd10dRoV = (Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	 dd10dRoV=(Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));
	
	}
	
	return dd10dRoV;
}

double ConcreteMcftNonLinear5::c1dd11dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd11dRoV;
	

if(e1 > fcr/Ec) {

		dd11dRoV = (Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV*(cott/2. - tan(theta)/2.) - (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*(-(fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*(cott/2. - tan(theta)/2.)*((sin(2*theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*
                 pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/(2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

		dd11dRoV =(Es*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      ((Ec*cott)/2. + Es*RoV*(cott/2. - tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*((Ec*cott)/2. - (fcu*pow(nE,2)*tan(theta)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*(cott/2. - tan(theta)/2.)*((sin(2*theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*(-(exy*pow(sect,2))/2. + pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex - (exy*tan(theta))/2.)*pow((ex - (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE),2)) + 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*((exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex - (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex - (exy*tan(theta))/2.)/ecu,nE)))));

}

if(exy<0) 
	dd11dRoV = -1*dd11dRoV;
	

		return dd11dRoV;
}


double ConcreteMcftNonLinear5::c2dd00dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd00dRoV;

if(e1 > fcr/Ec) {

	dd00dRoV = (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	dd00dRoV =  (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));
}

	
	return dd00dRoV;
}

double ConcreteMcftNonLinear5::c2dd01dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd01dRoV;

if(e1 > fcr/Ec) {

	dd01dRoV = (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*(-cott/2. + tan(theta)/2.)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
             fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. + 
        cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))) ;

	} else {

		dd01dRoV =(Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      (-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*(-cott/2. + tan(theta)/2.)*(-(exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
         (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) + 
        (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
        (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
             (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));


}
	return dd01dRoV;
}

double ConcreteMcftNonLinear5::c2dd10dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd10dRoV;
	if(e1 > fcr/Ec) {

	 dd10dRoV = (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV - (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr)/(sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	dd10dRoV= (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Ec + Es*RoV - (sin(2*theta)*tan(theta)*(Ec + (fcu*pow(nE,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - (fcu*nE)/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));

	}
	
	return dd10dRoV;

}

double ConcreteMcftNonLinear5::c2dd11dRoV(double ex, double exy, double theta, double Ec, double nE, double fcu, double ecu, double e1, double fcr, double Es, double RoV)
{
	
	double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);
	double dd11dRoV;

if(e1 > fcr/Ec) {

	dd11dRoV = (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (Es*RoV*(-cott/2. + tan(theta)/2.) + (5*sqrt(5.0)*fcr*cott)/
         (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
           pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
        (sin(2*theta)*tan(theta)*((fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) + 
             (5*sqrt(5.0)*fcr*cott)/
              (2.*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2.)*
      ((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    pow(Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))),2) - 
   (Es*(-cott/2. + tan(theta)/2.)*((sin(2*theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*
                pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/(2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
             (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                  2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
              (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
                pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. + 
        cos(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))))))/
    (Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
       (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
         pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2)) - 
      (sin(2*theta)*tan(theta)*((exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))) - 
           (5*sqrt(5.0)*fcr*(pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
                2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))/
            (sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
              pow(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))),2))))/2. - 
      (pow(sect,2)*sin(2*theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
           fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))))/2. - 
      cos(2*theta)*tan(theta)*(-((fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))) + 
         fcr/(1 + 10*sqrt(5.0)*sqrt(pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))))));

	} else {

	dd11dRoV = (Es*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
        2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)))*
      (-(Ec*cott)/2. + Es*RoV*(-cott/2. + tan(theta)/2.) - 
        (sin(2*theta)*tan(theta)*(-(Ec*cott)/2. + (fcu*pow(nE,2)*tan(theta)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (fcu*nE*tan(theta))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2.)*
      ((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    pow(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))),2) - 
   (Es*(-cott/2. + tan(theta)/2.)*((sin(2*theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
             2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
             (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
              (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
             (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. + 
        cos(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE))))))/
    (Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
      2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
      Es*RoV*((exy*pow(sect,2))/2. + pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
         2*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2))) - 
      (sin(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*pow(sect,2))/2. + 2*ex*pow(sect,2)*tan(theta)) - 
           2*Ec*cott*pow(csct,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) + 
           (exy*fcu*pow(nE,2)*pow(sect,2)*(ex + (exy*tan(theta))/2.)*pow((ex + (exy*tan(theta))/2.)/ecu,-1 + nE))/
            (2.*pow(ecu,2)*pow(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE),2)) - 
           (exy*fcu*nE*pow(sect,2))/(2.*ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      (pow(sect,2)*sin(2*theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
           (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))))/2. - 
      cos(2*theta)*tan(theta)*(Ec*pow(cott,2)*(-(exy*tan(theta))/2. + ex*pow(tan(theta),2)) - 
         (fcu*nE*(ex + (exy*tan(theta))/2.))/(ecu*(-1 + nE + pow((ex + (exy*tan(theta))/2.)/ecu,nE)))));

}
if(exy<0) 
	dd11dRoV = -1*dd11dRoV;

	return dd11dRoV;

}

	/*double cott = 1/tan(theta);
	double sect = 1/cos(theta);
	double csct = 1/sin(theta);*/
