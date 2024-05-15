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
                                                                        
// $Revision: 1.00 $
// $Date: 2021-Jan-14 12:15:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/HystereticPoly.cpp,v $
                                                                        
// Written: Salvatore Sessa Mail: salvatore.sessa2@unina.it
// Created: 01-2021
// Revision: A
//
// Description: This file contains the class implementation for 
// HystereticPoly. 
//
// What: "@(#) HystereticPoly.C, revA"


#include <HystereticPoly.h>
#include <Vector.h>
#include <Matrix.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <math.h>
#include <float.h>

#include <elementAPI.h>
#include <MaterialResponse.h>



void *
OPS_HystereticPoly()
{
	UniaxialMaterial *theMaterial = 0;
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs < 6) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: uniaxialMaterial HystereticPoly tag? ka? kb? a? b1? b2? <tol?>" << endln;
		return 0;
	}


	int iData[1];
	double dData[6];
	dData[5] = 1.0e-20;		// default tolerance


	int numData = 1;
	if (OPS_GetIntInput(&numData, iData) != 0) {
		opserr << "WARNING invalid tag for uniaxialMaterial HystereticPoly" << endln;
		return 0;
	}

	numData = numArgs - 1;
	if (OPS_GetDoubleInput(&numData, dData) != 0) {
		opserr << "Invalid data for uniaxial HystereticPoly " << iData[0] << endln;
		return 0;
	}
	if (dData[0] <= 0.0) {
		opserr << "uniaxialMaterial HystereticPoly ka must be positive" << endln;
		return 0;
	}
	if (dData[1] >= dData[0]) {
		opserr << "uniaxialMaterial HystereticPoly kb must be < ka" << endln;
		return 0;
	}
	if (dData[2] <= 0.0) {
		opserr << "uniaxialMaterial HystereticPoly a must be positive and <> 1" << endln;
		return 0;
	}
	if (dData[2] == 1.0) {
		opserr << "uniaxialMaterial HystereticPoly a must be positive and <> 1" << endln;
		return 0;
	}


	// Parsing was successful, allocate the material
	theMaterial = new HystereticPoly(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);
	if (theMaterial == 0) {
		opserr << "WARNING could not create uniaxialMaterial of type HystereticPoly\n";
		return 0;
	}

	return theMaterial;

}

HystereticPoly::HystereticPoly(int tag, double K1, double K2, double A, double C, double D, double TOL)
   :UniaxialMaterial(tag, MAT_TAG_HystereticPoly),
   k1(K1), k2(K2), a(A), c(C), d(D), tol(TOL)
{
   // Sets all history and state variables to initial values

   // History variables

   sc = 1.0;
   st = 1.0;
	
	// Parameters:

	uo = 0.5*(pow((k1-k2)/tol, (1/a))-1);
	Fbar = 0.5*(k1 - k2)*(pow(1+2*uo, (1-a))-1) / (1 - a);

	k12 = k1 - k2;
	Uoa = pow(1 + 2 * uo, 1 - a);
	am1 = 1 / (1 - a);
	Uo = 1 + 2 * uo;
	a1 = 1 - a;

	uj = Uo - pow(a1 / k12*(  k12*am1*Uoa / st - Fbar), am1);


   


   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = k2 + k12 * pow(1-uj+2*uo, -1*a);

   
   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Ctangent;

   
// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////
}

HystereticPoly::HystereticPoly():UniaxialMaterial(0, MAT_TAG_HystereticPoly),
 k1(0.0), k2(0.0), a(0.0), c(0.0), d(0.0), tol(0.0)
{

// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
// AddingSensitivity:END //////////////////////////////////////

}

HystereticPoly::~HystereticPoly ()
{
// AddingSensitivity:BEGIN /////////////////////////////////////
	if (SHVs != 0) 
		delete SHVs;
// AddingSensitivity:END //////////////////////////////////////
}

int HystereticPoly::setTrialStrain (double strain, double strainRate)
{
   // Determine change in strain from last converged state

   Tstrain = strain;
   dStrain = strain - Cstrain;
   st = signum(dStrain);	



   
   
   
   uj = Cstrain + st*Uo - st*pow(st*a1 / k12*(Cstress - c*pow(Cstrain , 3) - d*pow(Cstrain, 5) - k2*Cstrain - st*Fbar + k12*am1*Uoa / st ),am1);
   
   
   Tstress = c*pow(Tstrain, 3) + d*pow(Tstrain, 5) + k2*Tstrain + k12*(pow(1 + st*Tstrain - st*uj + 2 * uo, 1 - a) / st / (1 - a) - am1*Uoa/st ) + st*Fbar;

   
   Ttangent = 3 * c*pow(Tstrain, 2) + 5 * d*pow(Tstrain, 4) + k2 + k12*pow(1 + st*Tstrain - st*uj + 2 * uo, -1 * a);
   

   return 0;
}



double HystereticPoly::getStrain ()
{
   return Tstrain;
}

double HystereticPoly::getStress ()
{
   return Tstress;
}

double HystereticPoly::getTangent ()
{
   return Ttangent;
}

int HystereticPoly::commitState ()
{
   // State variables
   Cstrain = Tstrain;
   Cstress = Tstress;
   Ctangent = Ttangent;

   return 0;
}

int HystereticPoly::revertToLastCommit ()
{

   // State variables
   Tstrain = Cstrain;
   Tstress = Cstress;
   Ttangent = Ctangent;

   return 0;
}

int HystereticPoly::revertToStart ()
{
	uj = Uo - pow(a1 / k12*(k12*am1*Uoa / st - Fbar), am1);


   // State variables
   Cstrain = 0.0;
   Cstress = 0.0;
   Ctangent = k2 + k12 * pow(1 - uj + 2 * uo, -1 * a);


   Tstrain = 0.0;
   Tstress = 0.0;
   Ttangent = Ctangent;

   

// AddingSensitivity:BEGIN /////////////////////////////////
	if (SHVs != 0) 
		SHVs->Zero();
// AddingSensitivity:END //////////////////////////////////

   return 0;
}

UniaxialMaterial* HystereticPoly::getCopy ()
{
   HystereticPoly* theCopy = new HystereticPoly(this->getTag(), k1, k2, a, c, d, tol);

   // State variables
   theCopy->Cstrain = Cstrain;
   theCopy->Cstress = Cstress;
   theCopy->Ctangent = Ctangent;

   theCopy->Tstrain = Tstrain;
   theCopy->Tstress = Tstress;
   theCopy->Ttangent = Ttangent;

   return theCopy;
}

int HystereticPoly::sendSelf (int commitTag, Channel& theChannel)
{
   int res = 0;
   static Vector data(10);
   data(0) = this->getTag();

   // Material properties
   data(1) = k1;
   data(2) = k2;
   data(3) = a;
   data(4) = c;
   data(5) = d;
   
   // State variables from last converged state
   data(6) = Cstrain;
   data(7) = Cstress;
   data(8) = Ctangent;

   // New variables
   data(9) = tol;

   // Data is only sent after convergence, so no trial variables
   // need to be sent through data vector

   res = theChannel.sendVector(this->getDbTag(), commitTag, data);
   if (res < 0) 
      opserr << "HystereticPoly::sendSelf() - failed to send data\n";

   return res;
}

int HystereticPoly::recvSelf (int commitTag, Channel& theChannel,
                                FEM_ObjectBroker& theBroker)
{
   int res = 0;
   static Vector data(10);
   res = theChannel.recvVector(this->getDbTag(), commitTag, data);
  
   if (res < 0) {
      opserr << "HystereticPoly::recvSelf() - failed to receive data\n";
      this->setTag(0);      
   }
   else {
      this->setTag(int(data(0)));

   k1 = data(1);
   k2 = data(2);
   a = data(3);
   c = data(4);
   d = data(5);
   
   // State variables from last converged state
   Cstrain = data(6);
   Cstress = data(7);
   Ctangent = data(8);

   // New variables
   tol = data(9);

      // Copy converged state values into trial values
      Tstrain = Cstrain;
      Tstress = Cstress;
      Ttangent = Ctangent;
   }
    
   return res;
}


void HystereticPoly::Print (OPS_Stream& s, int flag)
{
   s << "HystereticPoly tag: " << this->getTag() << endln;
   /*s << "  ka: " << k1 << " ";
   s << "  kb: " << k2 << " ";
   s << "  a:  " << a << " ";
   s << "  b1: " << c << " ";
   s << "  b2: " << d << " ";
   s << " tol:" << tol << " ";*/
   s << " strain: " << this->getStrain() << endln;
   s << " stress: " << this->getStress() << endln;
   s << " tangent: " << this->getTangent() << endln;
}




// AddingSensitivity:BEGIN ///////////////////////////////////
int
HystereticPoly::setParameter(const char **argv, int argc, Parameter &param)
{

  if (strcmp(argv[0],"ka") == 0 )
    return param.addObject(1, this);
  
  if (strcmp(argv[0],"kb") == 0)
    return param.addObject(2, this);
  
  if (strcmp(argv[0],"a") == 0)
    return param.addObject(3, this);
  
  if (strcmp(argv[0],"b1") == 0)
    return param.addObject(4, this);
  
  if (strcmp(argv[0],"b2") == 0)
    return param.addObject(5, this);

  if (strcmp(argv[0], "tol") == 0)
	  return param.addObject(6, this);
  

  return -1;
}



int
HystereticPoly::updateParameter(int parameterID, Information &info)
{
	switch (parameterID) {
	case -1:
		return -1;
	case 1:
		this->k1 = info.theDouble;
		break;
	case 2:
		this->k2 = info.theDouble;
		break;
	case 3:
		this->a = info.theDouble;
		break;
	case 4:
		this->c = info.theDouble;
		break;
	case 5:
		this->d = info.theDouble;
		break;
	default:
		return -1;
	}

	Ttangent = k1;          // Initial stiffness

	return 0;
}




int
HystereticPoly::activateParameter(int passedParameterID)
{
	parameterID = passedParameterID;

	return 0;
}

double 
HystereticPoly::getStrainSensitivity(int gradNumber) {
	if (SHVs == 0) {
		opserr << "warning:HystereticPoly::getStrainsSensitivity, SHVs =0 " << endln;
		return 0.0;
	}
	double  Tsensitivity = (*SHVs)(0, (gradNumber));
	return Tsensitivity;

}

double
HystereticPoly::getStressSensitivity(int gradIndex, bool conditional)
{
	double gradient = 0.0;
	
	// Initialize return value
	


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = 0.0;
	if (SHVs != 0) {
		Duc = (*SHVs)(0,gradIndex);
		Dfc = (*SHVs)(1,gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	Dk1 = 0.0;
	Dk2 = 0.0;
	Da = 0.0;
	Dc = 0.0;
	Dd = 0.0;
	Dtol = 0.0;

	if (parameterID == 1) {
		Dk1 = 1.0;
	}
	else if (parameterID == 2) {
		Dk2 = 1.0;
	}
	else if (parameterID == 3) {
		Da = 1.0;
	}
	else if (parameterID == 4) {
		Dc = 1.0;
	}
	else if (parameterID == 5) {
		Dd = 1.0;
	}
	else if (parameterID == 6) {
		Dtol = 1.0;
	}


	// This is the Tstress derivative with respect to Tstrain

	double Delta1 = st*(1 - a) / (k1 - k2);
	double Delta3 = (k1 - k2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a);
	double Delta2 = Cstress - c*pow(Cstrain, 3) - d*pow(Cstrain, 5) - k2*Cstrain - st*Fbar + Delta3;
	uj = Cstrain + st*(1 + 2 * uo) - st*pow(Delta1*Delta2, 1 / (1 - a));
	double Delta4 = pow(1 + st*Tstrain - st*uj + 2 * uo, 1 - a) / st / (1 - a);
	double Delta5 = (pow(1 + 2 * uo, 1 - a) ) / (1 - a);
	Tstress = c*pow(Tstrain, 3) + d*pow(Tstrain, 5) + k2*Tstrain + (k1 - k2)*(Delta4 - Delta5) + st*Fbar;

	double Duo = 0.5*uo / a * (tol / (k1 - k2) * ((Dk1 - Dk2) / tol - (k1 - k2) / pow(tol, 2)*Dtol) - Da / a* log((k1 - k2) / tol));
	double DFbar = 0.5*(Dk1 - Dk2) * (pow(1 + 2 * uo, 1 - a) - 1) / (1 - a) +
		0.5*(k1 - k2) * (pow(1 + 2 * uo, 1 - a) / (1 - a) * (2*Duo*(1-a)/(1+2*uo) - Da*log(1+2*uo) ) + Da*(pow(1+2*uo,1-a)-1)/pow(1-a,2) );

	double DD1 = (st*Da * (k1 - k2) + st*(1 - a)*(Dk1 - Dk2)) / pow(k1 - k2, 2);
	double DD3 =	(Dk1 - Dk2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a) + (k1 - k2)*pow(1 + 2 * uo, 1 - a) / pow(st, 2) / pow(1 - a, 2)*Da +
					Delta3 * (2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo));

	double DD2 = Dfc - Dc*pow(Cstrain, 3) - Dd * pow(Cstrain, 5) - Dk2*Cstrain - (3 * c*pow(Cstrain, 2) + 5 * d*pow(Cstrain, 4) + k2)*Duc - st*DFbar + DD3;

	double Duj = Duc + 2 * st*Duo - st*pow(Delta1*Delta2, 1 / (1 - a)) * (Da*log(Delta1*Delta2) / pow(1 - a, 2) + (DD1*Delta2 + Delta1*DD2) / (Delta1*Delta2*(1 - a)));

	double DD4 = pow(1 + st*Tstrain - st*uj + 2 * uo, 1 - a) / pow(st, 2) / pow(1 - a, 2)*Da + Delta4 * (1 - a)*(st*Dut - st*Duj + 2 * Duo) / (1 + st*Tstrain - st*uj + 2 * uo);
	double DD5 = (pow(1 + 2 * uo, 1 - a)) / pow(1 - a, 2)*Da + pow(1 + 2 * uo, 1 - a) / (1 - a)*(2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo));

	
	gradient = Dc*pow(Tstrain, 3) + Dd*pow(Tstrain, 5) + Dk2*Tstrain + (3 * c*pow(Tstrain, 2) + 5 * d*pow(Tstrain, 4) + k2)*Dut +
		st*DFbar + (Dk1 - Dk2)*(Delta4 - Delta5) + (k1 - k2)*(DD4 - DD5);
	
	
	return gradient;
}




double
HystereticPoly::getInitialTangentSensitivity(int gradIndex)
{
	
	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = 0.0;


	// Assign values to parameter derivatives (depending on what's random)
	Dk1 = 0.0;
	Dk2 = 0.0;
	Da = 0.0;
	Dc = 0.0;
	Dd = 0.0;
	Dtol = 0.0;

	if (parameterID == 1) {
		Dk1 = 1.0;
	}
	else if (parameterID == 2) {
		Dk2 = 1.0;
	}
	else if (parameterID == 3) {
		Da = 1.0;
	}
	else if (parameterID == 4) {
		Dc = 1.0;
	}
	else if (parameterID == 5) {
		Dd = 1.0;
	}
	else if (parameterID == 6) {
		Dtol = 1.0;
	}


	double Duo = 0.5*uo / a * (tol / (k1 - k2) * ((Dk1 - Dk2) / tol - (k1 - k2) / pow(tol, 2)*Dtol) - Da / a* log((k1 - k2) / tol));
	double DFbar = 0.5*(Dk1 - Dk2) * (pow(1 + 2 * uo, 1 - a) - 1) / (1 - a) +
		0.5*(k1 - k2) * (pow(1 + 2 * uo, 1 - a) / (1 - a) * (2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo)) + Da*(pow(1 + 2 * uo, 1 - a) - 1) / pow(1 - a, 2));
	double Delta1 = st*(1 - a) / (k1 - k2);
	double Delta3 = (k1 - k2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a);
	double Delta2 = -1* st*Fbar + Delta3;
	double DD1 = (st*Da * (k1 - k2) + st*(1 - a)*(Dk1 - Dk2)) / pow(k1 - k2, 2);
	double DD3 = (Dk1 - Dk2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a) + (k1 - k2)*pow(1 + 2 * uo, 1 - a) / pow(st, 2) / pow(1 - a, 2)*Da +
		Delta3 * (2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo));

	double DD2 = Dfc -  k2*Duc - st*DFbar + DD3;
	double Duj = Duc + 2 * st*Duo - st*pow(Delta1*Delta2, 1 / (1 - a)) * (Da*log(Delta1*Delta2) / pow(1 - a, 2) + (DD1*Delta2 + Delta1*DD2) / (Delta1*Delta2*(1 - a)));

	gradient = (st * (Dk1 - Dk2)*(1 - a) - st*Da*(k1 - k2))*pow(1 - st*uj + 2 * uo, -1 * a) -
		st*(k1 - k2)*(1 - a)*pow(1 - st*uj + 2 * uo, -1 * a) * (Da * log(1 - st*uj + 2 * uo) + a / (1 - st*uj + 2 * uo)*(2 * Duo - st*Duj));



	return gradient;
}


int
HystereticPoly::commitSensitivity(double TstrainSensitivity, int gradIndex, int numGrads)
{
	if (SHVs == 0) {
		SHVs = new Matrix(2,numGrads);
	}


	// Initialize return value
	double gradient = 0.0;


	// Pick up sensitivity history variables
	Duc = 0.0;
	Dfc = 0.0;
	Dut = TstrainSensitivity;
	if (SHVs != 0) {
		Duc = (*SHVs)(0, gradIndex);
		Dfc = (*SHVs)(1, gradIndex);
	}


	// Assign values to parameter derivatives (depending on what's random)
	Dk1 = 0.0;
	Dk2 = 0.0;
	Da = 0.0;
	Dc = 0.0;
	Dd = 0.0;
	Dtol = 0.0;

	if (parameterID == 1) {
		Dk1 = 1.0;
	}
	else if (parameterID == 2) {
		Dk2 = 1.0;
	}
	else if (parameterID == 3) {
		Da = 1.0;
	}
	else if (parameterID == 4) {
		Dc = 1.0;
	}
	else if (parameterID == 5) {
		Dd = 1.0;
	}
	else if (parameterID == 6) {
		Dtol = 1.0;
	}


	// This is the Tstress derivative with respect to Tstrain

	double Delta1 = st*(1 - a) / (k1 - k2);
	double Delta3 = (k1 - k2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a);
	double Delta2 = Cstress - c*pow(Cstrain, 3) - d*pow(Cstrain, 5) - k2*Cstrain - st*Fbar + Delta3;
	uj = Cstrain + st*(1 + 2 * uo) - st*pow(Delta1*Delta2, 1 / (1 - a));
	double Delta4 = pow(1 + st*Tstrain - st*uj + 2 * uo, 1 - a) / st / (1 - a);
	double Delta5 = (pow(1 + 2 * uo, 1 - a)) / (1 - a);
	Tstress = c*pow(Tstrain, 3) + d*pow(Tstrain, 5) + k2*Tstrain + (k1 - k2)*(Delta4 - Delta5) + st*Fbar;

	double Duo = 0.5*uo / a * (tol / (k1 - k2) * ((Dk1 - Dk2) / tol - (k1 - k2) / pow(tol, 2)*Dtol) - Da / a* log((k1 - k2) / tol));
	double DFbar = 0.5*(Dk1 - Dk2) * (pow(1 + 2 * uo, 1 - a) - 1) / (1 - a) +
		0.5*(k1 - k2) * (pow(1 + 2 * uo, 1 - a) / (1 - a) * (2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo)) + Da*(pow(1 + 2 * uo, 1 - a) - 1) / pow(1 - a, 2));

	double DD1 = (st*Da * (k1 - k2) + st*(1 - a)*(Dk1 - Dk2)) / pow(k1 - k2, 2);
	double DD3 = (Dk1 - Dk2)*pow(1 + 2 * uo, 1 - a) / st / (1 - a) + (k1 - k2)*pow(1 + 2 * uo, 1 - a) / pow(st, 2) / pow(1 - a, 2)*Da +
		Delta3 * (2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo));

	double DD2 = Dfc - Dc*pow(Cstrain, 3) - Dd * pow(Cstrain, 5) - Dk2*Cstrain - (3 * c*pow(Cstrain, 2) + 5 * d*pow(Cstrain, 4) + k2)*Duc - st*DFbar + DD3;

	double Duj = Duc + 2 * st*Duo - st*pow(Delta1*Delta2, 1 / (1 - a)) * (Da*log(Delta1*Delta2) / pow(1 - a, 2) + (DD1*Delta2 + Delta1*DD2) / (Delta1*Delta2*(1 - a)));

	double DD4 = pow(1 + st*Tstrain - st*uj + 2 * uo, 1 - a) / pow(st, 2) / pow(1 - a, 2)*Da + Delta4 * (1 - a)*(st*Dut - st*Duj + 2 * Duo) / (1 + st*Tstrain - st*uj + 2 * uo);
	double DD5 = (pow(1 + 2 * uo, 1 - a)) / pow(1 - a, 2)*Da + pow(1 + 2 * uo, 1 - a) / (1 - a)*(2 * Duo*(1 - a) / (1 + 2 * uo) - Da*log(1 + 2 * uo));


	gradient = Dc*pow(Tstrain, 3) + Dd*pow(Tstrain, 5) + Dk2*Tstrain + (3 * c*pow(Tstrain, 2) + 5 * d*pow(Tstrain, 4) + k2)*Dut +
		st*DFbar + (Dk1 - Dk2)*(Delta4 - Delta5) + (k1 - k2)*(DD4 - DD5);

	// Commit 
	(*SHVs)(0,gradIndex) = TstrainSensitivity;
	(*SHVs)(1,gradIndex) = gradient;
	
	return 0;
}

// AddingSensitivity:END /////////////////////////////////////////////


double 
HystereticPoly::signum(double value)
{
	if (value > 0.0) {
		return 1.0;
	}
	else if (value == 0){
		return 1.0;
	}
	else {
		return -1.0;
	}
}
