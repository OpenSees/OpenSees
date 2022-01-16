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

// $Source$

// Written: HyungSuk Shin
// Created: Sept. 2004
//
// Description: This file contains the class definition for 
// pyUCLA.  pyUCLA provides the abstraction
// for a one-dimensional p-y contact material. 
// This material is based on model of E.Taciroglu, C.Rha, J.Wallace, UCLA. 

#include <pyUCLA.h>
#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <Matrix.h>
#include <Information.h>
#include <float.h>


#include <elementAPI.h>
#define OPS_Export 

static int num_pyUCLA = 0;

OPS_Export void * OPS_ADD_RUNTIME_VPV(OPS_pyUCLA)
{
  if (num_pyUCLA == 0) {
    num_pyUCLA++;
    //OPS_Error("pyUCLAMaterial unaxial material - Written by H.Shin, P.Arduino, U.Washington\n based on model of E.Taciroglu, C.Rha, J.Wallace, UCLA", 1);
    opserr << "pyUCLAMaterial unaxial material - Written by H.Shin, P.Arduino, U.Washington\n based on model of E.Taciroglu, C.Rha, J.Wallace, UCLA\n";
  }

  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  if (OPS_GetNumRemainingInputArgs() != 5) {
    opserr << "Invalid #args,  want: uniaxialMaterial pyUCLA tag? soilType? pult? y50? Cd? " << endln;
    return 0;
  }
    
  int    iData[2];
  double dData[3];

  int numData = 2;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag or soilType uniaxialMaterial pyUCLAMaterial" << endln;
    return 0;
  }

  numData = 3;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid pyData data for material uniaxial pyUCLA " << iData[0] << endln;
    return 0;	
  }

  // Parsing was successful, allocate the material
  theMaterial = new pyUCLA(iData[0], iData[1], dData[0], dData[1], dData[2]);

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type pyUCLAMaterial\n";
    return 0;
  }

  return theMaterial;
}


pyUCLA::pyUCLA(int tag, int type, double Pu,
	       double yc, double CdRatio)
:UniaxialMaterial(tag, MAT_TAG_pyUCLA),
 soilType(type), pult(Pu), y50(yc), Cd(CdRatio)
{
	// Initialize variables
    this->revertToStart();
}

pyUCLA::pyUCLA()
:UniaxialMaterial(0,0),
  soilType(1), pult(0.0), y50(0.0), Cd(0.0)
{
	// Initialize variables
	this->revertToStart();
}

pyUCLA::~pyUCLA()
{

}

int 
pyUCLA::setTrialStrain (double strain, double strainRate)
{
	pult50		= 0.5 * pult;
	n			= 1.0/3.0;

	E			= pult50/ y50;
	limitStress = 1.0;
	//limitStress = 0.0001;
	
	epsilonY	= pow(pult50/(E * pow(y50, n)), 1.0/(1-n));
	theta		= (pult50 * n) / pow(y50,n);
	

	Ed			= 1.0e6;
	dragStress	= pult * Cd;

	Tstrain		= strain;
	//Tstress1	= 0.0;
	//Tstress2	= 0.0;

	//--- (1) py (no tension) ==============================================
	int signLimitStress1 = (limitStress < 0) ? -1 : 1;
	int plumin1 = -signLimitStress1;


	if (strain == 0.0) {
		Tstress1 = 0.0;
		Ttangent1 = E;
		TplasticStrain1 = CplasticStrain1;
		Thardening1 = Chardening1;	
	}
	else {

		Tstress1  = E * (strain - CplasticStrain1);
		Ttangent1 = E ;
		double f1 = Tstress1 * plumin1 - Chardening1;		
		Tstrain1  = strain * plumin1 - epsilonY;


		//--- hyperelastic below epsilonY
		if (Tstrain1 < 1.0e-16 && Chardening1 == 0.0) {
			Tstress1 = E * strain;
			Ttangent1 = E ; 
			TplasticStrain1 = 0.0;
			Thardening1 = 0.0;

		}
		//--- elastic step
		else if (f1 < 1.0e-16) {
			//Tstress1 = Tstress1;
			//Ttangent1 = E ; 
			TplasticStrain1 = CplasticStrain1;
			Thardening1 = Chardening1;
		}
		//--- plastic step
		else {
			int signStrain1 = (strain < 0) ? -1 : 1;
			Tstress1 = pult50 * signStrain1 * pow( (fabs(strain)/y50) , n) ;
			Ttangent1 = theta *  pow(fabs(strain),n-1); 
			TplasticStrain1 = strain - (Tstress1 / E);
			Thardening1 = fabs(Tstress1);
		}
	}


	//--- (2) py (no compression) =======================================
	int signLimitStress2 = (-limitStress < 0) ? -1 : 1;
	int plumin2 = -signLimitStress2;

	if (strain == 0.0) {
		Tstress2 = 0.0;
		Ttangent2 = E;
		TplasticStrain2 = CplasticStrain2;
		Thardening2 = Chardening2;	
	}
	else {

		Tstress2  = E * (strain - CplasticStrain2);
		Ttangent2 = E ;
		double f2 = Tstress2 * plumin2 - Chardening2;		
		Tstrain2  = strain * plumin2 - epsilonY;


		//--- hyperelastic below epsilonY
		if (Tstrain2 < 1.0e-16 && Chardening2 == 0.0) {
			Tstress2 = E * strain;
			Ttangent2 = E ; 
			TplasticStrain2 = 0.0;
			Thardening2 = 0.0;

		}
		//--- elastic step
		else if (f2 < 1.0e-16) {
			//Tstress2 = Tstress2;
			//Ttangent2 = E ; 
			TplasticStrain2 = CplasticStrain2;
			Thardening2 = Chardening2;
		}
		//--- plastic step
		else {
			int signStrain2 = (strain < 0) ? -1 : 1;
			Tstress2 = pult50 * signStrain2 * pow( (fabs(strain)/y50) , n) ;
			Ttangent2 = theta *  pow(fabs(strain),n-1); 
			TplasticStrain2 = strain - (Tstress2 / E);
			Thardening2 = fabs(Tstress2);
		}
	}



	//--- (3) drag element =================================================
    Tstress3 = Ed * (strain - CplasticStrain3);   // Elastic trial stress
    double f3 = fabs(Tstress3) - dragStress;		// Compute yield criterion

	if (f3 <= 1.0e-16)
    {
		TplasticStrain3 = CplasticStrain3;			// Update plastic strain
		//Tstress3 = Tstress3;
		Ttangent3 = Ed;								// Set trial tangent
    }
    else {		
		double dGamma = f3 / (Ed);					// Compute consistency parameter
		int signTstress3 = (Tstress3 < 0) ? -1 : 1;
		TplasticStrain3 = CplasticStrain3 + dGamma * signTstress3;// Update plastic strain
		Tstress3 = (1.0 - dGamma * Ed / fabs(Tstress3) ) * Tstress3;
		Ttangent3 = 0.0;

	}
 


	//--- projecting the obtained stresses and tangents ==================
/*	
    opserr << "BeforeTstress1: " << Tstress1 << "\n";	
	opserr << "BeforeTstress2: " << Tstress2 << "\n";	
	opserr << "BeforeTstress3: " << Tstress3 << "\n";	
	opserr << "BeforeTtangent1: " << Ttangent1 << "\n";	
	opserr << "BeforeTtangent2: " << Ttangent2 << "\n";
*/	
	
	projectStressTangent();
	
/*	
	opserr << "plumin1: " << plumin1 << "\n";	
	opserr << "plumin2: " << plumin2 << "\n";
	opserr << "Tstress1: " << Tstress1 << "\n";	
	opserr << "Tstress2: " << Tstress2 << "\n";	
	opserr << "Tstress3: " << Tstress3 << "\n";	
	opserr << "Ttangent1: " << Ttangent1 << "\n";	
	opserr << "Ttangent2: " << Ttangent2 << "\n";
*/

	//--- sum the resultance from each component ===========================
	Tstress = Tstress1+Tstress2+Tstress3;
	Ttangent = Ttangent1+Ttangent2+Ttangent3;
/*
	opserr << "Tstress: " << Tstress << "\n";	
	opserr << "Ttangent: " << Ttangent << "\n";	
	opserr << "==========================" << "\n";
*/
    return 0;
}



void 
pyUCLA::projectStressTangent()
{

	double beta = log(2.0)/(2.0*limitStress);
	Tstress1 = Tstress1 - 1.0/(2.0*beta)*log(0.5* (exp(2.0*beta*Tstress1)+1) );	
	Tstress2 = Tstress2 - 1.0/(2.0*beta)*log(0.5* (exp(2.0*beta*Tstress2)+1) );
	Ttangent1  = 1.0 / (pow(2,Tstress1/limitStress)+1) * Ttangent1;
	Ttangent2  = 1.0 / (pow(2,Tstress2/limitStress)+1) * Ttangent2;

	if (Tstress1 > 1.0e10)
	{
		Tstress1 = limitStress;
	}

	if (Tstress2 > 1.0e10)
	{
		Tstress2 = limitStress;
	}
}


double 
pyUCLA::getStress(void)
{
    return Tstress;
}

double 
pyUCLA::getTangent(void)
{
    return Ttangent;
}

double 
pyUCLA::getStrain(void)
{
    return Tstrain;
}

int 
pyUCLA::commitState(void)
{
    // Commit trial history variables
    CplasticStrain1 = TplasticStrain1;
    CplasticStrain2 = TplasticStrain2;
    CplasticStrain3 = TplasticStrain3;
    Chardening1		= Thardening1;
    Chardening2		= Thardening2;    
    return 0;
}

int 
pyUCLA::revertToLastCommit(void)
{
	return 0;
}

int 
pyUCLA::revertToStart(void)
{
    // Reset committed history variables
    CplasticStrain1 = 0.0;
	CplasticStrain2 = 0.0;
	CplasticStrain3 = 0.0;
    Chardening1		= 0.0;
    Chardening2		= 0.0;

    // Reset trial history variables
    TplasticStrain1 = 0.0;
    TplasticStrain2 = 0.0;
    TplasticStrain3 = 0.0;
    Thardening1		= 0.0;
    Thardening2		= 0.0;

	// Initialize state variables
	Tstrain			= 0.0;
	Tstress			= 0.0;
	Ttangent		= E;
	//Ttangent		= 0.0;

    Tstress1		= 0.0;	
	Tstress2		= 0.0;		
	Tstress3		= 0.0;	
	Tstrain1		= 0.0;	
	Tstrain2		= 0.0;		// Trial strain 2
	Tstrain3		= 0.0;		// Trial strain 3
	Ttangent1		= 0.0;		// Trial tangent 1
    Ttangent2		= 0.0;		// Trial tangent 2
    Ttangent3		= 0.0;		// Trial tangent 3

    return 0;
}

UniaxialMaterial *
pyUCLA::getCopy(void)
{
    //pyUCLA *theCopy =
	//new pyUCLA(this->getTag(), soilType, pult, y50, Cd);
	pyUCLA *theCopy;
	theCopy = new pyUCLA();
	*theCopy = *this;

    // Copy committed history variables
    theCopy->CplasticStrain1 = CplasticStrain1;
	theCopy->CplasticStrain2 = CplasticStrain2;
	theCopy->CplasticStrain3 = CplasticStrain3;
    theCopy->Chardening1 = Chardening1;
    theCopy->Chardening2 = Chardening2;

    // Copy trial history variables
    theCopy->TplasticStrain1 = TplasticStrain1;
    theCopy->TplasticStrain2 = TplasticStrain2;
    theCopy->TplasticStrain3 = TplasticStrain3;
    theCopy->Thardening1 = Thardening1;
    theCopy->Thardening2 = Thardening2;

    // Copy trial state variables
    theCopy->Tstrain = Tstrain;
    theCopy->Tstress = Tstress;
    theCopy->Ttangent = Ttangent;

/*	theCopy->Tstrain1 = Tstrain1;
    theCopy->Tstress1 = Tstress1;
    theCopy->Ttangent1 = Ttangent1;

	theCopy->Tstrain2 = Tstrain2;
    theCopy->Tstress2 = Tstress2;
    theCopy->Ttangent2 = Ttangent2;

	theCopy->Tstrain3 = Tstrain3;
    theCopy->Tstress3 = Tstress3;
    theCopy->Ttangent3 = Ttangent3;
*/

    return theCopy;
}

int 
pyUCLA::sendSelf(int cTag, Channel &theChannel)
{
  int res = 0;
  
  static Vector data(13);
  
  data(0)  = this->getTag();
  data(1)  = soilType;
  data(2)  = pult;
  data(3)  = y50;
  data(4)  = Cd;
  data(5)  = CplasticStrain1;
  data(6)  = CplasticStrain2;
  data(7)  = CplasticStrain3;
  data(8)  = Chardening1;
  data(9)  = Chardening2;
  data(10) = Tstrain;
  data(11) = Tstress;
  data(12) = Ttangent;

  /*data(13) = Tstrain1;
  data(14) = Tstrain2;
  data(15) = Tstrain3;
  data(16) = Tstress1;
  data(17) = Tstress2;
  data(18) = Tstress3;
  data(19) = Ttangent1;
  data(20) = Ttangent2;  
  data(21) = Ttangent3;
  */

  
  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "pyUCLA::sendSelf() - failed to send data\n";

  return res;
}

int 
pyUCLA::recvSelf(int cTag, Channel &theChannel, 
			       FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  static Vector data(13);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  
  if (res < 0) {
      opserr << "pyUCLA::recvSelf() - failed to receive data\n";
      E = 0; 
      this->setTag(0);      
  }
  else {
    this->setTag((int)data(0));
    soilType		= data(1);
    pult			= data(2);
    y50				= data(3);
    Cd				= data(4);
    CplasticStrain1 = data(5);
    CplasticStrain2 = data(6);
    CplasticStrain3 = data(7);
    Chardening1		= data(8);
    Chardening2		= data(9);
    Tstrain			= data(10);
    Tstress			= data(11);
    Ttangent		= data(12);

	/*Tstrain1		= data(13);
	Tstrain2		= data(14);
	Tstrain3		= data(15);
	Tstress1		= data(16);
	Tstress2		= data(17);
	Tstress3		= data(18);
	Ttangent1		= data(19);
	Ttangent2		= data(20);  
	Ttangent3		= data(21);
	*/
  }
    
  return res;
}

void 
pyUCLA::Print(OPS_Stream &s, int flag)
{
    s << "pyUCLA, tag: " << this->getTag() << endln;
    s << "  SoilType: " << soilType << endln;
    s << "  Pult: " << pult << endln;
    s << "  Y50: " << y50 << endln;
    s << "  Cd: " << Cd << endln;
}
