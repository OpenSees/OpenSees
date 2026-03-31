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
// $Date: 2010-09-16 00:04:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/OriginCentered.cpp,v $

//
// Description: This file contains the class definition for 
// OriginCentered. 

#include <math.h>

#include <stdlib.h>
#include <OriginCentered.h>
#include <OPS_Globals.h>
#include <float.h>
#include <Channel.h>

#include <elementAPI.h>
#include <OPS_Globals.h>


void *
OPS_OriginCentered(void)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;

  int    iData[1];
  double dData[12];
  int numData = 1;

  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial OriginCentered tag" << endln;
    return 0;
  }

  numData = OPS_GetNumRemainingInputArgs();

  if (numData != 6) {
    opserr << "Invalid #args, want: uniaxialMaterial OriginCentered " << iData[0] << 
      " f1? e1? f2? e2? f3? e3?>>" << endln;
    return 0;
  }

    if (OPS_GetDoubleInput(&numData, dData) != 0) {
      opserr << "Invalid arggs: uniaxialMaterial OriginCentered " << iData[0] << 
	" fy? E? b? <R0? cR1? cR2? <a1? a2? a3? a4?>>" << endln;
      return 0;
    }

    // Parsing was successful, allocate the material
    theMaterial = new OriginCentered(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5]);    

  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial of type OriginCentered Material\n";
    return 0;
  }

  return theMaterial;
}

OriginCentered::OriginCentered(int tag,
	    double _f1, double _e1, double _f2,
	    double _e2, double _f3, double _e3):
UniaxialMaterial(tag, MAT_TAG_OriginCentered),
	f1(_f1), e1(_e1), f2(_f2), e2(_e2), f3(_f3), e3(_e3)
{
	E1 = (f1)/(e1);
	E2 = (f2-f1)/(e2-e1);
	E3 = (f3-f2)/(e3-e2);

	this->revertToStart();
}

OriginCentered::OriginCentered(void):
UniaxialMaterial(0, MAT_TAG_OriginCentered)
{
	f1 = 0.0;
	e1 = 0.0;
	f2 = 0.0;
	e2 = 0.0; 
	f3 = 0.0;
	e3 = 0.0;

	E1 = 0.0;
	E2 = 0.0;
	E3 = 0.0;

	this->revertToStart();
}

OriginCentered::~OriginCentered(void)
{
	// Does nothing
}

UniaxialMaterial*
	OriginCentered::getCopy(void)
{
	OriginCentered *theCopy = new OriginCentered(this->getTag(), f1, e1, f2, e2, f3, e3);

	return theCopy;
}

double
	OriginCentered::getInitialTangent(void)
{
	return E1;
}

int
	OriginCentered::setTrialStrain(double trialStrain, double strainRate)
{
	Teps = trialStrain;
	double deps = (Teps - Ceps);

	if (deps < 0) {
		if (Teps > 0) {
			Ttan = Csig/Ceps;
			Tsig = Ttan*Teps;
		} else if (Teps > CepsMin) {
			Ttan = CsigMin/CepsMin;
			Tsig = Ttan*Teps;
		} else if (Teps > -e1)  {
			Ttan = E1;
			Tsig = E1*Teps;
		} else if (Teps > -e2)  {
			Ttan = E2;
			Tsig = E2*(Teps+e1)-f1;
		} else if (Teps > -e3)  {
			Ttan = E3;
			Tsig = E3*(Teps+e2)-f2;
		} else {
			Ttan = 0;
			Tsig  = -f3;
		}
	} else if (deps > 0) {
		if (Teps < 0) {
			Ttan = Csig/Ceps;
			Tsig = Ttan*Teps;
		} else if (Teps < CepsMax) {
			Ttan = CsigMax/CepsMax;
			Tsig = Ttan*Teps;
		} else if (Teps < e1)  {
			Ttan = E1;
			Tsig = E1*Teps;
		} else if (Teps < e2)  {
			Ttan = E2;
			Tsig = E2*(Teps-e1)+f1;
		} else if (Teps < e3)  {
			Ttan = E3;
			Tsig = E3*(Teps-e2)+f2;
		} else {
			Ttan = 0;
			Tsig  = f3;
		}

	} else {
		Ttan = Ctan;
		Tsig = Csig;
	}

	if (Teps > TepsMax) {
		TepsMax = Teps;
		TsigMax = Tsig;
	} else if (Teps < TepsMin) {
		TepsMin = Teps;
		TsigMin = Tsig;
	}

	return 0;
}



double 
	OriginCentered::getStrain(void)
{
	return Teps;
}

double 
	OriginCentered::getStress(void)
{
	return Tsig;
}

double 
	OriginCentered::getTangent(void)
{
	return Ttan;
}

int 
	OriginCentered::commitState(void)
{
    CepsMax = TepsMax;
    CepsMin = TepsMin;
    CsigMax = TsigMax;
    CsigMin = TsigMin;
    
	Csig = Tsig;
	Ceps = Teps;
	Ctan = Ttan;

	return 0;
}

int 
	OriginCentered::revertToLastCommit(void)
{
    TepsMax = CepsMax;
    TepsMin = CepsMin;
    TsigMax = CsigMax;
    TsigMin = CsigMin;
    
	Tsig = Csig;
	Teps = Ceps;
	Ttan = Ctan;

	return 0;
}

int 
	OriginCentered::revertToStart(void)
{	
    CepsMax = 0.0;
    CepsMin = 0.0;
    CsigMax = 0.0;
    CsigMin = 0.0;
    
	Csig = 0.0;
	Ceps = 0.0;
	Ctan = E1;

	this->revertToLastCommit();
	return 0;
}

int 
	OriginCentered::sendSelf(int commitTag, Channel &theChannel)
{
	static Vector data(14);
	data(0) = f1;
	data(1) = e1;
	data(2) = f2;
	data(3) = e2;
	data(4) = f3;
	data(5) = e3;
	data(6) = CepsMax;
	data(7) = CepsMin;
	data(8) = CsigMax;
	data(9) = CsigMin;
	data(10) = Csig;
	data(11) = Ceps;
	data(12) = Ctan;
	data(13) = this->getTag();

	if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "OriginCentered::sendSelf() - failed to sendSelf\n";
		return -1;
	}
	return 0;
}

int 
	OriginCentered::recvSelf(int commitTag, Channel &theChannel, 
	FEM_ObjectBroker &theBroker)
{
	static Vector data(14);

	if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
		opserr << "OriginCentered::recvSelf() - failed to recvSelf\n";
		return -1;
	}

	f1 = data(0);
	e1 = data(1);
	f2 = data(2);
	e2 = data(3);
	f3 = data(4);
	e3 = data(5);

	E1 = (f1)/(e1);
	E2 = (f2-f1)/(e2-e1);
	E3 = (f3-f2)/(e3-e2);

	CepsMax = data(6);
	CepsMin = data(7);
	CsigMax = data(8);
	CsigMin = data(9);
	Csig = data(10);
	Ceps = data(11);
	Ctan = data(12);
	this->setTag(int(data(13)));

	this->revertToLastCommit();

	return 0;
}

void 
	OriginCentered::Print(OPS_Stream &s, int flag)
{
	s << "OriginCentered:(strain, stress, tangent) " << Ceps << " " << Csig << " " << Ctan << endln;
}
