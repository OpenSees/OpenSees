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
                                                                        
// $Revision: 1.0 $
// $Date: 2012/05/01 01:00:00 $
// $Source: /usr/local/cvs/OpenSees/PACKAGES/NewElement/cpp/RotationShearCurve.cpp,v $
                                                                        
// Written: MRL
//
// Description: This file contains the implementation for the rotation based shear curve.
//
// What: "RotationShearCurve.cpp Rev A"
#include <stdlib.h>

#include "RotationShearCurve.h"
#include <G3Globals.h>
#include <math.h>
#include <ElementResponse.h>
#include <Element.h>
#include <Node.h>
#include <Domain.h>
#include <Vector.h>
#include <Parameter.h>
#include <Channel.h>
#include <float.h>
#include <DummyStream.h>
#include <elementAPI.h>


static int shearCurveCount = 0;

#ifndef max
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b) (((a)<(b)) ? (a) : (b))
#endif

void *
OPS_RotationShearCurve(void)
{
  if (shearCurveCount == 0) {
    opserr << "RotationShearCurve limit curve - Written by MRL UT Austin Copyright 2012 -  Use at your Own Peril \n";
    shearCurveCount++;
  }

  int argc = OPS_GetNumRemainingInputArgs();

  if (!(argc == 9 || argc == 23)) { 
    opserr << "WARNING RotationShearCurve -- insufficient arguments\n";
    opserr << "For direct input of shear curve parameters and degrading slope want:\n\n";
    opserr << "limitCurve RotationShearCurve crvTag? eleTag? \n";
    opserr << "ndI? ndJ? rotAxis? Vn? Vr? Kdeg? rotLim? \n" << endln;
    
    opserr << "OR for calibrated shear curve and degrading slope want:\n\n";
    opserr << "limitCurve RotationShearCurve crvTag? eleTag?\n";
    opserr << "ndI? ndJ? rotAxis? Vn? Vr? Kdeg? defType?\n";
    opserr << "b? d? h? L? st? As? Acc? ld? db? rhot? f'c?\n";
    opserr << "fy? fyt? delta?\n" << endln;
    
    return 0;
  }

  int iTagData[2];	  
  int iNodeData[3];	  
  double dKdegData[3];  
  double dRotLimData[1];
  int iTypeData[1];	  
  double dPropData[14]; 
  
  int numData;
  numData = 2;
  if (OPS_GetIntInput(&numData, iTagData) != 0) {
    opserr << "WARNING RotationShearCurve -- invalid crvTag? eleTag?\n" << endln;
    return 0;
  }
  int eleTag = iTagData[1];
  Domain *theDomain = 0;
  theDomain = OPS_GetDomain();
  if (theDomain == 0) {
    opserr << "WARNING RotationShearCurve -- Pointer to Domain was not returned\n" << endln;
    return 0;
  }
  Element *theElement = 0;
  theElement = theDomain->getElement(eleTag);
  if (theElement == 0) {
    opserr << "WARNING RotationShearCurve -- Element with tag " << iTagData[1] << " does not exist for shear curve tag " << iTagData[0] << endln << endln;
    return 0;
  }
  numData = 3;
  if (OPS_GetIntInput(&numData, iNodeData) != 0) {
    opserr << "WARNING RotationShearCurve -- invalid ndI? ndJ? rotAxis?\n" << endln;
    return 0;
  }
  Node *theNodeI = 0;
  theNodeI = theDomain->getNode(iNodeData[0]);
  if (theNodeI == 0) {
    opserr << "WARNING RotationShearCurve -- Node with tag " << iNodeData[0] << " does not exist for shear curve tag " << iTagData[0] << endln << endln;
    return 0;
  }	
  Node *theNodeJ = 0;
  theNodeJ = theDomain->getNode(iNodeData[1]);
  if (theNodeJ == 0) {
    opserr << "WARNING RotationShearCurve -- Node with tag " << iNodeData[1] << " does not exist for shear curve tag " << iTagData[0] << endln << endln;
    return 0;
  }	
  if (iNodeData[2] < 3 || iNodeData[2] > 6) {
    opserr << "WARNING RotationShearCurve -- rotAxis is invalid\n";
    opserr << "rotAxis = 3 -- Rotation about z-axis - 2D\n";
    opserr << "rotAxis = 4 -- Rotation about x-axis - 3D\n";
    opserr << "rotAxis = 5 -- Rotation about y-axis - 3D\n";
    opserr << "rotAxis = 6 -- Rotation about z-axis - 3D\n" << endln;
    return 0;
  }
  if (argc == 9) {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dKdegData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid Vn? Vr? Kdeg?\n" << endln;
      return 0;
    }
    if (dKdegData[0] != -1 && !(dKdegData[0] > 0)) {
      opserr << "WARNING RotationShearCurve --  Vn input is invalid\n";
      opserr << "Vn = -1 -- Shear critical limit is not used\n";
      opserr << "Vn > 0 -- Shear critical limit is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[1] < -1) {
      opserr << "WARNING RotationShearCurve -- Vr input is invalid\n";
      opserr << "Vr = -1 -- Residual shear strength = 0.2*(maximum shear at failure)\n";
      opserr << "-1 < Vr < 0 -- Residual shear strength = Vr*(maximum shear at failure)\n";
      opserr << "Vr >= 0 -- Residual shear strength is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[2] >= 0) {
      opserr << "WARNING RotationShearCurve -- Kdeg input is invalid\n";
      opserr << "The degrading slope must be less than zero\n" << endln;
      return 0;
    }
    numData = 1;
    if (OPS_GetDoubleInput(&numData, dRotLimData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid rotLim?\n" << endln;
      return 0;
    }
    if (dRotLimData[0] <= 0) {
      opserr << "WARNING RotationShearCurve -- rotLim input must be greater than zero\n" << endln;
      return 0;
    }
    
    LimitCurve *theCurve = 0;
    theCurve = new RotationShearCurve(iTagData[0], iTagData[1],
				      iNodeData[0], iNodeData[1], iNodeData[2], dKdegData[0], dKdegData[1], dKdegData[2], dRotLimData[0], 0,
				      0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
				      0.0, 0.0, 0.0, theDomain, theElement, theNodeI, theNodeJ);
    if (theCurve == 0) {
      opserr << "WARNING RotationShearCurve -- could not create limitCurve with constructor " << iTagData[0] << endln << endln;
      return 0;
    }
    return theCurve;
  } 
  else {
    numData = 3;
    if (OPS_GetDoubleInput(&numData, dKdegData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid Vn? Vr? Kdeg?\n" << endln;
      return 0;
    }
    if (dKdegData[0] != -1 && !(dKdegData[0] >= 0)) {
      opserr << "WARNING RotationShearCurve --  Vn input is invalid\n";
      opserr << "Vn = -1 -- Shear critical limit is not used\n";
      opserr << "Vn = 0 -- Shear critical limit is calculated using ASCE 41 Eq. 6-4\n";
      opserr << "Vn > 0 -- Shear critical limit is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[1] < -1) {
      opserr << "WARNING RotationShearCurve -- Vr input is invalid\n";
      opserr << "Vr = -1 -- Residual shear strength from regression\n";
      opserr << "-1 < Vr < 0 -- Residual shear strength = Vr*(maximum shear at failure)\n";
      opserr << "Vr >= 0 -- Residual shear strength is the input value\n" << endln;
      return 0;
    }
    if (dKdegData[2] > 0) {
      opserr << "WARNING RotationShearCurve -- Kdeg input is invalid\n";
      opserr << "Kdeg = 0 -- Degrading slope calculated by regressions\n";
      opserr << "Kdeg < 0 -- Degrading slope is the input value\n" << endln;
      return 0;
    }
    numData = 1;
    if (OPS_GetIntInput(&numData, iTypeData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid defType?\n" << endln;
      return 0;
    }
    if (iTypeData[0] > 5 || iTypeData[0] <= 0) {
      opserr << "WARNING RotationShearCurve -- invalid defType input?\n" << endln;
      opserr << "1 -- Flexure-Shear capacity based on theta_f rotation capacity\n";
      opserr << "2 -- Flexure-Shear capacity based on theta_total rotation capacity\n";
      opserr << "3 -- Flexure-Shear capacity based on theta_flexural rotation capacity\n";
      opserr << "4 -- Flexure-Shear capacity based on theta_total-plastic rotation capacity\n";
      opserr << "5 -- Flexure-Shear capacity based on theta_flexural-plastic rotation capacity\n" << endln;
      return 0;
    }
    numData = 14;
    if (OPS_GetDoubleInput(&numData, dPropData) != 0) {
      opserr << "WARNING RotationShearCurve -- invalid b? d? h? L? st? As? Acc? ld? db? rhot? f'c? fy? fyt? delta?\n" << endln;
      return 0;
    }
    LimitCurve *theCurve = 0;
    theCurve = new RotationShearCurve(iTagData[0], iTagData[1],
				      iNodeData[0], iNodeData[1], iNodeData[2], 
				      dKdegData[0], dKdegData[1], dKdegData[2], 0.0, iTypeData[0],
				      fabs(dPropData[0]), fabs(dPropData[1]), fabs(dPropData[2]), fabs(dPropData[3]), 
				      fabs(dPropData[4]), fabs(dPropData[5]), fabs(dPropData[6]), fabs(dPropData[7]), fabs(dPropData[8]), 
				      fabs(dPropData[9]), fabs(dPropData[10]), fabs(dPropData[11]), fabs(dPropData[12]), 
				      dPropData[13], theDomain, theElement, theNodeI, theNodeJ);
    if (theCurve == 0) {
      opserr << "WARNING RotationShearCurve -- could not create limitCurve with constructor " << iTagData[0] << endln << endln;
      return 0;
    }
    return theCurve;
  }
}

RotationShearCurve::RotationShearCurve(int crvTag, int eTag,
	int ni, int nj, int axis, double Vnom, double Vres, double Kd, double rLim, int dfTyp,
	double B, double D, double H, double Lgth, double S, double AreaS, double AreaCC, double Ldev, double DiaBar, double rhoTrans, double Fc,
	double Fy, double Fytrans, double Delta, Domain *theDom, Element *theEle, Node *theNdI, Node *theNdJ)
	:LimitCurve(crvTag, LIMCRV_TAG_RotShearCurve), eleTag(eTag),
	ndI(ni), ndJ(nj), rotAxis(axis), Vn(Vnom), Vr(Vres), Kdeg(Kd), rotLim(rLim), defType(dfTyp),
	b(B), d(D), h(H), L(Lgth), st(S), As(AreaS), Acc(AreaCC), ld(Ldev), db(DiaBar), rhot(rhoTrans), fc(Fc),
	fy(Fy), fyt(Fytrans), delta(Delta), theDomain(theDom), theElement(theEle), theNodeI(theNdI), theNodeJ(theNdJ), curveTag(crvTag)
{
	if(revertToStart() != 0)
	{
		opserr << "FATAL ERROR RotationShearCurve -- could not initialize variables\n" << endln;
		exit(-1);
	}
}

RotationShearCurve::RotationShearCurve()
	:LimitCurve(0, LIMCRV_TAG_RotShearCurve), eleTag(0),
	ndI(0), ndJ(0), rotAxis(0), Vn(0.0), Vr(0.0), Kdeg(0.0), rotLim(0.0), defType(0),
	b(0.0), d(0.0), h(0.0), L(0.0), st(0.0), As(0.0), Acc(0.0), ld(0.0), db(0.0), rhot(0.0), fc(0.0),
	fy(0.0), fyt(0.0), delta(0.0), theDomain(0), theElement(0), theNodeI(0), theNodeJ(0)
{
	if(revertToStart() != 0)
	{
		opserr << "FATAL ERROR RotationShearCurve -- could not initialize variables\n" << endln;
		exit(-1);
	}
}

RotationShearCurve::~RotationShearCurve()
{
	//Nothing to Delete
}

LimitCurve*
RotationShearCurve::getCopy(void)
{
  RotationShearCurve *theCopy = new RotationShearCurve(this->getTag(), eleTag,
						       ndI, ndJ, rotAxis, Vn, Vr, Kdeg, rotLim, defType,
						       b, d, h, L, st, As, Acc, ld, db, rhot, fc,
						       fy, fyt, delta, theDomain, theElement, theNodeI, theNodeJ);
  
  theCopy->forceVec = forceVec;
  theCopy->thetaMin = thetaMin;
  theCopy->P = P;
  theCopy->M = M;
  theCopy->stateFlag = stateFlag;
  
  return theCopy;
}

int
RotationShearCurve::checkElementState(double springForce)
{	
	double shearForce = fabs(springForce);	
	getElemForces();	
	const Vector &dispI = theNodeI->getTrialDisp();
	const Vector &dispJ = theNodeJ->getTrialDisp();
	double rotDef = fabs(dispJ(rotAxis-1) - dispI(rotAxis-1));

	if (stateFlag == 0)
	{
		if (Vn == 0.0) {
			double ShearForceLimit = findCritLimit(shearForce, M);
			if(shearForce >= ShearForceLimit) {
				stateFlag = 1;
				setDegSlope(shearForce);
			}
		}
		else if (Vn > 0.0) {
			if (shearForce >= Vn) {
				stateFlag = 1;
				setDegSlope(shearForce);
			}
		}
		if (defType == 0) {
			if (rotDef >= rotLim) {
				stateFlag = 1;
				setDegSlope(shearForce);
			}
		} else {
			double shearRotLimit = findLimit(shearForce);
			if ((rotDef >= shearRotLimit) && (rotDef >= thetaMin)) {
				stateFlag = 1;
				setDegSlope(shearForce);
			}
		}
	}
	else 
	{
		stateFlag = 2;
	}
	return stateFlag;
}

double
RotationShearCurve::getDegSlope(void)
{
	return Kdeg;
}


double
RotationShearCurve::getResForce(void)
{
	return Vr;  
}

double
RotationShearCurve::getUnbalanceForce(void)
{
	//Do nothing for this class
	return 0.0;
}

int
RotationShearCurve::sendSelf(int commitTag, Channel &theChannel)
{
	return -1;
}

int
RotationShearCurve::recvSelf(int commitTag, Channel &theChannel, 
			FEM_ObjectBroker &theBroker)
{
	return -1;
}
    
void
RotationShearCurve::Print(OPS_Stream &s, int flag)
{
	s << "Shear Limit Curve, tag: " << this->getTag() << endln;
	s << "eleTag: " << eleTag << endln;
	s << "thetaMin: " << thetaMin << endln;
	s << "P: " << P << endln;
	s << "M: " << M << endln;
	s << "stateFlag: " << stateFlag << endln;
	s << "ndI: " << ndI << endln;
	s << "ndJ: " << ndJ << endln;
	s << "rotAxis: " << rotAxis << endln;
	s << "Vn: " << Vn << endln;
	s << "Vr: " << Vr << endln;
	s << "Kdeg: " << Kdeg << endln;
	s << "rotLim: " << rotLim << endln;
	s << "defType: " << defType << endln;
	s << "b: " << b << endln;
	s << "d: " << d << endln;
	s << "h: " << h << endln;
	s << "L: " << L << endln;
	s << "st: " << st << endln;
	s << "As: " << As << endln;
	s << "Acc: " << Acc << endln;
	s << "ld: " << ld << endln;
	s << "db: " << db << endln;
	s << "rhot: " << rhot << endln;	
	s << "fc: " << fc << endln;
	s << "fy: " << fy << endln;
	s << "fyt: " << fyt << endln;
	s << "delta: " << delta << endln;
	s << endln;
}

void
RotationShearCurve::setDegSlope(double V)
{	
  if (Vr == -1) {
    if ((st == 0.0) && (d == 0.0)) {
      Vr = V*0.2;
    }
    else {
      double ResShearRatio = 0.362283 + -0.170283*(st/d);
      ResShearRatio = max(ResShearRatio,0.0);
      Vr = V*ResShearRatio;
    }
  }
  else if(Vr > -1 && Vr < 0.0)
    Vr = fabs(V*Vr);
  else if(Vr >= 0.0)
    Vr = Vr;
  else
    {
      opserr << "FATAL ERROR RotationShearCurve -- Vr input is not implemented\n" << endln;
      exit(-1);
    }
  
  if (Kdeg == 0.0) {
    double Ag = b*h;
    double v = V/(b*d);
    double ResDriftRatio = -0.158370 - 15.437656*rhot - 0.009391*(ld/db) + 0.697682*(Acc/Ag) + 0.582667*(fy*As/(fc*Ag));
    ResDriftRatio = max(ResDriftRatio, 0.02);
    Kdeg = -V/(ResDriftRatio*L);
  } 
  else if (Kdeg < 0.0) {
    Kdeg = Kdeg;
  }
  else {
    opserr << "FATAL ERROR RotationShearCurve -- Regression Kdeg input is not implemented\n" << endln;
    exit(-1);
  }
}

double
RotationShearCurve::findLimit(double V) 
{
  double theta = 0.0;
  double Ag=b*h;
  double v=V/(b*d);
  
  if (defType == 0) {
    theta = rotLim;
    thetaMin = 0.0;
  }
  else if(defType == 1) {
    theta = 0.026515 - 0.033432*(P/(Ag*fc)) - 0.009963*(st/d)+delta;
    thetaMin = 0.006;
  }
  else if(defType == 2) {
    theta = 0.044-0.017*(st/d)-0.021*(P/(Ag*fc))-0.002*(v/(sqrt(fc*1000)/1000))+delta;
    thetaMin = 0.009;
  }
  else if(defType == 3) {
    theta = 0.45*(0.044-0.017*(st/d)-0.021*(P/(Ag*fc))-0.002*(v/(sqrt(fc*1000)/1000)))+delta;
    thetaMin = 0.00405;
  }	
  else if(defType == 4) {
    theta = 0.032-0.014*(st/d)-0.017*(P/(Ag*fc))-0.0016*(v/(sqrt(fc*1000)/1000))+delta;
    thetaMin = 0.0;
  }
  else if(defType == 5) {
    theta = 0.45*(0.032-0.014*(st/d)-0.017*(P/(Ag*fc))-0.0016*(v/(sqrt(fc*1000)/1000)))+delta;
    thetaMin = 0.0;
  }
  return theta;
}

double
RotationShearCurve::findCritLimit(double Vu, double Mu) 
{
	double rhow = As/(b*d);
	double Ag = b*h;
	double Mratio = Mu/(Vu*d);
	double Nu = P*1000;
	double Vc = 0.8*Ag*(6*sqrt(fc*1000)/Mratio*sqrt(1+Nu/(6*sqrt(fc*1000)*Ag)))/1000; //(kips)
	double Av = rhot*st*b;
	double Vs = Av*fyt*d/st;
	double V = Vc + Vs; 
	return V;
}

int RotationShearCurve::revertToStart(void)
{
	forceVec = 0;
	thetaMin = 0.0;
	P = 0.0;
	M = 0.0;
	stateFlag = 0;
	return 0;
}

void
RotationShearCurve::getElemForces(void)
{
	int trash;		
	const char *forceType2[1] = {"localForce"}; 
	Response *theForces = 0;	
	DummyStream dummy;	
	theForces =  theElement->setResponse(forceType2, 1, dummy); 
	trash = theForces->getResponse();		
	Information &theInfo = theForces->getInformation();
	Vector *forceVec = (theInfo.theVector);	
	if (forceVec == 0) {
		opserr << "FATAL ERROR RotationShearCurve -- unable to assign force vector\n" << endln;
		exit(-1);
	}
	P = fabs((*forceVec)(0));
	M = fabs((*forceVec)(2));
}
