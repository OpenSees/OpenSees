// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureDependMultiYield.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// PressureDependMultiYield.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include "PressureDependMultiYield.h"
#include <Information.h>

int PressureDependMultiYield::loadStage = 0;
double PressureDependMultiYield::AtmoPress = 0.;
Matrix PressureDependMultiYield::theTangent = Matrix(6,6);
T2Vector PressureDependMultiYield::trialStrain = Vector(6);
T2Vector PressureDependMultiYield::subStrainRate = Vector(6);
Vector PressureDependMultiYield::workV = Vector(3);
Matrix PressureDependMultiYield::workM = Matrix(3,3);

const Vector zeroVector(6);
const	double pi = 3.14159265358979;

PressureDependMultiYield::PressureDependMultiYield (int tag, int nd, 
															double refShearModul,
    			  		              double refBulkModul, double frictionAng,
								              double peakShearStra, double refPress, 
															double cohesi, 	double pressDependCoe,
															int numberOfYieldSurf, 
															double phaseTransformAng, 
                              double contractionParam1,
                              double contractionParam2,
                              double dilationParam1,
                              double dilationParam2,
															double volLimit,
                              double liquefactionParam1,
                              double liquefactionParam2,
                              double liquefactionParam3,
                              double liquefactionParam4,
															double atm)
 : NDMaterial(tag,MAT_TAG_PressureDependMultiYield), currentStress(zeroVector),
   trialStress(zeroVector), currentStrain(zeroVector), strainRate(zeroVector),
   reversalStress(zeroVector), PPZPivot(zeroVector), PPZCenter(zeroVector), 
	 lockStress(zeroVector), reversalStressCommitted(zeroVector), 
	 PPZPivotCommitted(zeroVector), PPZCenterCommitted(zeroVector),
	 lockStressCommitted(zeroVector)
{
	if (nd !=2 && nd !=3) {
		cerr << "FATAL:PressureDependMultiYield:: dimension error" << endl;
    cerr << "Dimension has to be 2 or 3, you give nd= " << nd << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refShearModul <= 0) {
		cerr << "FATAL:PressureDependMultiYield:: refShearModulus <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refBulkModul <= 0) {
		cerr << "FATAL:PressureDependMultiYield:: refBulkModulus <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
  if (frictionAng <= 0.) {
  	cerr << "FATAL:PressureDependMultiYield:: frictionAngle <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (frictionAng >= 90.) {
  	cerr << "FATAL:PressureDependMultiYield:: frictionAngle >= 90" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (phaseTransformAng <= 0.) {
  	cerr << "FATAL:PressureDependMultiYield:: phaseTransformAng <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (phaseTransformAng > frictionAng) {
  	cerr << "WARNING:PressureDependMultiYield:: phaseTransformAng > frictionAng" << endl;
		cerr << "Will set phaseTransformAng = frictionAng." <<endl;
    phaseTransformAng = frictionAng;
  }
	if (cohesi < 0) {
		cerr << "WARNING:PressureDependMultiYield:: cohesion < 0" << endl;
    cerr << "Will reset cohesion to zero." << endl;
	  cohesi = 0.;
  }
	if (peakShearStra <= 0) {
		cerr << "FATAL:PressureDependMultiYield:: peakShearStra <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refPress <= 0) {
		cerr << "FATAL:PressureDependMultiYield:: refPress <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (pressDependCoe < 0) {
		cerr << "WARNING:PressureDependMultiYield:: pressDependCoe < 0" << endl;
    cerr << "Will reset pressDependCoe to zero." << endl;
	  pressDependCoe = 0.;
  }
	if (numberOfYieldSurf <= 0) {
		cerr << "WARNING:PressureDependMultiYield:: numberOfSurfaces <= 0" << endl;
		cerr << "Will use 10 yield surfaces." << endl;
	  numberOfYieldSurf = 10;
  }
	if (volLimit < 0) {
		cerr << "WARNING:PressureDependMultiYield:: volLimit < 0" << endl;
		cerr << "Will reset volLimit to zero." << endl;
		volLimit = 0.;
	}
	if (atm <= 0) {
		cerr << "FATAL:PressureDependMultiYield:: atm <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }

	ndm = nd;
	loadStage = 0;   //default
  refShearModulus = refShearModul;
	refBulkModulus = refBulkModul;
	frictionAngle = frictionAng;
	peakShearStrain = peakShearStra;
	refPressure = -refPress;  //compression is negative
	cohesion = cohesi;
	pressDependCoeff = pressDependCoe;
	numOfSurfaces = numberOfYieldSurf;
	phaseTransfAngle = phaseTransformAng;
  contractParam1 = contractionParam1;
  contractParam2 = contractionParam2;
  dilateParam1 = dilationParam1;
  dilateParam2 = dilationParam2;
	volumeLimit = volLimit;
  liquefyParam1 = liquefactionParam1;
  liquefyParam2 = liquefactionParam2;
  liquefyParam3 = liquefactionParam3;
  liquefyParam4 = liquefactionParam4;
	AtmoPress = atm;

	e2p = committedActiveSurf = activeSurfaceNum = 0; 
  onPPZCommitted = onPPZ = -1 ; 
	PPZSizeCommitted = PPZSize = 0.;
	pressureDCommitted = pressureD = 0.;
	cumuDilateStrainOctaCommitted    = cumuDilateStrainOcta = 0.;
  maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta = 0.;
	cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta = 0.;
	prePPZStrainOctaCommitted = prePPZStrainOcta = 0.;
	oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta = 0.;

	theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 

	setUpSurfaces();  // residualPress and stressRatioPT are calculated inside.
}
   

PressureDependMultiYield::PressureDependMultiYield () 
 : NDMaterial(0,MAT_TAG_PressureDependMultiYield), 
   currentStress(zeroVector), trialStress(zeroVector), currentStrain(zeroVector), 
	 strainRate(zeroVector), reversalStress(zeroVector), PPZPivot(zeroVector), 
	 PPZCenter(zeroVector), lockStress(zeroVector), reversalStressCommitted(zeroVector), 
	 PPZPivotCommitted(zeroVector), PPZCenterCommitted(zeroVector),
	 lockStressCommitted(zeroVector)
{
	ndm = 3;
  refShearModulus = refBulkModulus = frictionAngle = 0.;
	peakShearStrain = refPressure = cohesion = pressDependCoeff = 0.;
	numOfSurfaces = 1;
	phaseTransfAngle = contractParam1 = contractParam2 = 0.;
  dilateParam1 = dilateParam2 = volumeLimit = liquefyParam1 = 0.;
  liquefyParam2 = liquefyParam3 = liquefyParam4 = 0.;

	theSurfaces = new MultiYieldSurface[1];
  committedSurfaces = new MultiYieldSurface[1];
  activeSurfaceNum = committedActiveSurf = 0; 
}


PressureDependMultiYield::PressureDependMultiYield (const PressureDependMultiYield & a)
 : NDMaterial(a.getTag(),MAT_TAG_PressureDependMultiYield), 
   currentStress(a.currentStress), trialStress(a.trialStress), 
	 currentStrain(a.currentStrain), strainRate(a.strainRate), 
   reversalStress(a.reversalStress), PPZPivot(a.PPZPivot), PPZCenter(a.PPZCenter), 
	 lockStress(a.lockStress), reversalStressCommitted(a.reversalStressCommitted), 
	 PPZPivotCommitted(a.PPZPivotCommitted), 
	 PPZCenterCommitted(a.PPZCenterCommitted),
	 lockStressCommitted(a.lockStressCommitted)
{
	ndm = a.ndm;
  refShearModulus = a.refShearModulus;
	refBulkModulus = a.refBulkModulus;
	frictionAngle = a.frictionAngle;
	peakShearStrain = a.peakShearStrain;
	refPressure = a.refPressure;
	cohesion = a.cohesion;
	pressDependCoeff = a.pressDependCoeff;
	numOfSurfaces = a.numOfSurfaces;
	phaseTransfAngle = a.phaseTransfAngle;
  contractParam1 = a.contractParam1;
  contractParam2 = a.contractParam2;
  dilateParam1 = a.dilateParam1;
  dilateParam2 = a.dilateParam2;
	volumeLimit = a.volumeLimit;
  liquefyParam1 = a.liquefyParam1;
  liquefyParam2 = a.liquefyParam2;
  liquefyParam3 = a.liquefyParam3;
  liquefyParam4 = a.liquefyParam4;

	e2p = a.e2p;
	residualPress = a.residualPress;
	stressRatioPT = a.stressRatioPT;
	strainPTOcta = a.strainPTOcta;
	modulusFactor = a.modulusFactor;
	activeSurfaceNum = a.activeSurfaceNum;
	committedActiveSurf = a.committedActiveSurf;
	pressureDCommitted     = a.pressureDCommitted;
  onPPZCommitted = a.onPPZCommitted; 
	PPZSizeCommitted      = a.PPZSizeCommitted;
	cumuDilateStrainOctaCommitted    = a.cumuDilateStrainOctaCommitted;
  maxCumuDilateStrainOctaCommitted = a.maxCumuDilateStrainOctaCommitted;
	cumuTranslateStrainOctaCommitted = a.cumuTranslateStrainOctaCommitted;
	prePPZStrainOctaCommitted        = a.prePPZStrainOctaCommitted;
	oppoPrePPZStrainOctaCommitted    = a.oppoPrePPZStrainOctaCommitted;
	pressureD     = a.pressureD;
  onPPZ = a.onPPZ; 
	PPZSize      = a.PPZSize;
	cumuDilateStrainOcta    = a.cumuDilateStrainOcta;
  maxCumuDilateStrainOcta = a.maxCumuDilateStrainOcta;
	cumuTranslateStrainOcta = a.cumuTranslateStrainOcta;
	prePPZStrainOcta        = a.prePPZStrainOcta;
	oppoPrePPZStrainOcta    = a.oppoPrePPZStrainOcta;

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];  
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];  
    theSurfaces[i] = a.theSurfaces[i];  
	}
}


PressureDependMultiYield::~PressureDependMultiYield ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
	if (committedSurfaces != 0) delete [] committedSurfaces;
}


void PressureDependMultiYield::elast2Plast(void)
{
	if (loadStage == 0 || e2p == 1) return;
	e2p = 1;

	if (currentStress.volume() > 0.) {
  	cerr << "WARNING:PressureDependMultiYield::elast2Plast(): material in tension." << endl;
    currentStress = T2Vector(currentStress.deviator(),0);
	}
  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
     if (committedActiveSurf == numOfSurfaces) {
        cerr <<"WARNING:PressureDependMultiYield::elast2Plast(): stress out of failure surface"<<endl;
				deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
        initSurfaceUpdate();
				//initStrainUpdate();
				return;
     }
  } 

  committedActiveSurf--;
	initSurfaceUpdate();
	//initStrainUpdate();
}


int PressureDependMultiYield::setTrialStrain (const Vector &strain)
{
	Vector temp(6);
	if (ndm==3 && strain.Size()==6) 
		temp = strain;
	else if (ndm==2 && strain.Size()==3) {
	  temp[0] = strain[0];
	  temp[1] = strain[1];
	  temp[3] = strain[2];
  }
	else {
		cerr << "Fatal:PressureDependMultiYield:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  strainRate = T2Vector(temp-currentStrain.t2Vector());
	return 0;
}


int PressureDependMultiYield::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureDependMultiYield::setTrialStrainIncr (const Vector &strain)
{
	Vector temp(6);
	if (ndm==3 && strain.Size()==6) 
		temp = strain;
	else if (ndm==2 && strain.Size()==3) {
	  temp[0] = strain[0];
	  temp[1] = strain[1];
	  temp[3] = strain[2];
  }
	else {
		cerr << "Fatal:PressureDependMultiYield:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  strainRate = T2Vector(temp);

	return 0;
}


int PressureDependMultiYield::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureDependMultiYield::getTangent (void)
{
	if (loadStage != 0 && e2p == 0) elast2Plast();

	if (loadStage==0) {  //linear elastic
  	for (int i=0;i<6;i++) 
	  	for (int j=0;j<6;j++) {
		  	theTangent(i,j) = 0.;
        if (i==j) theTangent(i,j) += 2.*refShearModulus;
			  if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
		}
	}
	else {
	  double coeff1, coeff2;
  	Vector devia(6);
  	double factor = getModulusFactor(currentStress);
  	double shearModulus = factor*refShearModulus;
  	double bulkModulus = factor*refBulkModulus;		
	
    if (loadStage!=0 && committedActiveSurf > 0) {
	  	T2Vector Q = getSurfaceNormal(currentStress);
	    devia = Q.deviator();
	    double volume = Q.volume();
	  	double Ho = 9.*bulkModulus*volume*volume + 2.*shearModulus*(devia && devia);
	    Vector devia = currentStress.deviator()-committedSurfaces[committedActiveSurf].center();
	    double plastModul = committedSurfaces[committedActiveSurf].modulus();
	    coeff1 = 9.*bulkModulus*bulkModulus*volume*volume/(Ho+plastModul);
	  	coeff2 = 4.*shearModulus*shearModulus/(Ho+plastModul);
		}
	  else coeff1 = coeff2 = 0.;

	  for (int i=0;i<6;i++) 
	  	for (int j=0;j<6;j++) {
		  	theTangent(i,j) = - coeff2*devia[i]*devia[j];
        if (i==j) theTangent(i,j) += 2.*shearModulus;
		  	if (i<3 && j<3) theTangent(i,j) += (bulkModulus - 2.*shearModulus/3. - coeff1);
			}
  }

	if (ndm==3) 
		return theTangent;
	else {
	  workM(0,0) = theTangent(0,0);
	  workM(0,1) = theTangent(0,1);
	  workM(0,2) = 0.;//theTangent(0,3);
	  workM(1,0) = theTangent(1,0);
	  workM(1,1) = theTangent(1,1);
	  workM(1,2) = 0.;//theTangent(1,3);
	  workM(2,0) = 0.;//theTangent(3,0);
	  workM(2,1) = 0.;//theTangent(3,1);
	  workM(2,2) = theTangent(3,3);
	
  	return workM;
	}
}


const Vector & PressureDependMultiYield::getStress (void)
{
  int i;
  if (loadStage != 0 && e2p == 0) 
    elast2Plast();

  if (loadStage==0) {  //linear elastic
    trialStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
    getTangent();
    Vector a = theTangent * trialStrain.t2Vector();
    trialStress = T2Vector(a);
  }
  else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
    pressureD = pressureDCommitted;
    reversalStress = reversalStressCommitted;
    onPPZ = onPPZCommitted;
    PPZSize = PPZSizeCommitted;
    cumuDilateStrainOcta = cumuDilateStrainOctaCommitted;
    maxCumuDilateStrainOcta = maxCumuDilateStrainOctaCommitted;
    cumuTranslateStrainOcta = cumuTranslateStrainOctaCommitted;
    prePPZStrainOcta = prePPZStrainOctaCommitted;
    oppoPrePPZStrainOcta = oppoPrePPZStrainOctaCommitted;
    PPZPivot = PPZPivotCommitted;
    PPZCenter = PPZCenterCommitted;
    lockStress = lockStressCommitted;

    subStrainRate = strainRate;
    setTrialStress(currentStress);
    if (isLoadReversal()) {
      updateInnerSurface();
      activeSurfaceNum = 0;
    }

    int numSubIncre = setSubStrainRate();

    for (i=0; i<numSubIncre; i++) {
      if (i==0)  
	setTrialStress(currentStress);
      else setTrialStress(trialStress); 
      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;
      trialStrain = T2Vector(currentStrain.t2Vector() 
			     + subStrainRate.t2Vector()*(i+1));
      int lock = stressCorrection(0);
      if(lock==0) updateActiveSurface();
    }
  }

	if (ndm==3)
    return trialStress.t2Vector();
	else {
    workV[0] = trialStress.t2Vector()[0];
    workV[1] = trialStress.t2Vector()[1];
    workV[2] = trialStress.t2Vector()[3];
    return workV;
  }
}


const Vector & PressureDependMultiYield::getStrain (void)
{
  return getCommittedStrain();
}


int PressureDependMultiYield::commitState (void)
{
	currentStress = trialStress;
	currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
	if (loadStage) {
    committedActiveSurf = activeSurfaceNum;
	  for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
		pressureDCommitted = pressureD;
    reversalStressCommitted = reversalStress;
	  onPPZCommitted = onPPZ;
	  PPZSizeCommitted = PPZSize;
	  cumuDilateStrainOctaCommitted = cumuDilateStrainOcta;
    maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta;
	  cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta;
	  prePPZStrainOctaCommitted = prePPZStrainOcta;
	  oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta;
	  PPZPivotCommitted = PPZPivot;
	  PPZCenterCommitted = PPZCenter;
	  lockStressCommitted = lockStress;
	}

	return 0;
}


int PressureDependMultiYield::revertToLastCommit (void)
{
	return 0;
}


NDMaterial * PressureDependMultiYield::getCopy (void)
{
  PressureDependMultiYield * copy = new PressureDependMultiYield(*this);
	return copy;
}


NDMaterial * PressureDependMultiYield::getCopy (const char *code)
{
	if (strcmp(code,"PressureDependMultiYield") == 0) {
     PressureDependMultiYield * copy = new PressureDependMultiYield(*this);
	   return copy;
	}

	return 0;
}


const char * PressureDependMultiYield::getType (void) const
{
  return "PressureDependMultiYield";
}


int PressureDependMultiYield::getOrder (void) const
{
	return (ndm == 2) ? 3 : 6;
}


int PressureDependMultiYield::updateParameter(int responseID, Information &info)
{
	loadStage = responseID;
	return 0;
}


int PressureDependMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
	// Need to implement
	return 0;
}


int PressureDependMultiYield::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	// Need to implement
	return 0;
}


int PressureDependMultiYield::getResponse (int responseID, Information &matInfo)
{
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = getCommittedStrain();
			return 0;
		case 3:
			if (matInfo.theMatrix != 0)
				*(matInfo.theMatrix) = getTangent();
			return 0;
		default:
			return -1;
	}
}


void PressureDependMultiYield::Print(ostream &s, int flag )
{
	s << "PressureDependMultiYield" << endl;
}


const Vector & PressureDependMultiYield::getCommittedStress (void)
{
	if (ndm==3)
    return currentStress.t2Vector();
	else {
    workV[0] = currentStress.t2Vector()[0];
    workV[1] = currentStress.t2Vector()[1];
    workV[2] = currentStress.t2Vector()[3];
    return workV;
  }
}


const Vector & PressureDependMultiYield::getCommittedStrain (void)
{
	if (ndm==3)
    return currentStrain.t2Vector();
	else {
    workV[0] = currentStrain.t2Vector()[0];
    workV[1] = currentStrain.t2Vector()[1];
    workV[2] = currentStrain.t2Vector()[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used 
void PressureDependMultiYield::setUpSurfaces (void)
{ 
  double refStrain, peakShear, coneHeight;

	double sinPhi = sin(frictionAngle * pi/180.);
	double Mnys = 6.*sinPhi/(3.-sinPhi);
  double sinPhiPT = sin(phaseTransfAngle * pi/180.);
	stressRatioPT = 6.*sinPhiPT/(3.-sinPhiPT);
	residualPress = 3.* cohesion / (sqrt(2.) * Mnys);

	// a small nonzero residualPress for numerical purpose only
	if (residualPress < 1.) residualPress = 1.; 
	coneHeight = - (refPressure - residualPress);
  peakShear = sqrt(2.) * coneHeight * Mnys / 3.; 
  refStrain = (peakShearStrain * peakShear) 
			        / (refShearModulus * peakShearStrain - peakShear);

  double stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double ratio1, ratio2;
	double stressInc = peakShear / numOfSurfaces;

  for (int ii=1; ii<=numOfSurfaces; ii++){
        stress1 = ii * stressInc; 
				stress2 = stress1 + stressInc;
				ratio1 = 3. * stress1 / sqrt(2.) / coneHeight;
        ratio2 = 3. * stress2 / sqrt(2.) / coneHeight;
        strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
        strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);

				if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
					double ratio = (ratio2 - stressRatioPT)/(ratio2 - ratio1);
					strainPTOcta = strain2 - ratio * (strain2 - strain1);
				}

        size = ratio1;
        elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);
        if ( (2.*refShearModulus - elasto_plast_modul) <= 0) 
					plast_modul = UP_LIMIT;
        else 
					plast_modul = (2.*refShearModulus * elasto_plast_modul)/
                        (2.*refShearModulus - elasto_plast_modul);
        if (plast_modul < 0) plast_modul = 0;
        if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
        if (ii==numOfSurfaces) plast_modul = 0;
        committedSurfaces[ii] = MultiYieldSurface(zeroVector,size,plast_modul);
  }  // ii   
}


double PressureDependMultiYield::yieldFunc(const T2Vector & stress, 
																	const MultiYieldSurface * surfaces, int surfaceNum)
{
	double coneHeight = stress.volume() - residualPress;
	Vector temp = stress.deviator() - surfaces[surfaceNum].center()*coneHeight;
	double sz = surfaces[surfaceNum].size()*coneHeight;

  return 3./2.*(temp && temp) - sz * sz;
}


void PressureDependMultiYield::deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
																			int surfaceNum)
{
	double diff = yieldFunc(stress, surfaces, surfaceNum);
	double coneHeight = stress.volume() - residualPress;

	if ( surfaceNum < numOfSurfaces && diff < 0. ) {
		double sz = -surfaces[surfaceNum].size()*coneHeight;
		double deviaSz = sqrt(sz*sz + diff);
    Vector devia = stress.deviator(); 
	  Vector temp = devia - surfaces[surfaceNum].center()*coneHeight;
		double coeff = (sz-deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
	  devia += temp * coeff;
	  stress = T2Vector(devia, stress.volume());
    deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
	}

	if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
    double sz = -surfaces[surfaceNum].size()*coneHeight;
    Vector newDevia = stress.deviator() * sz/sqrt(diff+sz*sz);
    stress = T2Vector(newDevia, stress.volume());
	}
}


void PressureDependMultiYield::initSurfaceUpdate(void)
{
	if (committedActiveSurf == 0) return; 

	double coneHeight = - (currentStress.volume() - residualPress);
	Vector devia = currentStress.deviator();
	double Ms = sqrt(3./2.*(devia && devia));
  Vector newCenter;

  if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
    newCenter = devia * (1. - committedSurfaces[committedActiveSurf].size()*coneHeight / Ms); 
    newCenter = newCenter / -coneHeight;
    committedSurfaces[committedActiveSurf].setCenter(newCenter);
  }

  for (int i=1; i<committedActiveSurf; i++) {
   	newCenter = devia * (1. - committedSurfaces[i].size()*coneHeight / Ms);
    newCenter = newCenter / -coneHeight;
    committedSurfaces[i].setCenter(newCenter); 
		theSurfaces[i] = committedSurfaces[i];
  }
  activeSurfaceNum = committedActiveSurf;
}


void PressureDependMultiYield::initStrainUpdate(void)
{
	// elastic strain state
  double stressRatio = currentStress.deviatorRatio(residualPress);
	double ratio = (-currentStress.volume()+residualPress)/(-refPressure+residualPress);
	ratio = pow(ratio, 1.-pressDependCoeff);
	modulusFactor = getModulusFactor(currentStress);
	double shearCoeff = 1./(2.*refShearModulus*modulusFactor);
	double bulkCoeff = 1./(3.*refBulkModulus*modulusFactor);
	currentStrain = currentStress.deviator()*shearCoeff
		              + currentStress.volume()*bulkCoeff;
	double octalStrain = currentStrain.octahedralShear(1);
	if (octalStrain <= LOW_LIMIT) octalStrain = LOW_LIMIT;

	// plastic strain state (scaled from elastic strain)
	double scale, PPZLimit;
	if (stressRatio >= stressRatioPT) {  //above PT
		onPPZ = 2;
		prePPZStrainOcta = strainPTOcta * ratio;
    PPZLimit = getPPZLimits(1,currentStress);
		scale = sqrt(prePPZStrainOcta+PPZLimit)/octalStrain;
	}
	else {  // below PT
		onPPZ = -1;
		prePPZStrainOcta = octalStrain;
		if (prePPZStrainOcta > strainPTOcta * ratio) prePPZStrainOcta=strainPTOcta*ratio;
		scale = sqrt(prePPZStrainOcta)/octalStrain;
	}
	currentStrain = T2Vector(currentStrain.deviator()*scale, currentStrain.volume());
	PPZPivot = currentStrain;
}


double PressureDependMultiYield::getModulusFactor(T2Vector & stress)
{
	double conHeig = stress.volume() - residualPress;
  double scale = conHeig / (refPressure-residualPress);
           
  return pow(scale, pressDependCoeff); 
}


void PressureDependMultiYield::setTrialStress(T2Vector & stress)
{
	modulusFactor = getModulusFactor(stress);
  Vector devia = stress.deviator() 
		             + subStrainRate.deviator()*2.*refShearModulus*modulusFactor;
	double volume = stress.volume() 
		             + subStrainRate.volume()*3.*refBulkModulus*modulusFactor;
	if (volume > 0.) volume = 0.;
  trialStress = T2Vector(devia, volume);
}


int PressureDependMultiYield::setSubStrainRate(void)
{
	if (activeSurfaceNum==numOfSurfaces) return 1;
  if (strainRate==T2Vector(zeroVector)) return 0;

	double elast_plast_modulus;
	double conHeig = -(currentStress.volume() - residualPress);
	double factor = getModulusFactor(currentStress);
	if (activeSurfaceNum==0) elast_plast_modulus = 2*refShearModulus*factor;
	else {
    double plast_modulus = theSurfaces[activeSurfaceNum].modulus()*factor;
		elast_plast_modulus = 2*refShearModulus*factor*plast_modulus 
			                    / (2*refShearModulus*factor+plast_modulus);
	}
  Vector incre = strainRate.deviator()*elast_plast_modulus;
	T2Vector increStress = T2Vector(incre,0);
  double singleCross = theSurfaces[numOfSurfaces].size()*conHeig / numOfSurfaces;
  double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	
	//int numOfSub1 = strainRate.octahedralShear() / 1.0e-4;
  //if (numOfSub1 > numOfSub) numOfSub = numOfSub1;

	incre = strainRate.t2Vector() / numOfSub;
  subStrainRate = T2Vector(incre);

	return numOfSub;
}


T2Vector PressureDependMultiYield::getContactStress(void)
{
	double conHeig = trialStress.volume() - residualPress;
	Vector center = theSurfaces[activeSurfaceNum].center(); 
  Vector devia = trialStress.deviator() - center*conHeig;
  double Ms = sqrt(3./2.*(devia && devia));
  devia = devia * theSurfaces[activeSurfaceNum].size()*(-conHeig) / Ms + center*conHeig;

  return T2Vector(devia,trialStress.volume()); 
}


int PressureDependMultiYield::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  T2Vector surfaceNormal = getSurfaceNormal(currentStress);
  if (((trialStress.t2Vector() - currentStress.t2Vector()) 
		&& surfaceNormal.t2Vector()) < 0) return 1;

  return 0;   
}


T2Vector PressureDependMultiYield::getSurfaceNormal(const T2Vector & stress)
{
	double conHeig = stress.volume() - residualPress;
  Vector devia = stress.deviator();
	Vector center = theSurfaces[activeSurfaceNum].center(); 
	double sz = theSurfaces[activeSurfaceNum].size();
	double volume = conHeig*((center && center) - 2./3.*sz*sz) - (devia && center);
  T2Vector Q = T2Vector((devia-center*conHeig)*3., volume);
	return T2Vector(Q.unitT2Vector());
}


double PressureDependMultiYield::getPlasticPotential(const T2Vector & contactStress,
																					const T2Vector & surfaceNormal)
{
	double plasticPotential, contractRule, unloadRule, dilateRule, shearLoading, temp;

  double contactRatio = contactStress.deviatorRatio(residualPress);
	temp = contactRatio/stressRatioPT;
	double factorPT = (temp*temp - 1)/(temp*temp + 1)/3.;
	double volume = contactStress.volume();
  contractRule = factorPT*contractParam1*exp(contractParam2*volume/AtmoPress);
  if (contractRule > 0.) contractRule = -contractRule;
	temp = currentStress.volume() - pressureD;
	if (temp >= 0.) unloadRule = 0.;
	else {
		double conHeiD = pressureD-residualPress;
		double temp1 = sqrt(3./2.)*currentStress.deviatorLength() 
			             + stressRatioPT*conHeiD;
		temp = temp1/(-temp);
		if (temp < theSurfaces[numOfSurfaces].size()) 
			temp = theSurfaces[numOfSurfaces].size();
		temp1 = (reversalStress.volume()-residualPress)/conHeiD;
    unloadRule = -sqrt(3./2.)*surfaceNormal.deviatorLength()*temp1/temp;
  }

  double currentRatio = currentStress.deviatorRatio(residualPress);
  double trialRatio = trialStress.deviatorRatio(residualPress);
	shearLoading = currentStress.deviator() && trialStress.deviator();

	if (factorPT < 0.) {  //below PT
		if (pressureD == 0.) plasticPotential = contractRule;
		else if (trialStress.volume() >= pressureD) {
			pressureD = 0.;
      plasticPotential = contractRule;
		}
		else if (trialRatio > currentRatio && shearLoading >= 0.)
      plasticPotential = contractRule;
		else  plasticPotential = unloadRule;
	}

	else {  //above PT
    if (trialRatio > currentRatio && shearLoading >= 0.) {  //dilation
			if (pressureD == 0.) pressureD = currentStress.volume();
			reversalStress = currentStress;
      updatePPZ(contactStress);  
			if (onPPZ==-1 || onPPZ==1) return LOCK_VALUE; 
			if (trialStrain.volume()*3. > volumeLimit)  
				dilateRule = 0;
			else
        dilateRule = factorPT*dilateParam1*exp(dilateParam2*cumuDilateStrainOcta);
			return dilateRule;
		}
		else if (pressureD == 0.) plasticPotential = contractRule;
		else if (trialStress.volume() >= pressureD) {
      pressureD = 0.;
      plasticPotential = contractRule;
		}
		else plasticPotential = unloadRule;
  }

	if (onPPZ > 0) onPPZ = 0;
  if (onPPZ != -1) PPZTranslation(contactStress);  
	return plasticPotential;
}


void PressureDependMultiYield::updatePPZ(const T2Vector & contactStress)
{
	T2Vector distance; 

	// PPZ inactive if liquefyParam1==0.
	/*if (liquefyParam1==0.) {
		if (onPPZ==2) {
      distance = T2Vector(trialStrain.t2Vector() - PPZPivot.t2Vector());
      cumuDilateStrainOcta = distance.octahedralShear(1);
		}
		else if (onPPZ != 2) {
			onPPZ = 2;
      PPZPivot = trialStrain;
      cumuDilateStrainOcta = 0.;
		}
		return;
	}*/

	// dilation: calc. cumulated dilative strain
  if (onPPZ==2) {
		PPZPivot = trialStrain;
    distance = T2Vector(PPZPivot.t2Vector() - PPZCenter.t2Vector());
		cumuDilateStrainOcta = distance.octahedralShear(1) - PPZSize;
		if (cumuDilateStrainOcta > maxCumuDilateStrainOcta) 
			maxCumuDilateStrainOcta = cumuDilateStrainOcta;
		return;
	}

	// calc. PPZ size.
  double PPZLimit = getPPZLimits(1,contactStress);
  if (onPPZ==-1 || onPPZ==0) {
		double volume = -contactStress.volume();
		oppoPrePPZStrainOcta = prePPZStrainOcta;
		double ratio = (volume+residualPress)/(-refPressure+residualPress);
		ratio = pow(ratio, 1.-pressDependCoeff);
    prePPZStrainOcta = ratio * strainPTOcta;
		if (oppoPrePPZStrainOcta == 0.) oppoPrePPZStrainOcta = prePPZStrainOcta;
	}
  PPZSize = PPZLimit 
		        + (prePPZStrainOcta+oppoPrePPZStrainOcta+maxCumuDilateStrainOcta)/2.;

	// calc. new PPZ center.
	if (onPPZ==0 || onPPZ==1) { 
    distance = T2Vector(PPZPivot.t2Vector() - PPZCenter.t2Vector());
		double coeff = PPZSize/distance.octahedralShear(1);
    PPZCenter = T2Vector(PPZPivot.t2Vector() - distance.t2Vector()*coeff);
	}

  distance = T2Vector(trialStrain.t2Vector() - PPZCenter.t2Vector());	
	if (distance.octahedralShear(1) > PPZSize) {  //outside PPZ
		onPPZ = 2;
    PPZPivot = trialStrain;
    cumuDilateStrainOcta = 0.;
		cumuTranslateStrainOcta = 0.;
		if (PPZLimit == 0.) maxCumuDilateStrainOcta = 0.;
	}
  else {  //inside PPZ
		if (onPPZ == 0 || onPPZ == 1) PPZTranslation(contactStress);
    if (onPPZ == -1 || onPPZ == 0) lockStress = contactStress;
		if (onPPZ == 0) onPPZ = 1;
	}
}


void PressureDependMultiYield::PPZTranslation(const T2Vector & contactStress)
{
	//if (liquefyParam1==0.) return;

	double PPZTranslationLimit = getPPZLimits(2,contactStress);

	T2Vector distance = T2Vector(PPZPivot.deviator()-PPZCenter.deviator());
	double temp = subStrainRate.deviator() && distance.unitDeviator();
	if (temp > 0.) {
		cumuTranslateStrainOcta += temp;
		if (cumuTranslateStrainOcta <= PPZTranslationLimit) { // PPZ translation
		  PPZPivot = T2Vector(PPZPivot.t2Vector() + distance.unitDeviator() * temp);
      PPZCenter = T2Vector(PPZCenter.t2Vector() + distance.unitDeviator() * temp);
		}
	}
}


double PressureDependMultiYield::getPPZLimits(int which, const T2Vector & contactStress)
{
  double PPZLimit, temp;
  double volume = -contactStress.volume();

  if (volume >= liquefyParam1) PPZLimit = 0.;
  else {
    temp = volume*pi/liquefyParam1/2.;
    PPZLimit = liquefyParam2 * pow(cos(temp), liquefyParam3);
  }
  
  if (which==1) 
    return PPZLimit;
  else if (which==2) 
    return liquefyParam4 * PPZLimit;
  else {
    cerr << "FATAL:PressureDependMultiYield::getPPZLimits: unknown argument value" << endl;
    g3ErrorHandler->fatal("");
    return 0.0;
  }
}


double PressureDependMultiYield::getLoadingFunc(const T2Vector & contactStress, 
																			 const T2Vector & surfaceNormal,
																			 double plasticPotential,
																			 int crossedSurface)
{
	double loadingFunc, limit;
	double modul = theSurfaces[activeSurfaceNum].modulus();
  double temp1 = 2. * refShearModulus * modulusFactor 
		             * (surfaceNormal.deviator() && surfaceNormal.deviator()) ;
  double temp2 = 9. * refBulkModulus * modulusFactor 
		             * surfaceNormal.volume() * plasticPotential ;
	
  //for the first crossing 
  double temp = temp1 + temp2 + modul * modulusFactor;
	if (activeSurfaceNum == numOfSurfaces) 
		limit = theSurfaces[activeSurfaceNum-1].modulus() * modulusFactor /2.;
	else limit = modul * modulusFactor /2.;
	if (temp < limit) {
    plasticPotential = (temp2 + limit - temp)/(9. * refBulkModulus * modulusFactor 
		                   * surfaceNormal.volume());
		temp = limit;
	}
  loadingFunc = (surfaceNormal.t2Vector() 
		             && (trialStress.deviator()-contactStress.deviator()))/temp;

	if (loadingFunc < 0.) loadingFunc = 0;

   //for more than one crossing 
  if(crossedSurface) {
		temp = (theSurfaces[activeSurfaceNum-1].modulus() - modul)
			     /theSurfaces[activeSurfaceNum-1].modulus();
		loadingFunc *= temp;
	}

	return loadingFunc;
}


int PressureDependMultiYield::stressCorrection(int crossedSurface)
{
	T2Vector contactStress = getContactStress();
	T2Vector surfaceNormal = getSurfaceNormal(contactStress);
	double plasticPotential = getPlasticPotential(contactStress,surfaceNormal);
	if (plasticPotential==LOCK_VALUE && (onPPZ == -1 || onPPZ == 1)) {
    trialStress = lockStress;
		return 1;
	}
  double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, 
		                                  plasticPotential, crossedSurface);
	double volume = trialStress.volume() 
		             - plasticPotential*3*refBulkModulus*modulusFactor*loadingFunc;
	if (volume > 0.) volume = 0.;
  Vector devia = trialStress.deviator() 
		         - surfaceNormal.deviator()*2*refShearModulus*modulusFactor*loadingFunc;
  trialStress = T2Vector(devia, volume);
  deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);
	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
    return stressCorrection(1);  //recursive call
	}
	return 0;
}


void PressureDependMultiYield::updateActiveSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return;

  double A, B, C, X;
	T2Vector direction;
  Vector t1, t2, temp;
	double conHeig = trialStress.volume() - residualPress;
	Vector center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	Vector outcenter = theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();

  t1 = trialStress.deviator() - center*conHeig;
  t2 = (center - outcenter)*conHeig;
  A = t1 && t1;
  B = 2. * (t1 && t2);
  C = (t2 && t2) - 2./3.* outsize * outsize * conHeig * conHeig;
  X = secondOrderEqn(A,B,C,0);

  if (X < 1.-LOW_LIMIT){
    cerr << "FATAL:PressureDependMultiYield::updateActiveSurface(): error in Direction of surface motion.\n" 
			   << "X-1= " << X-1 <<" A= "<<A<<" B= "<<B<<" C= "<<C<<endl; 
    g3ErrorHandler->fatal("");
  }

  temp = (t1 * X + center*conHeig) * (1. - size / outsize) 
		     - (center - outcenter * size / outsize) * conHeig;
	direction = T2Vector(temp);

	X = direction.deviatorLength();
	if (X < LOW_LIMIT) X = LOW_LIMIT;
  temp = direction.deviator()/X;  

  A = conHeig * conHeig;
  B = 2 * conHeig * (t1 && temp);
  C = (t1 && t1) - 2./3.* size * size * conHeig * conHeig;
	if (C < 1.e-13) C = 0.;
	if (C < 0.) {
    cerr << "FATAL:PressureDependMultiYield::updateActiveSurface(): error in surface motion." << endl; 
    g3ErrorHandler->fatal("");
	}
  X = secondOrderEqn(A,B,C,1);  

  center -= temp * X;
  theSurfaces[activeSurfaceNum].setCenter(center);
}      


void PressureDependMultiYield::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

	double conHeig = currentStress.volume() - residualPress;
	Vector devia = currentStress.deviator();
	Vector center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
  Vector newCenter;

	for (int i=1; i<activeSurfaceNum; i++) {
    newCenter = devia - (devia - center*conHeig) * theSurfaces[i].size() / size;
		newCenter = newCenter/conHeig;
    theSurfaces[i].setCenter(newCenter);
	}
}


int PressureDependMultiYield:: isCrossingNextSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return 0;  

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;
  
  return 0;
}
 
