// $Revision: 1.1 $
// $Date: 2000-12-19 03:35:02 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureIndependMultiYield.cpp,v $
                                                                        
// Written: ZHY
// Created: August 2000

//
// PressureIndependMultiYield.cpp
// -------------------
//
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>
#include <math.h>
#include <stdlib.h>
#include "PressureIndependMultiYield.h"
#include <Information.h>

int PressureIndependMultiYield::loadStage = 0;
Matrix PressureIndependMultiYield::theTangent = Matrix(6,6);
T2Vector PressureIndependMultiYield::subStrainRate = Vector(6);
Vector PressureIndependMultiYield::workV = Vector(3);
Matrix PressureIndependMultiYield::workM = Matrix(3,3);

const Vector zeroVector(6);

PressureIndependMultiYield::PressureIndependMultiYield (int tag, int nd, double refShearModul,
    			  		                  double refBulkModul, double frictionAng,
								                  double peakShearStra, double refPress, 
																	double cohesi, 	double pressDependCoe,
																	int numberOfYieldSurf)
 : NDMaterial(tag,MAT_TAG_PressureIndependMultiYield), currentStress(zeroVector),
   trialStress(zeroVector), currentStrain(zeroVector), strainRate(zeroVector)
{
	if (nd !=2 && nd !=3) {
		cerr << "FATAL:PressureIndependMultiYield:: dimension error" << endl;
    cerr << "Dimension has to be 2 or 3, you give nd= " << nd << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refShearModul <= 0) {
		cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refShearModulus <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refBulkModul <= 0) {
		cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refBulkModulus <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
  if (frictionAng < 0.) {
  	cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle < 0" << endl;
    cerr << "Will reset frictionAngle to zero." << endl;
	  frictionAng = 0.;
  }
  if (frictionAng == 0. && cohesi <= 0. ) {
  	cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle && cohesion <= 0." << endl;
	  g3ErrorHandler->fatal("");
  }
	if (cohesi <= 0) {
		cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: cohesion <= 0" << endl;
    cerr << "Will reset cohesion to zero." << endl;
	  cohesi = 0.;
  }
	if (peakShearStra <= 0) {
		cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: peakShearStra <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (refPress <= 0) {
		cerr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refPress <= 0" << endl;
	  g3ErrorHandler->fatal("");
  }
	if (pressDependCoe < 0) {
		cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: pressDependCoe < 0" << endl;
    cerr << "Will reset pressDependCoe to zero." << endl;
	  pressDependCoe = 0.;
  }
	if (numberOfYieldSurf <= 0) {
		cerr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: numberOfSurfaces <= 0" << endl;
		cerr << "Will use 10 yield surfaces." << endl;
	  numberOfYieldSurf = 10;
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

	e2p = 0;
	theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces+1]; 
  activeSurfaceNum = committedActiveSurf = 0; 

	setUpSurfaces();  // residualPress is calculated inside.
}
   

PressureIndependMultiYield::PressureIndependMultiYield () 
 : NDMaterial(0,MAT_TAG_PressureIndependMultiYield), 
   currentStress(zeroVector), trialStress(zeroVector), currentStrain(zeroVector), 
	 strainRate(zeroVector)
{
	ndm = 3;
  loadStage = 0;   
  refShearModulus = 0.;
	refBulkModulus = 0.;
	frictionAngle = 0.;
	peakShearStrain = 0.;
	refPressure = 0.;  //compression is negative
	cohesion = 0.;
	pressDependCoeff = 0.;
	numOfSurfaces = 1;
  residualPress = 0.;

	e2p = 0;
	theSurfaces = new MultiYieldSurface[1];
  committedSurfaces = new MultiYieldSurface[1];
  activeSurfaceNum = committedActiveSurf = 0; 
}


PressureIndependMultiYield::PressureIndependMultiYield (const PressureIndependMultiYield & a)
 : NDMaterial(a.getTag(),MAT_TAG_PressureIndependMultiYield), 
   currentStress(a.currentStress), trialStress(a.trialStress), 
	 currentStrain(a.currentStrain), strainRate(a.strainRate)
{
	ndm = a.ndm;
  loadStage = a.loadStage;  
  refShearModulus = a.refShearModulus;
	refBulkModulus = a.refBulkModulus;
	frictionAngle = a.frictionAngle;
	peakShearStrain = a.peakShearStrain;
	refPressure = a.refPressure;
	cohesion = a.cohesion;
	pressDependCoeff = a.pressDependCoeff;
	numOfSurfaces = a.numOfSurfaces;
	residualPress = a.residualPress;

	e2p = a.e2p;
	committedActiveSurf = a.committedActiveSurf;
  activeSurfaceNum = a.activeSurfaceNum; 

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];  
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];  
    theSurfaces[i] = a.theSurfaces[i];  
	}
}


PressureIndependMultiYield::~PressureIndependMultiYield ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
	if (committedSurfaces != 0) delete [] committedSurfaces;
}


void PressureIndependMultiYield::elast2Plast(void)
{
	if (loadStage == 0 || e2p == 1) return;
	e2p = 1;

	if (currentStress.volume() > 0. && frictionAngle > 0.) {
  	cerr << "WARNING:PressureIndependMultiYield::elast2Plast(): material in tension." << endl;
    currentStress = T2Vector(currentStress.deviator(),0);
	}

	paramScaling();  // scale surface parameters corresponding to initial confinement

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
     if (committedActiveSurf == numOfSurfaces) {
        cerr <<"WARNING:PressureIndependMultiYield::elast2Plast(): stress out of failure surface"<<endl;
				deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
        initSurfaceUpdate();
				return;
     }
  } 
  committedActiveSurf--;
	initSurfaceUpdate();
}


int PressureIndependMultiYield::setTrialStrain (const Vector &strain)
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
		cerr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  strainRate = T2Vector(temp-currentStrain.t2Vector());

	return 0;
}


int PressureIndependMultiYield::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain)
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
		cerr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endl;
		cerr << "But strain vector size is: " << strain.Size() << endl;
		g3ErrorHandler->fatal("");
	}

  strainRate = T2Vector(temp);
	return 0;
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureIndependMultiYield::getTangent (void)
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
	  double coeff;  Vector devia(6);

	  if (committedActiveSurf > 0) {
		  devia = currentStress.deviator()-committedSurfaces[committedActiveSurf].center();
	    double size = committedSurfaces[committedActiveSurf].size();
	    double plastModul = committedSurfaces[committedActiveSurf].modulus();
	    coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size;
		}
	  else coeff = 0.;

	  for (int i=0;i<6;i++) 
		  for (int j=0;j<6;j++) {
			  theTangent(i,j) = - coeff*devia[i]*devia[j];
        if (i==j) theTangent(i,j) += 2.*refShearModulus;
		  	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
			}
  }

	if (ndm==3) 
		return theTangent;
	else {
	  workM(0,0) = theTangent(0,0);
	  workM(0,1) = theTangent(0,1);
	  workM(0,2) = theTangent(0,3);
	  workM(1,0) = theTangent(1,0);
	  workM(1,1) = theTangent(1,1);
	  workM(1,2) = theTangent(1,3);
	  workM(2,0) = theTangent(3,0);
	  workM(2,1) = theTangent(3,1);
	  workM(2,2) = theTangent(3,3);
  	return workM;
	}
}


const Vector & PressureIndependMultiYield::getStress (void)
{
  int i;
	if (loadStage != 0 && e2p == 0) elast2Plast();

	if (loadStage==0) {  //linear elastic
		Vector trialStrain = currentStrain.t2Vector() + strainRate.t2Vector();
    getTangent();
    Vector a = theTangent * trialStrain;
		trialStress = T2Vector(a);
	}

	else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
		subStrainRate = strainRate;
    setTrialStress(currentStress);
  	if (isLoadReversal()) {
  		updateInnerSurface();
  		activeSurfaceNum = 0;
		}
  	int numSubIncre = setSubStrainRate();

  	for (i=0; i<numSubIncre; i++) {
  	  if (i==0)  setTrialStress(currentStress);
      else setTrialStress(trialStress);
  		if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;
  		stressCorrection(0);
      updateActiveSurface();
		}
	  //volume stress change
   	double volum = refBulkModulus*(strainRate.volume()*3.);
		volum += currentStress.volume();
		if (volum > 0) volum = 0.;
  	trialStress = T2Vector(trialStress.deviator(),volum);
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


const Vector & PressureIndependMultiYield::getStrain (void)
{
  return getCommittedStrain();
}


int PressureIndependMultiYield::commitState (void)
{
	currentStress = trialStress;
	currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
	if (loadStage) {
    committedActiveSurf = activeSurfaceNum;
	  for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
	}

	return 0;
}


int PressureIndependMultiYield::revertToLastCommit (void)
{
	return 0;
}


NDMaterial * PressureIndependMultiYield::getCopy (void)
{
  PressureIndependMultiYield * copy = new PressureIndependMultiYield(*this);
	return copy;
}


NDMaterial * PressureIndependMultiYield::getCopy (const char *code)
{
	if (strcmp(code,"PressureIndependMultiYield") == 0) {
     PressureIndependMultiYield * copy = new PressureIndependMultiYield(*this);
	   return copy;
	}

	return 0;
}


const char * PressureIndependMultiYield::getType (void) const
{
  return "PressureIndependMultiYield";
}


int PressureIndependMultiYield::getOrder (void) const
{
	return (ndm == 2) ? 3 : 6;
}


int PressureIndependMultiYield::updateParameter(int responseID, Information &info)
{
	loadStage = responseID;
	return 0;
}


int PressureIndependMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
	// Need to implement
	return 0;
}


int PressureIndependMultiYield::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)    
{
	// Need to implement
	return 0;
}


int PressureIndependMultiYield::getResponse (int responseID, Information &matInfo)
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


void PressureIndependMultiYield::Print(ostream &s, int flag )
{
	s << "PressureIndependMultiYield" << endl;
}


const Vector & PressureIndependMultiYield::getCommittedStress (void)
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


const Vector & PressureIndependMultiYield::getCommittedStrain (void)
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
void PressureIndependMultiYield::setUpSurfaces (void)
{ 
	double pi = 3.14159265358979;
  double refStrain, peakShear, coneHeight;

	if (frictionAngle > 0) {
		double sinPhi = sin(frictionAngle * pi/180.);
		double Mnys = 6.*sinPhi/(3.-sinPhi);
		residualPress = 3.* cohesion / (sqrt(2.) * Mnys);
		coneHeight = - (refPressure - residualPress);
    peakShear = sqrt(2.) * coneHeight * Mnys / 3.; 
    refStrain = (peakShearStrain * peakShear) 
			          / (refShearModulus * peakShearStrain - peakShear);
	}

	else if (frictionAngle == 0.) { // cohesion = peakShearStrength
    peakShear = cohesion;
    refStrain = (peakShearStrain * peakShear) 
			          / (refShearModulus * peakShearStrain - peakShear);
    residualPress = 0.;
	}

  double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double stressInc = peakShear / numOfSurfaces;

  for (int ii=1; ii<=numOfSurfaces; ii++){
        stress1 = ii * stressInc; 
				stress2 = stress1 + stressInc;
        strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
        strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);
        if (frictionAngle > 0.) size = 3. * stress1 / sqrt(2.) / coneHeight;
        else if (frictionAngle == 0.) size = 3. * stress1 / sqrt(2.);
 
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


double PressureIndependMultiYield::yieldFunc(const T2Vector & stress, 
																	const MultiYieldSurface * surfaces, int surfaceNum)
{
	Vector temp = stress.deviator() - surfaces[surfaceNum].center();
	double sz = surfaces[surfaceNum].size();
  return 3./2.*(temp && temp) - sz * sz;
}


void PressureIndependMultiYield::deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces, 
																			int surfaceNum)
{
	double diff = yieldFunc(stress, surfaces, surfaceNum);

	if ( surfaceNum < numOfSurfaces && diff < 0. ) {
		double sz = surfaces[surfaceNum].size();
		double deviaSz = sqrt(sz*sz + diff);
    Vector devia = stress.deviator(); 
	  Vector temp = devia - surfaces[surfaceNum].center();
		double coeff = (sz-deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
	  devia += temp * coeff;
	  stress = T2Vector(devia, stress.volume());

    deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
	}

	if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
    double sz = surfaces[surfaceNum].size();
    Vector newDevia = stress.deviator() * sz/sqrt(diff+sz*sz);
    stress = T2Vector(newDevia, stress.volume());
	}
}


void PressureIndependMultiYield::initSurfaceUpdate()
{
	if (activeSurfaceNum == 0) return; 

	Vector devia = currentStress.deviator();
	double Ms = sqrt(3./2.*(devia && devia));
  Vector newCenter;

  if (activeSurfaceNum < numOfSurfaces) { // failure surface can't move
    newCenter = devia * (1. - committedSurfaces[activeSurfaceNum].size() / Ms); 
    committedSurfaces[activeSurfaceNum].setCenter(newCenter);
  }

  for (int i=1; i<activeSurfaceNum; i++) {
   	newCenter = devia * (1. - committedSurfaces[i].size() / Ms);
    committedSurfaces[i].setCenter(newCenter); 
  }
}


void PressureIndependMultiYield::paramScaling(void)
{
	if (frictionAngle == 0.) return;

	double conHeig = - (currentStress.volume() - residualPress);
  double scale = -conHeig / (refPressure-residualPress);
           
  scale = pow(scale, pressDependCoeff); 
  refShearModulus *= scale;
   
  double plastModul, size;
	for (int i=1; i<=numOfSurfaces; i++) {
		 plastModul = committedSurfaces[i].modulus() * scale;
     size = committedSurfaces[i].size() * conHeig;
     committedSurfaces[i] =  MultiYieldSurface(zeroVector,size,plastModul);
	}
}


void PressureIndependMultiYield::setTrialStress(T2Vector & stress)
{
  Vector devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
  trialStress = T2Vector(devia, stress.volume());
}


int PressureIndependMultiYield::setSubStrainRate(void)
{
	if (activeSurfaceNum==numOfSurfaces) return 1;
  if (strainRate==T2Vector(zeroVector)) return 0;

	double elast_plast_modulus;
	if (activeSurfaceNum==0) elast_plast_modulus = 2*refShearModulus;
	else {
    double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
		elast_plast_modulus = 2*refShearModulus*plast_modulus / (2*refShearModulus+plast_modulus);
	}
  Vector incre = strainRate.deviator()*elast_plast_modulus;
	T2Vector increStress = T2Vector(incre,0);
  double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
  double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	incre = strainRate.t2Vector() / numOfSub;
  subStrainRate = T2Vector(incre);

	return numOfSub;
}


T2Vector PressureIndependMultiYield::getContactStress(void)
{
	Vector center = theSurfaces[activeSurfaceNum].center(); 
  Vector devia = trialStress.deviator() - center;
  double Ms = sqrt(3./2.*(devia && devia));
  devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;

  return T2Vector(devia,trialStress.volume()); 
}


int PressureIndependMultiYield::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  Vector surfaceNormal = getSurfaceNormal(currentStress);
  if(((trialStress.deviator() - currentStress.deviator()) && surfaceNormal) < 0) 
    return 1;

  return 0;   
}


Vector PressureIndependMultiYield::getSurfaceNormal(const T2Vector & stress)
{
  Vector Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();

  return Q / sqrt(Q && Q);
}


double PressureIndependMultiYield::getLoadingFunc(const T2Vector & contactStress, 
																			 const Vector & surfaceNormal,
																			 int crossedSurface)
{
	double loadingFunc;
  double temp1 = 2. * refShearModulus ;
  double temp2 = theSurfaces[activeSurfaceNum].modulus();

  //for crossing first surface
  double temp = temp1 + temp2;
  loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp;
  
   //for crossing more than one surface
  if(crossedSurface) {
		double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
		loadingFunc *= (temp3 - temp2)/temp3;
	}

	return loadingFunc;
}


void PressureIndependMultiYield::stressCorrection(int crossedSurface)
{
	T2Vector contactStress = getContactStress();
	Vector surfaceNormal = getSurfaceNormal(contactStress);
  double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, crossedSurface);
  Vector devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
  trialStress = T2Vector(devia, trialStress.volume());
  deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
    stressCorrection(1);  //recursive call
	}
}


void PressureIndependMultiYield::updateActiveSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return;

  double A, B, C, X;
	T2Vector direction;
  Vector t1, t2, temp;
	Vector center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	Vector outcenter = theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();

  t1 = trialStress.deviator() - center;
  t2 = center - outcenter;
  A = t1 && t1;
  B = 2. * (t1 && t2);
  C = (t2 && t2) - 2./3.* outsize * outsize;
  X = secondOrderEqn(A,B,C,0);

  if (X < 1. - LOW_LIMIT){
    cerr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in Direction of surface motion." 
			   << endl; 
    g3ErrorHandler->fatal("");
  }

  temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	direction = T2Vector(temp);

	X = direction.deviatorLength();
	if (X < LOW_LIMIT) X = LOW_LIMIT;
  temp = direction.deviator()/X;

  A = 1.;
  B = - 2 * (t1 && temp);
  C = (t1 && t1) - 2./3.* size * size;
	if (fabs(C) < LOW_LIMIT) C = 0.;
	if (C < 0.) {
    cerr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in surface motion." << endl; 
    g3ErrorHandler->fatal("");
	}
  X = secondOrderEqn(A,B,C,1);  

  center += temp * X;
  theSurfaces[activeSurfaceNum].setCenter(center);
}      


void PressureIndependMultiYield::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

	Vector devia = currentStress.deviator();
	Vector center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
  Vector newcenter;

	for (int i=1; i<activeSurfaceNum; i++) {
    newcenter = devia - (devia - center) * theSurfaces[i].size() / size;
    theSurfaces[i].setCenter(newcenter);
	}
}


int PressureIndependMultiYield:: isCrossingNextSurface(void)
{
  if (activeSurfaceNum == numOfSurfaces) return 0;  

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;
  
  return 0;
}
 
