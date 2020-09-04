// $Revision: 1.41 $
// $Date: 2009-10-07 20:14:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureIndependMultiYield.cpp,v $

// Written: ZHY
// Created: August 2000
// Last Modified: September 2009
//
// PressureIndependMultiYield.cpp
// -------------------
//
#include <math.h>
#include <stdlib.h>
#include <PressureIndependMultiYield.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <MultiYieldSurface.h>


Matrix PressureIndependMultiYield::theTangent(6,6);
T2Vector PressureIndependMultiYield::subStrainRate;
int PressureIndependMultiYield::matCount=0;
int* PressureIndependMultiYield::loadStagex=0;  //=0 if elastic; =1 if plastic
int* PressureIndependMultiYield::ndmx=0;        //num of dimensions (2 or 3)
double* PressureIndependMultiYield::rhox=0;
double* PressureIndependMultiYield::frictionAnglex=0;
double* PressureIndependMultiYield::peakShearStrainx=0;
double* PressureIndependMultiYield::refPressurex=0;
double* PressureIndependMultiYield::cohesionx=0;
double* PressureIndependMultiYield::pressDependCoeffx=0;
int*    PressureIndependMultiYield::numOfSurfacesx=0;
double* PressureIndependMultiYield::residualPressx=0;

void* OPS_PressureIndependMultiYield()
{
    const int numParam = 6;
    const int totParam = 10;

    int argc = OPS_GetNumRemainingInputArgs() + 2;

    char * arg[] = {"nd", "rho", "refShearModul", "refBulkModul",
		    "cohesi", "peakShearStra",
		    "frictionAng (=0)", "refPress (=100)", "pressDependCoe (=0.0)",
		    "numberOfYieldSurf (=20)"};
    if (argc < (3+numParam)) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: nDMaterial PressureIndependMultiYield tag? " << arg[0];
	opserr << "? "<< "\n";
	opserr << arg[1] << "? "<< arg[2] << "? "<< arg[3] << "? "<< "\n";
	opserr << arg[4] << "? "<< arg[5] << "? "<< arg[6] << "? "<< "\n";
	opserr << arg[7] << "? "<< arg[8] << "? "<< arg[9] << "? "<<"\n";
	return 0;
    }
    
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid PressureIndependMultiYield tag" << "\n";
	return 0;
    }

    int nd;
    if (OPS_GetIntInput(&numdata, &nd) < 0) {
	opserr << "WARNING invalid PressureIndependMultiYield nd" << "\n";
	return 0;
    }

    double param[8];
    param[5] = 0.0;
    param[6] = 100.;
    param[7] = 0.0;
    numdata = 8;
    if (OPS_GetDoubleInput(&numdata, &param[0]) < 0) {
	opserr << "WARNING invalid PressureIndependMultiYield double inputs" << "\n";
	return 0;
    }

    int numberOfYieldSurf = 20;
    numdata = 1;
    if (OPS_GetIntInput(&numdata, &numberOfYieldSurf) < 0) {
	opserr << "WARNING invalid PressureIndependMultiYield numberOfYieldSurf" << "\n";
	return 0;
    }

    static double * gredu = 0;
    // user defined yield surfaces
    if (numberOfYieldSurf < 0 && numberOfYieldSurf > -40) {
	numberOfYieldSurf = -int(numberOfYieldSurf);
	numdata = int(2*numberOfYieldSurf);
	gredu = new double[numdata];
	if (OPS_GetDoubleInput(&numdata, gredu) < 0) {
	    opserr << "WARNING invalid PressureIndependMultiYield double inputs" << "\n";
	    return 0;
	}
    }

    PressureIndependMultiYield * temp =
	new PressureIndependMultiYield (tag, nd, param[0], param[1], param[2],
					param[3], param[4], param[5], param[6],
					param[7], numberOfYieldSurf, gredu);
    if (gredu != 0) {
	delete [] gredu;
	gredu = 0;
    }

    return temp;
}

PressureIndependMultiYield::PressureIndependMultiYield (int tag, int nd,
							double r, double refShearModul,
							double refBulkModul,
							double cohesi, double peakShearStra,
							double frictionAng, double refPress, double pressDependCoe,
							int numberOfYieldSurf, double * gredu)
 : NDMaterial(tag,ND_TAG_PressureIndependMultiYield), currentStress(),
   trialStress(), currentStrain(), strainRate()
{
  if (nd !=2 && nd !=3) {
    opserr << "FATAL:PressureIndependMultiYield:: dimension error" << endln;
    opserr << "Dimension has to be 2 or 3, you give nd= " << nd << endln;
    exit(-1);
  }
  if (refShearModul <= 0) {
    opserr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refShearModulus <= 0" << endln;
    exit(-1);
  }
  if (refBulkModul <= 0) {
    opserr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refBulkModulus <= 0" << endln;
    exit(-1);
  }
  if (frictionAng < 0.) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle < 0" << endln;
    opserr << "Will reset frictionAngle to zero." << endln;
    frictionAng = 0.;
  }
  if (frictionAng == 0. && cohesi <= 0. ) {
    opserr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: frictionAngle && cohesion <= 0." << endln;
    exit(-1);
  }
  if (cohesi <= 0) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: cohesion <= 0" << endln;
    opserr << "Will reset cohesion to zero." << endln;
    cohesi = 0.;
  }
  if (peakShearStra <= 0) {
    opserr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: peakShearStra <= 0" << endln;
    exit(-1);
  }
  if (refPress <= 0) {
    opserr << "FATAL:PressureIndependMultiYield::PressureIndependMultiYield: refPress <= 0" << endln;
    exit(-1);
  }
  if (pressDependCoe < 0) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: pressDependCoe < 0" << endln;
    opserr << "Will reset pressDependCoe to zero." << endln;
    pressDependCoe = 0.;
  }
  if (pressDependCoe > 0 && frictionAng == 0) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: pressDependCoe > 0 while frictionAngle = 0" << endln;
    opserr << "Will reset pressDependCoe to zero." << endln;
    pressDependCoe = 0.;
  }
  if (numberOfYieldSurf <= 0) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: numberOfSurfaces <= 0" << endln;
    opserr << "Will use 10 yield surfaces." << endln;
    numberOfYieldSurf = 10;
  }
  if (numberOfYieldSurf > 100) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: numberOfSurfaces > 100" << endln;
    opserr << "Will use 100 yield surfaces." << endln;
    numberOfYieldSurf = 100;
  }
  if (r < 0) {
    opserr << "WARNING:PressureIndependMultiYield::PressureIndependMultiYield: mass density < 0" << endln;
    opserr << "Will use rho = 0." << endln;
    r = 0.;
  }

  int * temp1 = loadStagex;
  int * temp2 = ndmx;
  double * temp3 = rhox;
  double * temp6 = frictionAnglex;
  double * temp7 = peakShearStrainx;
  double * temp8 = refPressurex;
  double * temp9 = cohesionx;
  double * temp10 = pressDependCoeffx;
  int * temp11 = numOfSurfacesx;
  double * temp12 = residualPressx;
  
  int newCount = matCount+1;
  loadStagex = new int[newCount];
  ndmx = new int[newCount];
  rhox = new double[newCount];
  frictionAnglex = new double[newCount];
  peakShearStrainx = new double[newCount];
  refPressurex = new double[newCount];
  cohesionx = new double[newCount];
  pressDependCoeffx = new double[newCount];
  numOfSurfacesx = new int[newCount];
  residualPressx = new double[newCount];

  for (int i=0; i<matCount; i++) {
    loadStagex[i] = temp1[i];
    ndmx[i] = temp2[i];
    rhox[i] = temp3[i];
    frictionAnglex[i] = temp6[i];
    peakShearStrainx[i] = temp7[i];
    refPressurex[i] = temp8[i];
    cohesionx[i] = temp9[i];
    pressDependCoeffx[i] = temp10[i];
    numOfSurfacesx[i] = temp11[i];
    residualPressx[i] = temp12[i];
  }

  if (matCount > 0) {
    delete [] temp1; delete [] temp2; delete [] temp3;
    delete [] temp6; delete [] temp7; delete [] temp8;
    delete [] temp9; delete [] temp10; delete [] temp11;
    delete [] temp12;
  }

  ndmx[matCount] = nd;
  loadStagex[matCount] = 0;   //default
  refShearModulus = refShearModul;
  refBulkModulus = refBulkModul;
  frictionAnglex[matCount] = frictionAng;
  peakShearStrainx[matCount] = peakShearStra;
  refPressurex[matCount] = -refPress;  //compression is negative
  cohesionx[matCount] = cohesi;
  pressDependCoeffx[matCount] = pressDependCoe;
  numOfSurfacesx[matCount] = numberOfYieldSurf;
  rhox[matCount] = r;

  e2p = 0;
  matN = matCount;
  matCount=newCount;

  theSurfaces = new MultiYieldSurface[numberOfYieldSurf+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numberOfYieldSurf+1];
  activeSurfaceNum = committedActiveSurf = 0;

  mGredu = gredu;
  setUpSurfaces(gredu);  // residualPress is calculated inside.
}


PressureIndependMultiYield::PressureIndependMultiYield ()
 : NDMaterial(0,ND_TAG_PressureIndependMultiYield),
   currentStress(), trialStress(), currentStrain(),
  strainRate(), theSurfaces(0), committedSurfaces(0)
{
  //does nothing
}


PressureIndependMultiYield::PressureIndependMultiYield (const PressureIndependMultiYield & a)
 : NDMaterial(a.getTag(),ND_TAG_PressureIndependMultiYield),
   currentStress(a.currentStress), trialStress(a.trialStress),
  currentStrain(a.currentStrain), strainRate(a.strainRate)
{
  matN = a.matN;
  e2p = a.e2p;
  refShearModulus = a.refShearModulus;
  refBulkModulus = a.refBulkModulus;

  int numOfSurfaces = numOfSurfacesx[matN];

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
  int loadStage = loadStagex[matN];
  double frictionAngle = frictionAnglex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  if (loadStage != 1 || e2p == 1) return;
  e2p = 1;

  if (currentStress.volume() > 0. && frictionAngle > 0.) {
    //opserr << "WARNING:PressureIndependMultiYield::elast2Plast(): material in tension." << endln;
    currentStress.setData(currentStress.deviator(),0);
  }

  paramScaling();  // scale surface parameters corresponding to initial confinement

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
    if (committedActiveSurf == numOfSurfaces) {
      //opserr <<"WARNING:PressureIndependMultiYield::elast2Plast(): stress out of failure surface"<<endln;
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
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  static Vector temp(6);
  if (ndm==3 && strain.Size()==6)
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[2] = 0.0;
    temp[3] = strain[2];
    temp[4] = 0.0;
    temp[5] = 0.0;
  }
  else {
    opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
    exit(-1);
  }

  //strainRate.setData(temp-currentStrain.t2Vector(1),1);
  temp -= currentStrain.t2Vector(1);
  strainRate.setData(temp, 1);

  return 0;
}


int PressureIndependMultiYield::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain)
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  static Vector temp(6);
  if (ndm==3 && strain.Size()==6)
    temp = strain;
  else if (ndm==2 && strain.Size()==3) {
    temp[0] = strain[0];
    temp[1] = strain[1];
    temp[3] = strain[2];
  }
  else {
    opserr << "Fatal:D2PressDepMYS:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
    exit(-1);
  }

  strainRate.setData(temp,1);
  return 0;
}


int PressureIndependMultiYield::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureIndependMultiYield::getTangent (void)
{
  int loadStage = loadStagex[matN];
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  if (loadStage == 1 && e2p == 0) elast2Plast();

  if (loadStage!=1) {  //linear elastic
    for (int i=0;i<6;i++)
      for (int j=0;j<6;j++) {
	theTangent(i,j) = 0.;
	if (i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }

  }
  else {
    double coeff;
    static Vector devia(6);

    /*if (committedActiveSurf > 0) {
      //devia = currentStress.deviator()-committedSurfaces[committedActiveSurf].center();
      devia = currentStress.deviator();
      devia -= committedSurfaces[committedActiveSurf].center();

      double size = committedSurfaces[committedActiveSurf].size();
      double plastModul = committedSurfaces[committedActiveSurf].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size;
    }*/
    if (activeSurfaceNum > 0) {
      //devia = currentStress.deviator()-committedSurfaces[committedActiveSurf].center();
      devia = trialStress.deviator();
      devia -= theSurfaces[activeSurfaceNum].center();

      double size = theSurfaces[activeSurfaceNum].size();
      double plastModul = theSurfaces[activeSurfaceNum].modulus();
      coeff = 6.*refShearModulus*refShearModulus/(2.*refShearModulus+plastModul)/size/size;
    }

    else coeff = 0.;

    for (int i=0;i<6;i++)
      for (int j=0;j<6;j++) {
	theTangent(i,j) = - coeff*devia[i]*devia[j];
        if (i==j) theTangent(i,j) += refShearModulus;
        if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
	if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
      }
  }

  if (ndm==3)
    return theTangent;
  else {
    static Matrix workM(3,3);
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


const Matrix & PressureIndependMultiYield::getInitialTangent (void)
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  for (int i=0;i<6;i++)
    for (int j=0;j<6;j++) {
      theTangent(i,j) = 0.;
      if (i==j) theTangent(i,j) += refShearModulus;
      if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus;
      if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.);
    }

  if (ndm==3)
    return theTangent;
  else {
    static Matrix workM(3,3);
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
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  int i;
  if (loadStage == 1 && e2p == 0) elast2Plast();

  if (loadStage!=1) {  //linear elastic
    //trialStrain.setData(currentStrain.t2Vector() + strainRate.t2Vector());
    getTangent();
    static Vector a(6);
    a = currentStress.t2Vector();
	a.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
    trialStress.setData(a);
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
      if (i==0)
	setTrialStress(currentStress);
      else
	setTrialStress(trialStress);
      if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
      if (activeSurfaceNum==0) activeSurfaceNum++;
      stressCorrection(0);
      updateActiveSurface();
    }
    //volume stress change
    double volum = refBulkModulus*(strainRate.volume()*3.);
    volum += currentStress.volume();
     //if (volum > 0) volum = 0.;
    trialStress.setData(trialStress.deviator(),volum);
  }

  if (ndm==3)
    return trialStress.t2Vector();
  else {
    static Vector workV(3);
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
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  currentStress = trialStress;

  //currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
  static Vector temp(6);
  temp = currentStrain.t2Vector();
  temp += strainRate.t2Vector();
  currentStrain.setData(temp);
  temp.Zero();
  strainRate.setData(temp);

  if (loadStage==1) {
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
  if (strcmp(code,"PressureIndependMultiYield") == 0 || strcmp(code,"PlaneStrain") == 0
      || strcmp(code,"ThreeDimensional") == 0) {
    PressureIndependMultiYield * copy = new PressureIndependMultiYield(*this);
    return copy;
  }

  return 0;
}


const char * PressureIndependMultiYield::getType (void) const
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int PressureIndependMultiYield::getOrder (void) const
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  return (ndm == 2) ? 3 : 6;
}


int PressureIndependMultiYield::setParameter(const char **argv, int argc, Parameter &param)
{
  if (argc < 2)
    return -1;

  int theMaterialTag;
  theMaterialTag = atoi(argv[1]);

  // check for material tag
  if (theMaterialTag == this->getTag()) {

    if (strcmp(argv[0],"updateMaterialStage") == 0) {
      return param.addObject(1, this);
    } else if (strcmp(argv[0],"shearModulus") == 0) {
      return param.addObject(10, this);
	} else if (strcmp(argv[0],"bulkModulus") == 0) {
      return param.addObject(11, this);
    } else if (strcmp(argv[0],"frictionAngle") == 0) {
      return param.addObject(12, this);
    } else if (strcmp(argv[0],"cohesion") == 0) {
      return param.addObject(13, this);
    }
  }
  return -1;
}

int PressureIndependMultiYield::updateParameter(int responseID, Information &info)
{    
  if (responseID == 1) {
    loadStagex[matN] = info.theInt;
  } else if (responseID==10) {
    refShearModulus = info.theDouble;
  } else if (responseID==11) {
    refBulkModulus = info.theDouble;
  } else if (responseID==12) {
    frictionAnglex[matN] = info.theDouble;
    double *g = 0;
    setUpSurfaces(g);
    paramScaling();
    initSurfaceUpdate();
  } else if (responseID==13) {
    cohesionx[matN] = info.theDouble;
    double *g = 0;
    setUpSurfaces(g);
    paramScaling();
    initSurfaceUpdate();
  }

  // used by BBarFourNodeQuadUP element
  else if (responseID==20 && ndmx[matN] == 2)
		ndmx[matN] = 0;

  return 0;
}


int PressureIndependMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
  int loadStage = loadStagex[matN];
  int ndm = ndmx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  double rho = rhox[matN];
  double frictionAngle = frictionAnglex[matN];
  double peakShearStrain = peakShearStrainx[matN];
  double refPressure = refPressurex[matN];
  double cohesion = cohesionx[matN];
  double pressDependCoeff = pressDependCoeffx[matN];
  double residualPress = residualPressx[matN];

  int i, res = 0;

  static ID idData(6);
  idData(0) = this->getTag();
  idData(1) = numOfSurfaces;
  idData(2) = loadStage;
  idData(3) = ndm;
  idData(4) = matN;
  idData(5) = matCount;

  res += theChannel.sendID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PressureIndependMultiYield::sendSelf -- could not send ID\n";
    return res;
  }

  Vector data(24+numOfSurfaces*8);
  static Vector temp(6);
  data(0) = rho;
  data(1) = refShearModulus;
  data(2) = refBulkModulus;
  data(3) = frictionAngle;
  data(4) = peakShearStrain;
  data(5) = refPressure;
  data(6) = cohesion;
  data(7) = pressDependCoeff;
  data(8) = residualPress;
  data(9) = e2p;
  data(10) = committedActiveSurf;
  data(11) = activeSurfaceNum;

  temp = currentStress.t2Vector();
  for(i = 0; i < 6; i++) data(i+12) = temp[i];

  temp = currentStrain.t2Vector();
  for(i = 0; i < 6; i++) data(i+18) = temp[i];

  for(i = 0; i < numOfSurfaces; i++) {
    int k = 24 + i*8;
    data(k) = committedSurfaces[i+1].size();
    data(k+1) = committedSurfaces[i+1].modulus();
    temp = committedSurfaces[i+1].center();
    data(k+2) = temp(0);
    data(k+3) = temp(1);
    data(k+4) = temp(2);
    data(k+5) = temp(3);
    data(k+6) = temp(4);
    data(k+7) = temp(5);
  }

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureIndependMultiYield::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}


int PressureIndependMultiYield::recvSelf(int commitTag, Channel &theChannel,
					 FEM_ObjectBroker &theBroker)
{
  int i, res = 0;

  static ID idData(6);

  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PressureIndependMultiYield::recvSelf -- could not recv ID\n";
    return res;
  }

  this->setTag((int)idData(0));
  int numOfSurfaces = idData(1);
  int loadStage = idData(2);
  int ndm = idData(3);
  matN = idData(4);

  int matCountSendSide = idData(5);

  Vector data(24+idData(1)*8);
  static Vector temp(6);

  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureIndependMultiYield::recvSelf -- could not recv Vector\n";
    return res;
  }

  double rho = data(0);
  refShearModulus = data(1);
  refBulkModulus = data(2);
  double frictionAngle = data(3);
  double peakShearStrain = data(4);
  double refPressure = data(5);
  double cohesion = data(6);
  double pressDependCoeff = data(7);
  double residualPress = data(8);
  e2p = data(9);
  committedActiveSurf = data(10);
  activeSurfaceNum = data(11);

  for(i = 0; i < 6; i++) 
    temp[i] = data(i+12);
  currentStress.setData(temp);

  for(i = 0; i < 6; i++) 
    temp[i] = data(i+18);
  currentStrain.setData(temp);

  if (committedSurfaces != 0) {
    delete [] committedSurfaces;
    delete [] theSurfaces;
  }

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];

  for(i = 0; i < numOfSurfaces; i++) {
    int k = 24 + i*8;
    temp(0) = data(k+2);
    temp(1) = data(k+3);
    temp(2) = data(k+4);
    temp(3) = data(k+5);
    temp(4) = data(k+6);
    temp(5) = data(k+7);
    committedSurfaces[i+1].setData(temp, data(k), data(k+1));
  }

  int *temp1, *temp2, *temp11;
  double *temp3, *temp6, *temp7, *temp8, *temp9, *temp10, *temp12;

  if (matCountSendSide > matCount) {

	 temp1 = loadStagex;
	 temp2 = ndmx;
	 temp3 = rhox;
	 temp6 = frictionAnglex;
	 temp7 = peakShearStrainx;
	 temp8 = refPressurex;
	 temp9 = cohesionx;
	 temp10 = pressDependCoeffx;
	 temp11 = numOfSurfacesx;
	 temp12 = residualPressx;

	 loadStagex = new int[matCountSendSide];
	 ndmx = new int[matCountSendSide];
	 rhox = new double[matCountSendSide];
	 frictionAnglex = new double[matCountSendSide];
	 peakShearStrainx = new double[matCountSendSide];
	 refPressurex = new double[matCountSendSide];
	 cohesionx = new double[matCountSendSide];
	 pressDependCoeffx = new double[matCountSendSide];
	 numOfSurfacesx = new int[matCountSendSide];
	 residualPressx = new double[matCountSendSide];
	 
	 for (int i=0; i<matCount; i++) {
	   loadStagex[i] = temp1[i];
	   ndmx[i] = temp2[i];
	   rhox[i] = temp3[i];
	   frictionAnglex[i] = temp6[i];
	   peakShearStrainx[i] = temp7[i];
	   refPressurex[i] = temp8[i];
	   cohesionx[i] = temp9[i];
	   pressDependCoeffx[i] = temp10[i];
	   numOfSurfacesx[i] = temp11[i];
	   residualPressx[i] = temp12[i];
	 }
	 if( matCount > 0 ) {
	   delete [] temp1; delete [] temp2; delete [] temp3;
	   delete [] temp6; delete [] temp7; delete [] temp8;
	   delete [] temp9; delete [] temp10; delete [] temp11;
	   delete [] temp12;
	 }
	 matCount = matCountSendSide;
  }
  
  loadStagex[matN] = loadStage;
  ndmx[matN] = ndm;
  numOfSurfacesx[matN] = numOfSurfaces;
  rhox[matN] = rho;
  frictionAnglex[matN] = frictionAngle;
  peakShearStrainx[matN] = peakShearStrain;
  refPressurex[matN] = refPressure;
  cohesionx[matN] = cohesion;
  pressDependCoeffx[matN] = pressDependCoeff;
  residualPressx[matN] = residualPress;

  return res;
}


Response*
PressureIndependMultiYield::setResponse (const char **argv, int argc, OPS_Stream &output)
{
  // begin change by Alborz Ghofrani - UW --- get only 6 components of stress
  if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
	  if ((argc > 1) && (atoi(argv[1]) > 2) && (atoi(argv[1]) < 8)) 
		 return new MaterialResponse(this, 2 + atoi(argv[1]), this->getStressToRecord(atoi(argv[1])));
	  else
		 return new MaterialResponse(this, 1, this->getCommittedStress());
	// end change by Alborz Ghofrani - UW


  else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
    return new MaterialResponse(this, 2, this->getCommittedStrain());

  else if (strcmp(argv[0],"tangent") == 0)
    return new MaterialResponse(this, 3, this->getTangent());

  else if (strcmp(argv[0],"backbone") == 0) {
    int numOfSurfaces = numOfSurfacesx[matN];
    static Matrix curv(numOfSurfaces+1,(argc-1)*2);
    for (int i=1; i<argc; i++)
      curv(0,(i-1)*2) = atoi(argv[i]);
    return new MaterialResponse(this, 4, curv);
  }
  else
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
  case 4:
    if (matInfo.theMatrix != 0)
      getBackbone(*(matInfo.theMatrix));
    return 0;
	// begin change by Alborz Ghofrani UW --- get 6 components of stress
  case 5:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getStressToRecord(3);
    return 0;
  case 6:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getStressToRecord(4);
    return 0;
  case 7:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getStressToRecord(5);
    return 0;
  case 8:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getStressToRecord(6);
    return 0;
  case 9:
    if (matInfo.theVector != 0)
      *(matInfo.theVector) = getStressToRecord(7);
    return 0;
	// end change by Alborz Ghofrani UW
  default:
    return -1;
  }
}


void PressureIndependMultiYield::getBackbone (Matrix & bb)
{
  double residualPress = residualPressx[matN];
  double refPressure = refPressurex[matN];
  double pressDependCoeff =pressDependCoeffx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  double vol, conHeig, scale, factor, shearModulus, stress1,
		     stress2, strain1, strain2, plastModulus, elast_plast, gre;

	for (int k=0; k<bb.noCols()/2; k++) {
		vol = bb(0,k*2);
		if (vol<=0.) {
			opserr <<k<< "\nNDMaterial " <<this->getTag()
			  <<": invalid confinement for backbone recorder, " << vol << endln;
			continue;
		}
		conHeig = vol + residualPress;
		scale = -conHeig / (refPressure-residualPress);
		factor = pow(scale, pressDependCoeff);
		shearModulus = factor*refShearModulus;
		for (int i=1; i<=numOfSurfaces; i++) {
			if (i==1) {
				stress2 = committedSurfaces[i].size()*factor/sqrt(3.0);
				strain2 = stress2/shearModulus;
				bb(1,k*2) = strain2; bb(1,k*2+1) = shearModulus;
			} else {
				stress1 = stress2; strain1 = strain2;
				plastModulus = factor*committedSurfaces[i-1].modulus();
				elast_plast = 2*shearModulus*plastModulus/(2*shearModulus+plastModulus);
				stress2 = factor*committedSurfaces[i].size()/sqrt(3.0);
			  strain2 = 2*(stress2-stress1)/elast_plast + strain1;
				gre = stress2/strain2;
                bb(i,k*2) = strain2; bb(i,k*2+1) = gre;
			}
		}
	}

}

void PressureIndependMultiYield::Print(OPS_Stream &s, int flag )
{
  s << "PressureIndependMultiYield - loadStage: " <<  loadStagex[matN] << endln;
}


const Vector & PressureIndependMultiYield::getCommittedStress (void)
{
	int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;
	int numOfSurfaces = numOfSurfacesx[matN];

	double scale = sqrt(3./2.)*currentStress.deviatorLength()/committedSurfaces[numOfSurfaces].size();
	if (loadStagex[matN] != 1) scale = 0.;
	if (ndm==3) {
		static Vector temp7(7), temp6(6);
		temp6 = currentStress.t2Vector();
        temp7[0] = temp6[0];
        temp7[1] = temp6[1];
        temp7[2] = temp6[2];
        temp7[3] = temp6[3];
        temp7[4] = temp6[4];
        temp7[5] = temp6[5];
        temp7[6] = scale;
	    return temp7;
	}
    else {
        static Vector temp5(5), temp6(6);
		temp6 = currentStress.t2Vector();
        temp5[0] = temp6[0];
        temp5[1] = temp6[1];
        temp5[2] = temp6[2];
        temp5[3] = temp6[3];
        temp5[4] = scale;
        return temp5;
    }
}

// begin change by Alborz Ghofrani - UW --- get 6 components of stress
const Vector & PressureIndependMultiYield::getStressToRecord (int numOutput)
{
  int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;

  if (ndm==3) {
	static Vector temp7(7);
	temp7 = this->getCommittedStress();
	if (numOutput == 6)
	{
		static Vector temp6(6);
		temp6[0] = temp7[0];
		temp6[1] = temp7[1];
		temp6[2] = temp7[2];
		temp6[3] = temp7[3];
		temp6[4] = temp7[4];
		temp6[5] = temp7[5];
		return temp6;
	} else if (numOutput == 7) 
	{
		return temp7;
	} else {
		opserr << "Wrong number of stress components to record!" << endln;
		return temp7;
	}
  }

  else {
    static Vector temp5(5);
	temp5 = this->getCommittedStress();
	if (numOutput == 3)
	{
		static Vector temp3(3);
		temp3[0] = temp5[0];
		temp3[1] = temp5[1];
		temp3[2] = temp5[3];
		return temp3;
	} else if (numOutput == 4) 
	{
		static Vector temp4(4);
		temp4[0] = temp5[0];
		temp4[1] = temp5[1];
		temp4[2] = temp5[2];
		temp4[3] = temp5[3];
		return temp4;
	} else if (numOutput == 5) 
	{
		return temp5;
	} else {
		opserr << "Wrong number of stress components to record!" << endln;
		return temp5;
	}
  }
}
// end change by Alborz Ghofrani - UW 

const Vector & PressureIndependMultiYield::getCommittedStrain (void)
{
	int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;

  if (ndm==3)
    return currentStrain.t2Vector(1);
  else {
    static Vector workV(3), temp6(6);
		temp6 = currentStrain.t2Vector(1);
    workV[0] = temp6[0];
    workV[1] = temp6[1];
    workV[2] = temp6[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used
void PressureIndependMultiYield::setUpSurfaces (double * gredu)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff =pressDependCoeffx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
	double peakShearStrain = peakShearStrainx[matN];

	double  stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double pi = 3.14159265358979;
	double refStrain, peakShear, coneHeight;

	if (gredu == 0) {  //automatic generation of surfaces
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
			peakShear = 2*sqrt(2.)*cohesion/3;
			refStrain = (peakShearStrain * peakShear)
				/ (refShearModulus * peakShearStrain - peakShear);
			residualPress = 0.;
		}

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

			static Vector temp(6);
			committedSurfaces[ii] = MultiYieldSurface(temp,size,plast_modul);
		}  // ii
	}
	else {  //user defined surfaces
		if (frictionAngle > 0) {   // ignore user defined frictionAngle
			int ii = 2*(numOfSurfaces-1);
			double tmax = refShearModulus*gredu[ii]*gredu[ii+1];
			double Mnys = -(sqrt(3.) * tmax - 2. * cohesion) / refPressure;
			if (Mnys <= 0) {   // also ignore user defined cohesion
				cohesion = sqrt(3.)/2 * tmax;
				frictionAngle = 0.;
				coneHeight = 1.;
				residualPress = 0.;
			} else {
				double sinPhi = 3*Mnys /(6+Mnys);
				if (sinPhi<0. || sinPhi>1.) {
					opserr <<"\nNDMaterial " <<this->getTag()<<": Invalid friction angle, please modify ref. pressure or G/Gmax curve."<<endln;
					exit(-1);
				}
				residualPress = 2. * cohesion / Mnys;
				if (residualPress < 0.01*refPressure) residualPress = 0.01*refPressure;
				coneHeight = - (refPressure - residualPress);
				frictionAngle = asin(sinPhi)*180/pi;
			}
		}  else if (frictionAngle == 0.) {   // ignore user defined cohesion
			int ii = 2*(numOfSurfaces-1);
			double tmax = refShearModulus*gredu[ii]*gredu[ii+1];
			cohesion = sqrt(3.)/2 * tmax;
			coneHeight = 1.;
			residualPress = 0.;
		}

		/*
		opserr << "\nNDMaterial " <<this->getTag()<<": Friction angle = "<<frictionAngle
			<<", Cohesion = "<<cohesion<<"\n"<<endln;
		*/


		if (frictionAngle == 0.) pressDependCoeff = 0.; // ignore user defined pressDependCoeff

		for (int i=1; i<numOfSurfaces; i++) {
			int ii = 2*(i-1);
			strain1 = gredu[ii];
			stress1 = refShearModulus*gredu[ii+1]*strain1;
			strain2 = gredu[ii+2];
			stress2 = refShearModulus*gredu[ii+3]*strain2;

			size = sqrt(3.) * stress1 / coneHeight;
			elasto_plast_modul = 2.*(stress2 - stress1)/(strain2 - strain1);
			if ( (2.*refShearModulus - elasto_plast_modul) <= 0)
				plast_modul = UP_LIMIT;
			else
				plast_modul = (2.*refShearModulus * elasto_plast_modul)/
				(2.*refShearModulus - elasto_plast_modul);
			if (plast_modul <= 0) {
				opserr << "\nNDMaterial " <<this->getTag()<<": Surface " << i
					<< " has plastic modulus < 0.\n Please modify G/Gmax curve.\n"<<endln;
				exit(-1);
			}
			if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;

			static Vector temp(6);
			committedSurfaces[i] = MultiYieldSurface(temp,size,plast_modul);

			if (i==(numOfSurfaces-1)) {
				plast_modul = 0;
				size = sqrt(3.) * stress2 / coneHeight;
				committedSurfaces[i+1] = MultiYieldSurface(temp,size,plast_modul);
			}
		}
	}

	residualPressx[matN] = residualPress;
	frictionAnglex[matN] = frictionAngle;
	cohesionx[matN] = cohesion;
}


double PressureIndependMultiYield::yieldFunc(const T2Vector & stress,
											 const MultiYieldSurface * surfaces, int surfaceNum)
{
	static Vector temp(6);
	//temp = stress.deviator() - surfaces[surfaceNum].center();
	temp = stress.deviator();
	temp -= surfaces[surfaceNum].center();

	double sz = surfaces[surfaceNum].size();
	return 3./2.*(temp && temp) - sz * sz;
}


void PressureIndependMultiYield::deviatorScaling(T2Vector & stress, const MultiYieldSurface * surfaces,
																			int surfaceNum, int count)
{
	count++;
	int numOfSurfaces = numOfSurfacesx[matN];

	double diff = yieldFunc(stress, surfaces, surfaceNum);

	if ( surfaceNum < numOfSurfaces && diff < 0. ) {
		double sz = surfaces[surfaceNum].size();
		double deviaSz = sqrt(sz*sz + diff);
		static Vector devia(6);
		devia = stress.deviator();
		static Vector temp(6);
		temp = devia - surfaces[surfaceNum].center();
		double coeff = (sz-deviaSz) / deviaSz;
		if (coeff < 1.e-13) coeff = 1.e-13;
		devia.addVector(1.0, temp, coeff);
		stress.setData(devia, stress.volume());
		deviatorScaling(stress, surfaces, surfaceNum, count);  // recursive call
	}

	if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
		double sz = surfaces[surfaceNum].size();
		static Vector newDevia(6);
		newDevia.addVector(0.0, stress.deviator(), sz/sqrt(diff+sz*sz));
		stress.setData(newDevia, stress.volume());
	}
}


void PressureIndependMultiYield::initSurfaceUpdate()
{
	if (committedActiveSurf == 0) return;

	int numOfSurfaces = numOfSurfacesx[matN];

	static Vector devia(6);
	devia = currentStress.deviator();
	double Ms = sqrt(3./2.*(devia && devia));
	static Vector newCenter(6);

	if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
		//newCenter = devia * (1. - committedSurfaces[activeSurfaceNum].size() / Ms);
		newCenter.addVector(0.0, devia, 1.0-committedSurfaces[committedActiveSurf].size()/Ms);
		committedSurfaces[committedActiveSurf].setCenter(newCenter);
	}

	for (int i=1; i<committedActiveSurf; i++) {
	  newCenter = devia * (1. - committedSurfaces[i].size() / Ms);
	  committedSurfaces[i].setCenter(newCenter);
	}
}


void PressureIndependMultiYield::paramScaling(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];

	if (frictionAngle == 0.) return;

	double conHeig = - (currentStress.volume() - residualPress);
	double scale = -conHeig / (refPressure-residualPress);

	scale = pow(scale, pressDependCoeff);
	refShearModulus *= scale;
   	refBulkModulus *= scale;

	double plastModul, size;
	static Vector temp(6);
	for (int i=1; i<=numOfSurfaces; i++) {
	  plastModul = committedSurfaces[i].modulus() * scale;
	  size = committedSurfaces[i].size() * conHeig;
	  committedSurfaces[i] =  MultiYieldSurface(temp,size,plastModul);
	}

}


void PressureIndependMultiYield::setTrialStress(T2Vector & stress)
{
  static Vector devia(6);
  //devia = stress.deviator() + subStrainRate.deviator()*2.*refShearModulus;
  devia = stress.deviator();
  devia.addVector(1.0, subStrainRate.deviator(), 2.*refShearModulus);

  trialStress.setData(devia, stress.volume());
}


int PressureIndependMultiYield::setSubStrainRate(void)
{
    int numOfSurfaces = numOfSurfacesx[matN];

	//if (activeSurfaceNum==numOfSurfaces) return 1;

	//if (strainRate==T2Vector()) return 0;
	if (strainRate.isZero()) return 0;

	double elast_plast_modulus;
	if (activeSurfaceNum==0)
	  elast_plast_modulus = 2*refShearModulus;
	else {
	  double plast_modulus = theSurfaces[activeSurfaceNum].modulus();
	  elast_plast_modulus = 2*refShearModulus*plast_modulus
	    / (2*refShearModulus+plast_modulus);
	}
	static Vector incre(6);
	//incre = strainRate.deviator()*elast_plast_modulus;
	incre.addVector(0.0, strainRate.deviator(),elast_plast_modulus);

	static T2Vector increStress;
	increStress.setData(incre, 0);
	double singleCross = theSurfaces[numOfSurfaces].size() / numOfSurfaces;
	double totalCross = 3.*increStress.octahedralShear() / sqrt(2.);
	int numOfSub = totalCross/singleCross + 1;
	if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;
	//incre = strainRate.t2Vector() / numOfSub;
	incre = strainRate.t2Vector();
	incre /= numOfSub;
	subStrainRate.setData(incre);

	return numOfSub;
}


void
PressureIndependMultiYield::getContactStress(T2Vector &contactStress)
{
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	static Vector devia(6);
	//devia = trialStress.deviator() - center;
	devia = trialStress.deviator();
	devia -= center;

	double Ms = sqrt(3./2.*(devia && devia));
	//devia = devia * theSurfaces[activeSurfaceNum].size() / Ms + center;
	devia *= theSurfaces[activeSurfaceNum].size() / Ms;
	devia += center;

	contactStress.setData(devia,trialStress.volume());
}


int PressureIndependMultiYield::isLoadReversal(void)
{
  if(activeSurfaceNum == 0) return 0;

  static Vector surfaceNormal(6);
  getSurfaceNormal(currentStress, surfaceNormal);

  //(((trialStress.deviator() - currentStress.deviator()) && surfaceNormal) < 0)
  // return 1;
  static Vector a(6);
  a = trialStress.deviator();
  a-= currentStress.deviator();
  if((a && surfaceNormal) < 0)
    return 1;

  return 0;
}


void
PressureIndependMultiYield::getSurfaceNormal(const T2Vector & stress, Vector &surfaceNormal)
{
  //Q = stress.deviator() - theSurfaces[activeSurfaceNum].center();
  // return Q / sqrt(Q && Q);

  surfaceNormal = stress.deviator();
  surfaceNormal -= theSurfaces[activeSurfaceNum].center();
  surfaceNormal /= sqrt(surfaceNormal && surfaceNormal);
}


double PressureIndependMultiYield::getLoadingFunc(const T2Vector & contactStress,
									 const Vector & surfaceNormal, int crossedSurface)
{
  double loadingFunc;
  double temp1 = 2. * refShearModulus ;
  double temp2 = theSurfaces[activeSurfaceNum].modulus();

  //for crossing first surface
  double temp = temp1 + temp2;
  //loadingFunc = (surfaceNormal && (trialStress.deviator()-contactStress.deviator()))/temp;
  static Vector tmp(6);
  tmp =trialStress.deviator();
  tmp -= contactStress.deviator();
  loadingFunc = (surfaceNormal && tmp)/temp;
   //for crossing more than one surface
  if(crossedSurface) {
    double temp3 = theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= (temp3 - temp2)/temp3;
  }

  return loadingFunc;
}


void PressureIndependMultiYield::stressCorrection(int crossedSurface)
{
	static T2Vector contactStress;
	this->getContactStress(contactStress);
	static Vector surfaceNormal(6);
	this->getSurfaceNormal(contactStress, surfaceNormal);
	double loadingFunc = getLoadingFunc(contactStress, surfaceNormal, crossedSurface);
	static Vector devia(6);

	//devia = trialStress.deviator() - surfaceNormal * 2 * refShearModulus * loadingFunc;
	devia.addVector(0.0, surfaceNormal, -2*refShearModulus*loadingFunc);
	devia += trialStress.deviator();

	trialStress.setData(devia, trialStress.volume());
	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		stressCorrection(1);  //recursive call
	}
}


void PressureIndependMultiYield::updateActiveSurface(void)
{
  int numOfSurfaces = numOfSurfacesx[matN];

  if (activeSurfaceNum == numOfSurfaces) return;

	double A, B, C, X;
	static T2Vector direction;
	static Vector t1(6);
	static Vector t2(6);
	static Vector temp(6);
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector outcenter(6);
	outcenter= theSurfaces[activeSurfaceNum+1].center();
	double outsize = theSurfaces[activeSurfaceNum+1].size();


	//t1 = trialStress.deviator() - center;
	//t2 = center - outcenter;
	t1 = trialStress.deviator();
	t1 -= center;
	t2 = center;
	t2 -= outcenter;

	A = t1 && t1;
	B = 2. * (t1 && t2);
	C = (t2 && t2) - 2./3.* outsize * outsize;
	X = secondOrderEqn(A,B,C,0);
	if ( fabs(X-1.) < LOW_LIMIT ) X = 1.;
	if (X < 1.){
	  opserr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in Direction of surface motion." 
	       << endln;
	  exit(-1);
	}

	//temp = (t1 * X + center) * (1. - size / outsize) - (center - outcenter * size / outsize);
	temp = center;
	temp.addVector(1.0, t1, X);
	temp *= (1.0 - size/outsize);
	t2 = center;
	t2.addVector(1.0, outcenter, -size/outsize);
	temp -= t2;

	direction.setData(temp);

	if (direction.deviatorLength() < LOW_LIMIT) return;

	temp = direction.deviator();
	A = temp && temp;
	B = - 2 * (t1 && temp);
	if (fabs(B) < LOW_LIMIT) B = 0.;
	C = (t1 && t1) - 2./3.* size * size;
	if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;

	// Added by Alborz Ghofrani UW to avoid solving quadratic equations with very small coefficients B and C
	if ( (fabs(B) < 1.0e-10) && (fabs(C) < 1.0e-10) ){
		return;
	}
	// End Alborz Ghofrani UW

	if (B > 0. || C < 0.) {
	  opserr << "FATAL:PressureIndependMultiYield::updateActiveSurface(): error in surface motion.\n"
	       << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endln;
	  exit(-1);
	}
	X = secondOrderEqn(A,B,C,1);

	//center += temp * X;
	center.addVector(1.0, temp, X);
	theSurfaces[activeSurfaceNum].setCenter(center);
}


void PressureIndependMultiYield::updateInnerSurface(void)
{
	if (activeSurfaceNum <= 1) return;

	static Vector devia(6);
	devia = currentStress.deviator();
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	static Vector newcenter(6);

	for (int i=1; i<activeSurfaceNum; i++) {
		//newcenter = devia - (devia - center) * theSurfaces[i].size() / size;
		newcenter = center;
		newcenter -= devia;
		newcenter *= theSurfaces[i].size()/size;
		newcenter += devia;

		theSurfaces[i].setCenter(newcenter);
	}
}


int PressureIndependMultiYield:: isCrossingNextSurface(void)
{
  int numOfSurfaces = numOfSurfacesx[matN];
  if (activeSurfaceNum == numOfSurfaces) return 0;

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;

  return 0;
}


