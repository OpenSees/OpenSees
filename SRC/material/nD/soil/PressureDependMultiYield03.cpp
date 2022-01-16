// $Revision: 1.16 $
// $Date: 2018-12-18 20:14:00 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/soil/PressureDependMultiYield03.cpp,v $

// Written: ZHY
// Created: May 2004
// Last Modified: December 2018
//
// PressureDependMultiYield03.cpp
// -------------------
//

#include <math.h>
#include <stdlib.h>
#include <PressureDependMultiYield03.h>
#include <MultiYieldSurface.h>

#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>

int PressureDependMultiYield03::matCount=0;
int* PressureDependMultiYield03::loadStagex = 0;  //=0 if elastic; =1 if plastic
int* PressureDependMultiYield03::ndmx=0;  //num of dimensions (2 or 3)
double* PressureDependMultiYield03::rhox=0;
double* PressureDependMultiYield03::refShearModulusx=0;
double* PressureDependMultiYield03::refBulkModulusx=0;
double* PressureDependMultiYield03::frictionAnglex=0;
double* PressureDependMultiYield03::peakShearStrainx=0;
double* PressureDependMultiYield03::refPressurex=0;
double* PressureDependMultiYield03::cohesionx=0;
double* PressureDependMultiYield03::pressDependCoeffx=0;
int*    PressureDependMultiYield03::numOfSurfacesx=0;
double* PressureDependMultiYield03::phaseTransfAnglex=0;
int* PressureDependMultiYield03::mTypex = 0;
double* PressureDependMultiYield03::contractParam1x=0;
double* PressureDependMultiYield03::contractParam2x=0;
double* PressureDependMultiYield03::contractParam3x=0;
double* PressureDependMultiYield03::contractParam4x = 0;
double* PressureDependMultiYield03::contractParam5x = 0;
double* PressureDependMultiYield03::dilateParam1x=0;
double* PressureDependMultiYield03::dilateParam2x=0;
double* PressureDependMultiYield03::dilateParam3x = 0;
double* PressureDependMultiYield03::liquefyParam1x=0;
double* PressureDependMultiYield03::liquefyParam2x=0;
double* PressureDependMultiYield03::einitx=0;    //initial void ratio
double* PressureDependMultiYield03::volLimit1x=0;
double* PressureDependMultiYield03::volLimit2x=0;
double* PressureDependMultiYield03::volLimit3x=0;
double* PressureDependMultiYield03::residualPressx=0;
double* PressureDependMultiYield03::stressRatioPTx=0;

double* PressureDependMultiYield03::Hvx=0;
double* PressureDependMultiYield03::Pvx=0;

double PressureDependMultiYield03::pAtm = 101.;

Matrix PressureDependMultiYield03::theTangent(6,6);
T2Vector PressureDependMultiYield03::trialStrain;
T2Vector PressureDependMultiYield03::subStrainRate;
Vector PressureDependMultiYield03::workV6(6);
T2Vector PressureDependMultiYield03::workT2V;
const	double pi = 3.14159265358979;

void * OPS_ADD_RUNTIME_VPV(OPS_PressureDependMultiYield03)
{
    const int numParam = 18;
    const int totParam = 23;
    int tag;
    double param[totParam];
    param[numParam] = 20;
    param[numParam + 1] = 1.0;
    param[numParam + 2] = 0.;
    param[numParam + 3] = 101.;
    param[numParam + 4] = 1.73;

    const char* arg[] = {"nd", "rho", "refShearModul", "refBulkModul", "frictionAng",
                         "peakShearStra", "refPress", "pressDependCoe", "phaseTransformAngle",
                         "mType", "ca", "cb", "cc", "cd", "ce", "da", "db", "dc",
                         "numberOfYieldSurf (=20)", "liquefactionParam1=1.0", "liquefactionParam2=0.0",
                         "Atmospheric pressure (=101)", "cohesi (=1.73)"};

    int argc = OPS_GetNumRemainingInputArgs() + 2;
    if (argc < (3 + numParam)) {  
      // 3 refers to "nDMaterial PressureDependMultiYield03  $tag"
        opserr << "WARNING insufficient arguments\n";
        opserr << "Want: nDMaterial PressureDependMultiYield03 tag? " << arg[0];
        opserr << "? "
               << "\n";
        opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? "
               << "\n";
        opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? "
               << "\n";
        opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? "
               << "\n";
        opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? "
               << "\n";
        opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? "
               << "\n";
        opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? "
               << "\n";
        opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? "
               << "\n";
        opserr << arg[22] << "? "
               << "\n";
        return 0;
    }

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
        opserr << "WARNING invalid PressureDependMultiYield03 tag"
               << "\n";
        return 0;
    }

    int in = 22;
    for (int i = 3; (i < argc && i < in); i++) {
        if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
            opserr << "WARNING invalid " << arg[i - 3] << "\n";
            opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
            return 0;
        }
    }

    static double* gredu = 0;

    // user defined yield surfaces
    if (param[numParam] < 0 && param[numParam] > -100) {
        param[numParam] = -int(param[numParam]);
        gredu = new double[int(2 * param[numParam])];

        for (int i = 0; i < 2 * param[numParam]; i++) {
            if (OPS_GetDoubleInput(&numdata, &gredu[i]) < 0) {
                opserr << "WARNING invalid " << arg[i - 3]
                       << "\n";
                opserr << "nDMaterial PressureIndependMultiYield03: "
                       << tag << "\n";
                return 0;
            }
        }
    }

    if (gredu != 0) {
        for (int i = in + int(2 * param[numParam]); i < argc; i++) {
            if (OPS_GetDoubleInput(&numdata, &param[i - 3 - int(2 * param[numParam])]) < 0) {
                opserr << "WARNING invalid "
                       << arg[i - 3 - int(2 * param[numParam])]
                       << "\n";
                opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
                return 0;
            }
        }
    } else {
        for (int i = in; i < argc; i++) {
            if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
                opserr << "WARNING invalid "
                       << arg[i - 3 - int(2 * param[numParam])]
                       << "\n";
                opserr << "nDMaterial PressureDependMultiYield03: " << tag << "\n";
                return 0;
            }
        }
    }

    PressureDependMultiYield03* temp =
        new PressureDependMultiYield03(tag,
                                       param[0], param[1], param[2],
                                       param[3], param[4], param[5],
                                       param[6], param[7], param[8],
                                       param[9], param[10], param[11],
                                       param[12], param[13], param[14],
                                       param[15], param[16], param[17],
                                       param[18], gredu,
                                       param[19], param[20], param[21],
                                       param[22]);

    if (gredu != 0) {
        delete[] gredu;
        gredu = 0;
    }

    return temp;
}

PressureDependMultiYield03::PressureDependMultiYield03 (int tag, int nd,
						    double r, double refShearModul,
						    double refBulkModul, double frictionAng,
						    double peakShearStra, double refPress,
						    double pressDependCoe,
						    double phaseTransformAng,
							int mType,
							double ca,
							double cb,
							double cc,
							double cd,
							double ce,
							double da,
							double db,
							double dc,
					        int   numberOfYieldSurf,
						    double * gredu,
					        double liquefactionParam1,
						    double liquefactionParam2,
						    double atm, double cohesi)
 : NDMaterial(tag,ND_TAG_PressureDependMultiYield03), currentStress(),
   trialStress(), updatedTrialStress(), currentStrain(), strainRate(),
   PPZPivot(), PPZCenter(), PPZPivotCommitted(), PPZCenterCommitted(),
   PivotStrainRate(6), PivotStrainRateCommitted(6), check(0)
{
  if (nd !=2 && nd !=3) {
    opserr << "FATAL:PressureDependMultiYield03:: dimension error" << endln;
    opserr << "Dimension has to be 2 or 3, you give nd= " << nd << endln;
   exit(-1);
  }
  if (refShearModul <= 0) {
    opserr << "FATAL:PressureDependMultiYield03:: refShearModulus <= 0" << endln;
   exit(-1);
  }
  if (refBulkModul <= 0) {
    opserr << "FATAL:PressureDependMultiYield03:: refBulkModulus <= 0" << endln;
   exit(-1);
  }
  if (frictionAng <= 0.) {
    opserr << "FATAL:PressureDependMultiYield03:: frictionAngle <= 0" << endln;
   exit(-1);
  }
  if (frictionAng >= 90.) {
    opserr << "FATAL:PressureDependMultiYield03:: frictionAngle >= 90" << endln;
   exit(-1);
  }
  if (phaseTransformAng <= 0.) {
    opserr << "FATAL:PressureDependMultiYield03:: phaseTransformAng "<< phaseTransformAng << "<= 0" << endln;
   exit(-1);
  }
  /*if (phaseTransformAng > frictionAng) {
    opserr << "WARNING:PressureDependMultiYield03:: phaseTransformAng > frictionAng" << endln;
    opserr << "Will set phaseTransformAng = frictionAng." <<endln;
    phaseTransformAng = frictionAng+0.1;
  }*/
  if (cohesi < 0) {
    opserr << "WARNING:PressureDependMultiYield03:: cohesion < 0" << endln;
    opserr << "Will reset cohesion to 0.3." << endln;
    cohesi = 0.3;
  }
  if (peakShearStra <= 0) {
    opserr << "FATAL:PressureDependMultiYield03:: peakShearStra <= 0" << endln;
   exit(-1);
  }
  if (refPress <= 0) {
    opserr << "FATAL:PressureDependMultiYield03:: refPress <= 0" << endln;
   exit(-1);
  }
  if (pressDependCoe < 0) {
    opserr << "WARNING:PressureDependMultiYield03:: pressDependCoe < 0" << endln;
    opserr << "Will reset pressDependCoe to 0.5." << endln;
    pressDependCoe = 0.5;
  }
  if (numberOfYieldSurf <= 0) {
    opserr << "WARNING:PressureDependMultiYield03:: numberOfSurfaces " <<numberOfYieldSurf << "<= 0" << endln;
    opserr << "Will use 10 yield surfaces." << endln;
    numberOfYieldSurf = 10;
  }
  if (numberOfYieldSurf > 100) {
    opserr << "WARNING:PressureDependMultiYield03::PressureDependMultiYield03: numberOfSurfaces > 100" << endln;
    //opserr << "Will use 100 yield surfaces." << endln;
    //numberOfYieldSurf = 100;
  }
  if (r < 0) {
    opserr << "FATAL:PressureDependMultiYield03:: rho <= 0" << endln;
   exit(-1);
  }

  double ei = 0.6;
  double cs1 = 0.9;
  double cs2 = 0.02;
  double cs3 = 0.7;
  double hv = 0.0;
  double pv = 1.0;

  if (matCount%20 == 0) {
     int * temp1 = loadStagex;
	 int * temp2 = ndmx;
	 double * temp3 = rhox;
     double * temp4 = refShearModulusx;
     double * temp5 = refBulkModulusx;
     double * temp6 = frictionAnglex;
     double * temp7 = peakShearStrainx;
     double * temp8 = refPressurex;
     double * temp9 = cohesionx;
     double * temp10 = pressDependCoeffx;
	 int * temp11 = numOfSurfacesx;
     double * temp12 = residualPressx;
     double * temp13 = phaseTransfAnglex;
	 int * temp14 = mTypex;
	 double * temp15 = contractParam1x;
     double * temp16 = contractParam2x;
     double * temp17 = contractParam3x;
	 double * temp18 = contractParam4x;
	 double * temp19 = contractParam5x;

     double * temp20 = dilateParam1x;
     double * temp21 = dilateParam2x;
	 double * temp22 = dilateParam3x;

     double * temp23 = liquefyParam1x;
     double * temp24 = liquefyParam2x;
	 double * temp25 = stressRatioPTx;
	 
     loadStagex = new int[matCount+20];
     ndmx = new int[matCount+20];
     rhox = new double[matCount+20];
     refShearModulusx = new double[matCount+20];
     refBulkModulusx = new double[matCount+20];
     frictionAnglex = new double[matCount+20];
     peakShearStrainx = new double[matCount+20];
     refPressurex = new double[matCount+20];
	 cohesionx = new double[matCount+20];
     pressDependCoeffx = new double[matCount+20];
     numOfSurfacesx = new int[matCount+20];
     residualPressx = new double[matCount+20];
     phaseTransfAnglex = new double[matCount+20];

	 mTypex = new int[matCount + 20];
	 contractParam1x = new double[matCount+20];
     contractParam2x = new double[matCount+20];
     contractParam3x = new double[matCount+20];
	 contractParam4x = new double[matCount + 20];
	 contractParam5x = new double[matCount + 20];

     dilateParam1x = new double[matCount+20];
     dilateParam2x = new double[matCount+20];
	 dilateParam3x = new double[matCount + 20];

	 liquefyParam1x = new double[matCount+20];
     liquefyParam2x = new double[matCount+20];

	 stressRatioPTx = new double[matCount + 20];

     einitx = new double[matCount+20];    //initial void ratio
     volLimit1x = new double[matCount+20];
     volLimit2x = new double[matCount+20];
     volLimit3x = new double[matCount+20];
	 Hvx = new double[matCount+20];
	 Pvx = new double[matCount+20];
	 
	 for (int i=0; i<matCount; i++) {
         loadStagex[i] = temp1[i];
		 ndmx[i] = temp2[i];
         rhox[i] = temp3[i];
	     refShearModulusx[i] = temp4[i];
         refBulkModulusx[i] = temp5[i];
         frictionAnglex[i] = temp6[i];
         peakShearStrainx[i] = temp7[i];
         refPressurex[i] = temp8[i];
         cohesionx[i] = temp9[i];
         pressDependCoeffx[i] = temp10[i];
         numOfSurfacesx[i] = temp11[i];
         residualPressx[i] = temp12[i];
         phaseTransfAnglex[i] = temp13[i];

		 mTypex[i] = temp14[i];
		 contractParam1x[i] = temp15[i];
         contractParam2x[i] = temp16[i];
         contractParam3x[i] = temp17[i];
		 contractParam4x[i] = temp18[i];
		 contractParam5x[i] = temp19[i];

         dilateParam1x[i] = temp20[i];
         dilateParam2x[i] = temp21[i];
		 dilateParam3x[i] = temp22[i];

		 liquefyParam1x[i] = temp23[i];
         liquefyParam2x[i] = temp24[i];

         stressRatioPTx[i] = temp25[i];

		 einitx[i] = ei;    //default initial void ratio
		 volLimit1x[i] = cs1;
		 volLimit2x[i] = cs2;
		 volLimit3x[i] = cs3;
		 
		 Hvx[i] = hv;
		 Pvx[i] = pv;
     }

	 if (matCount > 0) {
	     delete [] temp1; delete [] temp2; delete [] temp3; delete [] temp4;
	     delete [] temp5; delete [] temp6; delete [] temp7; delete [] temp8;
	     delete [] temp9; delete [] temp10; delete [] temp11; delete [] temp12;
	     delete [] temp13; delete [] temp14; 
		 delete [] temp15; delete [] temp16;
	     delete [] temp17; delete [] temp18; delete [] temp19; delete [] temp20;
	     delete [] temp21; delete [] temp22; delete [] temp23; delete [] temp24;
         delete [] temp25; 
     }
  }

  ndmx[matCount] = nd;
  loadStagex[matCount] = 0;   //default
  refShearModulusx[matCount] = refShearModul;
  refBulkModulusx[matCount] = refBulkModul;
  frictionAnglex[matCount] = frictionAng;
  peakShearStrainx[matCount] = peakShearStra;
  refPressurex[matCount] = -refPress;  //compression is negative
  cohesionx[matCount] = cohesi;
  pressDependCoeffx[matCount] = pressDependCoe;
  numOfSurfacesx[matCount] = numberOfYieldSurf;
  rhox[matCount] = r;
  phaseTransfAnglex[matCount] = phaseTransformAng;

  mTypex[matCount] = mType;
  contractParam1x[matCount] = ca;
  contractParam2x[matCount] = cb;
  contractParam3x[matCount] = cc;
  contractParam4x[matCount] = cd;
  contractParam5x[matCount] = ce;

  dilateParam1x[matCount] = da;
  dilateParam2x[matCount] = db;
  dilateParam3x[matCount] = dc;

  liquefyParam1x[matCount] = liquefactionParam1;
  liquefyParam2x[matCount] = liquefactionParam2;

  einitx[matCount] = ei;    //default initial void ratio
  volLimit1x[matCount] = cs1;
  volLimit2x[matCount] = cs2;
  volLimit3x[matCount] = cs3;

  Hvx[matCount] = hv;
  Pvx[matCount] = pv;
  
  residualPressx[matCount] =0.;
  stressRatioPTx[matCount] =0.;
  
  matN = matCount;
  matCount ++;
  pAtm = atm;

  int numOfSurfaces = numOfSurfacesx[matN];
  initPress = refPressurex[matN];

  e2p = committedActiveSurf = activeSurfaceNum = 0;
  onPPZCommitted = onPPZ = -1 ;
  PPZSizeCommitted = PPZSize = 0.;
  pressureDCommitted = pressureD = modulusFactor = 0.;
  cumuDilateStrainOctaCommitted    = cumuDilateStrainOcta = 0.;
  maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta = 0.;
  cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta = 0.;
  prePPZStrainOctaCommitted = prePPZStrainOcta = 0.;
  oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta = 0.;
  maxPress = 0.;
  damage = 0.;

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];

  mGredu = gredu;
  setUpSurfaces(gredu);  // residualPress and stressRatioPT are calculated inside.
}


PressureDependMultiYield03::PressureDependMultiYield03 ()
 : NDMaterial(0,ND_TAG_PressureDependMultiYield03),
   currentStress(), trialStress(), currentStrain(),
  strainRate(), PPZPivot(), PPZCenter(), PivotStrainRate(6), PivotStrainRateCommitted(6),
  PPZPivotCommitted(), PPZCenterCommitted(), theSurfaces(0), committedSurfaces(0)
{
  //does nothing
}


PressureDependMultiYield03::PressureDependMultiYield03 (const PressureDependMultiYield03 & a)
 : NDMaterial(a.getTag(),ND_TAG_PressureDependMultiYield03),
   currentStress(a.currentStress), trialStress(a.trialStress),
  currentStrain(a.currentStrain), strainRate(a.strainRate), check(0),
  PPZPivot(a.PPZPivot), PPZCenter(a.PPZCenter), updatedTrialStress(a.updatedTrialStress),
  PPZPivotCommitted(a.PPZPivotCommitted), PPZCenterCommitted(a.PPZCenterCommitted),
  PivotStrainRate(a.PivotStrainRate), PivotStrainRateCommitted(a.PivotStrainRateCommitted)
{
  matN = a.matN;

  int numOfSurfaces = numOfSurfacesx[matN];

  e2p = a.e2p;
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
  initPress = a.initPress;
  maxPress = a.maxPress;
  damage = a.damage;

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1];  //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];
  for(int i=1; i<=numOfSurfaces; i++) {
    committedSurfaces[i] = a.committedSurfaces[i];
    theSurfaces[i] = a.theSurfaces[i];
  }
}


PressureDependMultiYield03::~PressureDependMultiYield03 ()
{
  if (theSurfaces != 0) delete [] theSurfaces;
  if (committedSurfaces != 0) delete [] committedSurfaces;
}


void PressureDependMultiYield03::elast2Plast(void)
{
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  if (loadStage != 1 || e2p == 1) 
		return;

  e2p = 1;

  if (currentStress.volume() > 0.) {
    //opserr << "WARNING:PressureDependMultiYield03::elast2Plast(): material in tension." << endln;
    currentStress.setData(currentStress.deviator(),0);
  }

  // Active surface is 0, return
  if (currentStress.deviatorLength() == 0.) return;

//  this->initStrainUpdate();

  // Find active surface
  while (yieldFunc(currentStress, committedSurfaces, ++committedActiveSurf) > 0) {
    if (committedActiveSurf == numOfSurfaces) {
      //opserr <<"WARNING:PressureDependMultiYield03::elast2Plast(): stress out of failure surface"<<endln;
      deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
      initSurfaceUpdate();
      return;
    }
  }

  committedActiveSurf--;
  initSurfaceUpdate();
}


int PressureDependMultiYield03::setTrialStrain (const Vector &strain)
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  if (ndm==3 && strain.Size()==6)
    workV6 = strain;
  else if (ndm==2 && strain.Size()==3) {
    workV6[0] = strain[0];
    workV6[1] = strain[1];
    workV6[2] = 0.0;
    workV6[3] = strain[2];
    workV6[4] = 0.0;
    workV6[5] = 0.0;
  }
  else {
    opserr << "Fatal:PressureDependMultiYield03:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
   exit(-1);
  }

  //strainRate.setData(workV6-currentStrain.t2Vector(1),1);
  workV6 -= currentStrain.t2Vector(1);
  strainRate.setData(workV6, 1);

  return 0;
}


int PressureDependMultiYield03::setTrialStrain (const Vector &strain, const Vector &rate)
{
  return setTrialStrain (strain);
}


int PressureDependMultiYield03::setTrialStrainIncr (const Vector &strain)
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  if (ndm==3 && strain.Size()==6)
    workV6 = strain;
  else if (ndm==2 && strain.Size()==3) {
    workV6[0] = strain[0];
    workV6[1] = strain[1];
    workV6[2] = 0.0;
    workV6[3] = strain[2];
    workV6[4] = 0.0;
    workV6[5] = 0.0;
  }
  else {
    opserr << "Fatal:PressureDependMultiYield03:: Material dimension is: " << ndm << endln;
    opserr << "But strain vector size is: " << strain.Size() << endln;
   exit(-1);
  }

  strainRate.setData(workV6,1);
  return 0;
}


int PressureDependMultiYield03::setTrialStrainIncr (const Vector &strain, const Vector &rate)
{
  return setTrialStrainIncr(strain);
}


const Matrix & PressureDependMultiYield03::getTangent (void)
{
  int loadStage = loadStagex[matN];
  double refShearModulus = refShearModulusx[matN];
  double refBulkModulus = refBulkModulusx[matN];
  double pressDependCoeff = pressDependCoeffx[matN];
  double refPressure = refPressurex[matN];
  double residualPress = residualPressx[matN];
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  if (loadStage == 1 && e2p == 0) {
//	  opserr << "PDMY03::getTang() - 1\n";
      initPress = currentStress.volume();
//	   opserr << "PDMY03::getTang() - 2\n";
      elast2Plast();
//	   opserr << "PDMY03::getTang() - 3\n";
  }
  if (loadStage>=2 && initPress==refPressure) 
	  initPress = currentStress.volume();

  if (loadStage != 1 ) {  
	double factor;
	if (loadStage==0) 
		factor = 1.0;
	else {
		factor = (initPress-residualPress)/(refPressure-residualPress);
		if (factor <= 1.e-10) factor = 1.e-10;
		else factor = pow(factor, pressDependCoeff);
		factor = (1.e-10>factor) ? 1.e-10 : factor;

		// use zero Bulk during compression (J. Lu 7/12/18)
		if (loadStage == 3 && strainRate.volume() <= 0.)
		{
			refBulkModulus = 0.0;
		}
	}
    for (int i=0;i<6;i++)
      for (int j=0;j<6;j++) {
	    theTangent(i,j) = 0.;
        if (i==j) 
			theTangent(i,j) += refShearModulus*factor;
		if (i<3 && j<3 && i==j) 
			theTangent(i,j) += refShearModulus*factor;
		if (i<3 && j<3) 
			theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.)*factor;
	  }
  } else {
    double coeff1, coeff2, coeff3, coeff4;
    double factor = getModulusFactor(updatedTrialStress);
    double shearModulus = factor*refShearModulus;
    double bulkModulus = factor*refBulkModulus;

	// volumetric plasticity
	if (Hvx[matN] != 0. && trialStress.volume()<=maxPress
		&& strainRate.volume()<0. && loadStage == 1) {
	  double tp = fabs(trialStress.volume() - residualPress);
      bulkModulus = (bulkModulus*Hvx[matN]*pow(tp,Pvx[matN]))/(bulkModulus+Hvx[matN]*pow(tp,Pvx[matN]));
	}

    /*if (loadStage!=0 && committedActiveSurf > 0) {
      getSurfaceNormal(updatedTrialStress, workT2V);
      workV6 = workT2V.deviator();
      double volume = workT2V.volume();
      double Ho = 9.*bulkModulus*volume*volume + 2.*shearModulus*(workV6 && workV6);
      double plastModul = factor*committedSurfaces[committedActiveSurf].modulus();
      coeff1 = 9.*bulkModulus*bulkModulus*volume*volume/(Ho+plastModul);
      coeff2 = 4.*shearModulus*shearModulus/(Ho+plastModul);
      // non-symmetric stiffness
      getSurfaceNormal(updatedTrialStress, workT2V);
      workV6 = workT2V.deviator();
      double qq = workT2V.volume();
			double pp=getPlasticPotential(updatedTrialStress, workT2V);
      double Ho = 9.*bulkModulus*pp*qq + 2.*shearModulus*(workV6 && workV6);
      double plastModul = factor*committedSurfaces[committedActiveSurf].modulus();
      coeff1 = 9.*bulkModulus*bulkModulus*pp*qq/(Ho+plastModul);
      coeff2 = 4.*shearModulus*shearModulus/(Ho+plastModul);
			coeff3 = 6.*shearModulus*pp/(Ho+plastModul);
			coeff4 = 6.*shearModulus*qq/(Ho+plastModul);
    }*/
    if (loadStage!=0 && activeSurfaceNum > 0) {
	//	 opserr << "PDMY03::getTang() - 5\n";
      factor = getModulusFactor(trialStress);
      shearModulus = factor*refShearModulus;
      bulkModulus = factor*refBulkModulus;
      getSurfaceNormal(trialStress, workT2V);
      workV6 = workT2V.deviator();
      double volume = workT2V.volume();
      double Ho = 9.*bulkModulus*volume*volume + 2.*shearModulus*(workV6 && workV6);
      double plastModul = factor*theSurfaces[activeSurfaceNum].modulus();
      coeff1 = 9.*bulkModulus*bulkModulus*volume*volume/(Ho+plastModul);
      coeff2 = 4.*shearModulus*shearModulus/(Ho+plastModul);
	  // opserr << "PDMY03::getTang() - 6\n";
	}

    else {
		// opserr << "PDMY03::getTang() - 7\n";
      coeff1 = coeff2 = coeff3 = coeff4 = 0.;
      workV6.Zero();
	   //opserr << "PDMY03::getTang() - 8\n";
    }

    for (int i=0;i<6;i++)
      for (int j=0;j<6;j++) {
	      theTangent(i,j) = - coeff2*workV6[i]*workV6[j];
        if (i==j) theTangent(i,j) += shearModulus;
        if (i<3 && j<3 && i==j) theTangent(i,j) += shearModulus;
	      if (i<3 && j<3) theTangent(i,j) += (bulkModulus - 2.*shearModulus/3. - coeff1);
/* non-symmetric stiffness
				if (i<3) theTangent(i,j) -= coeff3 * workV6[j];
				if (j<3) theTangent(i,j) -= coeff4 * workV6[i];*/
      }
  }
 //  opserr << "PDMY03::getTang() - 9\n";
  if (ndm==3)
    return theTangent;
  else {
    static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = 0.;

    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = 0.;

    workM(2,0) = 0.;
    workM(2,1) = 0.;
    workM(2,2) = theTangent(3,3);

    /* non-symmetric stiffness
       workM(0,2) = theTangent(0,3);
       workM(1,2) = theTangent(1,3);
       workM(2,0) = theTangent(3,0);
       workM(2,1) = theTangent(3,1);*/
    return workM;
  }
   //opserr << "PDMY03::getTang() - DONE\n";
}


const Matrix & PressureDependMultiYield03::getInitialTangent (void)
{
  int loadStage = loadStagex[matN];
  double refShearModulus = refShearModulusx[matN];
  double refBulkModulus = refBulkModulusx[matN];
  double pressDependCoeff = pressDependCoeffx[matN];
  double refPressure = refPressurex[matN];
  double residualPress = residualPressx[matN];
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  if (loadStage == 1 && e2p == 0) {
      initPress = currentStress.volume();
      elast2Plast();
  }
  if (loadStage>=2 && initPress==refPressure) 
	  initPress = currentStress.volume();
  double factor;
  if (loadStage==0)
	  factor = 1.;
  else if (loadStage>=2) { 
		factor = (initPress-residualPress)/(refPressure-residualPress);
		if (factor <= 1.e-10) factor = 1.e-10;
		else factor = pow(factor, pressDependCoeff);
		factor = (1.e-10>factor) ? 1.e-10 : factor;

		// use zero Bulk during compression (J. Lu 7/12/18)
		if (loadStage == 3 && strainRate.volume() <= 0.)
		{
			refBulkModulus = 0.0;
		}
  }
  else if (loadStage==1)
	  factor = getModulusFactor(currentStress);

  for (int i=0;i<6;i++)
    for (int j=0;j<6;j++) {
      theTangent(i,j) = 0.;
      if (i==j) theTangent(i,j) += refShearModulus*factor;
      if (i<3 && j<3 && i==j) theTangent(i,j) += refShearModulus*factor;
      if (i<3 && j<3) theTangent(i,j) += (refBulkModulus - 2.*refShearModulus/3.)*factor;
    }

  if (ndm==3)
    return theTangent;
  else {
    static Matrix workM(3,3);
    workM(0,0) = theTangent(0,0);
    workM(0,1) = theTangent(0,1);
    workM(0,2) = 0.;

    workM(1,0) = theTangent(1,0);
    workM(1,1) = theTangent(1,1);
    workM(1,2) = 0.;

    workM(2,0) = 0.;
    workM(2,1) = 0.;
    workM(2,2) = theTangent(3,3);

/* non-symmetric stiffness
    workM(0,2) = theTangent(0,3);
    workM(1,2) = theTangent(1,3);
    workM(2,0) = theTangent(3,0);
    workM(2,1) = theTangent(3,1);*/
    return workM;
  }
}


const Vector & PressureDependMultiYield03::getStress (void)
{
//	opserr << "PDMY03-getStress() -1\n";
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 3;

  int i, is;
  if (loadStage == 1 && e2p == 0) {
      initPress = currentStress.volume();
      elast2Plast();
  }

  if (loadStage!=1) {  //linear elastic
    //trialStrain.setData(currentStrain.t2Vector() + strainRate.t2Vector());
    getTangent();
    workV6 = currentStress.t2Vector();
	workV6.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
    trialStress.setData(workV6);
  }
  else {
    for (i=1; i<=numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
    activeSurfaceNum = committedActiveSurf;
    pressureD = pressureDCommitted;
    onPPZ = onPPZCommitted;
    PPZSize = PPZSizeCommitted;
    cumuDilateStrainOcta = cumuDilateStrainOctaCommitted;
    maxCumuDilateStrainOcta = maxCumuDilateStrainOctaCommitted;
    cumuTranslateStrainOcta = cumuTranslateStrainOctaCommitted;
    prePPZStrainOcta = prePPZStrainOctaCommitted;
    oppoPrePPZStrainOcta = oppoPrePPZStrainOctaCommitted;
    PPZPivot = PPZPivotCommitted;
	PivotStrainRate = PivotStrainRateCommitted;
    PPZCenter = PPZCenterCommitted;

    subStrainRate = strainRate;
    setTrialStress(currentStress);
    if (activeSurfaceNum>0 && isLoadReversal(currentStress)) {
      updateInnerSurface();
      activeSurfaceNum = 0;
    }

    if (activeSurfaceNum==0 && !isCrossingNextSurface()) {
 	    workV6 = currentStrain.t2Vector();
	    workV6.addVector(1.0, strainRate.t2Vector(), 1.0);
	    trialStrain.setData(workV6);
	}
	else {
      int numSubIncre = setSubStrainRate();

      for (i=0; i<numSubIncre; i++) {
//      trialStrain.setData(currentStrain.t2Vector()
//			     + subStrainRate.t2Vector()*(i+1));
 	      workV6 = currentStrain.t2Vector();
	      workV6.addVector(1.0, subStrainRate.t2Vector(), (i+1));
	      trialStrain.setData(workV6);

		if (i==0)  {
			  updatedTrialStress = currentStress;
			  setTrialStress(currentStress);
              is = isLoadReversal(currentStress);
		}
		else {
            updatedTrialStress = trialStress;
			workT2V.setData(trialStress.t2Vector());
			setTrialStress(trialStress);
            is = isLoadReversal(workT2V);
		}

        if (activeSurfaceNum>0 && is) {
          updateInnerSurface();
          activeSurfaceNum = 0;
		}
        if (activeSurfaceNum==0 && !isCrossingNextSurface()) continue;
        if (activeSurfaceNum==0) activeSurfaceNum++;
        stressCorrection(0);
        updateActiveSurface();

	    double refBulkModulus = refBulkModulusx[matN];
		//modulusFactor was calculated in setTrialStress
        double B = refBulkModulus*modulusFactor;
		//double deltaD = 3.*subStrainRate.volume()
		//	         - (trialStress.volume()-updatedTrialStress.volume())/B;
		//if (deltaD<0) deltaD /=2 ;
        //pressureD += deltaD;
		pressureD += 3.*subStrainRate.volume()
			         - (trialStress.volume()-updatedTrialStress.volume())/B;
		if (pressureD < 0.) pressureD = 0.;
		//opserr<<i<<" "<<activeSurfaceNum<<" "<<is<<" "<<subStrainRate.t2Vector()[3]<<endln;
	  }
	}
  }
 //  opserr << "PDMY03::getStress() - DONE\n";
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


const Vector & PressureDependMultiYield03::getStrain (void)
{
  return getCommittedStrain();
}


int PressureDependMultiYield03::commitState (void)
{
  int loadStage = loadStagex[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  currentStress = trialStress;
  //currentStrain = T2Vector(currentStrain.t2Vector() + strainRate.t2Vector());
  workV6 = currentStrain.t2Vector();
  workV6 += strainRate.t2Vector();
  currentStrain.setData(workV6);

  workV6.Zero();
  strainRate.setData(workV6);

  if (loadStage==1) {
    committedActiveSurf = activeSurfaceNum;
    for (int i=1; i<=numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
    pressureDCommitted = pressureD;
    onPPZCommitted = onPPZ;
    PPZSizeCommitted = PPZSize;
    cumuDilateStrainOctaCommitted = cumuDilateStrainOcta;
    maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta;
    cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta;
    prePPZStrainOctaCommitted = prePPZStrainOcta;
    oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta;
    PPZPivotCommitted = PPZPivot;
	PivotStrainRateCommitted = PivotStrainRate;
    PPZCenterCommitted = PPZCenter;
	if (currentStress.volume() < maxPress) maxPress = currentStress.volume();
  }

  return 0;
}


int PressureDependMultiYield03::revertToLastCommit (void)
{
  return 0;
}


NDMaterial * PressureDependMultiYield03::getCopy (void)
{
  PressureDependMultiYield03 * copy = new PressureDependMultiYield03(*this);
  return copy;
}


NDMaterial * PressureDependMultiYield03::getCopy (const char *code)
{
  if (strcmp(code,"PressureDependMultiYield03") == 0 || strcmp(code,"PlaneStrain") == 0
      || strcmp(code,"ThreeDimensional") == 0) {
    PressureDependMultiYield03 * copy = new PressureDependMultiYield03(*this);
    return copy;
  }

  return 0;
}


const char * PressureDependMultiYield03::getType (void) const
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int PressureDependMultiYield03::getOrder (void) const
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  return (ndm == 2) ? 3 : 6;
}


int PressureDependMultiYield03::setParameter(const char **argv, int argc, Parameter &param)
{
  /*if (argc < 1)
    return -1;

  if (strcmp(argv[0],"updateMaterialStage") == 0) {
    if (argc < 2)
      return -1;
    int matTag = atoi(argv[1]);
    if (this->getTag() == matTag)
      return param.addObject(1, this);
    else
      return -1;
  }

  else if (strcmp(argv[0],"shearModulus") == 0)
    return param.addObject(10, this);

  else if (strcmp(argv[0],"bulkModulus") == 0)
    return param.addObject(11, this);

  return -1;*/

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

int PressureDependMultiYield03::updateParameter(int responseID, Information &info)
{
 
  if (responseID == 1) {
      loadStagex[matN] = info.theInt;
  } else if (responseID==10) {
    refShearModulusx[matN] = info.theDouble;
  } else if (responseID==11) {
    refBulkModulusx[matN] = info.theDouble;
  } else if (responseID==12) {
    frictionAnglex[matN] = info.theDouble;
    setUpSurfaces(mGredu);
    initSurfaceUpdate();
  } else if (responseID==13) {
    cohesionx[matN] = info.theDouble;
    setUpSurfaces(mGredu);
    initSurfaceUpdate();
  }

  // used by BBarFourNodeQuadUP element
  else if (responseID==20 && ndmx[matN] == 2)
		ndmx[matN] = 0;

  return 0;
}


int PressureDependMultiYield03::sendSelf(int commitTag, Channel &theChannel)
{
    int loadStage = loadStagex[matN];
    int ndm = ndmx[matN];
	double rho = rhox[matN];
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];
    double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];
    double phaseTransfAngle = phaseTransfAnglex[matN];
	double stressRatioPT = stressRatioPTx[matN];
	double contractParam1 = contractParam1x[matN];
	double contractParam2 = contractParam2x[matN];
    double dilateParam1 = dilateParam1x[matN];
    double dilateParam2 = dilateParam2x[matN];
	double liquefyParam1 = liquefyParam1x[matN];
	double liquefyParam2 = liquefyParam2x[matN];
	double dilateParam3 = dilateParam3x[matN];
	double einit = einitx[matN];
	double volLimit1 = volLimit1x[matN];
	double volLimit2 = volLimit2x[matN];
	double volLimit3 = volLimit3x[matN];

     double contractionParam3 = contractParam3x[matN];
     double hv = Hvx[matN];
     double Pv = Pvx[matN];

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
    opserr << "PressureDependMultiYield03::sendSelf -- could not send ID\n";
    return res;
  }

  Vector data(72+numOfSurfaces*8);
  data(0) = rho;
  data(1) = einit;
  data(2) = refShearModulus;
  data(3) = refBulkModulus;
  data(4) = frictionAngle;
  data(5) = peakShearStrain;
  data(6) = refPressure;
  data(7) = cohesion;
  data(8) = pressDependCoeff;
  data(9) = phaseTransfAngle;
  data(10) = contractParam1;
  data(11) = dilateParam1;
  data(12) = dilateParam2;
  data(13) = volLimit1;
  data(14) = volLimit2;
  data(15) = volLimit3;
  data(16) = pAtm;
  data(17) = liquefyParam1;
  data(18) = liquefyParam2;
  data(19) = dilateParam3;
  data(20) = residualPress;
  data(21) = stressRatioPT;
  data(22) = e2p;
  data(23) = committedActiveSurf;
  data(24) = strainPTOcta;
  data(25) = pressureDCommitted;
  data(26) = onPPZCommitted;
  data(27) = PPZSizeCommitted;
  data(28) = cumuDilateStrainOctaCommitted;
  data(29) = maxCumuDilateStrainOctaCommitted;
  data(30) = cumuTranslateStrainOctaCommitted;
  data(31) = prePPZStrainOctaCommitted;
  data(32) = oppoPrePPZStrainOctaCommitted;

  data(33) = initPress;
  data(34) = contractParam2;
  data(35) = contractionParam3;// = contractParam3x[matN];
  data(36) =  hv; //  = Hvx[matN];
  data(37) = Pv; //  = Pvx[matN];

  workV6 = currentStress.t2Vector();
  for(i = 0; i < 6; i++) data(i+38) = workV6[i];

  workV6 = currentStrain.t2Vector();
  for(i = 0; i < 6; i++) data(i+44) = workV6[i];

  workV6 = PPZPivotCommitted.t2Vector();
  for(i = 0; i < 6; i++) data(i+50) = workV6[i];

  workV6 = PPZCenterCommitted.t2Vector();
  for(i = 0; i < 6; i++) data(i+56) = workV6[i];

  for(i = 0; i < numOfSurfaces; i++) {
    int k = 62 + i*8;
    data(k) = committedSurfaces[i+1].size();
    data(k+1) = committedSurfaces[i+1].modulus();
    workV6 = committedSurfaces[i+1].center();
    data(k+2) = workV6(0);
    data(k+3) = workV6(1);
    data(k+4) = workV6(2);
    data(k+5) = workV6(3);
    data(k+6) = workV6(4);
    data(k+7) = workV6(5);
  }

  i = 69 + numOfSurfaces * 8;
  data(i+1) = mTypex[matN];
  data(i+2) = contractParam4x[matN];
  data(i+3) = contractParam5x[matN];

  res += theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureDependMultiYield03::sendSelf -- could not send Vector\n";
    return res;
  }

  return res;
}


int PressureDependMultiYield03::recvSelf(int commitTag, Channel &theChannel,
				       FEM_ObjectBroker &theBroker)
{
  int i, res = 0;

  static ID idData(6);
  res += theChannel.recvID(this->getDbTag(), commitTag, idData);
  if (res < 0) {
    opserr << "PressureDependMultiYield03::recvelf -- could not recv ID\n";

    return res;
  }

  this->setTag(idData(0));
  int numOfSurfaces = idData(1);
  int loadStage = idData(2);
  int ndm = idData(3);
  matN = idData(4);

  int otherMatCount = idData(5);

  Vector data(72+idData(1)*8);
  res += theChannel.recvVector(this->getDbTag(), commitTag, data);
  if (res < 0) {
    opserr << "PressureDependMultiYield03::recvSelf -- could not recv Vector\n";
    return res;
  }

  double rho = data(0);
  double einit = data(1);
  double refShearModulus = data(2);
  double refBulkModulus = data(3);
  double frictionAngle = data(4);
  double peakShearStrain = data(5);
  double refPressure = data(6);
  double cohesion = data(7);
  double pressDependCoeff = data(8);
  double phaseTransfAngle = data(9);
  double contractParam1 = data(10);
  double dilateParam1 = data(11);
  double dilateParam2 = data(12);
  double volLimit1 = data(13);
  double volLimit2 = data(14);
  double volLimit3 = data(15);
  pAtm = data(16);
  double liquefyParam1 = data(17);
  double liquefyParam2 = data(18);
  double dilateParam3 = data(19);
  double residualPress = data(20);
  double stressRatioPT = data(21);
  e2p = data(22);
  committedActiveSurf = data(23);
  strainPTOcta = data(24);
  pressureDCommitted = data(25);
  onPPZCommitted = data(26);
  PPZSizeCommitted = data(27);
  cumuDilateStrainOctaCommitted = data(28);
  maxCumuDilateStrainOctaCommitted = data(29);
  cumuTranslateStrainOctaCommitted = data(30);
  prePPZStrainOctaCommitted = data(31);
  oppoPrePPZStrainOctaCommitted = data(32);

  initPress = data(33);
  double contractParam2 = data(34);
  double contractParam3 =  data(35); //  = contractionParam3;// = contractParam3x[matN];
  double hv = data(36); // =  hv; //  = Hvx[matN];
  double Pv = data(37); // = Pv; //  = Pvx[matN];


  for(i = 0; i < 6; i++) workV6[i] = data(i+38);
  currentStress.setData(workV6);

  for(i = 0; i < 6; i++) workV6[i] = data(i+44);
  currentStrain.setData(workV6);

  for(i = 0; i < 6; i++) workV6[i] = data(i+50);
  PPZPivotCommitted.setData(workV6);

  for(i = 0; i < 6; i++) workV6[i] = data(i+56);
  PPZCenterCommitted.setData(workV6);

  if (committedSurfaces != 0) {
      delete [] committedSurfaces;
      delete [] theSurfaces;
  }

  theSurfaces = new MultiYieldSurface[numOfSurfaces+1]; //first surface not used
  committedSurfaces = new MultiYieldSurface[numOfSurfaces+1];

  for(i = 0; i < numOfSurfaces; i++) {
    int k = 62 + i*8;
    workV6(0) = data(k+2);
    workV6(1) = data(k+3);
    workV6(2) = data(k+4);
    workV6(3) = data(k+5);
    workV6(4) = data(k+6);
    workV6(5) = data(k+7);
    committedSurfaces[i+1].setData(workV6, data(k), data(k+1));
  }

  i = 69 + numOfSurfaces * 8;
  int mType = data(i+1);
  double contractParam4 = data(i+2);
  double contractParam5 = data(i+3);
  
  
  int *temp1, *temp2, *temp11;
  double *temp3, *temp4, *temp5, *temp6, *temp7, *temp8, *temp9, *temp10, *temp12;
  double *temp13, *temp14, *temp15, *temp16, *temp17, *temp18, *temp19, *temp20;
  double *temp14a, *temp14b;
  double *temp21, *temp22, *temp23, *temp24, *temp25, *temp26;
  int *temp27;
  double *temp28, *temp29, *temp30, *temp31;

  if (matCount < otherMatCount) {  // allocate memory if not enough

     temp1 = loadStagex;
	 temp2 = ndmx;
	 temp3 = rhox;
     temp4 = refShearModulusx;
     temp5 = refBulkModulusx;
     temp6 = frictionAnglex;
     temp7 = peakShearStrainx;
     temp8 = refPressurex;
     temp9 = cohesionx;
     temp10 = pressDependCoeffx;
	 temp11 = numOfSurfacesx;
     temp12 = residualPressx;
     temp13 = phaseTransfAnglex;
     temp14 = contractParam1x;
     temp14a = contractParam2x;
	 temp14b = contractParam3x;
     temp15 = dilateParam1x;
     temp16 = dilateParam2x;
     temp17 = liquefyParam1x;
     temp18 = liquefyParam2x;
     temp19 = dilateParam3x;
     temp20 = einitx;    //initial void ratio
     temp21 = volLimit1x;
     temp22 = volLimit2x;
     temp23 = volLimit3x;
     temp24 = stressRatioPTx;
	 temp25 = Hvx;
	 temp26 = Pvx;

	 temp27 = mTypex;
	 temp28 = contractParam4x;
	 temp29 = contractParam5x;
	 
     loadStagex = new int[otherMatCount];
     ndmx = new int[otherMatCount];
     rhox = new double[otherMatCount];
     refShearModulusx = new double[otherMatCount];
     refBulkModulusx = new double[otherMatCount];
     frictionAnglex = new double[otherMatCount];
     peakShearStrainx = new double[otherMatCount];
     refPressurex = new double[otherMatCount];
	 cohesionx = new double[otherMatCount];
     pressDependCoeffx = new double[otherMatCount];
     numOfSurfacesx = new int[otherMatCount];
     residualPressx = new double[otherMatCount];
     phaseTransfAnglex = new double[otherMatCount];
     contractParam1x = new double[otherMatCount];
     contractParam2x = new double[otherMatCount];
	 contractParam3x = new double[otherMatCount];
     dilateParam1x = new double[otherMatCount];
     dilateParam2x = new double[otherMatCount];
     liquefyParam1x = new double[otherMatCount];
     liquefyParam2x = new double[otherMatCount];
     dilateParam3x = new double[otherMatCount];
     einitx = new double[otherMatCount];    //initial void ratio
     volLimit1x = new double[otherMatCount];
     volLimit2x = new double[otherMatCount];
     volLimit3x = new double[otherMatCount];
     stressRatioPTx = new double[otherMatCount];
	 Hvx = new double[otherMatCount];
	 Pvx = new double[otherMatCount];

	 mTypex = new int[otherMatCount];
	 contractParam4x = new double[otherMatCount];
	 contractParam5x = new double[otherMatCount];
	 
     if( matCount > 0 ) {
		 for (int i=0; i<matCount; i++) {
			 loadStagex[i] = temp1[i];
			 ndmx[i] = temp2[i];
			 rhox[i] = temp3[i];
			 refShearModulusx[i] = temp4[i];
			 refBulkModulusx[i] = temp5[i];
			 frictionAnglex[i] = temp6[i];
			 peakShearStrainx[i] = temp7[i];
			 refPressurex[i] = temp8[i];
			 cohesionx[i] = temp9[i];
			 pressDependCoeffx[i] = temp10[i];
			 numOfSurfacesx[i] = temp11[i];
			 residualPressx[i] = temp12[i];
			 phaseTransfAnglex[i] = temp13[i];
			 contractParam1x[i] = temp14[i];
	         contractParam2x[i] = temp14a[i];
			 contractParam3x[i] = temp14b[i];
			 dilateParam1x[i] = temp15[i];
			 dilateParam2x[i] = temp16[i];
			 liquefyParam1x[i] = temp17[i];
			 liquefyParam2x[i] = temp18[i];
	         dilateParam3x[i] = temp19[i];
			 einitx[i] = temp20[i];    //initial void ratio
			 volLimit1x[i] = temp21[i];
			 volLimit2x[i] = temp22[i];
			 volLimit3x[i] = temp23[i];
			 stressRatioPTx[i] = temp24[i];
			 Hvx[i] = temp25[i];
			 Pvx[i] = temp26[i];
			 
			 mTypex[i] = temp27[i];
			 contractParam4x[i] = temp28[i];
			 contractParam5x[i] = temp29[i];
			 
		 }

	     delete [] temp1; delete [] temp2; delete [] temp3; delete [] temp4;
	     delete [] temp5; delete [] temp6; delete [] temp7; delete [] temp8;
	     delete [] temp9; delete [] temp10; delete [] temp11; delete [] temp12;
	     delete [] temp13; delete [] temp14; delete [] temp15; delete [] temp16;
	     delete [] temp14a;
	     delete [] temp17; delete [] temp18; delete [] temp19; delete [] temp20;
	     delete [] temp21; delete [] temp22; delete [] temp23; delete [] temp24;
         delete [] temp25; delete [] temp26;
		 delete[] temp27; delete[] temp28; delete[] temp29; delete[] temp30; delete[] temp31;
     }
	 matCount = otherMatCount;
  }

    loadStagex[matN] = loadStage;
    ndmx[matN] = ndm;
	rhox[matN] = rho;
    residualPressx[matN] = residualPress;
    numOfSurfacesx[matN] = numOfSurfaces;
    refPressurex[matN] = refPressure;
    pressDependCoeffx[matN] =pressDependCoeff;
    refShearModulusx[matN] = refShearModulus;
	refBulkModulusx[matN] = refBulkModulus;
    frictionAnglex[matN] = frictionAngle;
	cohesionx[matN] = cohesion;
    peakShearStrainx[matN] = peakShearStrain;
    phaseTransfAnglex[matN] = phaseTransfAngle;
	stressRatioPTx[matN] = stressRatioPT;
	contractParam1x[matN] = contractParam1;
	contractParam2x[matN] = contractParam2;
	contractParam3x[matN] = contractParam3;
    dilateParam1x[matN] = dilateParam1;
    dilateParam2x[matN] = dilateParam2;
	liquefyParam1x[matN] = liquefyParam1;
	liquefyParam2x[matN] = liquefyParam2;
	dilateParam3x[matN] = dilateParam3;
	einitx[matN] = einit;
	volLimit1x[matN] = volLimit1;
	volLimit2x[matN] = volLimit2;
	volLimit3x[matN] = volLimit3;

	mTypex[matN] = mType;
	contractParam4x[matN] = contractParam4;
	contractParam5x[matN] = contractParam5;

  return res;
}


Response*
PressureDependMultiYield03::setResponse (const char **argv, int argc, OPS_Stream &s)
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
	    Matrix curv(numOfSurfaces+1,(argc-1)*2);
		for (int i=1; i<argc; i++)
			curv(0,(i-1)*2) = atoi(argv[i]);
		return new MaterialResponse(this, 4, curv);
  }
	else
		return 0;
}


void PressureDependMultiYield03::getBackbone (Matrix & bb)
{
  double residualPress = residualPressx[matN];
  double refPressure = refPressurex[matN];
  double pressDependCoeff =pressDependCoeffx[matN];
  double refShearModulus = refShearModulusx[matN];
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
				stress2 = committedSurfaces[i].size()*conHeig/sqrt(3.0);
				strain2 = stress2/shearModulus;
				bb(1,k*2) = strain2; bb(1,k*2+1) = shearModulus;
			} else {
				stress1 = stress2; strain1 = strain2;
				plastModulus = factor*committedSurfaces[i-1].modulus();
				elast_plast = 2*shearModulus*plastModulus/(2*shearModulus+plastModulus);
				stress2 = committedSurfaces[i].size()*conHeig/sqrt(3.0);
			    strain2 = 2*(stress2-stress1)/elast_plast + strain1;
				gre = stress2/strain2;
                bb(i,k*2) = strain2; bb(i,k*2+1) = gre;
			}
		}
	}

}


int PressureDependMultiYield03::getResponse (int responseID, Information &matInfo)
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


void PressureDependMultiYield03::Print(OPS_Stream &s, int flag )
{
  s << "PressureDependMultiYield03" << endln;
}


const Vector & PressureDependMultiYield03::getCommittedStress (void)
{
	int ndm = ndmx[matN];
    if (ndmx[matN] == 0) ndm = 2;
	int numOfSurfaces = numOfSurfacesx[matN];
    double residualPress = residualPressx[matN];

	double scale = currentStress.deviatorRatio(residualPress)/committedSurfaces[numOfSurfaces].size();
	if (loadStagex[matN] != 1) scale = 0.;
  if (ndm==3) {
		static Vector temp7(7);
		workV6 = currentStress.t2Vector();
    temp7[0] = workV6[0];
    temp7[1] = workV6[1];
    temp7[2] = workV6[2];
    temp7[3] = workV6[3];
    temp7[4] = workV6[4];
    temp7[5] = workV6[5];
    temp7[6] = scale;
    /*temp7[7] = committedActiveSurf;
	temp7[8] = stressRatioPTx[matN];
	temp7[9] = currentStress.deviatorRatio(residualPressx[matN]);
    temp7[10] = pressureDCommitted;
    temp7[11] = cumuDilateStrainOctaCommitted;
    temp7[12] = maxCumuDilateStrainOctaCommitted;
    temp7[13] = cumuTranslateStrainOctaCommitted;
    temp7[14] = onPPZCommitted;
    temp7[15] = PPZSizeCommitted;*/
		return temp7;
	}

  else {
    static Vector temp5(5);
	workV6 = currentStress.t2Vector();
    temp5[0] = workV6[0];
    temp5[1] = workV6[1];
    temp5[2] = workV6[2];
    temp5[3] = workV6[3];
    temp5[4] = scale;
    /*temp5[5] = committedActiveSurf;
	temp5[6] = PPZCenterCommitted.deviator()[3];
	temp5[7] = PPZPivotCommitted.deviator()[3];
    temp5[8] = pressureDCommitted;
    temp5[9] = cumuDilateStrainOctaCommitted;
    temp5[10] = maxCumuDilateStrainOctaCommitted;
    temp5[11] = cumuTranslateStrainOctaCommitted;
    temp5[12] = onPPZCommitted;
    temp5[13] = PPZSizeCommitted;
	temp5[14] = PivotStrainRateCommitted[3];
    temp5[15] = initPress;*/
    return temp5;
  }
}

// begin change by Alborz Ghofrani - UW --- get 6 components of stress
const Vector &
PressureDependMultiYield03::getStressToRecord (int numOutput)
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

const Vector & PressureDependMultiYield03::getCommittedStrain (void)
{
  int ndm = ndmx[matN];
  if (ndmx[matN] == 0) ndm = 2;

  if (ndm==3)
    return currentStrain.t2Vector(1);
  else {
		static Vector workV(3);
		workV6 = currentStrain.t2Vector(1);
    workV[0] = workV6[0];
    workV[1] = workV6[1];
    workV[2] = workV6[3];
    return workV;
  }
}


// NOTE: surfaces[0] is not used
void PressureDependMultiYield03::setUpSurfaces (double * gredu)
{
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];
    double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
    double peakShearStrain = peakShearStrainx[matN];
    double phaseTransfAngle = phaseTransfAnglex[matN];
	double stressRatioPT = stressRatioPTx[matN];
	int mType = mTypex[matN];

    double refStrain, peakShear, coneHeight;
    double stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
    double ratio1, ratio2;

	if (gredu==0) {
	  double sinPhi = sin(frictionAngle * pi/180.);
	double Mnys = 6.*sinPhi/(3.-sinPhi); // Triaxial Compression (default)
    double sinPhiPT = sin(phaseTransfAngle * pi/180.);
    stressRatioPT = 6.*sinPhiPT/(3.-sinPhiPT);
	
	if (mType == 1) {
		// Triaxial Extension
		//opserr << "Triaxial Extension" << endln;
		Mnys = 6.*sinPhi/(3.+sinPhi);
		stressRatioPT = 6.*sinPhiPT/(3.+sinPhiPT);
	}
	else if (mType == 2) {
		// Direct Shear
		//opserr << "Direct Shear" << endln;
		Mnys = sqrt(3.)*sinPhi;
		stressRatioPT = sqrt(3.)*sinPhiPT;
	}
	else 
	{
		// Triaxial Compression
		//opserr << "Triaxial Compression" << endln;
	}	
    
	// tao = cohesion * sqrt(8)/3.
    residualPress = 2 * cohesion / Mnys;
    // a small nonzero residualPress for numerical purpose only
    if (residualPress < 0.0001*pAtm) residualPress = 0.0001*pAtm;
    coneHeight = - (refPressure - residualPress);
    peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
    refStrain = (peakShearStrain * peakShear)
		  	        / (refShearModulus * peakShearStrain - peakShear);

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
      workV6.Zero();
	  //opserr<<ii<<" "<<size<<" "<<plast_modul<<endln;
      committedSurfaces[ii] = MultiYieldSurface(workV6,size,plast_modul);
		}  // ii
	}
	else {  //user defined surfaces
		int ii = 2*(numOfSurfaces-1);
		double tmax = refShearModulus*gredu[ii]*gredu[ii+1];
		double Mnys = -(sqrt(3.) * tmax - 2.* cohesion) / refPressure;
        residualPress = 2 * cohesion / Mnys;
        if (residualPress < 0.0001*pAtm) residualPress = 0.0001*pAtm;
        coneHeight = - (refPressure - residualPress);

        double sinPhi = 3*Mnys /(6+Mnys);
		if (sinPhi<0. || sinPhi>1.) {
			opserr <<"\nNDMaterial " <<this->getTag()<<": Invalid friction angle, please modify ref. pressure or G/Gmax curve."<<endln;
            exit(-1);
		}

		frictionAngle = asin(sinPhi)*180/pi;
		opserr << "\nNDMaterial " <<this->getTag()<<": Friction angle is "<<frictionAngle<<"\n"<<endln;
        if (phaseTransfAngle > frictionAngle) {
			opserr << "\nNDMaterial " <<this->getTag()<<": phase Transformation Angle > friction Angle,"
				   << "will set phase Transformation Angle = friction Angle.\n" <<endln;
			phaseTransfAngle = frictionAngle;
		}
		double sinPhiPT = sin(phaseTransfAngle * pi/180.);
        stressRatioPT = 6.*sinPhiPT/(3.-sinPhiPT);

		for (int i=1; i<numOfSurfaces; i++) {
			int ii = 2*(i-1);
			strain1 = gredu[ii];
            stress1 = refShearModulus*gredu[ii+1]*strain1;
			strain2 = gredu[ii+2];
            stress2 = refShearModulus*gredu[ii+3]*strain2;

            ratio1 = sqrt(3.) * stress1 / coneHeight;
            ratio2 = sqrt(3.) * stress2 / coneHeight;
            if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
               double ratio = (ratio2 - stressRatioPT)/(ratio2 - ratio1);
			  // gamma_oct = sqrt(6)/3*gamma12
              strainPTOcta = sqrt(6.)/3 * (strain2 - ratio * (strain2 - strain1));
			}

            size = ratio1;
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

      workV6.Zero();
			//opserr<<size<<" "<<i<<" "<<plast_modul<<" "<<gredu[ii]<<" "<<gredu[ii+1]<<endln;
      committedSurfaces[i] = MultiYieldSurface(workV6,size,plast_modul);

	  if (i==(numOfSurfaces-1)) {
		plast_modul = 0;
		size = ratio2;
		//opserr<<size<<" "<<i+1<<" "<<plast_modul<<" "<<gredu[ii+2]<<" "<<gredu[ii+3]<<endln;
        committedSurfaces[i+1] = MultiYieldSurface(workV6,size,plast_modul);
	  }
	}
  }

  residualPressx[matN] = residualPress;
  frictionAnglex[matN] = frictionAngle;
  cohesionx[matN] = cohesion;
  phaseTransfAnglex[matN] = phaseTransfAngle;
  stressRatioPTx[matN] = stressRatioPT;
}


double PressureDependMultiYield03::yieldFunc(const T2Vector & stress,
					   const MultiYieldSurface * surfaces,
					   int surfaceNum)
{
  double residualPress = residualPressx[matN];

  double coneHeight = stress.volume() - residualPress;
  //workV6 = stress.deviator() - surfaces[surfaceNum].center()*coneHeight;
  workV6 = stress.deviator();
  workV6.addVector(1.0, surfaces[surfaceNum].center(), -coneHeight);

  double sz = surfaces[surfaceNum].size()*coneHeight;

  return 3./2.*(workV6 && workV6) - sz * sz;
}


void PressureDependMultiYield03::deviatorScaling(T2Vector & stress,
					       const MultiYieldSurface * surfaces,
					       int surfaceNum)
{
  double residualPress = residualPressx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  double diff = yieldFunc(stress, surfaces, surfaceNum);
  double coneHeight = stress.volume() - residualPress;

  if ( surfaceNum < numOfSurfaces && diff < 0. ) {
    double sz = -surfaces[surfaceNum].size()*coneHeight;
    double deviaSz = sqrt(sz*sz + diff);
    static Vector devia(6);
    devia = stress.deviator();
    workV6 = devia;
    workV6.addVector(1.0, surfaces[surfaceNum].center(), -coneHeight);
    double coeff = (sz-deviaSz) / deviaSz;
    if (coeff < 1.e-13) coeff = 1.e-13;
    //devia += workV6 * coeff;
    devia.addVector(1.0, workV6, coeff);
    stress.setData(devia, stress.volume());
    deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
  }

  if (surfaceNum==numOfSurfaces && fabs(diff) > LOW_LIMIT) {
    double sz = -surfaces[surfaceNum].size()*coneHeight;
    //workV6 = stress.deviator() * sz/sqrt(diff+sz*sz);
    workV6 = stress.deviator();
    workV6 *= sz/sqrt(diff+sz*sz);
    stress.setData(workV6, stress.volume());
  }
}


void PressureDependMultiYield03::initSurfaceUpdate(void)
{
  double residualPress = residualPressx[matN];
  int numOfSurfaces = numOfSurfacesx[matN];

  if (committedActiveSurf == 0) return;

  double coneHeight = - (currentStress.volume() - residualPress);
  static Vector devia(6);
  devia = currentStress.deviator();
  double Ms = sqrt(3./2.*(devia && devia));

  if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
    //workV6 = devia * (1. - committedSurfaces[committedActiveSurf].size()*coneHeight / Ms);
    workV6.addVector(0.0, devia, (1. - committedSurfaces[committedActiveSurf].size()*coneHeight / Ms));

    //workV6 = workV6 / -coneHeight;
    workV6 /= -coneHeight;
    committedSurfaces[committedActiveSurf].setCenter(workV6);
  }

  for (int i=1; i<committedActiveSurf; i++) {
    //workV6 = devia * (1. - committedSurfaces[i].size()*coneHeight / Ms);
    workV6.addVector(0.0, devia, (1. - committedSurfaces[i].size()*coneHeight / Ms));
    //workV6 = workV6 / -coneHeight;
    workV6 /= -coneHeight;
    committedSurfaces[i].setCenter(workV6);
    theSurfaces[i] = committedSurfaces[i];
  }
  activeSurfaceNum = committedActiveSurf;
}


void PressureDependMultiYield03::initStrainUpdate(void)
{
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];
    double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];
    double stressRatioPT = stressRatioPTx[matN];

  // elastic strain state
  double stressRatio = currentStress.deviatorRatio(residualPress);
  double ratio = (-currentStress.volume()+residualPress)/(-refPressure+residualPress);
  ratio = pow(ratio, 1.-pressDependCoeff);
  modulusFactor = getModulusFactor(currentStress);
  double shearCoeff = 1./(2.*refShearModulus*modulusFactor);
  double bulkCoeff = 1./(3.*refBulkModulus*modulusFactor);

  //currentStrain = currentStress.deviator()*shearCoeff
  //              + currentStress.volume()*bulkCoeff;

  // modified fmk as discussed with z.yang
  workV6.addVector(0.0, currentStress.deviator(), shearCoeff);
  currentStrain.setData(workV6, currentStress.volume()*bulkCoeff);

  double octalStrain = currentStrain.octahedralShear(1);
  if (octalStrain <= LOW_LIMIT) octalStrain = LOW_LIMIT;

  // plastic strain state (scaled from elastic strain)
  double scale, PPZLimit;
  if (stressRatio >= stressRatioPT) {  //above PT
    onPPZCommitted = 2;
    prePPZStrainOctaCommitted = strainPTOcta * ratio;
    PPZLimit = getPPZLimits(1,currentStress);
    scale = sqrt(prePPZStrainOctaCommitted+PPZLimit)/octalStrain;
  }
  else {  // below PT
    onPPZCommitted = -1;
    prePPZStrainOctaCommitted = octalStrain;
    if (prePPZStrainOctaCommitted > strainPTOcta * ratio)
		prePPZStrainOctaCommitted=strainPTOcta*ratio;
    scale = sqrt(prePPZStrainOctaCommitted)/octalStrain;
  }
  //currentStrain.setData(currentStrain.deviator()*scale, currentStrain.volume());
  workV6.addVector(0.0, currentStrain.deviator(), scale);
  currentStrain.setData(workV6, currentStrain.volume());
  PPZPivotCommitted = currentStrain;
}


double PressureDependMultiYield03::getModulusFactor(T2Vector & stress)
{
    double residualPress = residualPressx[matN];
    double refPressure = refPressurex[matN];
    double pressDependCoeff =pressDependCoeffx[matN];

  double conHeig = stress.volume() - residualPress;
  double scale = conHeig / (refPressure-residualPress);
  scale = pow(scale, pressDependCoeff);

  return (1.e-10>scale) ? 1.e-10 : scale;
}


void PressureDependMultiYield03::setTrialStress(T2Vector & stress)
{
    double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];

  modulusFactor = getModulusFactor(stress);
  //workV6 = stress.deviator()
  //	             + subStrainRate.deviator()*2.*refShearModulus*modulusFactor;
  workV6 = stress.deviator();
  workV6.addVector(1.0, subStrainRate.deviator(), 2*refShearModulus*modulusFactor);

  double B = refBulkModulus*modulusFactor;

  if (Hvx[matN] != 0. && trialStress.volume()<=maxPress
	  && subStrainRate.volume()<0. && loadStagex[matN] == 1) {
     double tp = fabs(trialStress.volume() - residualPressx[matN]);
     B = (B*Hvx[matN]*pow(tp,Pvx[matN]))/(B+Hvx[matN]*pow(tp,Pvx[matN]));
  }

  double volume = stress.volume() + subStrainRate.volume()*3.*B;

  if (volume > 0.) volume = 0.;
  trialStress.setData(workV6, volume);
}


int PressureDependMultiYield03::setSubStrainRate(void)
{
    double residualPress = residualPressx[matN];
    double refShearModulus = refShearModulusx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

  //if (activeSurfaceNum==numOfSurfaces) return 1;
  if (strainRate.isZero()) return 0;

  double elast_plast_modulus;
  double conHeig = -(currentStress.volume() - residualPress);
  double factor = getModulusFactor(currentStress);
  if (activeSurfaceNum==0)
    elast_plast_modulus = 2*refShearModulus*factor;
  else {
    double plast_modulus = theSurfaces[activeSurfaceNum].modulus()*factor;
    elast_plast_modulus = 2*refShearModulus*factor*plast_modulus
      / (2*refShearModulus*factor+plast_modulus);
  }
  //workV6 = strainRate.deviator()*elast_plast_modulus;
  workV6.addVector(0.0, strainRate.deviator(), elast_plast_modulus);
  workT2V.setData(workV6,0);

  double singleCross = theSurfaces[numOfSurfaces].size()*conHeig / numOfSurfaces;
  double totalCross = 3.*workT2V.octahedralShear() / sqrt(2.);
  int numOfSub = totalCross/singleCross + 1;
  if (numOfSub > numOfSurfaces) numOfSub = numOfSurfaces;

  int numOfSub1 = strainRate.octahedralShear(1) / 1.0e-5;
  int numOfSub2 = strainRate.volume() / 1.e-5;
  if (numOfSub1 > numOfSub) numOfSub = numOfSub1;
  if (numOfSub2 > numOfSub) numOfSub = numOfSub2;

  workV6.addVector(0.0, strainRate.t2Vector(), 1.0/numOfSub);

  subStrainRate.setData(workV6);

  return numOfSub;
}


void
PressureDependMultiYield03::getContactStress(T2Vector &contactStress)
{
    double residualPress = residualPressx[matN];

  double conHeig = trialStress.volume() - residualPress;
  static Vector center(6);
  center = theSurfaces[activeSurfaceNum].center();
  //workV6 = trialStress.deviator() - center*conHeig;
  workV6 = trialStress.deviator();
  workV6.addVector(1.0, center, -conHeig);
  double Ms = sqrt(3./2.*(workV6 && workV6));
  //workV6 = workV6 * theSurfaces[activeSurfaceNum].size()*(-conHeig) / Ms + center*conHeig;
  workV6.addVector(theSurfaces[activeSurfaceNum].size()*(-conHeig) / Ms, center, conHeig);
  //return T2Vector(workV6,trialStress.volume());
  contactStress.setData(workV6,trialStress.volume());
}


int PressureDependMultiYield03::isLoadReversal(const T2Vector & stress)
{
  if(activeSurfaceNum == 0) return 0;

  getSurfaceNormal(stress, workT2V);

  //if (((trialStress.t2Vector() - currentStress.t2Vector())
  //	&& workT2V.t2Vector()) < 0) return 1;
  workV6 = trialStress.t2Vector();
  workV6 -= currentStress.t2Vector();

  if ((workV6 && workT2V.t2Vector()) < 0) return 1;

  return 0;
}


void
PressureDependMultiYield03::getSurfaceNormal(const T2Vector & stress, T2Vector &normal)
{
    double residualPress = residualPressx[matN];

  double conHeig = stress.volume() - residualPress;
  workV6 = stress.deviator();
  static Vector center(6);
  center = theSurfaces[activeSurfaceNum].center();
  double sz = theSurfaces[activeSurfaceNum].size();
  double volume = conHeig*((center && center) - 2./3.*sz*sz) - (workV6 && center);
  //workT2V.setData((workV6-center*conHeig)*3., volume);
  workV6.addVector(1.0, center, -conHeig);
  workV6 *= 3.0;
  workT2V.setData(workV6, volume);

  normal.setData(workT2V.unitT2Vector());
}


double PressureDependMultiYield03::getPlasticPotential(const T2Vector & contactStress,
						     const T2Vector & surfaceNormal)
{
    double residualPress = residualPressx[matN];
    double stressRatioPT = stressRatioPTx[matN];
	double contractParam1 = contractParam1x[matN];
	double contractParam2 = contractParam2x[matN];
	double contractParam3 = contractParam3x[matN];
    double dilateParam1 = dilateParam1x[matN];
    double dilateParam2 = dilateParam2x[matN];

  double plasticPotential, contractRule, shearLoading, angle;

  double contactRatio = contactStress.deviatorRatio(residualPress);
  double factorPT = contactRatio/stressRatioPT;
  
  double C = ContractionFactorC();
  
  //double currentRatio = currentStress.deviatorRatio(residualPress);
  double currentRatio = updatedTrialStress.deviatorRatio(residualPress);
  double trialRatio = trialStress.deviatorRatio(residualPress);
  shearLoading = updatedTrialStress.deviator() && trialStress.deviator();
  //shearLoading = currentStress.deviator() && trialStress.deviator();

  if (factorPT >= 1. && trialRatio >= currentRatio && shearLoading >= 0.) {  //dilation
      updatePPZ(contactStress);
      if (onPPZ==1)
		  plasticPotential = 0.;
      else if (onPPZ==2) {
          factorPT -= 1.0;
		    	double dilateParam3 = dilateParam3x[matN];
		  double ppp=pow((fabs(contactStress.volume())+fabs(residualPress))/pAtm, -dilateParam3);
          plasticPotential = ppp*factorPT*(factorPT)*(dilateParam1+pow(cumuDilateStrainOcta,dilateParam2));
          if (plasticPotential < 0.) plasticPotential = -plasticPotential;
		  if (plasticPotential>5.0e4) plasticPotential = 5.0e4;
      }
	  else {
          opserr << "FATAL: Wrong onPPZ value: " << onPPZ << endln;
          exit(-1);
	  }
  } else {  //contraction
	  if (currentRatio==0.) angle = 1.0;
      else {
         workV6 = trialStress.deviator();
		 workV6 /= (fabs(trialStress.volume())+fabs(residualPress));
		 workV6	-= updatedTrialStress.deviator()/(fabs(updatedTrialStress.volume())+fabs(residualPress));
		 //workV6	-= currentStress.deviator()/(fabs(currentStress.volume())+fabs(residualPress));
		 //workV6.Normalize();
		 //angle = updatedTrialStress.unitDeviator() && workV6;
		 workT2V = T2Vector(workV6);
		 if (workT2V.deviatorLength() == 0.) angle = 1.0;
		 //angle = (currentStress.deviator() && workV6)/workT2V.deviatorLength()/currentStress.deviatorLength();
		 else angle = (updatedTrialStress.deviator() && workV6)/workT2V.deviatorLength()/updatedTrialStress.deviatorLength();
	  }
      factorPT = factorPT*angle - 1.0;

	  contractRule = pow((fabs(contactStress.volume())+fabs(residualPress))/pAtm, contractParam3);
      if (contractRule < 0.1) contractRule = 0.1;

	  //plasticPotential = factorPT*(contractParam1+pressureD*contractParam2)*contractRule;
	  //plasticPotential = factorPT*(contractParam1+cumuDilateStrainOcta*contractParam2)*contractRule;
	  //double dd = maxCumuDilateStrainOcta > cumuDilateStrainOcta ? maxCumuDilateStrainOcta : cumuDilateStrainOcta;
	  //plasticPotential = -(factorPT)*(factorPT)*(contractParam1+dd*contractParam2)*contractRule;
	  
	  plasticPotential = -(factorPT)*(factorPT)*(contractParam1+maxCumuDilateStrainOcta*contractParam2)*contractRule;
	  
	  // PDMY03
	  plasticPotential *= C;
	  
      if (plasticPotential > 0.) plasticPotential = -plasticPotential;
	  //if (contractRule<-5.0e4) contractRule = -5.0e4;
      if (onPPZ > 0) onPPZ = 0;
      if (onPPZ != -1) PPZTranslation(contactStress);
  }

  if (isCriticalState(contactStress)) plasticPotential = 0;
  return plasticPotential;
}

double PressureDependMultiYield03::ContractionFactorC(void)
{
	double residualPress = residualPressx[matN];
	double cd = contractParam4x[matN];
	double ce = contractParam5x[matN];

	// Start - by Arash K.
	static Vector stress0 = currentStress.t2Vector();
	double sig110 = stress0[0];
	double sig220 = stress0[1];
	double sig330 = stress0[2];
	double tau120 = stress0[3];
	double tau230 = stress0[4];
	double tau130 = stress0[5];
	double factorP0 = fabs(sig110+sig220+sig330)/3 + fabs(residualPress);
	double factorCSR0 = pow(tau120*tau120 + tau230*tau230 + tau130*tau130 , 0.5) / factorP0;

	Vector stress = currentStress.t2Vector();
	double sig11 = stress[0];
	double sig22 = stress[1];
	double sig33 = stress[2];
	double tau12 = stress[3];
	double tau23 = stress[4];
	double tau13 = stress[5];
	double factorCSR = pow(tau12*tau12 + tau23*tau23 + tau13*tau13 , 0.5) / factorP0;
	
	double C = (1.0 + pow(fabs(factorCSR-factorCSR0)*cd , 3.0)) * pow(1.0 + ce*factorCSR0,2.0);
	// End	
	
	return C;
}

int PressureDependMultiYield03::isCriticalState(const T2Vector & stress)
{
	double einit = einitx[matN];
	double volLimit1 = volLimit1x[matN];
	double volLimit2 = volLimit2x[matN];
	double volLimit3 = volLimit3x[matN];

  double vol = trialStrain.volume()*3.0;
	double etria = einit + vol + vol*einit;
	vol = currentStrain.volume()*3.0;
	double ecurr = einit + vol + vol*einit;

	double ecr1, ecr2;
	if (volLimit3 != 0.) {
		ecr1 = volLimit1 - volLimit2*pow(fabs(-stress.volume()/pAtm), volLimit3);
	  ecr2 = volLimit1 - volLimit2*pow(fabs(-updatedTrialStress.volume()/pAtm), volLimit3);
	} else {
		ecr1 = volLimit1 - volLimit2*log(fabs(-stress.volume()/pAtm));
	  ecr2 = volLimit1 - volLimit2*log(fabs(-updatedTrialStress.volume()/pAtm));
	}

	if (ecurr < ecr2 && etria < ecr1) return 0;
	if (ecurr > ecr2 && etria > ecr1) return 0;
	return 1;
}


void PressureDependMultiYield03::updatePPZ(const T2Vector & contactStress)
{
  double liquefyParam1 = liquefyParam1x[matN];
  double residualPress = residualPressx[matN];
  double refPressure = refPressurex[matN];
  double pressDependCoeff =pressDependCoeffx[matN];
  	double liquefyParam2 = liquefyParam2x[matN];

  // onPPZ=-1 may not be needed. can start with onPPZ=0  ****

  double temp = strainRate.deviator() && PivotStrainRateCommitted;
  check = strainRate.deviator()[3];

  if (onPPZ < 1) {
    damage = 0.0;
    if ((maxPress-currentStress.volume())/(maxPress-residualPress) > 0.)
	  damage = pow((maxPress-currentStress.volume())/(maxPress-residualPress),0.25);
  }

  // PPZ inactive if liquefyParam1==0.
  if (liquefyParam1==0. || (onPPZ < 1 && damage < 0.)) {
    if (onPPZ==2) {
       PPZPivot = trialStrain;
	   //PivotStrainRate = strainRate.deviator();
	   cumuDilateStrainOcta += subStrainRate.octahedralShear(1);
    }
    else if (onPPZ != 2) {
       onPPZ = 2;
       PPZPivot = trialStrain;
	   PivotStrainRate = strainRate.deviator();
       if (temp < 0.) cumuDilateStrainOcta = 0.;
	  //if (temp < 0.) maxCumuDilateStrainOcta = 0.;
    }
    return;
  }

  // dilation: calc. cumulated dilative strain
  if (onPPZ==2) {
    PPZPivot = trialStrain;
	//PivotStrainRate = strainRate.deviator();
	cumuDilateStrainOcta += subStrainRate.octahedralShear(1);

	double zzz = 0.;
	if (damage > zzz) zzz = damage;
	maxCumuDilateStrainOcta += zzz*liquefyParam1*subStrainRate.octahedralShear(1);
    return;
  }

  if (onPPZ==-1 || onPPZ==0) {

	  // moved to the opposite direction, update prePPZStrainOcta and
	  // oppoPrePPZStrainOcta (they are small values, may consider drop these
	  // parameters altogether).

	  if (temp < 0.) {
        double volume = -contactStress.volume();
        oppoPrePPZStrainOcta = prePPZStrainOcta;
        double ratio = (volume+residualPress)/(-refPressure+residualPress);
        ratio = pow(ratio, 1.-pressDependCoeff);
        prePPZStrainOcta = ratio * strainPTOcta;
        if (oppoPrePPZStrainOcta == 0.) oppoPrePPZStrainOcta = prePPZStrainOcta;
	  }
  }

  //double dd = maxCumuDilateStrainOcta > cumuDilateStrainOcta ? maxCumuDilateStrainOcta : cumuDilateStrainOcta;

  //PPZSize = (cumuTranslateStrainOcta+dd)/2.;
  PPZSize = (cumuTranslateStrainOcta+maxCumuDilateStrainOcta)/2.;
	 //(prePPZStrainOcta+oppoPrePPZStrainOcta+cumuTranslateStrainOcta+maxCumuDilateStrainOcta)/2.;

  // calc. new PPZ center.
  if (onPPZ==0 || (onPPZ==1 && temp < 0.0)) {
	// new center lies on the vector of PPZPivot-PPZCenter,
	// its distance from PPZPivot is (PPZSize-cumuTranslateStrainOcta).

    //workT2V.setData(PPZPivot.t2Vector() - PPZCenter.t2Vector());
    workV6 = PPZPivot.t2Vector();
    workV6.addVector(1.0, PPZCenter.t2Vector(), -1.);
    workT2V.setData(workV6);

    double coeff;
	if (workT2V.octahedralShear(1)==0.) coeff = 0.;
	else coeff = (PPZSize-cumuTranslateStrainOcta)/workT2V.octahedralShear(1);
    //PPZCenter.setData(PPZPivot.t2Vector() - workT2V.t2Vector()*coeff);
    workV6 = PPZPivot.t2Vector();
    workV6.addVector(1.0, workT2V.t2Vector(), -coeff);
    PPZCenter.setData(workV6);
  }

  //workT2V.setData(trialStrain.t2Vector() - PPZCenter.t2Vector());
  workV6 = trialStrain.t2Vector();
  workV6.addVector(1.0, PPZCenter.t2Vector(), -1.);
  workT2V.setData(workV6);

  //outside PPZ
  //if (workT2V.octahedralShear(1) > PPZSize && temp > 0. || PPZLimit==0.) {
  if (workT2V.octahedralShear(1) > PPZSize) {

	// rezero cumuDilateStrainOcta only when load reverses
	// (partial unloading does not erase cumuDilateStrainOcta.

    //if (temp < 0.) {
         cumuDilateStrainOcta = 0.;
         //if (PPZLimit == 0.) maxCumuDilateStrainOcta = 0.;
	//}
    onPPZ = 2;
    PPZPivot = trialStrain;
	PivotStrainRate = strainRate.deviator();
    cumuTranslateStrainOcta = 0.;
  }
  else {  //inside PPZ
    if (onPPZ == 0 || onPPZ == 1) PPZTranslation(contactStress);
    //if (onPPZ == 0) onPPZ = 1;
	if (onPPZ == -1 || onPPZ == 0) onPPZ = 1;
  }
}


void PressureDependMultiYield03::PPZTranslation(const T2Vector & contactStress)
{
	double liquefyParam1 = liquefyParam1x[matN];
  	double liquefyParam2 = liquefyParam2x[matN];
    double residualPress = residualPressx[matN];

	//cumuDilateStrainOcta -= subStrainRate.octahedralShear(1);
    //if (cumuDilateStrainOcta < 0.) cumuDilateStrainOcta = 0.;

    if (liquefyParam1==0.) return;

  //Amount of translation is proportional to the amount of previous unloading,
  //and is also limited by the amount of previous dilation (no translation
  //if there was no dilation), as damage is really caused by dilation.
  damage = 0.0;
  if ((maxPress-currentStress.volume())/(maxPress-residualPress) > 0.)
	damage = pow((maxPress-currentStress.volume())/(maxPress-residualPress),0.25);

  double zzz = 0.;
  if (damage > zzz) zzz = damage;

  double temp = strainRate.deviator() && PivotStrainRateCommitted;

  if (temp < 0.0) {  //update only when load reverses

  //workT2V.setData(trialStrain.deviator()-PPZPivot.deviator());
     workV6 = trialStrain.deviator();
     workV6 -= PPZPivot.deviator();
     workT2V.setData(workV6);

     temp = workT2V.octahedralShear(1);
     if (cumuTranslateStrainOcta < zzz*liquefyParam2*temp)
	     cumuTranslateStrainOcta = zzz*liquefyParam2*temp;
  //if (maxCumuDilateStrainOcta == 0.) temp = 0.; //PPZTransLimit;
	//\\// else temp = PPZTransLimit*cumuDilateStrainOcta/maxCumuDilateStrainOcta;
  //else temp = dilateParam3*cumuDilateStrainOcta/maxCumuDilateStrainOcta;
  //if (cumuTranslateStrainOcta > temp) cumuTranslateStrainOcta = temp;

  //\\// cumuTranslateStrainOcta = dilateParam3*cumuDilateStrainOcta;
  }
}


double PressureDependMultiYield03::getPPZLimits(int which, const T2Vector & contactStress)
{
	double liquefyParam1 = liquefyParam1x[matN];
	double liquefyParam2 = liquefyParam2x[matN];
	double dilateParam3 = dilateParam3x[matN];

  double PPZLimit, temp;
  double volume = -contactStress.volume();

  if (volume >= liquefyParam1) PPZLimit = 0.;
  else {
    temp = volume*pi/liquefyParam1/2.;
		// liquefyParam3 = 3.0 default
    PPZLimit = liquefyParam2 * pow(cos(temp), 3.);
	//\\//
    PPZLimit = 0.0;
  }

  if (which==1)
    return PPZLimit;
  else if (which==2)
    return dilateParam3 * PPZLimit;
  else {
    opserr << "FATAL:PressureDependMultiYield03::getPPZLimits: unknown argument value" << endln;
   exit(-1);
    return 0.0;
  }
}


double PressureDependMultiYield03::getLoadingFunc(const T2Vector & contactStress,
						const T2Vector & surfaceNormal,
						double * plasticPotential,
						int crossedSurface)
{
    int numOfSurfaces = numOfSurfacesx[matN];
    double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];

  double loadingFunc, limit;
  double modul = theSurfaces[activeSurfaceNum].modulus();
  double temp1 = 2. * refShearModulus * modulusFactor
    * (surfaceNormal.deviator() && surfaceNormal.deviator()) ;
  double temp2 = 9. * refBulkModulus * modulusFactor
    * surfaceNormal.volume() * (*plasticPotential) ;

  //for the first crossing
  double temp = temp1 + temp2 + modul * modulusFactor;
  if (activeSurfaceNum == numOfSurfaces)
    limit = theSurfaces[activeSurfaceNum-1].modulus() * modulusFactor /2.;
  else limit = modul * modulusFactor /2.;
  if (temp < limit) {
    (*plasticPotential) = (temp2 + limit - temp)/(9. * refBulkModulus * modulusFactor
					       * surfaceNormal.volume());
    temp = limit;
  }
  //loadingFunc = (surfaceNormal.t2Vector()
  //	             && (trialStress.deviator()-contactStress.deviator()))/temp;
  workV6 = trialStress.deviator();
  workV6 -= contactStress.deviator();
  loadingFunc = (surfaceNormal.t2Vector() && workV6) /temp;

  if (loadingFunc < 0.) loadingFunc = 0;

  //for more than one crossing
  if(crossedSurface) {
    temp = (theSurfaces[activeSurfaceNum-1].modulus() - modul)
      /theSurfaces[activeSurfaceNum-1].modulus();
    loadingFunc *= temp;
  }

  return loadingFunc;
}


int PressureDependMultiYield03::stressCorrection(int crossedSurface)
{
    double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];

  static T2Vector contactStress;
  getContactStress(contactStress);
  static T2Vector surfNormal;
  getSurfaceNormal(contactStress, surfNormal);
  double plasticPotential = getPlasticPotential(contactStress,surfNormal);
  double tVolume = trialStress.volume();
  double loadingFunc = getLoadingFunc(contactStress, surfNormal,
				      &plasticPotential, crossedSurface);
  double volume = tVolume - plasticPotential*3*refBulkModulus*modulusFactor*loadingFunc;

  //workV6 = trialStress.deviator()
  //	         - surfNormal.deviator()*2*refShearModulus*modulusFactor*loadingFunc;
  workV6 = trialStress.deviator();

  if (volume > 0. && volume != tVolume) {
		double coeff = tVolume / (tVolume - volume);
    coeff *= -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
	volume = 0.;
  } else if (volume > 0.) {
	volume = 0.;
  } else {
	double coeff = -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
  }
	/*
	  if (volume>0.)volume = 0.;
		double coeff = -2*refShearModulus*modulusFactor*loadingFunc;
    workV6.addVector(1.0, surfNormal.deviator(), coeff);
  */
  trialStress.setData(workV6, volume);
  deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

  if (isCrossingNextSurface()) {
    activeSurfaceNum++;
    return stressCorrection(1);  //recursive call
  }

  return 0;
}


void PressureDependMultiYield03::updateActiveSurface(void)
{
    double residualPress = residualPressx[matN];
    int numOfSurfaces = numOfSurfacesx[matN];

  if (activeSurfaceNum == numOfSurfaces) return;

  double A, B, C, X;
  static Vector t1(6);
  static Vector t2(6);
  static Vector center(6);
  static Vector outcenter(6);
  double conHeig = trialStress.volume() - residualPress;
  center = theSurfaces[activeSurfaceNum].center();
  double size = theSurfaces[activeSurfaceNum].size();
  outcenter = theSurfaces[activeSurfaceNum+1].center();
  double outsize = theSurfaces[activeSurfaceNum+1].size();

  //t1 = trialStress.deviator() - center*conHeig;
  //t2 = (center - outcenter)*conHeig;
  t1 = trialStress.deviator();
  t1.addVector(1.0, center, -conHeig);
  t2 = center;
  t2 -= outcenter;
  t2 *=conHeig;

  A = t1 && t1;
  B = 2. * (t1 && t2);
  C = (t2 && t2) - 2./3.* outsize * outsize * conHeig * conHeig;
  X = secondOrderEqn(A,B,C,0);
  if ( fabs(X-1.) < LOW_LIMIT ) X = 1.;
  if (X < 1.) return;

  if (X < 1.){
    //t2 = trialStress.deviator() - outcenter*conHeig;
    t2 = trialStress.deviator();
    t2.addVector(1.0, outcenter, -conHeig);

    double xx1 = (t2 && t2) - 2./3.* outsize * outsize * conHeig * conHeig;
    double xx2 = (t1 && t1) - 2./3.* size * size * conHeig * conHeig;
    opserr << "FATAL:PressureDependMultiYield03::updateActiveSurface(): error in Direction of surface motion." << endln;
    opserr << "X-1= " << X-1 <<" A= "<<A<<" B= "<<B<<" C= "<<C <<" M= "<<activeSurfaceNum<<" low_limit="<<LOW_LIMIT<< endln;
    opserr << "diff1= "<<xx1 <<" diff2= "<<xx2 <<" p= "<<conHeig<<" size= "<<size<<" outs= "<<outsize<<endln;
   exit(-1);
  }

  //workV6 = (t1 * X + center*conHeig) * (1. - size / outsize)
  //	     - (center - outcenter * size / outsize) * conHeig;

  workV6.addVector(0.0, t1, X);
  workV6.addVector(1.0, center, conHeig);
  workV6 *= (1.0 - size / outsize);
  t2 = center;
  t2.addVector(1.0, outcenter, -size/outsize);
  t2 *= conHeig;
  workV6 -= t2;

  workT2V.setData(workV6);
  if (workT2V.deviatorLength() < LOW_LIMIT) return;

  workV6 = workT2V.deviator();
  A = conHeig * conHeig * (workV6 && workV6);
  B = 2 * conHeig * (t1 && workV6);
  if (fabs(B) < LOW_LIMIT) B = 0.;
  C = (t1 && t1) - 2./3.* size * size * conHeig * conHeig;
  if ( fabs(C) < LOW_LIMIT || fabs(C)/(t1 && t1) < LOW_LIMIT ) return;
  if (B > 0. || C < 0.) {
    opserr << "FATAL:PressureDependMultiYield03::updateActiveSurface(): error in surface motion.\n"
	 << "A= " <<A <<" B= " <<B <<" C= "<<C <<" (t1&&t1)= "<<(t1&&t1) <<endln;
   exit(-1);
  }
  X = secondOrderEqn(A,B,C,1);

  center.addVector(1.0, workV6, -X);
  theSurfaces[activeSurfaceNum].setCenter(center);
}


void PressureDependMultiYield03::updateInnerSurface(void)
{
    double residualPress = residualPressx[matN];

	if (activeSurfaceNum <= 1) return;
	static Vector devia(6);
	static Vector center(6);

	double conHeig = currentStress.volume() - residualPress;
	devia = currentStress.deviator();
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();

	for (int i=1; i<activeSurfaceNum; i++) {
		workV6.addVector(0.0, center, conHeig);
		workV6 -= devia;
		workV6 *= theSurfaces[i].size()/size;
		workV6 += devia;

		//workV6 = devia - (devia - center*conHeig) * theSurfaces[i].size() / size;

		workV6 /= conHeig;
		theSurfaces[i].setCenter(workV6);
	}
}


int PressureDependMultiYield03:: isCrossingNextSurface(void)
{
    int numOfSurfaces = numOfSurfacesx[matN];

  if (activeSurfaceNum == numOfSurfaces) return 0;

  if(yieldFunc(trialStress, theSurfaces, activeSurfaceNum+1) > 0) return 1;

  return 0;
}

