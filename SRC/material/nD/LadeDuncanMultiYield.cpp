/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright for this material model is by Zhijian Qiu, UCSD      **
**                                                                    **
** Commercial use of this program without express permission of the   **
** Authors is strictly prohibited.                                    **
**                                                                    **
** Developed by:                                                      **
**   Zhijian Qiu (zhq009@eng.ucsd.edu)                                **
** Qiu, Z., Prabhakaran, A., & Elgamal, A. (2023).                    **
** A three-dimensional multi-surface plasticity soil model for        **
** seismically-induced liquefaction and earthquake loading            **
** applications. Acta Geotechnica, 1-24.                              **
** 										                              **
**                                                                    **
** Yang, Z. and Elgamal, A., 2008. Multi-surface cyclic plasticity    **
** sand model with Lode angle effect. Geotechnical and Geological     **
** Engineering, 26(3), pp.335-348.                                    **
**                                                                    **
** ****************************************************************** */

    
#include <math.h>
#include <stdlib.h>
#include <LadeDuncanMultiYield.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>
#include <string.h>
#include <elementAPI.h>
#include <complex>
#include <iostream>
using namespace std;
typedef complex<double> dcmplx;

int LadeDuncanMultiYield::matCount = 0;
int* LadeDuncanMultiYield::loadStagex = 0;  //=0 if elastic; =1 if plastic
int* LadeDuncanMultiYield::ndmx = 0;  //num of dimensions (2 or 3)
double* LadeDuncanMultiYield::rhox = 0;
double* LadeDuncanMultiYield::refShearModulusx = 0;
double* LadeDuncanMultiYield::refBulkModulusx = 0;
double* LadeDuncanMultiYield::frictionAnglex = 0;
double* LadeDuncanMultiYield::peakShearStrainx = 0;
double* LadeDuncanMultiYield::refPressurex = 0;
double* LadeDuncanMultiYield::cohesionx = 0;
double* LadeDuncanMultiYield::pressDependCoeffx = 0;
int*    LadeDuncanMultiYield::numOfSurfacesx = 0;
double* LadeDuncanMultiYield::phaseTransfAnglex = 0;
double* LadeDuncanMultiYield::contractParam1x = 0;
double* LadeDuncanMultiYield::contractParam2x = 0;
double* LadeDuncanMultiYield::contractParam3x = 0;
double* LadeDuncanMultiYield::dilateParam1x = 0;
double* LadeDuncanMultiYield::dilateParam2x = 0;
double* LadeDuncanMultiYield::liquefyParam1x = 0;
double* LadeDuncanMultiYield::liquefyParam2x = 0;
double* LadeDuncanMultiYield::dilateParam3x = 0;
double* LadeDuncanMultiYield::einitx = 0;    //initial void ratio
double* LadeDuncanMultiYield::volLimit1x = 0;
double* LadeDuncanMultiYield::volLimit2x = 0;
double* LadeDuncanMultiYield::volLimit3x = 0;
double* LadeDuncanMultiYield::residualPressx = 0;
double* LadeDuncanMultiYield::stressRatioPTx = 0;
double* LadeDuncanMultiYield::a1x = 0;
double* LadeDuncanMultiYield::Hvx = 0;
double* LadeDuncanMultiYield::Pvx = 0;
double* LadeDuncanMultiYield::dilateParam4x = 0;
double* LadeDuncanMultiYield::contractParam4x = 0;
double* LadeDuncanMultiYield::contractParam5x = 0;
double* LadeDuncanMultiYield::liquefyParam3x = 0;
double LadeDuncanMultiYield::pAtm = 101.;
Matrix LadeDuncanMultiYield::theTangent(6, 6); 
T2Vector LadeDuncanMultiYield::trialStrain;
T2Vector LadeDuncanMultiYield::subStrainRate;
Vector LadeDuncanMultiYield::workV6(6);
T2Vector LadeDuncanMultiYield::workT2V;

const	double PI = 3.14159265358979;

void* OPS_LadeDuncanMultiYield()
{
	const int numParam = 13;
	const int totParam = 30;
	int tag;
	double param[totParam];
	param[numParam] = 20;
	param[numParam + 1] = 5.0;
	param[numParam + 2] = 3.;
	param[numParam + 3] = 1.;
	param[numParam + 4] = 0.;
	param[numParam + 5] = 0.6;
	param[numParam + 6] = 0.9;
	param[numParam + 7] = 0.02;
	param[numParam + 8] = 0.7;
	param[numParam + 9] = 101.;
	param[numParam + 10] = 0.1;

	param[numParam + 11] = 1.0;
	param[numParam + 12] = 0.0;
	param[numParam + 13] = 0.0;
	param[numParam + 14] = 1.0;

	param[numParam + 15] = 0.;
	param[numParam + 16] = 1.;

	int argc = OPS_GetNumRemainingInputArgs() + 2;
	char * arg[] = { "nd", "rho", "refShearModul",
		"refBulkModul", "frictionAng",
		"peakShearStra", "refPress", "pressDependCoe",
		"phaseTransformAngle", "contractionParam1",
		"contractionParam3","dilationParam1","dilationParam3",
		"numberOfYieldSurf (=20)",
		"contractionParam2=5.0", "dilationParam2=3.0",
		"liquefactionParam1=1.0", "liquefactionParam2=0.0",
		"e (=0.6)", "volLimit1 (=0.9)", "volLimit2 (=0.02)",
		"volLimit3 (=0.7)", "Atmospheric pressure (=101)", "cohesi (=.1)",
		"dilationParam4(=1.0)", "contractionParam4(=0.0)", "contractionParam5(=0.0)", "liquefactionParam3(=1.0)",
		"Hv (=0)", "Pv (=1.)" };
	if (argc < (3 + numParam)) {
		opserr << "WARNING insufficient arguments\n";
		opserr << "Want: nDMaterial LadeDuncanMultiYield tag? " << arg[0];
		opserr << "? " << "\n";
		opserr << arg[1] << "? " << arg[2] << "? " << arg[3] << "? " << "\n";
		opserr << arg[4] << "? " << arg[5] << "? " << arg[6] << "? " << "\n";
		opserr << arg[7] << "? " << arg[8] << "? " << arg[9] << "? " << "\n";
		opserr << arg[10] << "? " << arg[11] << "? " << arg[12] << "? " << "\n";
		opserr << arg[13] << "? " << arg[14] << "? " << arg[15] << "? " << "\n";
		opserr << arg[16] << "? " << arg[17] << "? " << arg[18] << "? " << "\n";
		opserr << arg[19] << "? " << arg[20] << "? " << arg[21] << "? " << "\n";
		opserr << arg[22] << "? " << arg[23] << "? " << arg[24] << "? " << "\n";
		opserr << arg[25] << "? " << arg[26] << "? " << arg[27] << "? " << endln;
		return 0;
	}

	int numdata = 1;
	if (OPS_GetIntInput(&numdata, &tag) < 0) {
		opserr << "WARNING invalid LadeDuncanMultiYield tag" << "\n";
		return 0;
	}

	int in = 17;
	for (int i = 3; (i<argc && i<in); i++)
		if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
			opserr << "WARNING invalid " << arg[i - 3] << "\n";
			opserr << "nDMaterial LadeDuncanMultiYield: " << tag << "\n";
			return 0;
		}

	static double * gredu = 0;

	// user defined yield surfaces
	if (param[numParam] < 0 && param[numParam] > -100) {
		param[numParam] = -int(param[numParam]);
		gredu = new double[int(2 * param[numParam])];

		for (int i = 0; i<2 * param[numParam]; i++)
			if (OPS_GetDoubleInput(&numdata, &gredu[i]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial PressureIndependMultiYield: " << tag << "\n";
				return 0;
			}
	}

	if (gredu != 0) {
		for (int i = in + int(2 * param[numParam]); i<argc; i++)
			if (OPS_GetDoubleInput(&numdata, &param[i - 3 - int(2 * param[numParam])]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial LadeDuncanMultiYield: " << tag << "\n";
				return 0;
			}
	}
	else {
		for (int i = in; i<argc; i++)
			if (OPS_GetDoubleInput(&numdata, &param[i - 3]) < 0) {
				opserr << "WARNING invalid " << " double" << "\n";
				opserr << "nDMaterial LadeDuncanMultiYield: " << tag << "\n";
				return 0;
			}
	}

	LadeDuncanMultiYield * temp =
		new LadeDuncanMultiYield(tag, param[0], param[1], param[2],
			param[3], param[4], param[5],
			param[6], param[7], param[8],
			param[9], param[10], param[11],
			param[12], param[13], gredu, param[14],
			param[15], param[16], param[17],
			param[18], param[19], param[20], param[21],
			param[22], param[23], param[24], param[25], param[26], param[27], param[28], param[29]);

	if (gredu != 0) {
		delete[] gredu;
		gredu = 0;
	}

	return temp;
}

LadeDuncanMultiYield::LadeDuncanMultiYield(int tag, int nd,
	double r, double refShearModul,
	double refBulkModul, double frictionAng,
	double peakShearStra, double refPress,
	double pressDependCoe,
	double phaseTransformAng,
	double contractionParam1,
	double contractionParam3,
	double dilationParam1,
	double dilationParam3,
	int   numberOfYieldSurf,
	double * gredu,
	double contractionParam2,
	double dilationParam2,
	double liquefactionParam1,
	double liquefactionParam2,
	double ei,
	double volLim1, double volLim2, double volLim3,
	double atm, double cohesi, double dilationParam4, double contractionParam4, double contractionParam5, double liquefactionParam3,
	double hv, double pv)
	: NDMaterial(tag, ND_TAG_LadeDuncanMultiYield), currentStress(),
	trialStress(), updatedTrialStress(), currentStrain(), strainRate(),
	PPZPivot(), PPZCenter(), PPZPivotCommitted(), PPZCenterCommitted(),
	PivotStrainRate(6), PivotStrainRateCommitted(6)
{
	if (nd != 2 && nd != 3) {
		opserr << "FATAL:LadeDuncanMultiYield:: dimension error" << endln;
		opserr << "Dimension has to be 2 or 3, you give nd= " << nd << endln;
		exit(-1);
	}
	if (refShearModul <= 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: refShearModulus <= 0" << endln;
		exit(-1);
	}
	if (refBulkModul <= 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: refBulkModulus <= 0" << endln;
		exit(-1);
	}
	if (frictionAng <= 0.) {
		opserr << "FATAL:LadeDuncanMultiYield:: frictionAngle <= 0" << endln;
		exit(-1);
	}
	if (frictionAng >= 90.) {
		opserr << "FATAL:LadeDuncanMultiYield:: frictionAngle >= 90" << endln;
		exit(-1);
	}
	if (phaseTransformAng <= 0.) {
		opserr << "FATAL:LadeDuncanMultiYield:: phaseTransformAng " << phaseTransformAng << "<= 0" << endln;
		exit(-1);
	}

	if (cohesi < 0) {
		opserr << "WARNING:LadeDuncanMultiYield:: cohesion < 0" << endln;
		opserr << "Will reset cohesion to 0.3." << endln;
		cohesi = 0.3;
	}
	if (peakShearStra <= 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: peakShearStra <= 0" << endln;
		exit(-1);
	}
	if (refPress <= 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: refPress <= 0" << endln;
		exit(-1);
	}
	if (pressDependCoe < 0) {
		opserr << "WARNING:LadeDuncanMultiYield:: pressDependCoe < 0" << endln;
		opserr << "Will reset pressDependCoe to 0.5." << endln;
		pressDependCoe = 0.5;
	}
	if (numberOfYieldSurf <= 0) {
		opserr << "WARNING:LadeDuncanMultiYield:: numberOfSurfaces " << numberOfYieldSurf << "<= 0" << endln;
		opserr << "Will use 10 yield surfaces." << endln;
		numberOfYieldSurf = 10;
	}
	if (numberOfYieldSurf > 100) {
		opserr << "WARNING:LadeDuncanMultiYield::LadeDuncanMultiYield: numberOfSurfaces > 100" << endln;
		opserr << "Will use 100 yield surfaces." << endln;
		numberOfYieldSurf = 100;
	}
	if (volLim1 < 0) {
		opserr << "WARNING:LadeDuncanMultiYield:: volLim1 < 0" << endln;
		opserr << "Will reset volLimit to 0.8" << endln;
		volLim1 = 0.8;
	}
	if (r < 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: rho <= 0" << endln;
		exit(-1);
	}
	if (ei < 0) {
		opserr << "FATAL:LadeDuncanMultiYield:: e <= 0" << endln;
		exit(-1);
	}

	if (matCount % 20 == 0) {
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
		double * temp14 = contractParam1x;
		double * temp14a = contractParam2x;
		double * temp14b = contractParam3x;
		double * temp15 = dilateParam1x;
		double * temp16 = dilateParam2x;
		double * temp17 = liquefyParam1x;
		double * temp18 = liquefyParam2x;
		double * temp19 = dilateParam3x;
		double * temp20 = einitx;    //initial void ratio
		double * temp21 = volLimit1x;
		double * temp22 = volLimit2x;
		double * temp23 = volLimit3x;
		double * temp24 = stressRatioPTx;
		double * temp25 = Hvx;
		double * temp26 = Pvx;
		double * temp27 = a1x;
		double * temp28 = dilateParam4x;
		double * temp29 = contractParam4x;
		double * temp30 = contractParam5x;
		double * temp31 = liquefyParam3x;

		loadStagex = new int[matCount + 20];
		ndmx = new int[matCount + 20];
		rhox = new double[matCount + 20];
		refShearModulusx = new double[matCount + 20];
		refBulkModulusx = new double[matCount + 20];
		frictionAnglex = new double[matCount + 20];
		peakShearStrainx = new double[matCount + 20];
		refPressurex = new double[matCount + 20];
		cohesionx = new double[matCount + 20];
		pressDependCoeffx = new double[matCount + 20];
		numOfSurfacesx = new int[matCount + 20];
		residualPressx = new double[matCount + 20];
		phaseTransfAnglex = new double[matCount + 20];
		contractParam1x = new double[matCount + 20];
		contractParam2x = new double[matCount + 20];
		contractParam3x = new double[matCount + 20];
		dilateParam1x = new double[matCount + 20];
		dilateParam2x = new double[matCount + 20];
		liquefyParam1x = new double[matCount + 20];
		liquefyParam2x = new double[matCount + 20];
		dilateParam3x = new double[matCount + 20];
		einitx = new double[matCount + 20];    //initial void ratio
		volLimit1x = new double[matCount + 20];
		volLimit2x = new double[matCount + 20];
		volLimit3x = new double[matCount + 20];
		stressRatioPTx = new double[matCount + 20];
		Hvx = new double[matCount + 20];
		Pvx = new double[matCount + 20];
		a1x = new double[matCount + 20];
		dilateParam4x = new double[matCount + 20];
		contractParam4x = new double[matCount + 20];
		contractParam5x = new double[matCount + 20];
		liquefyParam3x = new double[matCount + 20];

		for (int i = 0; i<matCount; i++) {
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
			a1x[i] = temp27[i];
			dilateParam4x[i] = temp28[i];
			contractParam4x[i] = temp29[i];
			contractParam5x[i] = temp30[i];
			liquefyParam3x[i] = temp31[i];
		}

		if (matCount > 0) {
			delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
			delete[] temp5; delete[] temp6; delete[] temp7; delete[] temp8;
			delete[] temp9; delete[] temp10; delete[] temp11; delete[] temp12;
			delete[] temp13; delete[] temp14; delete[] temp14a; delete[] temp14b;
			delete[] temp15; delete[] temp16;
			delete[] temp17; delete[] temp18; delete[] temp19; delete[] temp20;
			delete[] temp21; delete[] temp22; delete[] temp23; delete[] temp24;
			delete[] temp25; delete[] temp26; delete[] temp27; delete[] temp28; delete[] temp29; delete[] temp30; delete[] temp31;
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
	contractParam1x[matCount] = contractionParam1;
	contractParam2x[matCount] = contractionParam2;
	contractParam3x[matCount] = contractionParam3;
	dilateParam1x[matCount] = dilationParam1;
	dilateParam2x[matCount] = dilationParam2;
	volLimit1x[matCount] = volLim1;
	volLimit2x[matCount] = volLim2;
	volLimit3x[matCount] = volLim3;
	liquefyParam1x[matCount] = liquefactionParam1;
	liquefyParam2x[matCount] = liquefactionParam2;
	dilateParam3x[matCount] = dilationParam3;
	einitx[matCount] = ei;
	Hvx[matCount] = hv;
	Pvx[matCount] = pv;
	dilateParam4x[matCount] = dilationParam4;
	contractParam4x[matCount] = contractionParam4;
	contractParam5x[matCount] = contractionParam5;
	liquefyParam3x[matCount] = liquefactionParam3;
	residualPressx[matCount] = 0.;
	stressRatioPTx[matCount] = 0.;
	matN = matCount;
	matCount++;
	pAtm = atm;
	stress0.Zero();
	stress0Initialized = false;
	int numOfSurfaces = numOfSurfacesx[matN];
	initPress = refPressurex[matN];

	e2p = committedActiveSurf = activeSurfaceNum = 0;
	onPPZCommitted = onPPZ = -1;
	PPZSizeCommitted = PPZSize = 0.;
	pressureDCommitted = pressureD = modulusFactor = 0.;
	cumuDilateStrainOctaCommitted = cumuDilateStrainOcta = 0.;
	maxCumuDilateStrainOctaCommitted = maxCumuDilateStrainOcta = 0.;
	cumuTranslateStrainOctaCommitted = cumuTranslateStrainOcta = 0.;
	prePPZStrainOctaCommitted = prePPZStrainOcta = 0.;
	oppoPrePPZStrainOctaCommitted = oppoPrePPZStrainOcta = 0.;
	maxPress = 0.;
	damage = 0.;
	Dilation_Rule = 0;
	Dilation_Rule1 = 0;
	theSurfaces = new MultiYieldSurface[numOfSurfaces + 1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces + 1];
	mGredu = gredu;
	setUpSurfaces(gredu);  // residualPress and stressRatioPT are calculated inside.
}


LadeDuncanMultiYield::LadeDuncanMultiYield()
	: NDMaterial(0, ND_TAG_LadeDuncanMultiYield),
	currentStress(), trialStress(), currentStrain(),
	strainRate(), PPZPivot(), PPZCenter(), PivotStrainRate(6), PivotStrainRateCommitted(6),
	PPZPivotCommitted(), PPZCenterCommitted(), theSurfaces(0), committedSurfaces(0)
{
	//does nothing
}


LadeDuncanMultiYield::LadeDuncanMultiYield(const LadeDuncanMultiYield & a)
	: NDMaterial(a.getTag(), ND_TAG_LadeDuncanMultiYield),
	currentStress(a.currentStress), trialStress(a.trialStress),
	currentStrain(a.currentStrain), strainRate(a.strainRate),
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
	pressureDCommitted = a.pressureDCommitted;
	onPPZCommitted = a.onPPZCommitted;
	PPZSizeCommitted = a.PPZSizeCommitted;
	cumuDilateStrainOctaCommitted = a.cumuDilateStrainOctaCommitted;
	maxCumuDilateStrainOctaCommitted = a.maxCumuDilateStrainOctaCommitted;
	cumuTranslateStrainOctaCommitted = a.cumuTranslateStrainOctaCommitted;
	prePPZStrainOctaCommitted = a.prePPZStrainOctaCommitted;
	oppoPrePPZStrainOctaCommitted = a.oppoPrePPZStrainOctaCommitted;
	pressureD = a.pressureD;
	onPPZ = a.onPPZ;
	PPZSize = a.PPZSize;
	cumuDilateStrainOcta = a.cumuDilateStrainOcta;
	maxCumuDilateStrainOcta = a.maxCumuDilateStrainOcta;
	cumuTranslateStrainOcta = a.cumuTranslateStrainOcta;
	prePPZStrainOcta = a.prePPZStrainOcta;
	oppoPrePPZStrainOcta = a.oppoPrePPZStrainOcta;
	initPress = a.initPress;
	maxPress = a.maxPress;
	damage = a.damage;
	Dilation_Rule = 0;
	Dilation_Rule1 = 0;
	stress0.Zero();
	stress0Initialized = false;
	theSurfaces = new MultiYieldSurface[numOfSurfaces + 1];  //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces + 1];
	for (int i = 1; i <= numOfSurfaces; i++) {
		committedSurfaces[i] = a.committedSurfaces[i];
		theSurfaces[i] = a.theSurfaces[i];
	}
}


LadeDuncanMultiYield::~LadeDuncanMultiYield()
{
	if (theSurfaces != 0) delete[] theSurfaces;
	if (committedSurfaces != 0) delete[] committedSurfaces;
}


void LadeDuncanMultiYield::elast2Plast(void)

{
	int loadStage = loadStagex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double residualPress = residualPressx[matN];
	stress0.Zero();
	if (loadStage != 1 || e2p == 1)
		return;

	e2p = 1;

	if (currentStress.volume() > 0.) {
		currentStress.setData(currentStress.deviator(), 0);
	}

	workV6 = currentStress.deviator();
	double conHeig = currentStress.volume() - residualPress;
	committedActiveSurf++;
	workV6.addVector(1.0, committedSurfaces[committedActiveSurf].center(), -conHeig);
	double rat = stressRatio(workV6, conHeig);

	// Find active surface
	while (rat > committedSurfaces[committedActiveSurf].size()) {
		committedActiveSurf++;
		workV6 = currentStress.deviator();
		workV6.addVector(1.0, committedSurfaces[committedActiveSurf].center(), -conHeig);
		rat = stressRatio(workV6, conHeig);

		if (committedActiveSurf == numOfSurfaces) {
			deviatorScaling(currentStress, committedSurfaces, numOfSurfaces);
			initSurfaceUpdate();
			return;
		}
	}

	committedActiveSurf--;
	initSurfaceUpdate();
}


int LadeDuncanMultiYield::setTrialStrain(const Vector &strain)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3 && strain.Size() == 6)
		workV6 = strain;
	else if (ndm == 2 && strain.Size() == 3) {
		workV6[0] = strain[0];
		workV6[1] = strain[1];
		workV6[2] = 0.0;
		workV6[3] = strain[2];
		workV6[4] = 0.0;
		workV6[5] = 0.0;
	}
	else {
		opserr << "Fatal:LadeDuncanMultiYield:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);
	}

	workV6 -= currentStrain.t2Vector(1);
	strainRate.setData(workV6, 1);

	return 0;
}


int LadeDuncanMultiYield::setTrialStrain(const Vector &strain, const Vector &rate)
{
	return setTrialStrain(strain);
}


int LadeDuncanMultiYield::setTrialStrainIncr(const Vector &strain)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3 && strain.Size() == 6)
		workV6 = strain;
	else if (ndm == 2 && strain.Size() == 3) {
		workV6[0] = strain[0];
		workV6[1] = strain[1];
		workV6[2] = 0.0;
		workV6[3] = strain[2];
		workV6[4] = 0.0;
		workV6[5] = 0.0;
	}
	else {
		opserr << "Fatal:LadeDuncanMultiYield:: Material dimension is: " << ndm << endln;
		opserr << "But strain vector size is: " << strain.Size() << endln;
		exit(-1);
	}

	strainRate.setData(workV6, 1);
	return 0;
}


int LadeDuncanMultiYield::setTrialStrainIncr(const Vector &strain, const Vector &rate)
{
	return setTrialStrainIncr(strain);
}


const Matrix & LadeDuncanMultiYield::getTangent(void)
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
	if (loadStage == 2 && initPress == refPressure)
		initPress = currentStress.volume();

	if (loadStage == 0 || loadStage == 2) {  //linear elastic
		double factor;
		if (loadStage == 0)
			factor = 1.0;
		else {
			factor = (initPress - residualPress) / (refPressure - residualPress);
			if (factor <= 1.e-10) factor = 1.e-10;
			else factor = pow(factor, pressDependCoeff);
			factor = (1.e-10>factor) ? 1.e-10 : factor;
		}
		for (int i = 0;i<6;i++)
			for (int j = 0;j<6;j++) {
				theTangent(i, j) = 0.;
				if (i == j)
					theTangent(i, j) += refShearModulus*factor;
				if (i<3 && j<3 && i == j)
					theTangent(i, j) += refShearModulus*factor;
				if (i<3 && j<3)
					theTangent(i, j) += (refBulkModulus - 2.*refShearModulus / 3.)*factor;
			}
	}
	else {
		double coeff1, coeff2, coeff3, coeff4;
		double factor = getModulusFactor(updatedTrialStress);
		double shearModulus = factor*refShearModulus;
		double bulkModulus = factor*refBulkModulus;

		// volumetric plasticity
		if (Hvx[matN] != 0. && trialStress.volume() <= maxPress
			&& strainRate.volume()<0. && loadStage == 1) {
			double tp = fabs(trialStress.volume() - residualPress);
			bulkModulus = (bulkModulus*Hvx[matN] * pow(tp, Pvx[matN])) / (bulkModulus + Hvx[matN] * pow(tp, Pvx[matN]));
		}

		if (loadStage != 0 && activeSurfaceNum > 0) {
			factor = getModulusFactor(trialStress);
			shearModulus = factor*refShearModulus;
			bulkModulus = factor*refBulkModulus;
			getSurfaceNormal(trialStress, workT2V);
			workV6 = workT2V.deviator();
			double volume = workT2V.volume();
			double Ho = 9.*bulkModulus*volume*volume + 2.*shearModulus*(workV6 && workV6);
			double plastModul = factor*theSurfaces[activeSurfaceNum].modulus();
			coeff1 = 9.*bulkModulus*bulkModulus*volume*volume / (Ho + plastModul);
			coeff2 = 4.*shearModulus*shearModulus / (Ho + plastModul);
		}

		else {
			coeff1 = coeff2 = coeff3 = coeff4 = 0.;
			workV6.Zero();
		}

		for (int i = 0;i<6;i++)
			for (int j = 0;j<6;j++) {
				theTangent(i, j) = -coeff2*workV6[i] * workV6[j];
				if (i == j) theTangent(i, j) += shearModulus;
				if (i<3 && j<3 && i == j) theTangent(i, j) += shearModulus;
				if (i<3 && j<3) theTangent(i, j) += (bulkModulus - 2.*shearModulus / 3. - coeff1);
			}
	}

	if (ndm == 3)
		return theTangent;
	else {
		static Matrix workM(3, 3);
		workM(0, 0) = theTangent(0, 0);
		workM(0, 1) = theTangent(0, 1);
		workM(0, 2) = 0.;

		workM(1, 0) = theTangent(1, 0);
		workM(1, 1) = theTangent(1, 1);
		workM(1, 2) = 0.;

		workM(2, 0) = 0.;
		workM(2, 1) = 0.;
		workM(2, 2) = theTangent(3, 3);

		return workM;
	}
}


const Matrix & LadeDuncanMultiYield::getInitialTangent(void)
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
	if (loadStage == 2 && initPress == refPressure)
		initPress = currentStress.volume();
	double factor;
	if (loadStage == 0)
		factor = 1.;
	else if (loadStage == 2) {
		factor = (initPress - residualPress) / (refPressure - residualPress);
		if (factor <= 1.e-10) factor = 1.e-10;
		else factor = pow(factor, pressDependCoeff);
		factor = (1.e-10>factor) ? 1.e-10 : factor;
	}
	else if (loadStage == 1)
		factor = getModulusFactor(currentStress);

	for (int i = 0;i<6;i++)
		for (int j = 0;j<6;j++) {
			theTangent(i, j) = 0.;
			if (i == j) theTangent(i, j) += refShearModulus*factor;
			if (i<3 && j<3 && i == j) theTangent(i, j) += refShearModulus*factor;
			if (i<3 && j<3) theTangent(i, j) += (refBulkModulus - 2.*refShearModulus / 3.)*factor;
		}

	if (ndm == 3)
		return theTangent;
	else {
		static Matrix workM(3, 3);
		workM(0, 0) = theTangent(0, 0);
		workM(0, 1) = theTangent(0, 1);
		workM(0, 2) = 0.;

		workM(1, 0) = theTangent(1, 0);
		workM(1, 1) = theTangent(1, 1);
		workM(1, 2) = 0.;

		workM(2, 0) = 0.;
		workM(2, 1) = 0.;
		workM(2, 2) = theTangent(3, 3);

		return workM;
	}
}


const Vector & LadeDuncanMultiYield::getStress(void)
{

	int loadStage = loadStagex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 3;

	int i, is;
	if (loadStage == 1 && e2p == 0) {
		initPress = currentStress.volume();
		elast2Plast();
	}

	if (loadStage != 1) {  //linear elastic
		getTangent();
		workV6 = currentStress.t2Vector();
		workV6.addMatrixVector(1.0, theTangent, strainRate.t2Vector(1), 1.0);
		trialStress.setData(workV6);
	}
	else {
		for (i = 1; i <= numOfSurfaces; i++) theSurfaces[i] = committedSurfaces[i];
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

		if (activeSurfaceNum == 0 && !isCrossingNextSurface()) {
			workV6 = currentStrain.t2Vector();
			workV6.addVector(1.0, strainRate.t2Vector(), 1.0);
			trialStrain.setData(workV6);
		}
		else {
			int numSubIncre = setSubStrainRate();

			for (i = 0; i<numSubIncre; i++) {
				workV6 = currentStrain.t2Vector();
				workV6.addVector(1.0, subStrainRate.t2Vector(), (i + 1));
				trialStrain.setData(workV6);

				if (i == 0) {
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
				if (activeSurfaceNum == 0 && !isCrossingNextSurface()) continue;
				if (activeSurfaceNum == 0) activeSurfaceNum++;
				stressCorrection(0);
				updateActiveSurface();

				double refBulkModulus = refBulkModulusx[matN];
				double B = refBulkModulus*modulusFactor;
				pressureD += 3.*subStrainRate.volume()
					- (trialStress.volume() - updatedTrialStress.volume()) / B;
				if (pressureD < 0.) pressureD = 0.;
			}
		}
	}
	if (ndm == 3)
		return trialStress.t2Vector();
	else {
		static Vector workV(3);
		workV[0] = trialStress.t2Vector()[0];
		workV[1] = trialStress.t2Vector()[1];
		workV[2] = trialStress.t2Vector()[3];
		return workV;
	}
}


const Vector & LadeDuncanMultiYield::getStrain(void)
{
	return getCommittedStrain();
}


int LadeDuncanMultiYield::commitState(void)
{
	int loadStage = loadStagex[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	currentStress = trialStress;
	workV6 = currentStrain.t2Vector();
	workV6 += strainRate.t2Vector();
	currentStrain.setData(workV6);

	workV6.Zero();
	strainRate.setData(workV6);

	if (loadStage == 1) {
		committedActiveSurf = activeSurfaceNum;
		for (int i = 1; i <= numOfSurfaces; i++) committedSurfaces[i] = theSurfaces[i];
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


int LadeDuncanMultiYield::revertToLastCommit(void)
{
	return 0;
}


NDMaterial * LadeDuncanMultiYield::getCopy(void)
{
	LadeDuncanMultiYield * copy = new LadeDuncanMultiYield(*this);
	return copy;
}


NDMaterial * LadeDuncanMultiYield::getCopy(const char *code)
{
	if (strcmp(code, "LadeDuncanMultiYield") == 0 || strcmp(code, "PlaneStrain") == 0
		|| strcmp(code, "ThreeDimensional") == 0) {
		LadeDuncanMultiYield * copy = new LadeDuncanMultiYield(*this);
		return copy;
	}

	return 0;
}


const char * LadeDuncanMultiYield::getType(void) const
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	return (ndm == 2) ? "PlaneStrain" : "ThreeDimensional";
}


int LadeDuncanMultiYield::getOrder(void) const
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	return (ndm == 2) ? 3 : 6;
}


int LadeDuncanMultiYield::setParameter(const char **argv, int argc, Parameter &param)
{
	
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	// check for material tag
	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0], "shearModulus") == 0) {
			return param.addObject(10, this);
		}
		else if (strcmp(argv[0], "bulkModulus") == 0) {
			return param.addObject(11, this);
		}
		else if (strcmp(argv[0], "frictionAngle") == 0) {
			return param.addObject(12, this);
		}
		else if (strcmp(argv[0], "cohesion") == 0) {
			return param.addObject(13, this);
		}
	}
	return -1;
}

int LadeDuncanMultiYield::updateParameter(int responseID, Information &info)
{

	if (responseID == 1) {
		loadStagex[matN] = info.theInt;
	}
	else if (responseID == 10) {
		refShearModulusx[matN] = info.theDouble;
	}
	else if (responseID == 11) {
		refBulkModulusx[matN] = info.theDouble;
	}
	else if (responseID == 12) {
		frictionAnglex[matN] = info.theDouble;
		setUpSurfaces(mGredu);
		initSurfaceUpdate();
	}
	else if (responseID == 13) {
		cohesionx[matN] = info.theDouble;
		setUpSurfaces(mGredu);
		initSurfaceUpdate();
	}

	// used by BBarFourNodeQuadUP element
	else if (responseID == 20 && ndmx[matN] == 2)
		ndmx[matN] = 0;

	return 0;
}


int LadeDuncanMultiYield::sendSelf(int commitTag, Channel &theChannel)
{
	
	int loadStage = loadStagex[matN];
	int ndm = ndmx[matN];
	double rho = rhox[matN];
	double residualPress = residualPressx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
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
		opserr << "LadeDuncanMultiYield::sendSelf -- could not send ID\n";
		return res;
	}

	Vector data(73 + numOfSurfaces * 8);
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
	data(36) = hv; //  = Hvx[matN];
	data(37) = Pv; //  = Pvx[matN];

	workV6 = currentStress.t2Vector();
	for (i = 0; i < 6; i++) data(i + 38) = workV6[i];

	workV6 = currentStrain.t2Vector();
	for (i = 0; i < 6; i++) data(i + 44) = workV6[i];

	workV6 = PPZPivotCommitted.t2Vector();
	for (i = 0; i < 6; i++) data(i + 50) = workV6[i];

	workV6 = PPZCenterCommitted.t2Vector();
	for (i = 0; i < 6; i++) data(i + 56) = workV6[i];

	for (i = 0; i < numOfSurfaces; i++) {
		int k = 62 + i * 8;
		data(k) = committedSurfaces[i + 1].size();
		data(k + 1) = committedSurfaces[i + 1].modulus();
		workV6 = committedSurfaces[i + 1].center();
		data(k + 2) = workV6(0);
		data(k + 3) = workV6(1);
		data(k + 4) = workV6(2);
		data(k + 5) = workV6(3);
		data(k + 6) = workV6(4);
		data(k + 7) = workV6(5);
	}
	i = 69 + numOfSurfaces * 8;
	data(i + 2) = dilateParam4x[matN];
	data(i + 3) = contractParam4x[matN];
	data(i + 4) = contractParam5x[matCount];
	data(i + 5) = liquefyParam3x[matN];

	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "LadeDuncanMultiYield::sendSelf -- could not send Vector\n";
		return res;
	}

	return res;
}


int LadeDuncanMultiYield::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int i, res = 0;

	static ID idData(6);
	res += theChannel.recvID(this->getDbTag(), commitTag, idData);
	if (res < 0) {
		opserr << "LadeDuncanMultiYield::recvelf -- could not recv ID\n";

		return res;
	}

	this->setTag(idData(0));
	int numOfSurfaces = idData(1);
	int loadStage = idData(2);
	int ndm = idData(3);
	matN = idData(4);

	int otherMatCount = idData(5);

	Vector data(69 + idData(1) * 8);
	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "LadeDuncanMultiYield::recvSelf -- could not recv Vector\n";
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
	double contractParam3 = data(35); //  = contractionParam3;// = contractParam3x[matN];
	double hv = data(36); // =  hv; //  = Hvx[matN];
	double Pv = data(37); // = Pv; //  = Pvx[matN];


	for (i = 0; i < 6; i++) workV6[i] = data(i + 38);
	currentStress.setData(workV6);

	for (i = 0; i < 6; i++) workV6[i] = data(i + 44);
	currentStrain.setData(workV6);

	for (i = 0; i < 6; i++) workV6[i] = data(i + 50);
	PPZPivotCommitted.setData(workV6);

	for (i = 0; i < 6; i++) workV6[i] = data(i + 56);
	PPZCenterCommitted.setData(workV6);

	if (committedSurfaces != 0) {
		delete[] committedSurfaces;
		delete[] theSurfaces;
	}

	theSurfaces = new MultiYieldSurface[numOfSurfaces + 1]; //first surface not used
	committedSurfaces = new MultiYieldSurface[numOfSurfaces + 1];

	for (i = 0; i < numOfSurfaces; i++) {
		int k = 62 + i * 8;
		workV6(0) = data(k + 2);
		workV6(1) = data(k + 3);
		workV6(2) = data(k + 4);
		workV6(3) = data(k + 5);
		workV6(4) = data(k + 6);
		workV6(5) = data(k + 7);
		committedSurfaces[i + 1].setData(workV6, data(k), data(k + 1));
	}
	i = 69 + numOfSurfaces * 8;
	int mType = data(i + 1);
	double dilateParam4 = data(i + 2);
	double contractParam4 = data(i + 3);
	double contractParam5 = data(i + 4);
	double liquefyParam3 = data(i + 5);

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
		temp28 = dilateParam4x;
		temp29 = contractParam4x;
		temp30 = contractParam5x;
		temp31 = liquefyParam3x;

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

		if (matCount > 0) {
			for (int i = 0; i<matCount; i++) {
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
				dilateParam4x[i] = temp28[i];
				contractParam4x[i] = temp29[i];
				contractParam5x[i] = temp30[i];
				liquefyParam3x[i] = temp31[i];
			}

			delete[] temp1; delete[] temp2; delete[] temp3; delete[] temp4;
			delete[] temp5; delete[] temp6; delete[] temp7; delete[] temp8;
			delete[] temp9; delete[] temp10; delete[] temp11; delete[] temp12;
			delete[] temp13; delete[] temp14; delete[] temp15; delete[] temp16;
			delete[] temp14a;
			delete[] temp17; delete[] temp18; delete[] temp19; delete[] temp20;
			delete[] temp21; delete[] temp22; delete[] temp23; delete[] temp24;
			delete[] temp25; delete[] temp26;delete[] temp28; delete[] temp29; delete[] temp30; delete[] temp31;
		}
		matCount = otherMatCount;
	}

	loadStagex[matN] = loadStage;
	ndmx[matN] = ndm;
	rhox[matN] = rho;
	residualPressx[matN] = residualPress;
	numOfSurfacesx[matN] = numOfSurfaces;
	refPressurex[matN] = refPressure;
	pressDependCoeffx[matN] = pressDependCoeff;
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
	dilateParam4x[matN] = dilateParam4;
	contractParam4x[matN] = contractParam4;
	contractParam5x[matN] = contractParam5;
	liquefyParam3x[matN] = liquefyParam3;

	return res;
}

const Vector & LadeDuncanMultiYield::getCommittedCenter(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3)
		return PPZCenter.t2Vector(1);
	else {
		static Vector workV(3);
		workV6 = PPZCenter.t2Vector(1);
		workV[0] = workV6[0];
		workV[1] = workV6[1];
		workV[2] = workV6[3];
		return workV;
	}
}

const double  &LadeDuncanMultiYield::getCommittedPPZSize(void)
{
	return PPZSize;
}

const Vector & LadeDuncanMultiYield::getCommittedPivot(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3)
		return PPZPivotCommitted.t2Vector(1);
	else {
		static Vector workV(3);
		workV6 = PPZPivotCommitted.t2Vector(1);
		workV[0] = workV6[0];
		workV[1] = workV6[1];
		workV[2] = workV6[3];
		return workV;
	}
}

void LadeDuncanMultiYield::getBackbone(Matrix & bb)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	double refShearModulus = refShearModulusx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	double vol, conHeig, scale, factor, shearModulus, stress1,
		stress2, strain1, strain2, plastModulus, elast_plast, gre;

	for (int k = 0; k<bb.noCols() / 2; k++) {
		vol = bb(0, k * 2);
		if (vol <= 0.) {
			opserr << k << "\nNDMaterial " << this->getTag()
				<< ": invalid confinement for backbone recorder, " << vol << endln;
			continue;
		}
		conHeig = vol + residualPress;
		scale = -conHeig / (refPressure - residualPress);
		factor = pow(scale, pressDependCoeff);
		shearModulus = factor*refShearModulus;

		for (int i = 1; i <= numOfSurfaces; i++) {
			if (i == 1) {
				stress2 = committedSurfaces[i].size()*conHeig / sqrt(3.0);
				strain2 = stress2 / shearModulus;
				bb(1, k * 2) = strain2; bb(1, k * 2 + 1) = shearModulus;
			}
			else {
				stress1 = stress2; strain1 = strain2;
				plastModulus = factor*committedSurfaces[i - 1].modulus();
				elast_plast = 2 * shearModulus*plastModulus / (2 * shearModulus + plastModulus);
				stress2 = committedSurfaces[i].size()*conHeig / sqrt(3.0);
				strain2 = 2 * (stress2 - stress1) / elast_plast + strain1;
				gre = stress2 / strain2;
				bb(i, k * 2) = strain2; bb(i, k * 2 + 1) = gre;
			}
		}
	}

}


int LadeDuncanMultiYield::getResponse(int responseID, Information &matInfo)
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
	case 12:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getCommittedCenter();
		return 0;
	case 14:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getCommittedPivot();
		return 0;
	case 13:
		matInfo.setDouble(this->getCommittedPPZSize());
		return 0;

	case 80:
		matInfo.setDouble(this->Dilation_Rule);
		return 0;
	case 90:
		matInfo.setDouble(this->Dilation_Rule1);
		return 0;
	default:
		return -1;
	}
}


Response*
LadeDuncanMultiYield::setResponse(const char **argv, int argc, OPS_Stream &s)
{
	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getCommittedStress());

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getCommittedStrain());

	else if (strcmp(argv[0], "tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());

	else if (strcmp(argv[0], "PPZPivot") == 0 || strcmp(argv[0], "PPZPivots") == 0)
		return new MaterialResponse(this, 14, this->getCommittedPivot());

	else if (strcmp(argv[0], "PPZCenter") == 0 || strcmp(argv[0], "PPZCenters") == 0)
		return new MaterialResponse(this, 12, this->getCommittedCenter());

	else if (strcmp(argv[0], "PPZSize") == 0)
		return new MaterialResponse(this, 13, this->getCommittedPPZSize());

	else if (strcmp(argv[0], "DD") == 0)
		return new MaterialResponse(this, 80, Dilation_Rule);
	else if (strcmp(argv[0], "SS") == 0)
		return new MaterialResponse(this, 90, Dilation_Rule1);

	else if (strcmp(argv[0], "backbone") == 0) {
		int numOfSurfaces = numOfSurfacesx[matN];
		Matrix curv(numOfSurfaces + 1, (argc - 1) * 2);
		for (int i = 1; i<argc; i++)
			curv(0, (i - 1) * 2) = atoi(argv[i]);
		return new MaterialResponse(this, 4, curv);
	}
	else
		return 0;
}


void LadeDuncanMultiYield::Print(OPS_Stream &s, int flag)
{
	s << "LadeDuncanMultiYield" << endln;
}


const Vector & LadeDuncanMultiYield::getCommittedStress(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;
	int numOfSurfaces = numOfSurfacesx[matN];
	double residualPress = residualPressx[matN];
	double conHeig = currentStress.volume() - residualPress;
	double scale = stressRatio(currentStress.deviator(), conHeig);

	scale /= committedSurfaces[numOfSurfaces].size();
	if (loadStagex[matN] != 1) scale = 0.;

	if (ndm == 3) {
		static Vector temp7(7);
		workV6 = currentStress.t2Vector();
		temp7[0] = workV6[0];
		temp7[1] = workV6[1];
		temp7[2] = workV6[2];
		temp7[3] = workV6[3];
		temp7[4] = workV6[4];
		temp7[5] = workV6[5];
		temp7[6] = scale;
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
		return temp5;
	}
}

const Vector & LadeDuncanMultiYield::getCommittedStrain(void)
{
	int ndm = ndmx[matN];
	if (ndmx[matN] == 0) ndm = 2;

	if (ndm == 3)
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
void LadeDuncanMultiYield::setUpSurfaces(double * gredu)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	double refShearModulus = refShearModulusx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];
	double frictionAngle = frictionAnglex[matN];
	double cohesion = cohesionx[matN];
	double peakShearStrain = peakShearStrainx[matN];
	double phaseTransfAngle = phaseTransfAnglex[matN];
	double stressRatioPT = stressRatioPTx[matN];
	double a1 = a1x[matN];

	double refStrain, peakShear, coneHeight;
	double stress1, stress2, strain1, strain2, size, elasto_plast_modul, plast_modul;
	double ratio1, ratio2;


	if (gredu == 0) {
		double sinPhi = sin(frictionAngle * PI / 180.);
		double temp = 4.*sinPhi / 3. / (3. - sinPhi);
		a1 = (temp*temp - temp*temp*temp) / 4.;
		double Mnys = 6.*sinPhi / (3. - sinPhi);
		double sinPhiPT = sin(phaseTransfAngle * PI / 180.);
		stressRatioPT = 6.*sinPhiPT / (3. - sinPhiPT);

		residualPress = 2 * cohesion / Mnys;
		if (residualPress < 0.0001*pAtm) residualPress = 0.0001*pAtm;

		coneHeight = -(refPressure - residualPress);
		peakShear = sqrt(2.) * coneHeight * Mnys / 3.;
		refStrain = (peakShearStrain * peakShear)
			/ (refShearModulus * peakShearStrain - peakShear);

		double stressInc = peakShear / numOfSurfaces;

		for (int ii = 1; ii <= numOfSurfaces; ii++) {
			stress1 = ii * stressInc;
			stress2 = stress1 + stressInc;

			ratio1 = 3. * stress1 / sqrt(2.) / coneHeight;
			ratio2 = 3. * stress2 / sqrt(2.) / coneHeight;

			strain1 = stress1 * refStrain / (refShearModulus * refStrain - stress1);
			strain2 = stress2 * refStrain / (refShearModulus * refStrain - stress2);

			if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
				double ratio = (ratio2 - stressRatioPT) / (ratio2 - ratio1);
				strainPTOcta = strain2 - ratio * (strain2 - strain1);
			}

			size = ratio1 / Mnys;
			elasto_plast_modul = 2.*(stress2 - stress1) / (strain2 - strain1);

			if ((2.*refShearModulus - elasto_plast_modul) <= 0)
				plast_modul = UP_LIMIT;
			else
				plast_modul = (2.*refShearModulus * elasto_plast_modul) /
				(2.*refShearModulus - elasto_plast_modul);

			if (plast_modul < 0) plast_modul = 0;
			if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;
			if (ii == numOfSurfaces) plast_modul = 0;

			workV6.Zero();
			committedSurfaces[ii] = MultiYieldSurface(workV6, size, plast_modul);
		}  // ii

		stressRatioPT /= Mnys;
	}

	else {  //user defined surfaces
		int ii = 2 * (numOfSurfaces - 1);
		double tmax = refShearModulus*gredu[ii] * gredu[ii + 1];
		double Mnys = -(sqrt(3.) * tmax - 2.* cohesion) / refPressure;
		residualPress = 2 * cohesion / Mnys;
		if (residualPress < 0.0001*pAtm) residualPress = 0.0001*pAtm;
		coneHeight = -(refPressure - residualPress);
		double sinPhi = 3 * Mnys / (6 + Mnys);

		if (sinPhi<0. || sinPhi>1.) {
			opserr << "\nNDMaterial " << this->getTag() << ": Invalid friction angle, please modify ref. pressure or G/Gmax curve." << endln;
			exit(-1);
		}

		frictionAngle = asin(sinPhi) * 180 / PI;
		double temp = 4.*sinPhi / 3. / (3. - sinPhi);
		a1 = (temp*temp - temp*temp*temp) / 4.;

		opserr << "\nNDMaterial " << this->getTag() << ": Friction angle is " << frictionAngle << "\n" << endln;

		if (phaseTransfAngle > frictionAngle) {
			opserr << "\nNDMaterial " << this->getTag() << ": phase Transformation Angle > friction Angle,"
				<< "will set phase Transformation Angle = friction Angle.\n" << endln;
			phaseTransfAngle = frictionAngle;
		}

		double sinPhiPT = sin(phaseTransfAngle * PI / 180.);
		stressRatioPT = 6.*sinPhiPT / (3. - sinPhiPT);

		for (int i = 1; i<numOfSurfaces; i++) {
			int ii = 2 * (i - 1);
			strain1 = gredu[ii];
			stress1 = refShearModulus*gredu[ii + 1] * strain1;
			strain2 = gredu[ii + 2];
			stress2 = refShearModulus*gredu[ii + 3] * strain2;

			ratio1 = sqrt(3.) * stress1 / coneHeight;
			ratio2 = sqrt(3.) * stress2 / coneHeight;
			if (ratio1 <= stressRatioPT && ratio2 >= stressRatioPT) {
				double ratio = (ratio2 - stressRatioPT) / (ratio2 - ratio1);
				strainPTOcta = sqrt(6.) / 3 * (strain2 - ratio * (strain2 - strain1));
			}

			size = ratio1 / Mnys;
			elasto_plast_modul = 2.*(stress2 - stress1) / (strain2 - strain1);

			if ((2.*refShearModulus - elasto_plast_modul) <= 0)
				plast_modul = UP_LIMIT;
			else
				plast_modul = (2.*refShearModulus * elasto_plast_modul) /
				(2.*refShearModulus - elasto_plast_modul);

			if (plast_modul <= 0) {
				opserr << "\nNDMaterial " << this->getTag() << ": Surface " << i
					<< " has plastic modulus < 0.\n Please modify G/Gmax curve.\n" << endln;
				exit(-1);
			}

			if (plast_modul > UP_LIMIT) plast_modul = UP_LIMIT;

			workV6.Zero();
			committedSurfaces[i] = MultiYieldSurface(workV6, size, plast_modul);

			if (i == (numOfSurfaces - 1)) {
				plast_modul = 0;
				size = ratio2 / Mnys;
				committedSurfaces[i + 1] = MultiYieldSurface(workV6, size, plast_modul);
			}
		}

		stressRatioPT /= Mnys;
	}

	residualPressx[matN] = residualPress;
	frictionAnglex[matN] = frictionAngle;
	cohesionx[matN] = cohesion;
	phaseTransfAnglex[matN] = phaseTransfAngle;
	stressRatioPTx[matN] = stressRatioPT;
	a1x[matN] = a1;
}


void LadeDuncanMultiYield::deviatorScaling(T2Vector & stress,
	const MultiYieldSurface * surfaces,
	int surfaceNum)
{
	double residualPress = residualPressx[matN];

	int numOfSurfaces = numOfSurfacesx[matN];
	double coneHeight = stress.volume() - residualPress;

	workV6 = stress.deviator();
	workV6.addVector(1.0, surfaces[surfaceNum].center(), -coneHeight);
	double rat = stressRatio(workV6, coneHeight);
	double sz = surfaces[surfaceNum].size();

	double coeff = sz / rat - 1.;

	// make sure stress point is NOT inside the surface
	if (surfaceNum < numOfSurfaces && rat < sz) {
		if (coeff < 1.e-13) coeff = 1.e-13;
		workV6 *= (coeff + 1);
		workV6.addVector(1.0, surfaces[surfaceNum].center(), coneHeight);
		stress.setData(workV6, stress.volume());
		deviatorScaling(stress, surfaces, surfaceNum);  // recursive call
	}


	// make sure stress point is on last surface if it's active
	if (surfaceNum == numOfSurfaces && fabs(coeff) > LOW_LIMIT) {
		workV6 *= (coeff + 1);
		stress.setData(workV6, stress.volume());
	}
}


void LadeDuncanMultiYield::initSurfaceUpdate(void)
{
	double residualPress = residualPressx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	if (committedActiveSurf == 0) return;

	double coneHeight = currentStress.volume() - residualPress;
	workV6 = currentStress.deviator();
	double rat = stressRatio(workV6, coneHeight);

	if (committedActiveSurf < numOfSurfaces) { // failure surface can't move
		workV6 *= 1. - committedSurfaces[committedActiveSurf].size() / rat;
		workV6 /= coneHeight;
		committedSurfaces[committedActiveSurf].setCenter(workV6);
		theSurfaces[committedActiveSurf] = committedSurfaces[committedActiveSurf];
	}

	for (int i = 1; i<committedActiveSurf; i++) {
		workV6 = currentStress.deviator();
		workV6 *= 1. - committedSurfaces[i].size() / rat;
		workV6 /= coneHeight;
		committedSurfaces[i].setCenter(workV6);
		theSurfaces[i] = committedSurfaces[i];
	}


	activeSurfaceNum = committedActiveSurf;
}


void LadeDuncanMultiYield::initStrainUpdate(void)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];
	double stressRatioPT = stressRatioPTx[matN];

	// elastic strain state
	double conHeig = currentStress.volume() - residualPress;
	double rat = stressRatio(currentStress.deviator(), conHeig);
	double ratio = (-currentStress.volume() + residualPress) / (-refPressure + residualPress);
	ratio = pow(ratio, 1. - pressDependCoeff);
	modulusFactor = getModulusFactor(currentStress);
	double shearCoeff = 1. / (2.*refShearModulus*modulusFactor);
	double bulkCoeff = 1. / (3.*refBulkModulus*modulusFactor);

	workV6.addVector(0.0, currentStress.deviator(), shearCoeff);
	currentStrain.setData(workV6, currentStress.volume()*bulkCoeff);

	double octalStrain = currentStrain.octahedralShear(1);
	if (octalStrain <= LOW_LIMIT) octalStrain = LOW_LIMIT;

	// plastic strain state (scaled from elastic strain)
	double scale, PPZLimit;

	if (rat >= stressRatioPT) {  //above PT
		onPPZCommitted = 2;
		prePPZStrainOctaCommitted = strainPTOcta * ratio;
		PPZLimit = getPPZLimits(1, currentStress);
		scale = sqrt(prePPZStrainOctaCommitted + PPZLimit) / octalStrain;
	}
	else {  // below PT
		onPPZCommitted = -1;
		prePPZStrainOctaCommitted = octalStrain;
		if (prePPZStrainOctaCommitted > strainPTOcta * ratio)
			prePPZStrainOctaCommitted = strainPTOcta*ratio;
		scale = sqrt(prePPZStrainOctaCommitted) / octalStrain;
	}

	workV6.addVector(0.0, currentStrain.deviator(), scale);
	currentStrain.setData(workV6, currentStrain.volume());
	PPZPivotCommitted = currentStrain;
}


double LadeDuncanMultiYield::getModulusFactor(T2Vector & stress)
{
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];

	double conHeig = stress.volume() - residualPress;
	double scale = conHeig / (refPressure - residualPress);
	scale = pow(scale, pressDependCoeff);

	return (1.e-10>scale) ? 1.e-10 : scale;
}


void LadeDuncanMultiYield::setTrialStress(T2Vector & stress)
{
	double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];

	modulusFactor = getModulusFactor(stress);
	workV6 = stress.deviator();
	workV6.addVector(1.0, subStrainRate.deviator(), 2 * refShearModulus*modulusFactor);

	double B = refBulkModulus*modulusFactor;

	if (Hvx[matN] != 0. && trialStress.volume() <= maxPress
		&& subStrainRate.volume()<0. && loadStagex[matN] == 1) {
		double tp = fabs(trialStress.volume() - residualPressx[matN]);
		B = (B*Hvx[matN] * pow(tp, Pvx[matN])) / (B + Hvx[matN] * pow(tp, Pvx[matN]));
	}

	double volume = stress.volume() + subStrainRate.volume()*3.*B;
	if (volume > 0.) volume = 0.;
	trialStress.setData(workV6, volume);
}


int LadeDuncanMultiYield::setSubStrainRate(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	int numOfSub = strainRate.octahedralShear(1) / 1.0e-4 + 1;
	int numOfSub1 = strainRate.volume() / 1.0e-4 + 1;

	if (numOfSub1 > numOfSub) numOfSub = numOfSub1;
	if (numOfSub > numOfSurfaces) numOfSub = 2*numOfSurfaces;

	workV6.addVector(0.0, strainRate.t2Vector(), 1.0 / numOfSub);
	subStrainRate.setData(workV6);
	return numOfSub;
}


void
LadeDuncanMultiYield::getContactStress(T2Vector &contactStress)
{
	double residualPress = residualPressx[matN];
	double conHeig = trialStress.volume() - residualPress;
	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	workV6 = trialStress.deviator();
	workV6.addVector(1.0, center, -conHeig);
	double rat = stressRatio(workV6, conHeig);
	workV6.addVector(theSurfaces[activeSurfaceNum].size() / rat, center, conHeig);
	contactStress.setData(workV6, trialStress.volume());
}


int LadeDuncanMultiYield::isLoadReversal(const T2Vector & stress)
{
	if (activeSurfaceNum == 0) return 0;

	getSurfaceNormal(stress, workT2V);
	workV6 = trialStress.t2Vector();
	workV6 -= currentStress.t2Vector();

	if ((workV6 && workT2V.t2Vector()) < 0) return 1;

	return 0;
}


void
LadeDuncanMultiYield::getSurfaceNormal(const T2Vector & stress, T2Vector &normal)
{

	double a1 = a1x[matN];
	double residualPress = residualPressx[matN];
	double conHeig = stress.volume() - residualPress;
	Vector sds(6);
	Vector sds_full(6);

	Vector unitVector(6);
	unitVector.Zero();
	for (int i = 0; i<3; i++) {
		unitVector(i) = 1.0;
	}

	T2Vector norm;
	  
	workV6 = stress.deviator();

	static Vector center(6);
	center = theSurfaces[activeSurfaceNum].center();
	double sz = theSurfaces[activeSurfaceNum].size();


	workV6.addVector(1.0, center, -conHeig);
	double ddspa = workV6 && workV6;
	if (ddspa == 0.) {
		opserr << "In CALQ: zero vector length, no normal exists" << endln;
		workV6.Zero();
		normal.setData(workV6, 0.);
		return;
	}

	sds(0) = workV6(0)*workV6(0) + workV6(3)*workV6(3) + workV6(5)*workV6(5);
	sds(1) = workV6(1)*workV6(1) + workV6(3)*workV6(3) + workV6(4)*workV6(4);
	sds(2) = workV6(2)*workV6(2) + workV6(4)*workV6(4) + workV6(5)*workV6(5);
	sds(3) = workV6(0)*workV6(3) + workV6(1)*workV6(3) + workV6(4)*workV6(5);
	sds(4) = workV6(1)*workV6(4) + workV6(2)*workV6(4) + workV6(3)*workV6(5);
	sds(5) = workV6(2)*workV6(5) + workV6(0)*workV6(5) + workV6(3)*workV6(4);

	double dsda = sds && center;
	double ddsa = workV6 && center;
	double aux = (-sz - 2.)*ddspa / 6. - dsda / 3. + sz*conHeig*ddsa / 3.
		+ 3.*a1*pow(sz, 3.)*pow(3.*conHeig, 2.);
	sds.addVector(1.0, workV6, -conHeig*sz);

	sds_full = sds;
	sds_full.addVector(1.0, unitVector, aux);
	norm.setData(sds_full, 0);    // Zero is correct
	normal.setData(norm.unitT2Vector());
	return;
}

double LadeDuncanMultiYield::getPlasticPotential(const T2Vector & contactStress, const T2Vector & surfaceNormal)
{
	double residualPress = residualPressx[matN];
	double stressRatioPT = stressRatioPTx[matN];
	double contractParam1 = contractParam1x[matN];
	double contractParam2 = contractParam2x[matN];
	double contractParam3 = contractParam3x[matN];
	double dilateParam1 = dilateParam1x[matN];
	double dilateParam2 = dilateParam2x[matN];
	double C = ContractionFactorC();

	Dilation_Rule1 = 0;
	double plasticPotential, contractRule, shearLoading, angle;
	double dilateParam4 = dilateParam4x[matN];
	double contractParam4 = contractParam4x[matN];
	double contractParam5 = contractParam5x[matN];

	double conHeiC = contactStress.volume() - residualPress;
	double conHeiS = updatedTrialStress.volume() - residualPress;
	double conHeiT = trialStress.volume() - residualPress;
	double ratc = stressRatio(contactStress.deviator(), conHeiC);
	double rats = stressRatio(updatedTrialStress.deviator(), conHeiS);
	double ratt = stressRatio(trialStress.deviator(), conHeiT);

	double factorPT = ratc / stressRatioPT;
	shearLoading = trialStress.deviator() && updatedTrialStress.deviator();
	double factorPT1 = factorPT;
	
	if (factorPT >= 1.&& ratt >= rats && shearLoading >= 0.) {  //dilation
		updatePPZ(contactStress);
		if (onPPZ == 1)
			plasticPotential = 0.;
		else if (onPPZ == 2) {
			double dilateParam3 = dilateParam3x[matN];
			double ppp = pow((fabs(contactStress.volume()) + fabs(residualPress)) / pAtm, -dilateParam3);
			plasticPotential = ppp*(dilateParam1 + pow(cumuDilateStrainOcta, dilateParam2))*pow(C, dilateParam4);
			if (plasticPotential < 0.) plasticPotential = -plasticPotential;
			if (plasticPotential>5.0e4) plasticPotential = 5.0e4;
		}
		else {
			std::cout << "FATAL: Wrong onPPZ value: " << onPPZ << std::endl;
			exit(-1);
		}
	}

	else {  //contraction
			contractRule = pow((fabs(contactStress.volume()) + fabs(residualPress)) / pAtm, contractParam3);
			if (contractRule < 0.1) contractRule = 0.1;
			plasticPotential = -contractParam4 *pow(C, contractParam5)*(contractParam1 + maxCumuDilateStrainOcta*contractParam2)*contractRule;
			if (plasticPotential > 0.) plasticPotential = -plasticPotential;
			if (onPPZ > 0) onPPZ = 0;
			if (onPPZ != -1) PPZTranslation(contactStress);
			if (isCriticalState(contactStress)) plasticPotential = 0;
	}

	Dilation_Rule = plasticPotential;
	if (isCriticalState(contactStress)) plasticPotential = 0;
	return plasticPotential;
}


double LadeDuncanMultiYield::ContractionFactorC(void)
{
	double residualPress = residualPressx[matN];
	Vector stress = currentStress.t2Vector();
	double sig11 = stress[0];
	double sig22 = stress[1];
	double sig33 = stress[2];
	double tau12 = stress[3];
	double tau23 = stress[4];
	double tau13 = stress[5];
	double factorP = fabs(sig11 + sig22 + sig33) / 3 + fabs(residualPress);
	double factorCSR = pow(((sig11 - sig22) / 2.0)*((sig11 - sig22) / 2.0) + ((sig22 - sig33) / 2.0)*((sig22 - sig33) / 2.0) + ((sig11 - sig33) / 2.0)*((sig11 - sig33) / 2.0) + 6 * tau12*tau12 + 6 * tau23*tau23 + 6 * tau13*tau13, 0.5) / factorP / 3.0;
	return factorCSR;
}


int LadeDuncanMultiYield::isCriticalState(const T2Vector & stress)
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
		ecr1 = volLimit1 - volLimit2*pow(fabs(-stress.volume() / pAtm), volLimit3);
		ecr2 = volLimit1 - volLimit2*pow(fabs(-updatedTrialStress.volume() / pAtm), volLimit3);
	}
	else {
		ecr1 = volLimit1 - volLimit2*log(fabs(-stress.volume() / pAtm));
		ecr2 = volLimit1 - volLimit2*log(fabs(-updatedTrialStress.volume() / pAtm));
	}

	if (ecurr < ecr2 && etria < ecr1) return 0;
	if (ecurr > ecr2 && etria > ecr1) return 0;

	return 1;
}


void LadeDuncanMultiYield::updatePPZ(const T2Vector & contactStress)
{
	double liquefyParam1 = liquefyParam1x[matN];
	double residualPress = residualPressx[matN];
	double refPressure = refPressurex[matN];
	double pressDependCoeff = pressDependCoeffx[matN];
	double liquefyParam2 = liquefyParam2x[matN];
	double qq = ContractionFactorC();
	// onPPZ=-1 may not be needed. can start with onPPZ=0  ****

	double temp;
	temp = strainRate.deviator() && PivotStrainRateCommitted;
	double liquefyParam3 = liquefyParam3x[matN];
	if (onPPZ < 1) {
		damage = 0.0;
		if ((maxPress - currentStress.volume()) / (maxPress - residualPress) > 0.)
			damage = pow((maxPress - currentStress.volume()) / (maxPress - residualPress), 0.25);
	}

	// PPZ inactive if liquefyParam1==0.
	if (liquefyParam1 == 0. || (onPPZ < 1 && damage < 0.)) {
		if (onPPZ == 2) {
			PPZPivot = trialStrain;
			cumuDilateStrainOcta += subStrainRate.octahedralShear(1);
		}
		else if (onPPZ != 2) {
			onPPZ = 2;
			PPZPivot = trialStrain;
			PivotStrainRate = strainRate.deviator();
			if (temp < 0.) cumuDilateStrainOcta = 0.;
		}
		return;
	}

	// dilation: calc. cumulated dilative strain

	if (onPPZ == 2) {
		PPZPivot = trialStrain;
		cumuDilateStrainOcta += subStrainRate.octahedralShear(1);
		double zzz = 0.;
		if (damage > zzz) zzz = damage;
		maxCumuDilateStrainOcta += zzz*liquefyParam1*subStrainRate.octahedralShear(1)*pow(qq, liquefyParam3);
		return;
	}


	if (onPPZ == -1 || onPPZ == 0) {
		if (temp < 0.) {
			double volume = -contactStress.volume();
			oppoPrePPZStrainOcta = prePPZStrainOcta;
			double ratio = (volume + residualPress) / (-refPressure + residualPress);
			ratio = pow(ratio, 1. - pressDependCoeff);
			prePPZStrainOcta = ratio * strainPTOcta;
			if (oppoPrePPZStrainOcta == 0.) oppoPrePPZStrainOcta = prePPZStrainOcta;
		}
	}
	PPZSize = (cumuTranslateStrainOcta + maxCumuDilateStrainOcta) / 2.;

	// calc. new PPZ center.
	if (onPPZ == 0 || (onPPZ == 1 && temp < 0.0)) {

		workV6 = PPZPivot.t2Vector();
		workV6.addVector(1.0, PPZCenter.t2Vector(), -1.);
		workT2V.setData(workV6);

		double coeff;
		if (workT2V.octahedralShear(1) == 0.) coeff = 0.;
		else coeff = (PPZSize - cumuTranslateStrainOcta) / workT2V.octahedralShear(1);
		workV6 = PPZPivot.t2Vector();
		workV6.addVector(1.0, workT2V.t2Vector(), -coeff);
		PPZCenter.setData(workV6);
	}

	workV6 = trialStrain.t2Vector();
	workV6.addVector(1.0, PPZCenter.t2Vector(), -1.);
	workT2V.setData(workV6);


		//opserr << "trialStrain is " << trialStrain.t2Vector()[0] << "  " << trialStrain.t2Vector()[1] << "  " << trialStrain.t2Vector()[2] << "  " << "  " << trialStrain.t2Vector()[3] << "  " << trialStrain.t2Vector()[4] << "  " << trialStrain.t2Vector()[5] << endln;
		//opserr << "PPZCenter is " << PPZCenter.t2Vector()[0] << "  " << PPZCenter.t2Vector()[1] << "  " << PPZCenter.t2Vector()[2] << "  " << "  " << PPZCenter.t2Vector()[3] << "  " << PPZCenter.t2Vector()[4] << "  " << PPZCenter.t2Vector()[5] << endln;
		//opserr << "trialStrain-PPZCenter is " << workT2V.t2Vector()[0] << "  " << workT2V.t2Vector()[1] << "  " << workT2V.t2Vector()[2] << "  " << "  " << workT2V.t2Vector()[3] << "  " << workT2V.t2Vector()[4] << "  " << workT2V.t2Vector()[5] << endln;
		//opserr << "octahedralShear is " << workT2V.octahedralShear(1) << "   " << "PPZSize  is " << PPZSize << endln;
		//opserr << "===============================================" << endln;

	//outside PPZ
	if (workT2V.octahedralShear(1) > PPZSize) {


		Dilation_Rule1 = 1;
		cumuDilateStrainOcta = 0.;
		onPPZ = 2;
		PPZPivot = trialStrain;
		PivotStrainRate = strainRate.deviator();
		cumuTranslateStrainOcta = 0.;
	}
	else {  //inside PPZ
		if (onPPZ == 0 || onPPZ == 1) PPZTranslation(contactStress);
		if (onPPZ == -1 || onPPZ == 0) onPPZ = 1;
	}
}


void LadeDuncanMultiYield::PPZTranslation(const T2Vector & contactStress)
{
	double liquefyParam1 = liquefyParam1x[matN];
	double liquefyParam2 = liquefyParam2x[matN];
	double residualPress = residualPressx[matN];

	if (liquefyParam1 == 0.) return;
	damage = 0.0;
	if ((maxPress - currentStress.volume()) / (maxPress - residualPress) > 0.)
		damage = pow((maxPress - currentStress.volume()) / (maxPress - residualPress), 0.25);

	double zzz = 0.;
	if (damage > zzz) zzz = damage;

	double temp = strainRate.deviator() && PivotStrainRateCommitted;

	if (temp < 0.0) {  //update only when load reverses
		workV6 = trialStrain.deviator();
		workV6 -= PPZPivot.deviator();
		workT2V.setData(workV6);

		temp = workT2V.octahedralShear(1);
		if (cumuTranslateStrainOcta < zzz*liquefyParam2*temp)
			cumuTranslateStrainOcta = zzz*liquefyParam2*temp;
	}
}


double LadeDuncanMultiYield::getPPZLimits(int which, const T2Vector & contactStress)
{
	double liquefyParam1 = liquefyParam1x[matN];
	double liquefyParam2 = liquefyParam2x[matN];
	double dilateParam3 = dilateParam3x[matN];

	double PPZLimit, temp;
	double volume = -contactStress.volume();

	if (volume >= liquefyParam1) PPZLimit = 0.;
	else {
		temp = volume*PI / liquefyParam1 / 2.;
		PPZLimit = liquefyParam2 * pow(cos(temp), 3.);
		PPZLimit = 0.0;
	}

	if (which == 1)
		return PPZLimit;
	else if (which == 2)
		return dilateParam3 * PPZLimit;
	else {
		opserr << "FATAL:LadeDuncanMultiYield::getPPZLimits: unknown argument value" << endln;
		exit(-1);
		return 0.0;
	}
}


double LadeDuncanMultiYield::getLoadingFunc(const T2Vector & contactStress,
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
		* (surfaceNormal.deviator() && surfaceNormal.deviator());
	double temp2 = 9. * refBulkModulus * modulusFactor
		* surfaceNormal.volume() * (*plasticPotential);

	//for the first crossing
	double temp = temp1 + temp2 + modul * modulusFactor;
	if (activeSurfaceNum == numOfSurfaces)
		limit = theSurfaces[activeSurfaceNum - 1].modulus() * modulusFactor / 2.;
	else limit = modul * modulusFactor / 2.;
	if (temp < limit) {
		(*plasticPotential) = (temp2 + limit - temp) / (9. * refBulkModulus * modulusFactor
			* surfaceNormal.volume());
		temp = limit;
	}

	workV6 = trialStress.deviator();
	workV6 -= contactStress.deviator();
	loadingFunc = (surfaceNormal.t2Vector() && workV6) / temp;

	if (loadingFunc < 0.) loadingFunc = 0;

	//for more than one crossing
	if (crossedSurface) {
		temp = (theSurfaces[activeSurfaceNum - 1].modulus() - modul)
			/ theSurfaces[activeSurfaceNum - 1].modulus();
		loadingFunc *= temp;
	}

	return loadingFunc;
}


int LadeDuncanMultiYield::stressCorrection(int crossedSurface)
{
	double refShearModulus = refShearModulusx[matN];
	double refBulkModulus = refBulkModulusx[matN];

	static T2Vector contactStress;
	getContactStress(contactStress);
	static T2Vector surfNormal;
	getSurfaceNormal(contactStress, surfNormal);
	double plasticPotential = getPlasticPotential(contactStress, surfNormal);
	double tVolume = trialStress.volume();
	double loadingFunc = getLoadingFunc(contactStress, surfNormal,
		&plasticPotential, crossedSurface);
	double volume = tVolume - plasticPotential * 3 * refBulkModulus*modulusFactor*loadingFunc;

	workV6 = trialStress.deviator();

	if (volume > 0. && volume != tVolume) {
		double coeff = tVolume / (tVolume - volume);
		coeff *= -2 * refShearModulus*modulusFactor*loadingFunc;
		workV6.addVector(1.0, surfNormal.deviator(), coeff);
		volume = 0.;
	}
	else if (volume > 0.) {
		volume = 0.;
	}
	else {
		double coeff = -2 * refShearModulus*modulusFactor*loadingFunc;
		workV6.addVector(1.0, surfNormal.deviator(), coeff);
	} 

	trialStress.setData(workV6, volume);
	deviatorScaling(trialStress, theSurfaces, activeSurfaceNum);

	if (isCrossingNextSurface()) {
		activeSurfaceNum++;
		return stressCorrection(1);  //recursive call
	}

	return 0;
}

void LadeDuncanMultiYield::updateActiveSurface(void)
{
	double residualPress = residualPressx[matN];
	int numOfSurfaces = numOfSurfacesx[matN];

	double a1 = a1x[matN];

	if (activeSurfaceNum == numOfSurfaces) return;

	double bound = 1.e-5;
	double A, B, C, D, X;

	static Vector t1(6);
	static Vector t2(6);
	static Vector center(6);
	static Vector outcenter(6);

	// Step 1: Find ST==================================
	double conHeig = trialStress.volume() - residualPress;
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();
	outcenter = theSurfaces[activeSurfaceNum + 1].center();
	double outsize = theSurfaces[activeSurfaceNum + 1].size();

	t1 = trialStress.deviator();
	t1.addVector(1.0, center, -conHeig);
	double rat = stressRatio(t1, conHeig);

	t2 = center;
	t2 -= outcenter;
	t2 *= conHeig;

	A = TRIPRO(t1, t1, t1, 1);
	B = TRIPRO(t2, t1, t1, 0);
	C = TRIPRO(t2, t2, t1, 0);
	D = TRIPRO(t2, t2, t2, 1);
	A = A / 3.;
	D = D / 3.;

	double t11 = t1 && t1;
	double t12 = t1 && t2;
	double t22 = t2 && t2;
	double ppm = conHeig * outsize;
	B = B - ppm*t11 / 2.;
	C = C - ppm*t12;
	D = D - ppm*t22 / 2. + a1*pow(3.*ppm, 3);

	if (abs(A) >= bound) {
		X = thirdOrderEqn(A, B, C, D, 0); 
		t1 *= X;
		t1.addVector(1.0, center, conHeig);  //ST
		
		//t1.addVector(1.0, outcenter, -conHeig); This should match the backbone
		//rat = stressRatio(t1, conHeig);         This should match the backbone
	}
	else {
		t1 = trialStress.deviator();
		rat = stressRatio(t1, conHeig);
	}
	// Step 1: Find ST==================================


	// Step 2: Find direction U ==================================
	t2 = t1;
	t2 *= (1. - size / outsize);
	t2.addVector(1., center, -conHeig);
	t2.addVector(1., outcenter, conHeig*size / outsize);

	double ddv = t2 && t2;
	if (ddv<LOW_LIMIT) ddv = LOW_LIMIT;
	t2 /= sqrt(ddv);    /// This is moving direction
	// Step 2: Find direction U ==================================

	

	//Step 3: Move active surface ==================================

	t1 = trialStress.deviator();
	t1.addVector(1., center, -conHeig);
	D = A;
	A = TRIPRO(t2, t2, t2, 1);
	B = TRIPRO(t2, t2, t1, 0);
	C = TRIPRO(t2, t1, t1, 0);

	t12 = t1 && t2;
	t11 = t1 && t1;
	t22 = t2 && t2;

	ppm = conHeig * size;
	A = -A / 3.;
	B = B - t22*ppm / 2.;
	C = -C + ppm*t12;
	D = D - ppm*t11 / 2. + a1*pow(3.*ppm, 3);

	if (abs(A) > bound) {
		X = thirdOrderEqn(A, B, C, D, 1);
		t1 = center;
		t1.addVector(1., t2, X / conHeig);
		theSurfaces[activeSurfaceNum].setCenter(t1);
		t1.addVector(-conHeig, trialStress.deviator(), 1.);
		ppm = stressRatio(t1, conHeig);
	} 
	else
		ppm = size + 99.99;
	//Step 3: Move active surface ==================================

	return;
}



void LadeDuncanMultiYield::updateInnerSurface(void)
{
	double residualPress = residualPressx[matN];

	if (activeSurfaceNum <= 1) return;
	static Vector devia(6);
	static Vector center(6);

	double conHeig = currentStress.volume() - residualPress;
	devia = currentStress.deviator();
	center = theSurfaces[activeSurfaceNum].center();
	double size = theSurfaces[activeSurfaceNum].size();

	for (int i = 1; i<activeSurfaceNum; i++) {
		workV6.addVector(0.0, center, conHeig);
		workV6 -= devia;
		workV6 *= theSurfaces[i].size() / size;
		workV6 += devia;

		//workV6 = devia - (devia - center*conHeig) * theSurfaces[i].size() / size;

		workV6 /= conHeig;
		theSurfaces[i].setCenter(workV6);
	}
}


int LadeDuncanMultiYield::isCrossingNextSurface(void)
{
	int numOfSurfaces = numOfSurfacesx[matN];
	double residualPress = residualPressx[matN];
	double yieldsize = theSurfaces[activeSurfaceNum + 1].size();
	if (activeSurfaceNum == numOfSurfaces) return 0;

	workV6 = trialStress.deviator();
	double conHeig = trialStress.volume() - residualPress;
	workV6.addVector(1.0, theSurfaces[activeSurfaceNum + 1].center(), -conHeig);
	double rat = stressRatio(workV6, conHeig);

	if (rat > theSurfaces[activeSurfaceNum + 1].size())
		return 1;
	return 0;
}


 
double LadeDuncanMultiYield::stressRatio(const Vector & tdev, double tdil)
{
	double a1 = a1x[matN];
	double residualPress = residualPressx[matN];
	double as, cs, ds, ratio;
	static Vector Qiu;
	Qiu = tdev;
	cs = tdev && tdev;
	ds = TRIPRO(Qiu, Qiu, Qiu, 1);
	as = a1 * pow(3.*tdil, 3);
	cs = -tdil * cs / 2.;
	ds = ds / 3.;
	ratio = thirdOrderEqn(as, 0., cs, ds, 2);
	double yieldfunc = as*pow(ratio, 3) + cs*ratio + ds;
	if (abs(yieldfunc) - 1e-3 > 0) {
		ratio = thirdOrderEqn(as, 0., cs, ds, 2);
	}

	return ratio;
}



double
LadeDuncanMultiYield::thirdOrderEqn(double A, double B, double C, double D, double I) { // result = Aij*Bjk*Cki
	double DLIMIT = 1.e-6;
	double A1 = B / A;
	double B1 = C / A;
	double C1 = D / A;
	double TEMP = A1 / 3.0;
	double P, Q, X1, X2, X3, X;
	dcmplx AA, BB, QQ, Y1, Y2, Y3, zhijian, zhijian1;
	if (A1 == 0) {
		P = B1;
		Q = C1;
	}
	else {
		P = -A1*A1 / 3.0 + B1;
		Q = 2 * pow((A1 / 3.0), 3) - A1*B1 / 3 + C1;
	}
	if (P == 0 && Q == 0) {
		X1 = -TEMP;
		X2 = -TEMP;
		X3 = -TEMP;
	}

	else {
		QQ = pow((P / 3), 3) + pow((Q / 2.0), 2);
		AA = pow(-Q / 2.0 + sqrt(QQ), 1.0 / 3.0);
		BB = pow(-Q / 2.0 - sqrt(QQ), 1.0 / 3.0);
		Y1 = AA + BB;
		zhijian = dcmplx(0, 1)*dcmplx(0, 1);
		zhijian1 = dcmplx(0, 1)*AA;
		QQ = dcmplx(0, 1) * (AA - BB)*sqrt(3.0) / 2.0;
		Y2 = -Y1 / 2.0 + QQ;
		Y3 = -Y1 / 2.0 - QQ;
		X1 = real(Y1) - TEMP;
		X2 = real(Y2) - TEMP;
		X3 = real(Y3) - TEMP;
	}

	if (I == 0) { // FIND THE DIRECTION OF ACTIVE SURFACE MOTION

		static Vector t0(6);
		static Vector t1(6);
		static Vector t2(6);
		static Vector center(6);
		static Vector outcenter(6);
		double residualPress = residualPressx[matN];
		int numOfSurfaces = numOfSurfacesx[matN];
		double conHeig = trialStress.volume() - residualPress;
		center = theSurfaces[activeSurfaceNum].center();
		double size = theSurfaces[activeSurfaceNum].size();
		outcenter = theSurfaces[activeSurfaceNum + 1].center();
		double outsize = theSurfaces[activeSurfaceNum + 1].size();
		t0 = trialStress.deviator();
		t0.addVector(1.0, center, -conHeig);

		if (X1 > 0) {
			X = X1;
		}

		else if (X2 > 0) {
			X = X2;
		}
		else if (X3 > 0) {
			X = X3;
		}
	}

	else if (I == 1) { // ! FIND THE AMOUNT OF ACTIVE SURFACE MOTION
		if (X2 < 0)
			X = X3;
		else if (X3 < 0)
			X = X2;
		else if (X2 < X3)
			X = X2;
		else
			X = X3;
	}
	else if (I == 2) {// ! STRESS RATIO
		X = X1;
	}

	return X;

};

double
LadeDuncanMultiYield::TRIPRO(Vector &A, Vector &B, Vector &C, int I) { // result = Aij*Bjk*Cki

	if ((A.Size() != 6) || (B.Size() != 6) || (C.Size() != 6)) {
		opserr << "Fatal: CapPlasticity::TRIPRO() size does not match! " << endln;
		exit(-1);
	}
	double result = 0.0;

	if (I == 1) {
		result = 3.*(A(0)*A(1)*A(2) - A(0)*A(4)*A(4) - A(1)*A(5)*A(5) - A(2)*A(3)*A(3) + 2.*A(3)*A(4)*A(5));
	}
	else

		result = A(0)*(B(0)*C(0) + B(3)*C(3) + B(5)*C(5)) + A(1)*(B(1)*C(1) + B(3)*C(3) + B(4)*C(4)) + A(2)*(B(2)*C(2) + B(4)*C(4) + B(5)*C(5)) + A(3)*(B(3)*C(0) + B(1)*C(3) + B(4)*C(5)) + A(3)*(B(0)*C(3) + B(3)*C(1) + B(5)*C(4)) + A(4)*(B(5)*C(3) + B(4)*C(1) + B(2)*C(4)) + A(4)*(B(3)*C(5) + B(1)*C(4) + B(4)*C(2)) + A(5)*(B(5)*C(0) + B(4)*C(3) + B(2)*C(5)) + A(5)*(B(0)*C(5) + B(3)*C(4) + B(5)*C(2));

	return result;

};

