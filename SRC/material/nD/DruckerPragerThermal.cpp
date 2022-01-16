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

// $Revision: 1.2 $
// $Date: 2010-02-04 20:50:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DruckerPragerThermal.cpp,v $

// Based on DruckerPrager written by Kathryn Petek, Peter Mackenzie-Helnwein, and Pedro Arduino
// Modified by Jian Jiang, Liming Jiang [http://openseesforfire.github.io]
//
// Description: This file contains the implementation for the DruckerPragerThermal class.
// In its original version this file was called DruckerPragerTensionCutoff.
//				
//		This implementation of Drucker-Prager allows for non-associative flow 
//		through the use of the rho_bar parameter.  It also includes linear 
//		kinematic hardening and linear and nonlinear isotropic hardening.
//

#include <DruckerPragerThermal.h>
#include <DruckerPrager3DThermal.h>
#include <DruckerPragerPlaneStrain.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <string.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

const double DruckerPragerThermal::one3 = 1.0 / 3.0;
const double DruckerPragerThermal::two3 = 2.0 / 3.0;
const double DruckerPragerThermal::root23 = sqrt(2.0 / 3.0);

#include <elementAPI.h>
#define OPS_Export 

static int numDruckerPragerMaterials = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_DruckerPragerMaterialThermal)
{
	if (numDruckerPragerMaterials == 0) {
		numDruckerPragerMaterials++;
		opserr << "DruckerPragerThermal nDmaterial - Written: K.Petek, P.Mackenzie-Helnwein, P.Arduino, U.Washington\n";
	}

	// Pointer to a uniaxial material that will be returned
	NDMaterial *theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 12) {
		opserr << "Want: nDMaterial DruckerPragerThermal tag? K? G? sigma_y? rho? rho_bar? Kinf? Ko? delta1? delta2? H? theta? <massDensity? atm?>" << endln;
		return 0;
	}

	int tag;
	double dData[14];

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid nDMaterial DruckerPragerThermal material  tag" << endln;
		return 0;
	}
	if (numArgs == 12) {
		numData = 11;
	}
	else if (numArgs == 13) {
		numData = 12;
	}
	else {
		numData = 13;
	}

	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial DruckerPragerThermal material  with tag: " << tag << endln;
		return 0;
	}

	if (numArgs == 12) {
		theMaterial = new DruckerPragerThermal(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
			dData[6], dData[7], dData[8], dData[9], dData[10]);
	}
	else if (numArgs == 13) {
		theMaterial = new DruckerPragerThermal(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
			dData[6], dData[7], dData[8], dData[9], dData[10], dData[11]);
	}
	else {
		theMaterial = new DruckerPragerThermal(tag, 0, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
			dData[6], dData[7], dData[8], dData[9], dData[10], dData[11], dData[12]);
	}

	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory for nDMaterial DruckerPragerThermal material  with tag: " << tag << endln;
	}

	return theMaterial;
}



//full constructor
DruckerPragerThermal::DruckerPragerThermal(int tag, int classTag, double bulk, double shear, double s_y, double r,
	double r_bar, double Kinfinity, double Kinit, double d1,
	double d2, double H, double t, double mDen, double atm)
	: NDMaterial(tag, classTag),
	mEpsilon(6),
	mEpsilon_n_p(6),
	mEpsilon_n1_p(6),
	mSigma(6),
	mBeta_n(6),
	mBeta_n1(6),
	mCe(6, 6),
	mCep(6, 6),
	mI1(6),
	mIIvol(6, 6),
	mIIdev(6, 6),
	mState(5)
{
	massDen = mDen;
	mKref = bulk;
	mGref = shear;
	mPatm = atm;
	mK = bulk;
	mG = shear;
	msigma_y = s_y;
	mrho = r;
	mrho_bar = r_bar;
	mKinf = Kinfinity;
	mKo = Kinit;
	mdelta1 = d1;
	mdelta2 = d2;
	mHard = H;
	mtheta = t;

	mKref0 = bulk;
	mGref0 = shear;
	msigma_y0 = s_y;
	ThermalElongation = 0.0;

	if (mrho == 0.0) {
		mTo = 1e10;
	}
	else {
		mTo = root23*msigma_y / mrho;
	}
	// set the elastic flag
	//  0 = elastic+no param update; 1 = elastic+param update; 2 = elastoplastic+no param update (default)
	mElastFlag = 2;

	// Use these values to deactivate yield surface 1 - Create Pure Tension Cutoff
	//msigma_y = 1e10;
	//mTo      = 100;

	this->initialize();
}

//null constructor
DruckerPragerThermal::DruckerPragerThermal()
	: NDMaterial(),
	mEpsilon(6),
	mEpsilon_n_p(6),
	mEpsilon_n1_p(6),
	mSigma(6),
	mBeta_n(6),
	mBeta_n1(6),
	mCe(6, 6),
	mCep(6, 6),
	mI1(6),
	mIIvol(6, 6),
	mIIdev(6, 6),
	mState(5)
{
	massDen = 0.0;
	mKref = 0.0;
	mGref = 0.0;
	mPatm = 101.0;
	mK = 0.0;
	mG = 0.0;
	msigma_y = 1e+10;
	mrho = 0.0;
	mrho_bar = 0.0;
	mKinf = 0.0;
	mKo = 0.0;
	mdelta1 = 0.0;
	mdelta2 = 0.0;
	mHard = 0.0;
	mtheta = 0.0;
	mTo = 0.0;

	mKref0 = 0.0;
	mGref0 = 0.0;
	msigma_y0 = 0.0;
	ThermalElongation = 0.0;

	mElastFlag = 2;

	this->initialize();
}

//destructor
DruckerPragerThermal::~DruckerPragerThermal()
{
}

//zero internal variables
void DruckerPragerThermal::initialize()
{
	mEpsilon.Zero();
	mEpsilon_n_p.Zero();
	mEpsilon_n1_p.Zero();

	mSigma.Zero();
	mBeta_n.Zero();
	mBeta_n1.Zero();

	mAlpha1_n = 0.0;
	mAlpha1_n1 = 0.0;
	mAlpha2_n = 0.0;
	mAlpha2_n1 = 0.0;
	mFlag = 1;

	mHprime = (1 - mtheta)*mHard;

	// 2nd order Identity Tensor
	mI1.Zero();
	mI1(0) = 1;
	mI1(1) = 1;
	mI1(2) = 1;

	// 4th order Volumetric Tensor
	// IIvol = I1 tensor I1
	mIIvol.Zero();
	mIIvol(0, 0) = 1;
	mIIvol(0, 1) = 1;
	mIIvol(0, 2) = 1;
	mIIvol(1, 0) = 1;
	mIIvol(1, 1) = 1;
	mIIvol(1, 2) = 1;
	mIIvol(2, 0) = 1;
	mIIvol(2, 1) = 1;
	mIIvol(2, 2) = 1;

	// 4th order Deviatoric Tensor
	// Note:  this is the contravariant form!
	//        useable for s^a = 2G * IIdev^ab * epsilon_b
	// (Need a different form for s^a = IIdev ^a_b * sigma^a)
	mIIdev.Zero();
	mIIdev(0, 0) = two3;
	mIIdev(0, 1) = -one3;
	mIIdev(0, 2) = -one3;
	mIIdev(1, 0) = -one3;
	mIIdev(1, 1) = two3;
	mIIdev(1, 2) = -one3;
	mIIdev(2, 0) = -one3;
	mIIdev(2, 1) = -one3;
	mIIdev(2, 2) = two3;
	mIIdev(3, 3) = 0.5;
	mIIdev(4, 4) = 0.5;
	mIIdev(5, 5) = 0.5;

	mCe = mK * mIIvol + 2 * mG*mIIdev;
	mCep = mCe;
	mState.Zero();
}


NDMaterial * DruckerPragerThermal::getCopy(const char *type)
{
	if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
		DruckerPragerPlaneStrain *clone;
		clone = new DruckerPragerPlaneStrain(this->getTag(), mK, mG, msigma_y, mrho, mrho_bar, mKinf, mKo,
			mdelta1, mdelta2, mHard, mtheta, massDen, mPatm);
		return clone;
	}
	else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
		DruckerPrager3DThermal *clone;
		clone = new DruckerPrager3DThermal(this->getTag(), mK, mG, msigma_y, mrho, mrho_bar, mKinf, mKo,
			mdelta1, mdelta2, mHard, mtheta, massDen, mPatm);
		return clone;
	}
	else {
		opserr << "DruckerPragerThermal::getCopy failed to get copy: " << type << endln;
		return 0;
	}
}

int DruckerPragerThermal::commitState(void)
{
	mEpsilon_n_p = mEpsilon_n1_p;
	mAlpha1_n = mAlpha1_n1;
	mAlpha2_n = mAlpha2_n1;
	mBeta_n = mBeta_n1;

	return 0;
}

int DruckerPragerThermal::revertToLastCommit(void)
{
	return 0;
}

int DruckerPragerThermal::revertToStart(void)
{
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	}
	else {
		// normal call for revertToStart (not initialStateAnalysis)
		this->initialize();
	}

	return 0;
}

NDMaterial*
DruckerPragerThermal::getCopy(void)
{
	opserr << "DruckerPragerThermal::getCopy -- subclass responsibility\n";
	exit(-1);
	return 0;
}

const char*
DruckerPragerThermal::getType(void) const
{
	opserr << "DruckerPragerThermal::getType -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int
DruckerPragerThermal::getOrder(void) const
{
	opserr << "DruckerPragerThermal::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}


double
DruckerPragerThermal::setThermalTangentAndElongation(double &tempT, double&ET, double&Elong)
{
	double TempT = tempT + 20;
	double feT;//elastic compressive strength
	double ftT;//tensile strength
	double ee0; // elastic strain corresponding to fe0
	double eeT;
	double ec0; //strain at fc0
	double fc0;//compressive strength at ambient
	double ecT;//strain at fcT;
	double fcT;// compressive strength at elevated temperature
	double mRatio;//rho=feT/ftT
	double H_H0; //=H/H0
	fc0 = root23*msigma_y / (sqrt(5.0 / 9) + mrho / 3);
	//fc0 = 30e6;
	ec0 = 0.0025;
	H_H0 = 1.71;
	ET = fc0 / ec0;
	//ee0 = 0.0002;
	// EN 1992 pt 1-2-1. concrete at elevated temperatures
	//Case  JJ ET=2*fcT/ec0
	/*
	if (TempT <= 20) {

	ET = fc0/ec0;

	}
	else if (TempT <= 100) {

	mKref = mKref0/(0.0025 + (0.004-0.0025)*(TempT - 20)/80)*ec0;
	mGref = mGref0/(0.0025 + (0.004-0.0025)*(TempT - 20)/80)*ec0;;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0;
	fcT = fc0;
	ecT = (0.0025 + (0.004-0.0025)*(TempT - 20)/80);
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;//mHard is positive
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT;


	}
	else if (TempT <= 200) {


	mKref = mKref0*(1 - (TempT - 100)*0.05/100)/(0.0040 + (0.0055-0.0040)*(TempT - 100)/100)*0.0025;
	mGref = mGref0*(1 - (TempT - 100)*0.05/100)/(0.0040 + (0.0055-0.0040)*(TempT - 100)/100)*0.0025;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(1 - (TempT - 100)*0.05/100);
	fcT = fc0*(1 - (TempT - 100)*0.05/100);
	ecT = (0.0040 + (0.0055-0.0040)*(TempT - 100)/100);
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT*(1-(TempT-100)/500);
	}
	else if (TempT <= 300) {


	mKref = mKref0*(0.95 - (TempT - 200)*0.1/100)/(0.0055 + (0.0070-0.0055)*(TempT - 200)/100)*0.0025;
	mGref = mGref0*(0.95 - (TempT - 200)*0.1/100)/(0.0055 + (0.0070-0.0055)*(TempT - 200)/100)*0.0025;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.95 - (TempT - 200)*0.1/100);
	fcT = fc0*(0.95 - (TempT - 200)*0.1/100);
	ecT = (0.0055 + (0.0070-0.0055)*(TempT - 200)/100);
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT*(1-(TempT-100)/500);

	}
	else if (TempT <= 400) {

	mKref = mKref0*(0.85 - (TempT - 300)*0.1/100)/(0.0070 + (0.0100-0.0070)*(TempT - 300)/100)*0.0025;
	mGref = mGref0*(0.85 - (TempT - 300)*0.1/100)/(0.0070 + (0.0100-0.0070)*(TempT - 300)/100)*0.0025;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.85 - (TempT - 300)*0.1/100);
	fcT = fc0*(0.85 - (TempT - 300)*0.1/100);
	ecT = (0.0070 + (0.0100-0.0070)*(TempT - 300)/100);
	ET = fcT/ecT;
	ee0 = 0.33*fcT/ET;
	//mHard = 0.66*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT*(1-(TempT-100)/500);
	}
	else if (TempT <= 500) {

	mKref = mKref0*(0.75 - (TempT - 400)*0.15/100)/(0.0100 + (0.0150-0.0100)*(TempT - 400)/100)*0.0025;
	mGref = mGref0*(0.75 - (TempT - 400)*0.15/100)/(0.0100 + (0.0150-0.0100)*(TempT - 400)/100)*0.0025;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.75 - (TempT - 400)*0.15/100);
	fcT = fc0*(0.75 - (TempT - 400)*0.15/100);
	ecT = (0.0100 + (0.0150-0.0100)*(TempT - 400)/100);
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT*(1-(TempT-100)/500);
	}
	else if (TempT <= 600) {

	mKref = mKref0*(0.60 - (TempT - 500)*0.15/100)/(0.0150 + (0.0250-0.0150)*(TempT - 500)/100)*0.0025;
	mGref = mGref0*(0.60 - (TempT - 500)*0.15/100)/(0.0150 + (0.0250-0.0150)*(TempT - 500)/100)*0.0025;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.60 - (TempT - 500)*0.15/100);
	fcT = fc0*(0.60 - (TempT - 500)*0.15/100);
	ecT = (0.0150 + (0.0250-0.0150)*(TempT - 500)/100);
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= SigT*(1-(TempT-100)/500);
	if(mTo<1e-3)
	mTo = 1e-2*SigT;
	}
	else if (TempT <= 700) {

	mKref = mKref0*(0.45 - (TempT - 600)*0.15/100)/10;
	mGref = mGref0*(0.45 - (TempT - 600)*0.15/100)/10;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.45 - (TempT - 600)*0.15/100);
	fcT = fc0*(0.45 - (TempT - 600)*0.15/100);
	ecT = 0.025;
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo= 1e-2*SigT;
	}
	else if (TempT <= 800) {

	mKref = mKref0*(0.30 - (TempT - 700)*0.15/100)/10;
	mGref = mGref0*(0.30 - (TempT - 700)*0.15/100)/10;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.30 - (TempT - 700)*0.15/100);
	fcT = fc0*(0.30 - (TempT - 700)*0.15/100);
	ecT = 0.025;
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo = 1e-2*SigT;
	}
	else if (TempT <= 900) {

	mKref = mKref0*(0.15 - (TempT - 800)*0.07/100)/10;
	mGref = mGref0*(0.15 - (TempT - 800)*0.07/100)/10;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.15 - (TempT - 800)*0.07/100);
	fcT = fc0*(0.15 - (TempT - 800)*0.07/100);
	ecT = 0.025;
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	//mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo = 1e-2*SigT;
	}
	else if (TempT <= 1000) {

	mKref = mKref0*(0.08 - (TempT - 900)*0.04/100)/10;
	mGref = mGref0*(0.08 - (TempT - 900)*0.04/100)/10;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.08 - (TempT - 900)*0.04/100);
	fcT = fc0*(0.08 - (TempT - 900)*0.04/100);
	ecT = 0.025;
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	// mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo = 1e-2*SigT;
	}
	else if (TempT <= 1100) {

	mKref = mKref0*(0.04 - (TempT - 1000)*0.03/100)/10;
	mGref = mGref0*(0.04 - (TempT - 1000)*0.03/100)/10;
	mK = mKref;
	mG = mGref;
	msigma_y = msigma_y0*(0.04 - (TempT - 1000)*0.03/100);
	fcT = fc0*(0.04 - (TempT - 1000)*0.03/100);
	ecT = 0.025;
	ET = fcT/ecT;
	ee0 = 1/3*fcT/ET;
	// mHard = 2/3*fcT/(ecT-ee0)/2.8;
	mHard = 0.01*ET;//mHard is negative
	mTo = 1e-2*SigT;
	}


	else  {
	opserr << "the temperature " <<TempT<<"is invalid\n";
	}

	*/
	//ET = 1.5*fc0/ec0;  //Modificaiton
	// Calculate thermal elongation 
	if (TempT <= 0) {
		ThermalElongation = 0.0;
	}
	else if (TempT <= 700) {
		ThermalElongation = -1.8e-4 + 9e-6*TempT + 2.3e-11*TempT*TempT*TempT;

	}
	else if (TempT <= 1200) {
		ThermalElongation = 14e-3;
		//ThermalElongation = 12e-6*(tempT);
	}
	else {
		opserr << "the temperature is invalid\n";
	}

	//ET = E;  
	//ET = 3.84e10;
	ThermalElongation = 12e-6*(tempT);
	Elong = ThermalElongation;
	//TempAndElong(0) = TempT-20;
	//  TempAndElong(1) = ThermalElongation ;
	//this->initialize();
	//this->updateElasticParam( );
	this->plastic_integrator();
	return 0;
}

//--------------------Plasticity-------------------------------------
//plasticity integration routine
void DruckerPragerThermal::plastic_integrator()
{
	bool okay;		// boolean variable to ensure satisfaction of multisurface kuhn tucker conditions
	double f1;
	double f2;
	double norm_eta;
	double Invariant_1;
	double Invariant_ep;
	double norm_ep;
	double norm_dev_ep;
	Vector epsilon_e(6);
	Vector s(6);
	Vector eta(6);
	Vector dev_ep(6);
	Vector Jact(2);

	double fTOL;
	double gTOL;
	fTOL = 0.0;
	gTOL = -1.0e-10;

	double NormCep;

	double alpha1;			// hardening parameter for DP surface
	double alpha2;			// hardening parameter for tension cut-off
	Vector n(6);			// normal to the yield surface in strain space
	Vector R(2);			// residual vector
	Vector gamma(2);		// vector of consistency parameters
	Vector dgamma(2);		// incremental vector of consistency parameters
	Matrix g(2, 2);			// jacobian of the corner region (return map)
	Matrix g_contra(2, 2);	// inverse of jacobian of the corner region

							// set trial state:

							// epsilon_n1_p_trial = ..._n1_p  = ..._n_p
	mEpsilon_n1_p = mEpsilon_n_p;

	// alpha1_n+1_trial
	mAlpha1_n1 = mAlpha1_n;
	// alpha2_n+1_trial
	mAlpha2_n1 = mAlpha2_n;

	// beta_n+1_trial
	mBeta_n1 = mBeta_n;

	// epsilon_elastic = epsilon_n+1 - epsilon_n_p
	epsilon_e = mEpsilon - mEpsilon_n1_p;

	// trial stress
	mSigma = mCe*epsilon_e;

	// deviator stress tensor: s = 2G * IIdev * epsilon_e
	//I1_trial
	Invariant_1 = (mSigma(0) + mSigma(1) + mSigma(2));

	// s_n+1_trial
	s = mSigma - (Invariant_1 / 3.0)*mI1;

	//eta_trial = s_n+1_trial - beta_n;
	eta = s - mBeta_n;

	// compute yield function value (contravariant norm)
	norm_eta = sqrt(eta(0)*eta(0) + eta(1)*eta(1) + eta(2)*eta(2) + 2 * (eta(3)*eta(3) + eta(4)*eta(4) + eta(5)*eta(5)));

	// f1_n+1_trial
	f1 = norm_eta + mrho*Invariant_1 - root23*Kiso(mAlpha1_n1);

	// f2_n+1_trial
	f2 = Invariant_1 - T(mAlpha2_n1);

	// update elastic bulk and shear moduli 
	this->updateElasticParam();

	// check trial state
	int count = 1;
	if ((f1 <= fTOL) && (f2 <= fTOL) || mElastFlag < 2) {

		okay = true;
		// trial state = elastic state - don't need to do any updates.
		mCep = mCe;
		count = 0;

		// set state variables for recorders
		Invariant_ep = mEpsilon_n1_p(0) + mEpsilon_n1_p(1) + mEpsilon_n1_p(2);

		norm_ep = sqrt(mEpsilon_n1_p(0)*mEpsilon_n1_p(0) + mEpsilon_n1_p(1)*mEpsilon_n1_p(1) + mEpsilon_n1_p(2)*mEpsilon_n1_p(2)
			+ 0.5*(mEpsilon_n1_p(3)*mEpsilon_n1_p(3) + mEpsilon_n1_p(4)*mEpsilon_n1_p(4) + mEpsilon_n1_p(5)*mEpsilon_n1_p(5)));

		dev_ep = mEpsilon_n1_p - one3*Invariant_ep*mI1;

		norm_dev_ep = sqrt(dev_ep(0)*dev_ep(0) + dev_ep(1)*dev_ep(1) + dev_ep(2)*dev_ep(2)
			+ 0.5*(dev_ep(3)*dev_ep(3) + dev_ep(4)*dev_ep(4) + dev_ep(5)*dev_ep(5)));

		mState(0) = Invariant_1;
		mState(1) = norm_eta;
		mState(2) = Invariant_ep;
		mState(3) = norm_dev_ep;
		mState(4) = norm_ep;
		return;
	}
	else {
		// plastic correction required
		okay = false;

		// determine number of active surfaces.  size & fill Jact
		if ((f1 > fTOL) && (f2 <= fTOL)) {
			// f1 surface only
			Jact(0) = 1;
			Jact(1) = 0;
		}
		else if ((f1 <= fTOL) && (f2 > fTOL)) {
			// f2 surface only
			Jact(0) = 0;
			Jact(1) = 1;
		}
		else if ((f1 > fTOL) && (f2 > fTOL)) {
			// both surfaces active
			Jact(0) = 1;
			Jact(1) = 1;
		}
	}

	//-----------------MultiSurface Placity Return Map--------------------------------------
	while (!okay) {

		alpha1 = mAlpha1_n;
		alpha2 = mAlpha2_n;

		//  n = eta / norm_eta;  (contravaraint)
		if (norm_eta < 1.0e-13) {
			n.Zero();
		}
		else {
			n = eta / norm_eta;
		}

		// initialize R, gamma1, gamma2, dgamma1, dgamma2 = 0
		R.Zero();
		gamma.Zero();
		dgamma.Zero();
		// initialize g such that det(g) = 1
		g(0, 0) = 1;
		g(1, 1) = 1;
		g(1, 0) = 0;
		g(0, 1) = 0;

		// Newton procedure to compute nonlinear gamma1 and gamma2
		//initialize terms
		for (int i = 0; i < 2; i++) {
			if (Jact(i) == 1) {
				R(0) = norm_eta - (2 * mG + two3*mHprime)*gamma(0) + mrho*Invariant_1
					- 9 * mK*mrho*mrho_bar*gamma(0) - 9 * mK*mrho*gamma(1) - root23*Kiso(alpha1);
				g(0, 0) = -2 * mG - two3*(mHprime + Kisoprime(alpha1)) - 9 * mK*mrho*mrho_bar;
			}
			else if (Jact(i) == 2) {
				R(1) = Invariant_1 - 9 * mK*mrho_bar*gamma(0) - 9 * mK*gamma(1) - T(alpha2);
				g(1, 1) = -9 * mK + mdelta2*T(alpha2);
			}
		}
		if (Jact(0) == 1 && Jact(1) == 1) {
			g(0, 1) = -9 * mK*mrho;
			g(1, 0) = mrho_bar*(-9 * mK + mdelta2*T(alpha2));
		}
		g.Invert(g_contra);

		// iteration counter
		int m = 0;

		//iterate
		while ((fabs(R.Norm()) > 1e-10) && (m < 10)) {

			dgamma = -1 * g_contra * R;
			gamma += dgamma;

			//update alpha1 and alpha2
			alpha1 = mAlpha1_n + root23*gamma(0);
			alpha2 = mAlpha2_n + mrho_bar*gamma(0) + gamma(1);

			// reset g & R matrices
			g(0, 0) = 1;
			g(1, 1) = 1;
			g(1, 0) = 0;
			g(0, 1) = 0;
			R.Zero();
			for (int i = 0; i < 2; i++) {
				if (Jact(i) == 1) {
					R(0) = norm_eta - (2 * mG + two3*mHprime)*gamma(0) + mrho*Invariant_1
						- 9 * mK*mrho*mrho_bar*gamma(0) - 9 * mK*mrho*gamma(1) - root23*Kiso(alpha1);
					g(0, 0) = -2 * mG - two3*(mHprime + Kisoprime(alpha1)) - 9 * mK*mrho*mrho_bar;
				}
				else if (Jact(i) == 2) {
					R(1) = Invariant_1 - 9 * mK*mrho_bar*gamma(0) - 9 * mK*gamma(1) - T(alpha2);
					g(1, 1) = -9 * mK + mdelta2*T(alpha2);
				}
			}
			if (Jact(0) == 1 && Jact(1) == 1) {
				g(0, 1) = -9 * mK*mrho;
				g(1, 0) = mrho_bar*(-9 * mK + mdelta2*T(alpha2));
			}
			g.Invert(g_contra);

			m++;
		}

		// check maintain Kuhn-Tucker conditions
		f1 = norm_eta - (2 * mG + two3*mHprime)*gamma(0) + mrho*Invariant_1
			- 9 * mK*mrho*mrho_bar*gamma(0) - 9 * mK*mrho*gamma(1) - root23*Kiso(alpha1);
		f2 = Invariant_1 - 9 * mK*mrho_bar*gamma(0) - 9 * mK*gamma(1) - T(alpha2);

		if (count > 100) {
			okay = true;
			break;
		}

		// check active surfaces
		if ((Jact(0) == 1) && (Jact(1) == 0)) {
			// f2 may be > or < f2_tr because of softening of f2 related to alpha1
			if (f2 >= fTOL) {
				// okay = false;
				Jact(0) = 1;
				Jact(1) = 1;
				count += 1;

			}
			else {
				okay = true;

			}
		}
		else if ((Jact(0) == 0) && (Jact(1) == 1)) {
			// f1 will always be less than f1_tr
			okay = true;
		}
		else if ((Jact(0) == 1) && (Jact(1) == 1)) {
			if ((gamma(0) <= gTOL) && (gamma(1) > gTOL)) {
				// okay = false;
				Jact(0) = 0;
				Jact(1) = 1;
				count += 1;
			}
			else if ((gamma(0) > gTOL) && (gamma(1) <= gTOL)) {
				// okay = false;
				Jact(0) = 1;
				Jact(1) = 0;
				count += 1;
			}
			else if ((gamma(0) > gTOL) && (gamma(1) > gTOL)) {
				okay = true;
			}
		}

		if ((count > 3) && (!okay)) {
			Jact(0) = 1;
			Jact(1) = 1;
			count += 100;
		}

		if (count > 3) {
			opserr << "Jact = " << Jact;
			opserr << "count = " << count << endln;
		}

	} // end of while(!okay) loop


	  //update everything and exit!

	Vector b1(6);
	Vector b2(6);
	Vector n_covar(6);
	Vector temp1(6);
	Vector temp2(6);

	// update alpha1 and alpha2
	mAlpha1_n1 = alpha1;
	mAlpha2_n1 = alpha2;

	//update epsilon_n1_p
	//first calculate n_covar
	// n_a = G_ab * n^b = covariant
	n_covar(0) = n(0);
	n_covar(1) = n(1);
	n_covar(2) = n(2);
	n_covar(3) = 2 * n(3);
	n_covar(4) = 2 * n(4);
	n_covar(5) = 2 * n(5);
	mEpsilon_n1_p = mEpsilon_n_p + (mrho_bar*gamma(0) + gamma(1))*mI1 + gamma(0)*n_covar;


	Invariant_ep = mEpsilon_n1_p(0) + mEpsilon_n1_p(1) + mEpsilon_n1_p(2);

	norm_ep = sqrt(mEpsilon_n1_p(0)*mEpsilon_n1_p(0) + mEpsilon_n1_p(1)*mEpsilon_n1_p(1) + mEpsilon_n1_p(2)*mEpsilon_n1_p(2)
		+ 0.5*(mEpsilon_n1_p(3)*mEpsilon_n1_p(3) + mEpsilon_n1_p(4)*mEpsilon_n1_p(4) + mEpsilon_n1_p(5)*mEpsilon_n1_p(5)));

	dev_ep = mEpsilon_n1_p - one3*Invariant_ep*mI1;

	norm_dev_ep = sqrt(dev_ep(0)*dev_ep(0) + dev_ep(1)*dev_ep(1) + dev_ep(2)*dev_ep(2)
		+ 0.5*(dev_ep(3)*dev_ep(3) + dev_ep(4)*dev_ep(4) + dev_ep(5)*dev_ep(5)));

	// update sigma
	mSigma -= (3 * mK*mrho_bar*gamma(0) + 3 * mK*gamma(1))*mI1 + 2 * mG*gamma(0)*n;

	s -= 2 * mG*gamma(0) * n;
	Invariant_1 -= 9 * mK*mrho_bar*gamma(0) + 9 * mK*gamma(1);
	//mSigma        = s + Invariant_1/3.0 * mI1;

	//update beta_n1
	mBeta_n1 = mBeta_n - (two3*mHprime*gamma(0))*n;

	//eta_n+1 = s_n+1 - beta_n+1;
	eta = s - mBeta_n1;
	norm_eta = sqrt(eta(0)*eta(0) + eta(1)*eta(1) + eta(2)*eta(2) + 2 * (eta(3)*eta(3) + eta(4)*eta(4) + eta(5)*eta(5)));

	// update Cep
	// note: Cep is contravariant
	if ((Jact(0) == 1) && (Jact(1) == 0)) {
		b1 = 2 * mG*n + 3 * mK*mrho*mI1;
		b2.Zero();
	}
	else if ((Jact(0) == 0) && (Jact(1) == 1)) {
		b1.Zero();
		b2 = 3 * mK*mI1;
	}
	else if ((Jact(0) == 1) && (Jact(1) == 1)) {
		b1 = 2 * mG*n + 3 * mK*mrho*mI1;
		b2 = 3 * mK*mI1;
	}

	temp1 = g_contra(0, 0)*b1 + g_contra(0, 1)*b2;
	temp2 = mrho_bar*temp1 + g_contra(1, 0)*b1 + g_contra(1, 1)*b2;

	NormCep = 0.0;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			mCep(i, j) = mCe(i, j)
				+ 3 * mK * mI1(i)*temp2(j)
				+ 2 * mG * n(i)*temp1(j)
				- 4 * mG*mG / norm_eta*gamma(0) * (mIIdev(i, j) - n(i)*n(j));
			NormCep += mCep(i, j)*mCep(i, j);
		}
	}

	if (NormCep < 1e-10) {
		mCep = 1.0e-3 * mCe;
		opserr << "NormCep = " << NormCep << endln;
	}

	mState(0) = Invariant_1;
	mState(1) = norm_eta;
	mState(2) = Invariant_ep;
	mState(3) = norm_dev_ep;
	mState(4) = norm_ep;

	return;
}

int DruckerPragerThermal::updateElasticParam()
{
	double Sigma_mean = 0.0;
	if (mElastFlag == 1 && mFlag == 1) {
		Sigma_mean = -one3*(mSigma(0) + mSigma(1) + mSigma(2));
		if (Sigma_mean < 0.0) Sigma_mean = 0.0;  // prevents modulus update for cases where tension exists 
		mK = mKref * pow(1 + (Sigma_mean / mPatm), 0.5);
		mG = mGref * pow(1 + (Sigma_mean / mPatm), 0.5);
		mCe = mK * mIIvol + 2 * mG*mIIdev;
		mFlag = 0;
		//opserr << "Plastic Integrator -->" << "K = " << mK  << "  G =" << mG << endln;
	}
	else if (mElastFlag != 1) {
		mFlag = 1;
	}

	return 0;
}

double DruckerPragerThermal::Kiso(double alpha1)
{
	return msigma_y + mtheta * mHard * alpha1 + (mKinf - mKo) * (1 - exp(-mdelta1 * alpha1));
}


double DruckerPragerThermal::Kisoprime(double alpha1)
{
	return mtheta * mHard + (mKinf - mKo) * mdelta1*  exp(-mdelta1 * alpha1);
}

double DruckerPragerThermal::T(double alpha2)
{
	return mTo * exp(-mdelta2 * alpha2);
}


double DruckerPragerThermal::deltaH(double dGamma)
{
	return mHprime * root23 * dGamma;
}


Vector DruckerPragerThermal::getState()
{
	return mState;
}


Response*
DruckerPragerThermal::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	Response *theResponse = 0;
	const char *matType = this->getType();

	output.tag("NdMaterialOutput");
	output.attr("matType", this->getClassType());
	output.attr("matTag", this->getTag());

	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int DruckerPragerThermal::getResponse(int responseID, Information &matInfo)
{
	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;
	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;
	case 3:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getState();
		return 0;
	default:
		return -1;
	}
}

int DruckerPragerThermal::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 2)
		return -1;

	int theMaterialTag = atoi(argv[1]);

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
	}

	return -1;
}

int
DruckerPragerThermal::updateParameter(int responseID, Information &info)
{
	// updateMaterialStage called
	if (responseID == 1) {
		mElastFlag = (int)info.theDouble;
	}
	// materialState called
	else if (responseID == 5) {
		mElastFlag = (int)info.theDouble;
	}
	// frictionalStrength called
	else if (responseID == 7) {
		mrho = info.theDouble;
		// update tension cutoff
		if (mrho == 0.0) {
			mTo = 1e10;
		}
		else {
			mTo = root23*msigma_y / mrho;
		}
	}
	// nonassociativeTerm called
	else if (responseID == 8) {
		mrho_bar = info.theDouble;
	}
	// cohesiveIntercept called
	else if (responseID == 9) {
		msigma_y = info.theDouble;
		// update tension cutoff
		if (mrho == 0.0) {
			mTo = 1e10;
		}
		else {
			mTo = root23*msigma_y / mrho;
		}
	}
	// shearModulus called
	else if (responseID == 10) {
		mG = info.theDouble;
		mCe = mK*mIIvol + 2 * mG*mIIdev;
	}
	// bulkModulus called
	else if (responseID == 11) {
		mK = info.theDouble;
		mCe = mK*mIIvol + 2 * mG*mIIdev;
	}

	return 0;
}

int DruckerPragerThermal::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	// place data in a vector
	static Vector data(45);
	data(0) = this->getTag();
	data(1) = mKref;
	data(2) = mGref;
	data(3) = mK;
	data(4) = mG;
	data(5) = msigma_y;
	data(6) = mrho;
	data(7) = mrho_bar;
	data(8) = mKinf;
	data(9) = mKo;
	data(10) = mdelta1;
	data(11) = mdelta2;
	data(12) = mHard;
	data(13) = mtheta;
	data(14) = massDen;
	data(15) = mPatm;
	data(16) = mTo;
	data(17) = mHprime;
	data(18) = mAlpha1_n;
	data(19) = mAlpha2_n;
	data(20) = mElastFlag;
	data(21) = mFlag;

	data(22) = mEpsilon(0);
	data(23) = mEpsilon(1);
	data(24) = mEpsilon(2);
	data(25) = mEpsilon(3);
	data(26) = mEpsilon(4);
	data(27) = mEpsilon(5);

	data(28) = mEpsilon_n_p(0);
	data(29) = mEpsilon_n_p(1);
	data(30) = mEpsilon_n_p(2);
	data(31) = mEpsilon_n_p(3);
	data(32) = mEpsilon_n_p(4);
	data(33) = mEpsilon_n_p(5);

	data(34) = mBeta_n(0);
	data(35) = mBeta_n(1);
	data(36) = mBeta_n(2);
	data(37) = mBeta_n(3);
	data(38) = mBeta_n(4);
	data(39) = mBeta_n(5);

	data(40) = mState(0);
	data(41) = mState(1);
	data(42) = mState(2);
	data(43) = mState(3);
	data(44) = mState(4);

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DruckerPragerThermal::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int DruckerPragerThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	// receive data
	static Vector data(45);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: DruckerPragerThermal::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}

	// set member variables
	this->setTag((int)data(0));
	mKref = data(1);
	mGref = data(2);
	mK = data(3);
	mG = data(4);
	msigma_y = data(5);
	mrho = data(6);
	mrho_bar = data(7);
	mKinf = data(8);
	mKo = data(9);
	mdelta1 = data(10);
	mdelta2 = data(11);
	mHard = data(12);
	mtheta = data(13);
	massDen = data(14);
	mPatm = data(15);
	mTo = data(16);
	mHprime = data(17);
	mAlpha1_n = data(18);
	mAlpha2_n = data(19);
	mElastFlag = (int)data(20);
	mFlag = (int)data(21);

	mEpsilon(0) = data(22);
	mEpsilon(1) = data(23);
	mEpsilon(2) = data(24);
	mEpsilon(3) = data(25);
	mEpsilon(4) = data(26);
	mEpsilon(5) = data(27);

	mEpsilon_n_p(0) = data(28);
	mEpsilon_n_p(1) = data(29);
	mEpsilon_n_p(2) = data(30);
	mEpsilon_n_p(3) = data(31);
	mEpsilon_n_p(4) = data(32);
	mEpsilon_n_p(5) = data(33);

	mBeta_n(0) = data(34);
	mBeta_n(1) = data(35);
	mBeta_n(2) = data(36);
	mBeta_n(3) = data(37);
	mBeta_n(4) = data(38);
	mBeta_n(5) = data(39);

	mState(0) = data(40);
	mState(1) = data(41);
	mState(2) = data(42);
	mState(3) = data(43);
	mState(4) = data(44);

	mCe = mK*mIIvol + 2 * mG*mIIdev;
	mCep = mCe;

	return 0;
}

void DruckerPragerThermal::Print(OPS_Stream &s, int flag)
{
	s << "DruckerPragerThermal" << endln;
}

