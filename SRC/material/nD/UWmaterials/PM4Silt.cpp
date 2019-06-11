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

// Computational Geomechanics Group
// Author: Long Chen, Pedro Arduino
// University of Washington
// Date:      Sep 2018
// Last Modified: May 2019

// Description: This file contains the implementation for the PM4Silt class.
// PM4Silt(Version 1): A Silt Plasticity Model For Earthquake Engineering Applications
// by R.W.Boulanger and K.Ziotopoulou
// Jan 2018

#include <PM4Silt.h>
#include <MaterialResponse.h>

// #include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <algorithm>
#define fmax std::max
#define fmin std::min
#endif

#define INT_ModifiedEuler 1
#define INT_ForwardEuler  2
#define INT_RungeKutta4   3
#define INT_MAXSTR_FE     4
#define INT_MAXSTR_ME     5

const double		PM4Silt::root12 = sqrt(1.0 / 2.0);
const double		PM4Silt::one3 = 1.0 / 3.0;
const double		PM4Silt::two3 = 2.0 / 3.0;
const double		PM4Silt::root23 = sqrt(2.0 / 3.0);
const double		PM4Silt::small = 1e-10;
const double		PM4Silt::maxStrainInc = 1e-6;
const bool  		PM4Silt::debugFlag = false;
const char unsigned	PM4Silt::mMaxSubStep = 10;
char  unsigned		PM4Silt::me2p = 1;

Vector 			PM4Silt::mI1(3);
Matrix  		PM4Silt::mIIco(3, 3);
Matrix 			PM4Silt::mIIcon(3, 3);
Matrix 			PM4Silt::mIImix(3, 3);
Matrix 			PM4Silt::mIIvol(3, 3);
Matrix 			PM4Silt::mIIdevCon(3, 3);
Matrix 			PM4Silt::mIIdevMix(3, 3);
Matrix 			PM4Silt::mIIdevCo(3, 3);
PM4Silt::initTensors PM4Silt::initTensorOps;

static int numPM4SiltMaterials = 0;

void *
OPS_PM4SiltMaterial(void)
{
	if (numPM4SiltMaterials == 0) {
		numPM4SiltMaterials++;
		opserr << "PM4Silt nDmaterial - Written: L.Chen, P.Arduino, U.Washington\n";
	}

	NDMaterial *theMaterial = 0;

	int numArgs = OPS_GetNumRemainingInputArgs();

	if (numArgs < 6) {
		opserr << "Want: nDMaterial PM4Silt tag? Su? Su_rate? G0? hpo? rho?" << endln;
		return 0;
	}

	int tag;
	double dData[5];
	double oData[24];

	oData[0] = 1.0;      // Fsu
	oData[1] = 101.3;    // p_atm
	oData[2] = 0.3;      // nu
	oData[3] = 0.75;     // nG
	oData[4] = 0.5;      // h0
	oData[5] = 0.9;      // e_init
	oData[6] = 0.06;     // Lambda
	oData[7] = 32.0;     // phi_cv
	oData[8] = 0.8;      // nb_wet
	oData[9] = 0.5;      // nb_dry
	oData[10] = 0.3;     // nd
	oData[11] = 0.8;     // Ado
	oData[12] = -1;      // ru_max
	oData[13] = -1;      // z_max
	oData[14] = 100.0;   // Cz
	oData[15] = -1;      // Ce
	oData[16] = 3.0;     // Cgd
	oData[17] = 4.0;     // Ckaf
	oData[18] = 0.01;    // m
	oData[19] = 2.0;	 // CG_consol
	oData[20] = 1;       // integration scheme
	oData[21] = 0;       // tangent type
	oData[22] = 1.0e-7;	 // TolF
	oData[23] = 1.0e-10;	// TolR

	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING invalid nDMaterial PM4Silt material tag" << endln;
		return 0;
	}

	numData = 5;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial PM4Silt material  with tag: " << tag << endln;
		return 0;
	}

	numData = numArgs - 6;
	if (numData != 0)
		if (OPS_GetDouble(&numData, oData) != 0) {
			opserr << "WARNING invalid material data for nDMaterial PM4Silt material  with tag: " << tag << endln;
			return 0;
		}

	theMaterial = new PM4Silt(tag, ND_TAG_PM4Silt, dData[0], dData[1], dData[2], dData[3], dData[4], oData[0], oData[1], oData[2],
		oData[3], oData[4], oData[5], oData[6], oData[7], oData[8], oData[9], oData[10], oData[11], oData[12], oData[13], oData[14],
		oData[15], oData[16], oData[17], oData[18], oData[19], (int)oData[20], (int)oData[21], oData[22], oData[23]);


	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory for nDMaterial PM4Silt material with tag: " << tag << endln;
	}

	return theMaterial;
}

// full constructor
PM4Silt::PM4Silt(int tag, int classTag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nu, double nG, double h0,
	double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd, double Ado, double ru_max, double z_max, double cz, double ce,
	double Cgd, double Ckaf, double m, double CG_consol, int integrationScheme, int tangentType, double TolF, double TolR) : NDMaterial(tag, classTag),
	mEpsilon(3),
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_b(3),
	mSigma_n(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
	mAlpha_in_p_n(3),
	mAlpha_in_true(3),
	mAlpha_in_true_n(3),
	mAlpha_in_max(3),
	mAlpha_in_max_n(3),
	mAlpha_in_min(3),
	mAlpha_in_min_n(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mFabric_in_n(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3),
	mTracker(3)
{
	m_Su = Su;
	m_Su_rate = Su_rate;
	m_G0 = G0;
	m_hpo = hpo;
	massDen = mDen;
	m_Fsu = Fsu;
	m_P_atm = P_atm;
	m_nu = nu;
	if (m_nu < 0.0) {
		m_nu = 0.3;
	}
	else if (m_nu >= 0.5) {
		opserr << "Warning, Poisson's ratio is larger than 0.5, using 0.49 instead. \n";
		m_nu = 0.49;
	}
	m_nG = nG;
	if (m_nG < 0.0) m_nG = 0.75;
	m_h0 = h0;
	m_e_init = einit;
	m_lambda = lambda;
	if (phi_cv < 0.0) phi_cv = 32.0;
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nbwet = nbwet;
	m_nbdry = nbdry;
	m_nd = nd;
	m_Ado = Ado;
	m_ru_max = ru_max;
	m_z_max = z_max;
	m_cz = cz;
	m_ce = ce;
	m_Cgd = Cgd;
	m_Ckaf = Ckaf;
	m_m = m;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = CG_consol;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;
	mIter = 0;

	initialize();
}

// full constructor
PM4Silt::PM4Silt(int tag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nu, double nG, double h0,
	double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd, double Ado, double ru_max, double z_max, double cz,
	double ce, double Cgd, double Ckaf, double m, double CG_consol, int integrationScheme, int tangentType, double TolF, double TolR)
	: NDMaterial(tag, ND_TAG_PM4Silt),
	mEpsilon(3),
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_n(3),
	mSigma_b(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
	mAlpha_in_p_n(3),
	mAlpha_in_true(3),
	mAlpha_in_true_n(3),
	mAlpha_in_max(3),
	mAlpha_in_max_n(3),
	mAlpha_in_min(3),
	mAlpha_in_min_n(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mFabric_in_n(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3),
	mTracker(3)
{
	m_Su = Su;
	m_Su_rate = Su_rate;
	m_G0 = G0;
	m_hpo = hpo;
	massDen = 1.5;
	m_Fsu = Fsu;
	m_P_atm = P_atm;
	m_nu = nu;
	if (m_nu < 0.0) {
		m_nu = 0.3;
	}
	else if (m_nu >= 0.5) {
		opserr << "Warning, Poisson's ratio is larger than 0.5, using 0.49 instead. \n";
		m_nu = 0.49;
	}
	m_nG = nG;
	if (m_nG < 0.0) m_nG = 0.75;
	m_h0 = h0;
	m_e_init = einit;
	m_lambda = lambda;
	if (phi_cv < 0.0) phi_cv = 32.0;
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nbwet = nbwet;
	m_nbdry = nbdry;
	m_nd = nd;
	m_Ado = Ado;
	m_ru_max = ru_max;
	m_z_max = z_max;
	m_cz = cz;
	m_ce = ce;
	m_Cgd = Cgd;
	m_Ckaf = Ckaf;
	m_m = m;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = CG_consol;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;
	mIter = 0;

	initialize();
}

// null constructor
PM4Silt::PM4Silt()
	: NDMaterial(),
	mEpsilon(3),
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_n(3),
	mSigma_b(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
	mAlpha_in_p_n(3),
	mAlpha_in_true(3),
	mAlpha_in_true_n(3),
	mAlpha_in_max(3),
	mAlpha_in_max_n(3),
	mAlpha_in_min(3),
	mAlpha_in_min_n(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mFabric_in_n(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3),
	mTracker(3)
{
	m_Su = 0.0;
	m_Su_rate = -1.0;
	m_G0 = 0.0;
	m_hpo = 0.0;
	massDen = 1.4;
	m_Fsu = 0.0;
	m_P_atm = 101.3;
	m_nu = 0.3;
	m_nG = 0.75;
	m_h0 = 0.5;
	m_e_init = 0.9;
	m_lambda = 0.06;
	m_Mc = 0.0;
	m_nbwet = 0.8;
	m_nbdry = 0.5;
	m_nd = 0.3;
	m_Ado = 0.8;
	m_ru_max = -1.0;
	m_z_max = -1.0;
	m_cz = 100.0;
	m_ce = -1.0;
	m_Cgd = 3.0;
	m_Ckaf = 4.0;
	m_m = 0.01;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = 2.0;
	mScheme = 1;
	mTangType = 0;
	mTolF = 1.0e-7;
	mTolR = 1.0e-7;

	mIter = 0;

	this->initialize();
}

// destructor
PM4Silt::~PM4Silt()
{
}

NDMaterial*
PM4Silt::getCopy(const char *type)
{
	if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
		PM4Silt *clone;
		double phi_cv = asin(m_Mc / 2.0) * 180.0 / 3.14159265359;
		clone = new PM4Silt(this->getTag(), m_Su, m_Su_rate, m_G0, m_hpo, massDen, m_Fsu, m_P_atm, m_nu, m_nG, m_h0,
			m_e_init, m_lambda, phi_cv, m_nbwet, m_nbdry, m_nd, m_Ado, m_ru_max, m_z_max, m_cz, m_ce, m_Cgd, m_Ckaf, m_m,
			m_CG_consol, mScheme, mTangType, mTolF, mTolR);
		return clone;
	}
	else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
		opserr << "This is a 2D model and it's not compatible with " << type << endln;
		return 0;
	}
	else {
		opserr << "PM4Silt::getCopy failed to get copy: " << type << endln;
		return 0;
	}
}

int
PM4Silt::commitState(void)
{
	Vector n(3), R(3), dFabric(3);

	mAlpha_in_n = mAlpha_in;
	mAlpha_n = mAlpha;
	mAlpha_in_p_n = mAlpha_in_p;
	mAlpha_in_true_n = mAlpha_in_true;
	mAlpha_in_max_n = mAlpha_in_max;
	mAlpha_in_min_n = mAlpha_in_min;
	mSigma_n = mSigma;
	mEpsilon_n = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	dFabric = mFabric - mFabric_n;
	// update cumulated fabric
	mzcum = mzcum + sqrt(DoubleDot2_2_Contr(dFabric, dFabric) / 2.0);
	mzpeak = fmax(sqrt(DoubleDot2_2_Contr(mFabric, mFabric) / 2.0), mzpeak);
	mFabric_n = mFabric;
	mFabric_in_n = mFabric_in;
	mDGamma_n = mDGamma;
	mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);

	this->GetElasticModuli(mSigma, mK, mG, mMcur, mzcum);
	mCe = GetStiffness(mK, mG);
	mCep = GetElastoPlasticTangent(mSigma_n, mCe, R, n, mKp);
	mCep_Consistent = mCe;
	return 0;
}

int PM4Silt::revertToLastCommit(void)
{
	// need to be added
	return 0;
}

int PM4Silt::revertToStart(void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	}
	else {
		// normal call for revertToStart (not initialStateAnalysis)
		this->initialize(mSigma);
	}

	return 0;
}

NDMaterial*
PM4Silt::getCopy(void)
{
	opserr << "PM4Silt::getCopy -- subclass responsibility\n";
	exit(-1);
	return 0;
}

const char*
PM4Silt::getType(void) const
{
	opserr << "PM4Silt::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int
PM4Silt::getOrder(void) const
{
	opserr << "PM4Silt::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}


Response*
PM4Silt::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else if (strcmp(argv[0], "alpha") == 0 || strcmp(argv[0], "backstressratio") == 0)
		return new MaterialResponse(this, 4, this->getAlpha());
	else if (strcmp(argv[0], "fabric") == 0)
		return new MaterialResponse(this, 5, this->getFabric());
	else if (strcmp(argv[0], "alpha_in") == 0 || strcmp(argv[0], "alphain") == 0)
		return new MaterialResponse(this, 6, this->getAlpha_in());
	else if (strcmp(argv[0], "trackers") == 0 || strcmp(argv[0], "tracker") == 0)
		return new MaterialResponse(this, 7, this->getTracker());
	else
		return 0;
}

int
PM4Silt::getResponse(int responseID, Information &matInfo)
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
	case 4:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getAlpha();
		return 0;
	case 5:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getFabric();
		return 0;
	case 6:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getAlpha_in();
		return 0;
	case 7:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getTracker();
		return 0;
	default:
		return -1;
	}
}

int
PM4Silt::sendSelf(int commitTag, Channel &theChannel)
{

	int res = 0;
	static Vector data(105);

	data(0) = this->getTag();

	data(1) = m_Su;
	data(2) = m_Su_rate;
	data(3) = m_G0;
	data(4) = m_hpo;
	data(5) = massDen;
	data(6) = m_Fsu;
	data(7) = m_P_atm;
	data(8) = m_nG;
	data(9) = m_h0;
	data(10) = m_e_init;
	data(11) = m_lambda;
	data(12) = m_nbwet;
	data(13) = m_nbdry;
	data(14) = m_nd;
	data(15) = m_Ado;
	data(16) = m_ru_max;
	data(17) = m_z_max;
	data(18) = m_cz;
	data(19) = m_ce;
	data(20) = m_Mc;
	data(21) = m_Cgd;
	data(22) = m_Ckaf;
	data(23) = m_nu;
	data(24) = m_m;
	data(25) = m_CG_consol;
	data(26) = m_FirstCall;
	data(27) = m_PostShake;

	data(28) = mTolF;
	data(29) = mTolR;
	data(30) = mScheme;
	data(31) = mTangType;
	data(32) = m_Pmin;
	data(35) = m_pzpFlag;
	data(36) = me2p;

	data(37) = mDGamma;
	data(38) = mDGamma_n;
	data(39) = me0;
	data(40) = mpcs;
	data(41) = mK;
	data(42) = mG;
	data(43) = mVoidRatio;
	data(44) = mKp;
	data(45) = mzcum;
	data(46) = mzpeak;
	data(47) = mpzp;
	data(48) = mzxp;
	data(49) = mMb;
	data(50) = mMb_max;
	data(51) = mC_MB;
	data(52) = mMd;
	data(53) = mMcur;

	data(54) = mEpsilon(0);		  data(57) = mEpsilon_n(0);	    data(60) = mSigma(0);	data(63) = mSigma_n(0);   data(66) = mSigma_b(0);
	data(55) = mEpsilon(1);		  data(58) = mEpsilon_n(1);	    data(61) = mSigma(1);	data(64) = mSigma_n(1);	  data(67) = mSigma_b(1);
	data(56) = mEpsilon(2);		  data(59) = mEpsilon_n(2);	    data(62) = mSigma(2);	data(65) = mSigma_n(2);	  data(68) = mSigma_b(2);

	data(69) = mEpsilonE(0);	  data(72) = mEpsilonE_n(0);	data(75) = mAlpha(0);	data(78) = mAlpha_n(0);   data(81) = mAlpha_in_n(0);
	data(70) = mEpsilonE(1);	  data(73) = mEpsilonE_n(1);	data(76) = mAlpha(1);	data(79) = mAlpha_n(1);	  data(82) = mAlpha_in_n(1);
	data(71) = mEpsilonE(2);	  data(74) = mEpsilonE_n(2);	data(77) = mAlpha(2);	data(80) = mAlpha_n(2);	  data(83) = mAlpha_in_n(2);

	data(84) = mAlpha_in_p_n(0);  data(87) = mAlpha_in_true_n(0);    data(90) = mAlpha_in_max_n(0);      data(93) = mAlpha_in_min_n(0);
	data(85) = mAlpha_in_p_n(1);  data(88) = mAlpha_in_true_n(1);    data(91) = mAlpha_in_max_n(1);      data(94) = mAlpha_in_min_n(1);
	data(86) = mAlpha_in_p_n(2);  data(89) = mAlpha_in_true_n(2);    data(92) = mAlpha_in_max_n(2);      data(95) = mAlpha_in_min_n(2);

	data(96) = mFabric(0);		data(99) = mFabric_n(0);	 data(102) = mFabric_in_n(0);
	data(97) = mFabric(1);		data(100) = mFabric_n(1);	 data(103) = mFabric_in_n(1);
	data(98) = mFabric(2);		data(101) = mFabric_n(2);	 data(104) = mFabric_in_n(2);

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: PM4Silt::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int
PM4Silt::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(105);

	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: PM4Silt::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}
	this->setTag((int)data(0));

	m_Su = data(1);
	m_Su_rate = data(2);
	m_G0 = data(3);
	m_hpo = data(4);
	massDen = data(5);
	m_Fsu = data(6);
	m_P_atm = data(7);
	m_nG = data(8);
	m_h0 = data(9);
	m_e_init = data(10);
	m_lambda = data(11);
	m_nbwet = data(12);
	m_nbdry = data(13);
	m_nd = data(14);
	m_Ado = data(15);
	m_ru_max = data(16);
	m_z_max = data(17);
	m_cz = data(18);
	m_ce = data(19);
	m_Mc = data(20);
	m_Cgd = data(21);
	m_Ckaf = data(22);
	m_nu = data(23);
	m_m = data(24);
	m_CG_consol = data(25);
	m_FirstCall = data(26);
	m_PostShake = data(27);

	mTolF = data(28);
	mTolR = data(29);
	mScheme = data(30);
	mTangType = data(31);
	m_Pmin = data(32);
	m_pzpFlag = data(35);
	me2p = data(36);

	mDGamma = data(37);
	mDGamma_n = data(38);
	me0 = data(39);
	mpcs = data(40);
	mK = data(41);
	mG = data(42);
	mVoidRatio = data(43);
	mKp = data(44);
	mzcum = data(45);
	mzpeak = data(46);
	mpzp = data(47);
	mzxp = data(48);
	mMb = data(49);
	mMb_max = data(50);
	mC_MB = data(51);
	mMd = data(52);
	mMcur = data(53);

	mEpsilon(0) = data(54); 		mEpsilon_n(0) = data(57); 	    mSigma(0) = data(60); 	  mSigma_n(0) = data(63);       mSigma_b(0) = data(66);
	mEpsilon(1) = data(55); 		mEpsilon_n(1) = data(58); 	    mSigma(1) = data(61); 	  mSigma_n(1) = data(64); 	     mSigma_b(1) = data(67);
	mEpsilon(2) = data(56); 		mEpsilon_n(2) = data(59); 	    mSigma(2) = data(62); 	  mSigma_n(2) = data(65); 	     mSigma_b(2) = data(68);

	mEpsilonE(0) = data(69); 	    mEpsilonE_n(0) = data(72); 	  mAlpha(0) = data(75); 	   mAlpha_n(0) = data(78);     mAlpha_in_n(0) = data(81);
	mEpsilonE(1) = data(70); 	    mEpsilonE_n(1) = data(73); 	  mAlpha(1) = data(76); 	   mAlpha_n(1) = data(79);	    mAlpha_in_n(1) = data(82);
	mEpsilonE(2) = data(71); 	    mEpsilonE_n(2) = data(74); 	  mAlpha(2) = data(77); 	   mAlpha_n(2) = data(80);	    mAlpha_in_n(2) = data(83);

	mAlpha_in_p_n(0) = data(84);    mAlpha_in_true_n(0) = data(87);       mAlpha_in_max_n(0) = data(90);       mAlpha_in_min_n(0) = data(93);
	mAlpha_in_p_n(1) = data(85);    mAlpha_in_true_n(1) = data(88);       mAlpha_in_max_n(1) = data(91);       mAlpha_in_min_n(1) = data(94);
	mAlpha_in_p_n(2) = data(86);    mAlpha_in_true_n(2) = data(89);       mAlpha_in_max_n(2) = data(92);       mAlpha_in_min_n(2) = data(95);

	mFabric(0) = data(96);		  mFabric_n(0) = data(99);  	   mFabric_in_n(0) = data(102);
	mFabric(1) = data(97);		  mFabric_n(1) = data(100); 	   mFabric_in_n(1) = data(103);
	mFabric(2) = data(98);		  mFabric_n(2) = data(101); 	   mFabric_in_n(2) = data(104);

	return 0;
}

void PM4Silt::Print(OPS_Stream &s, int flag)
{
	s << "PM4Silt Material, tag: " << this->getTag() << endln;
	s << "Type: " << this->getType() << endln;
}

int
PM4Silt::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {     // enforce elastic/elastoplastic response
			opserr << this->getTag() << " update Material Stage\n";
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0], "materialState") == 0) {     // enforce elastic/elastoplastic response
			return param.addObject(5, this);
		}
		else if (strcmp(argv[0], "IntegrationScheme") == 0) { // change integration scheme (Explicit/Implicit)
			return param.addObject(2, this);
		}
		else if ((strcmp(argv[0], "refShearModulus") == 0) ||
			(strcmp(argv[0], "ShearModulus") == 0)) {         // change G0
			return param.addObject(6, this);
		}
		else if (strcmp(argv[0], "poissonRatio") == 0) {      // change nu
			return param.addObject(7, this);
		}
		else if (strcmp(argv[0], "FirstCall") == 0) {       // update first call, remove fabric
			return param.addObject(8, this);
		}
		else if (strcmp(argv[0], "voidRatio") == 0) {        // change e_init
			return param.addObject(9, this);
		}
		else if (strcmp(argv[0], "PostShake") == 0) {        // activate post shaking reconsolidation 
			return param.addObject(13, this);
		}
		else if (strcmp(argv[0], "Su_factor") == 0) {        // Undraned shear strength reduction factor
			return param.addObject(14, this);
		}
	}
	return -1;
}

int
PM4Silt::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		me2p = info.theInt;
	}
	// called materialState in tcl file
	else if (responseID == 5) {
		me2p = (int)info.theDouble;
	}
	// called update IntegrationScheme
	else if (responseID == 2) {
		mScheme = (int)info.theDouble;
	}
	// called update refShearModulus
	else if (responseID == 6) {
		m_G0 = info.theDouble;
	}
	// called update poissonRatio
	else if (responseID == 7) {
		m_nu = info.theDouble;
	}
	//called update first call
	else if (responseID == 8) {
		m_FirstCall = 0;
		initialize(mSigma_n);
		opserr << this->getTag() << " initialize" << endln;
	}
	// called update voidRatio
	else if (responseID == 9) {
		double eps_v = GetTrace(mEpsilon);
		m_e_init = (info.theDouble + eps_v) / (1 - eps_v);
	}
	// called PostShake
	else if (responseID == 13) {
		m_PostShake = 1;
		// mElastFlag = 1;
		GetElasticModuli(mSigma, mK, mG, mMcur, mzcum);
		opserr << this->getTag() << " activate post shaking reconsolidation" << endln;
	}
	// update undrained shear strength reduction factor Fsu
	else if (responseID == 14) {
		m_Fsu = info.theDouble;
	}
	else {
		return -1;
	}

	return 0;
}

int
PM4Silt::initialize(Vector initStress)
{
	double p0;
	p0 = 0.5 * GetTrace(initStress);
	// check secondary input parameters, use default values if negative values are assigned
	if (m_Fsu <= 0.0) m_Fsu = 1.0;
	if (m_h0 < 0.0) m_h0 = 0.5;
	if (m_e_init < 0.0) m_e_init = 0.9;
	if (m_lambda < 0.0) m_lambda = 0.06;
	if (m_nbwet < 0.0) m_nbwet = 0.8;
	if (m_nbdry < 0.0) m_nbdry = 0.5;
	if (m_nd < 0.0) m_nd = 0.3;
	if (m_Ado < 0.0) m_Ado = 0.8;
	if (m_cz < 0.0) m_cz = 100.0;
	if (m_Cgd < 0.0) m_Cgd = 3.0;
	if (m_Ckaf < 0.0) m_Ckaf = 4.0;
	if (m_m < 0.0) m_m = 0.01;
	if (m_CG_consol < 0.0) m_CG_consol = 2.0;
	// set phi_max as 60 degree
	mMb_max = 2 * sin(60.0 / 180.0 * 3.14159265359);
	mC_MB = 1 / (pow(mMb_max / m_Mc, 1 / m_nbdry) - 1.0);
	//
	if (p0 < m_P_atm / 200.0) {
		//stress tensile, increase residualP to prevent negative initial p
		if (debugFlag)
			opserr << "Warning, initial p is small. \n";
		p0 = m_P_atm / 200.0;
		mSigma_n = p0 * mI1;
		mSigma_b = initStress - mSigma_n;
		mAlpha.Zero();
		mAlpha_n.Zero();
	}
	else {
		mSigma_n = initStress;
		mSigma_b.Zero();
		mAlpha_n = GetDevPart(initStress) / p0;
	}
	if (m_Su <= 0.0) {
		m_Su = m_Su_rate * initStress(1);
	}
	else {
		// if both Su and Su_rate are given, only Su will be used
		m_Su_rate = m_Su / initStress(1);
	}
	mpcs = 2 * m_Su / m_Mc;
	if (m_ru_max < 0.0)
		// minimum p'
		m_Pmin = fmin(p0, mpcs / 8.0);
	else {
		m_ru_max = fmin(0.99, m_ru_max);
		m_Pmin = (1 - m_ru_max) * p0 / 2.0;
	}
	m_Pmin = fmax(m_Pmin, m_P_atm / 200.0);
	// initialize zmax
	if (m_z_max < 0) {
		if (m_Su_rate <= 0.25)
			m_z_max = 10.0;
		else if (m_Su_rate > 0.50)
			m_z_max = 20.0;
		else
			m_z_max = 40 * m_Su_rate;
	}
	// initialize Ce
	if (m_ce < 0.0)
		m_ce = fmin(1.3, 0.5 + 1.2 * Macauley(m_Su_rate - 0.25));
	// positioning the critical state line
	me0 = m_e_init + m_lambda * log(101.3 * 2 * m_Su / m_Mc / m_P_atm);
	double ksi = GetKsi(m_e_init, p0);
	mMd = fmin(1.4142136, m_Mc * exp(m_nd * ksi / m_lambda));
	if (ksi < 0) {
		// dense of critical
		mMb = m_Mc * pow((1 + mC_MB) / (p0 / mpcs + mC_MB), m_nbdry);
	}
	else {
		//loose of critical
		mMb = m_Mc * exp(-1.0 * m_nbwet * ksi / m_lambda);
	}

	// check if initial stresses are inside bounding/dilatancy surface 
	double Mcut = fmax(mMb, mMd);
	double Mfin = sqrt(2.0) * GetNorm_Contr(GetDevPart(mSigma_n));
	Mfin = Mfin / p0;
	if (Mfin > Mcut)
	{
		Vector r = (mSigma_n - p0 * mI1) / p0 * Mcut / Mfin;
		mSigma_n = p0 * mI1 + r * p0;
		mSigma_b = initStress - mSigma_n;
		mAlpha_n = r * (Mcut - m_m) / Mcut;
	}
	mzcum = 0.0;
	GetElasticModuli(mSigma_n, mK, mG, mMcur, mzcum);
	mCe = mCep = mCep_Consistent = GetStiffness(mK, mG);
	mKp = 100 * mG;
	mAlpha = mAlpha_n;
	mAlpha_in_n = mAlpha_n;
	mAlpha_in_p = mAlpha_n;
	mAlpha_in_p_n.Zero();
	mAlpha_in_true = mAlpha_n;
	mAlpha_in_true_n = mAlpha_n;
	mAlpha_in_max = mAlpha_n;
	mAlpha_in_max_n = mAlpha_n;
	mAlpha_in_min = mAlpha_n;
	mAlpha_in_min_n = mAlpha_n;
	mFabric.Zero();
	mFabric_in.Zero();
	mFabric_in_n.Zero();
	mFabric_n.Zero();
	mzpeak = m_z_max / 100000.0;
	mpzp = fmax(p0, m_Pmin) / 100.0;
	mzxp = 0.0;
	m_pzpFlag = true;
	mTracker.Zero();

	return 0;
}

// Initialize PM4Silt Material
int
PM4Silt::initialize()
{
	// set Initial parameters with p = p_atm
	Vector mSig(3);
	m_Pmin = m_P_atm / 200.0;
	mSig(0) = m_P_atm;
	mSig(1) = m_P_atm;
	mSig(2) = 0.0;

	mzcum = 0.0;
	mzpeak = m_z_max / 100000.0;
	GetElasticModuli(mSig, mK, mG);
	mCe = mCep = mCep_Consistent = GetStiffness(mK, mG);

	return 0;
}

int
PM4Silt::setTrialStrain(const Vector &strain_from_element) {
	// mEpsilon = -1.0 * strain_from_element;   // -1.0 is for geotechnical sign convention
	mEpsilon = strain_from_element; mEpsilon *= (-1.0);
	integrate();
	return 0;
}

// unused trial strain functions
int
PM4Silt::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

//send back the state parameters to the recorders
const Vector
PM4Silt::getState()
{
	Vector result(16);
	result.Assemble(mEpsilonE, 0, 1.0);
	result.Assemble(mAlpha, 3, 1.0);
	result.Assemble(mFabric, 6, 1.0);
	result.Assemble(mAlpha_in, 9, 1.0);
	result(12) = mVoidRatio;
	result(13) = mDGamma;
	result(14) = mG;
	result(15) = mKp;

	return result;
}
//send back alpha tensor
const Vector
PM4Silt::getAlpha()
{
	return mAlpha_n;
}
//send back fabric tensor
const Vector
PM4Silt::getFabric()
{
	return mFabric_n;
}
//send back alpha_in tensor
const Vector
PM4Silt::getAlpha_in()
{
	return mAlpha_in_n;
}
//send back internal parameter for tracking
const Vector
PM4Silt::getTracker()
{
	return mTracker;
}
//send back Kp
double
PM4Silt::getKp()
{
	return mKp;
}
//send back shear modulus
double
PM4Silt::getG()
{
	return mG;
}

//send back previous alpha_in tensor
const Vector
PM4Silt::getAlpha_in_p()
{
	return mAlpha_in_p_n;
}
//send back previous L
double
PM4Silt::getDGamma()
{
	return mDGamma;
}
/*************************************************************/
const Matrix&
PM4Silt::getTangent() {
	if (mTangType == 0)
		return mCe;
	else if (mTangType == 1)
		return mCep;
	else return mCep_Consistent;
}
/*************************************************************/
const Matrix &
PM4Silt::getInitialTangent() {
	return mCe;
}
/*************************************************************/
const Vector &
PM4Silt::getStress() {
	// mSigma_r = -1.0 * (mSigma + mSigma_b);
	mSigma_r = mSigma; mSigma_r += mSigma_b; mSigma_r *= (-1.0);
	return  mSigma_r;  // -1.0 is for geotechnical sign convention
}
/*************************************************************/
const Vector &
PM4Silt::getStrain() {
	// mEpsilon_r = -1.0 * mEpsilon;   // -1.0 is for geotechnical sign convention
	mEpsilon_r = mEpsilon; mEpsilon_r *= (-1.0);
	return mEpsilon_r;
}
/*************************************************************/
const Vector &
PM4Silt::getElasticStrain() {
	// mEpsilonE_r = -1.0 * mEpsilonE;   // -1.0 is for geotechnical sign convention
	mEpsilonE_r = mEpsilon; mEpsilonE_r *= (-1.0);
	return mEpsilonE_r;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void PM4Silt::integrate()
{
	mAlpha = mAlpha_n;
	mAlpha_in = mAlpha_in_n;
	mAlpha_in_true = mAlpha_in_true_n;
	mAlpha_in_p = mAlpha_in_p_n;
	mAlpha_in_max = mAlpha_in_max_n;
	mAlpha_in_min = mAlpha_in_min_n;
	mFabric = mFabric_n;
	mFabric_in = mFabric_in_n;

	Vector n_tr(3), tmp0(3), tmp1(3), mAlpha_mAlpha_in_true(3);
	// n_tr = GetNormalToYield(mSigma_n + mCe*(mEpsilon - mEpsilon_n), mAlpha);
	tmp0 += mSigma_n; tmp1 = mEpsilon; tmp1 -= mEpsilon_n;
	tmp0 += (mCe * tmp1);
	n_tr = GetNormalToYield(tmp0, mAlpha);
	// n_tr = GetNormalToYield(mSigma_n, mAlpha);

	// if ((DoubleDot2_2_Contr(mAlpha - mAlpha_in_true, n_tr) < 0.0) && me2p) {
	mAlpha_mAlpha_in_true = mAlpha; mAlpha_mAlpha_in_true -= mAlpha_in_true;
	if ((DoubleDot2_2_Contr(mAlpha_mAlpha_in_true, n_tr) < 0.0) && me2p) {
		mAlpha_in_p = mAlpha_in;
		mAlpha_in_true = mAlpha;
		mFabric_in = mFabric;
		// This is a loading reversal
		// update pzp
		double p = 0.5 * GetTrace(mSigma_n);
		p = (p <= m_Pmin) ? (m_Pmin) : p;
		double zxpTemp = GetNorm_Contr(mFabric_n) * p;
		if (((zxpTemp > mzxp) && (p > mpzp)) || m_pzpFlag) {
			mzxp = zxpTemp;
			mpzp = p;
			m_pzpFlag = false;
		}
		// track initial back-stress ratio history 
		for (int ii = 0; ii < 3; ii++) {
			if (mAlpha_in(ii) > 0.0)
				// minimum positive value
				mAlpha_in_min(ii) = fmin(mAlpha_in_min(ii), mAlpha(ii));
			else
				// maximum negative value
				mAlpha_in_max(ii) = fmax(mAlpha_in_max(ii), mAlpha(ii));
		}
		if (mAlpha(2) * mAlpha_in_p(2) > 0) {
			for (int ii = 0; ii < 3; ii++) {
				if (n_tr(ii) > 0.0)
					// positive loading direction
					mAlpha_in(ii) = fmax(0.0, mAlpha_in_min(ii));
				else
					// negative loading direction
					mAlpha_in(ii) = fmin(0.0, mAlpha_in_max(ii));
			}
		}
		else
			mAlpha_in = mAlpha;
	}

	// Force elastic response
	if (me2p == 0) {
		elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha,
			mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent);
	}
	// ElastoPlastic response
	else {
		// explicit schemes
		explicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in,
			mAlpha_in_p, mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG,
			mK, mCe, mCep, mCep_Consistent);
	}

}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Elastic Integrator
/*************************************************************/
void PM4Silt::elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
	double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	Vector dStrain(3);

	// calculate elastic response
	// dStrain = NextStrain - CurStrain;
	dStrain = NextStrain; dStrain -= CurStrain;
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	// NextElasticStrain = CurElasticStrain + dStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;
	aCep_Consistent = aCep = aC = GetStiffness(K, G);
	// NextStress = CurStress + DoubleDot4_2(aC, dStrain);
	NextStress = CurStress; NextStress += DoubleDot4_2(aC, dStrain);
	double p = 0.5 * GetTrace(NextStress);
	if (p > m_Pmin) {
		NextAlpha = GetDevPart(NextStress) / p;
	}

}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Explicit Integrator
/*************************************************************/
void PM4Silt::explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// function pointer to the integration scheme
	void (PM4Silt::*exp_int) (const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
		const Vector&, const Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&, double&, double&,
		Matrix&, Matrix&, Matrix&);

	switch (mScheme) {
	case INT_ForwardEuler:	// Forward Euler
		exp_int = &PM4Silt::ForwardEuler;
		break;

	case INT_ModifiedEuler:	// Modified Euler with error control
		exp_int = &PM4Silt::ModifiedEuler;
		break;

	case INT_RungeKutta4:  // 4th order Ruge-Kutta
		exp_int = &PM4Silt::RungeKutta4;
		break;
	case INT_MAXSTR_FE:   // Forward Euler constraining maximum strain increment
	case INT_MAXSTR_ME:	 // Modified Euler constraining maximum strain increment
		exp_int = &PM4Silt::MaxStrainInc;
		break;
	default:
		exp_int = &PM4Silt::MaxStrainInc;
		break;
	}

	double elasticRatio, f, fn, dVolStrain;
	Vector dStrain(3), dSigma(3), dDevStrain(3), n(3), tmp(3), dElasStrain(3);

	// NextElasticStrain = CurElasticStrain + NextStrain - CurStrain;
	// dVolStrain = GetTrace(NextStrain - CurStrain);
	// dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
	dStrain = NextStrain; dStrain -= CurStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;
	dVolStrain = GetTrace(dStrain);
	dDevStrain += mI1; dDevStrain *= (-1.0 * dVolStrain / 3.0); dDevStrain += dStrain;

	aC = GetStiffness(K, G);
	// dSigma = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
	tmp = ToContraviant(dDevStrain); tmp *= (2 * mG);
	dSigma = mI1; dSigma *= (mK * dVolStrain); dSigma += tmp;
	// NextStress = CurStress + dSigma;
	NextStress = CurStress; NextStress += dSigma;

	f = GetF(NextStress, CurAlpha);

	fn = GetF(CurStress, CurAlpha);

	n = GetNormalToYield(NextStress, CurAlpha);

	if (f <= mTolF)
	{
		// This is a pure elastic loading/unloading
		NextAlpha = CurAlpha;
		NextFabric = CurFabric;
		NextL = 0;
		aCep_Consistent = aCep = aC;
		// Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, 
		// NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, G, K , aC, aCep, aCep_Consistent);

		return;

	}
	else if (fn < -mTolF) {
		// This is a transition from elastic to plastic
		elasticRatio = IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, 0.0, 1.0);
		// dSigma = DoubleDot4_2(aC, elasticRatio*(dStrain));
		dElasStrain = dStrain; dElasStrain *= elasticRatio;
		dSigma = DoubleDot4_2(aC, dElasStrain);
		// (this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
		// 	CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio,
		// 	G, K, aC, aCep, aCep_Consistent);
		(this->*exp_int)(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain,
			CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio,
			G, K, aC, aCep, aCep_Consistent);

		return;
	}
	else if (fabs(fn) < mTolF) {
		if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha), dSigma) / (GetNorm_Contr(dSigma) == 0 ? 1.0 : GetNorm_Contr(dSigma)) > (-sqrt(mTolF))) {
			// This is a pure plastic step
			(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha,
				NextFabric, NextL, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

			return;
		}
		else {
			// This is an elastic unloding followed by plastic loading
			elasticRatio = IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha);
			// dSigma = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
			dElasStrain = dStrain; dElasStrain *= elasticRatio;
			dSigma = DoubleDot4_2(aC, dElasStrain);
			(this->*exp_int)(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain,
				CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio,
				G, K, aC, aCep, aCep_Consistent);

			return;
		}
	}
	else {
		// This is an illegal stress state! This shouldn't happen.
		if (debugFlag) opserr << "PM4Silt : Encountered an illegal stress state! Tag: " << this->getTag() << endln;
		if (debugFlag) opserr << "                  f = " << GetF(CurStress, CurAlpha) << endln;
		(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha,
			NextFabric, NextL, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
		return;
	}
}

// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Forward-Euler Integrator
/*************************************************************/
void PM4Silt::ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	double CurVoidRatio, Cka, h, p, dVolStrain, D, AlphaAlphaBDotN;
	Vector n(3), R(3), alphaD(3), dPStrain(3), b(3), dDevStrain(3), r(3), dStrain(3);
	Vector dSigma(3), dAlpha(3), dFabric(3);

	this->GetElasticModuli(NextStress, K, G, mMcur, mzcum);
	CurVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	p = 0.5 * GetTrace(CurStress);
	p = p < m_Pmin ? m_Pmin : p;
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	dStrain = NextStrain; dStrain -= CurStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;
	// using NextStress instead of CurStress to get correct n
	GetStateDependent(NextStress, CurAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum,
		mzpeak, mpzp, mMcur, CurVoidRatio, n, D, R, mKp, alphaD, Cka, h, b, AlphaAlphaBDotN);
	// dVolStrain = GetTrace(NextStrain - CurStrain);
	dVolStrain = GetTrace(dStrain);
	// dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
	dDevStrain = mI1;
	dDevStrain *= (-1.0 * dVolStrain / 3.0);
	dDevStrain += dStrain;
	// r = GetDevPart(CurStress) / p;
	r = GetDevPart(NextStress);  r /= p;
	double temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
	// if (temp4 < 0.0) {
	// 	mKp = -0.5 * (2 * G - K* D *DoubleDot2_2_Contr(n, r));
	// 	temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
	// 	h = 1.5 * mKp / (p * AlphaAlphaBDotN);
	// }
	if (fabs(temp4) < small) {
		// Neutral loading
		dSigma.Zero();
		dAlpha.Zero();
		dFabric.Zero();
		// dPStrain = dDevStrain + dVolStrain * mI1;
		dPStrain = dStrain;
	}
	else {
		NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4;
		mDGamma = NextL;
		if (NextL < 0) {
			if (debugFlag) {
				opserr << "NextL is smaller than 0\n";
				opserr << "NextL = " << NextL << endln;
			}
			dSigma = 2 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1;
			// dAlpha = GetDevPart(NextStress + dSigma) / (0.5 * GetTrace(NextStress + dSigma))
			// 	- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress));
			dAlpha.Zero();
			dFabric.Zero();
			dPStrain.Zero();
		}
		else {
			// dSigma = 2.0*mG*mIIcon*dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
			// 	(2.0 * mG * n + mK * D * mI1);
			Vector tmp0(n), tmp1(mI1), tmp2(mI1);
			tmp0 *= (2.0 * G);
			tmp1 *= (K * D); tmp1 += tmp0; tmp1 *= (-Macauley(NextL));
			tmp2 *= (K * dVolStrain);
			dSigma = ToContraviant(dDevStrain); dSigma *= (2.0 * G);
			dSigma += tmp2; dSigma += tmp1;
			// update fabric
			if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
				// dFabric = m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
				dFabric = n;
				dFabric *= m_z_max;
				dFabric += CurFabric;
				dFabric *= (-1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL) * MacauleyIndex(-D));
			}
			// update alpha
			// dAlpha = two3 * NextL * h * b;
			// dPStrain = NextL * mIIco * R;
			dPStrain = ToCovariant(R);  dPStrain *= NextL;
			dAlpha = b; dAlpha *= (two3 * NextL * h);
		}
	}
	// NextFabric = CurFabric + dFabric;
	NextFabric = CurFabric; NextFabric += dFabric;
	// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
	NextElasticStrain -= dPStrain;
	// NextStress = CurStress + dSigma;
	NextStress = CurStress;  NextStress += dSigma;
	// NextAlpha = CurAlpha + dAlpha;
	NextAlpha = CurAlpha;  NextAlpha += dAlpha;
	Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
	// Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, R, n, r);
	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Integrator Constraining Maximum Strain Increment
/*************************************************************/
void PM4Silt::MaxStrainInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// function pointer to the integration scheme
	void (PM4Silt::*exp_int) (const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
		const Vector&, const Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&, double&, double&,
		Matrix&, Matrix&, Matrix&);

	switch (mScheme)
	{
	case INT_MAXSTR_FE:
		exp_int = &PM4Silt::ForwardEuler;
		break;
	case INT_MAXSTR_ME:
		exp_int = &PM4Silt::ModifiedEuler;
		break;
	default:
		exp_int = &PM4Silt::ModifiedEuler;
		break;
	}
	Vector StrainInc(3); StrainInc = NextStrain;  StrainInc -= CurStrain;
	double maxInc = StrainInc(0);

	for (int ii = 1; ii < 3; ii++)
		if (fabs(StrainInc(ii)) > fabs(maxInc))
			maxInc = StrainInc(ii);

	if (fabs(maxInc) > maxStrainInc) {
		int numSteps = (int)floor(fabs(maxInc) / maxStrainInc) + 1;
		// StrainInc = (NextStrain - CurStrain) / (double)numSteps;
		StrainInc = NextStrain; StrainInc -= CurStrain; StrainInc /= (double)numSteps;

		Vector cStress(3), cStrain(3), cAlpha(3), cFabric(3), cAlpha_in(3), cAlpha_in_p(3), cEStrain(3);
		Vector nStrain(3);
		Matrix nCe(3, 3), nCep(3, 3), nCepC(3, 3);
		double nL, nVoidRatio, nG, nK;

		// create temporary variables
		cStress = CurStress; cStrain = CurStrain; cAlpha = CurAlpha; cFabric = CurFabric;
		cAlpha_in = alpha_in; cAlpha_in_p = alpha_in_p; cEStrain = CurElasticStrain;


		for (int ii = 1; ii <= numSteps; ii++)
		{
			nStrain = cStrain + StrainInc;

			(this->*exp_int)(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, cAlpha_in_p, nStrain, NextElasticStrain, NextStress, NextAlpha,
				NextFabric, nL, nVoidRatio, nG, nK, nCe, nCep, nCepC);

			cStress = NextStress; cStrain = nStrain; cEStrain = NextElasticStrain;  cAlpha = NextAlpha; cFabric = NextFabric;
		}

	}
	else {
		(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, alpha_in_p, NextStrain, NextElasticStrain, NextStress, NextAlpha,
			NextFabric, NextL, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
	}
	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Modified-Euler Integrator
/*************************************************************/
void PM4Silt::ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	double dVolStrain, p, Cka, temp4, curStepError, q, stressNorm, h, D, AlphaAlphaBDotN;
	Vector n(3), R1(3), R2(3), alphaD(3), dDevStrain(3), r(3), b(3), tmp0(3), tmp1(3), tmp2(3), alphaD_NextAlpha(3);
	Vector nStress(3), nAlpha(3), nFabric(3);
	Vector dSigma1(3), dSigma2(3), dAlpha1(3), dAlpha2(3), dFabric1(3), dFabric2(3), dPStrain1(3), dPStrain2(3);
	double T = 0.0, dT = 1.0, dT_min = 1e-4, TolE = 1e-5;

	// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextElasticStrain = CurElasticStrain; NextElasticStrain += NextStrain; NextElasticStrain -= CurStrain;
	NextStress = CurStress;
	NextAlpha = CurAlpha;
	NextFabric = CurFabric;

	this->GetElasticModuli(NextStress, K, G, mMcur, mzcum);

	p = 0.5 * GetTrace(CurStress);
	if (p < m_Pmin / 5.0)
	{
		if (debugFlag)
			opserr << "Tag = " << this->getTag() << " : p < pmin / 5, should not happen" << endln;
		NextStress = GetDevPart(NextStress) + m_Pmin / 5.0 * mI1;
	}
	while (T < 1.0)
	{
		//NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + T*(NextStrain - CurStrain));
		tmp0 = NextStrain; tmp0 -= CurStrain; tmp0 *= T; tmp0 += CurStrain;
		NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(tmp0);
		//dVolStrain = dT * GetTrace(NextStrain - CurStrain);
		tmp0 = NextStrain; tmp0 -= CurStrain;
		dVolStrain = dT * GetTrace(tmp0);
		// dDevStrain = dT * (NextStrain - CurStrain)-dVolStrain / 3.0 * mI1;
		dDevStrain = mI1;
		dDevStrain *= (-1.0 * dVolStrain / 3.0);
		tmp0 *= dT;
		dDevStrain += tmp0;

		p = 0.5 * GetTrace(NextStress);
		// Calc Delta 1
		GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, NextFabric, mFabric_in, G, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R1, mKp, alphaD, Cka, h, b, AlphaAlphaBDotN);

		// r += GetDevPart(NextStress) / p;
		r = GetDevPart(NextStress);  r /= p;

		temp4 = mKp + 2 * G - K * D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma1.Zero();
			dAlpha1.Zero();
			dFabric1.Zero();
			// dPStrain1 = dDevStrain + dVolStrain * mI1;
			dPStrain1 = tmp0;
		}
		else {
			NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4;
			if (NextL < 0) {
				if (debugFlag) {
					opserr << "1 NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				// dSigma1 = 2 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1;
				tmp2 = mI1; tmp2 *= (K * dVolStrain);
				dSigma1 = ToContraviant(dDevStrain); dSigma1 *= (2.0 * G);
				dSigma1 += tmp2;
				// dAlpha1 = 2.0*(GetDevPart(NextStress + dSigma1) / GetTrace(NextStress + dSigma1) - GetDevPart(NextStress) / GetTrace(NextStress));
				dAlpha1.Zero();
				dFabric1.Zero();
				dPStrain1.Zero();
			}
			else {
				// dSigma1 = 2.0 * G * mIIcon * dDevStrain + K*dVolStrain*mI1 - Macauley(NextL) * (2.0 * G * n + K * D * mI1);
				tmp0 = n; tmp0 *= (2.0 * G);
				tmp1 = mI1; tmp1 *= (K * D); tmp1 += tmp0; tmp1 *= (-Macauley(NextL));
				tmp2 = mI1; tmp2 *= (K * dVolStrain);
				dSigma1 = ToContraviant(dDevStrain); dSigma1 *= (2.0 * G);
				dSigma1 += tmp2; dSigma1 += tmp1;
				// update fabric
				// if (DoubleDot2_2_Contr(alphaD - NextAlpha, n) < 0.0) {
				alphaD_NextAlpha = alphaD; alphaD_NextAlpha -= NextAlpha;
				if (DoubleDot2_2_Contr(alphaD_NextAlpha, n) < 0.0) {
					// dFabric1 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + NextFabric);
					dFabric1 = n;
					dFabric1 *= m_z_max;
					dFabric1 += NextFabric;
					dFabric1 *= (-1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL) * MacauleyIndex(-D));
				}
				// dPStrain1 = NextL * ToCovariant(R1);
				// dAlpha1 = two3 * NextL * h * b;
				dPStrain1 = ToCovariant(R1);  dPStrain1 *= NextL;
				dAlpha1 = b; dAlpha1 *= (two3 * NextL * h);
			}
		}
		//Calc Delta 2
		// p = 0.5 * GetTrace(NextStress + dSigma1);
		tmp0 = NextStress; tmp0 += dSigma1;   // tmp0 is NextStress + dSigma1 until line 1511
		p = 0.5 * GetTrace(tmp0);
		if (p < 0) {
			if (dT == dT_min) {
				if (debugFlag)
					opserr << "Delta 1: p < 0";
				NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
				NextStress = CurStress;
				NextAlpha = CurAlpha;
				NextFabric = CurFabric;
				return;
			}
			dT = fmax(0.1*dT, dT_min);
			continue;
		}

		// GetStateDependent(NextStress + dSigma1, NextAlpha + dAlpha1, alpha_in, alpha_in_p, NextFabric + dFabric1, mFabric_in, G, mzcum
		// 	, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R2, mKp, alphaD, Cka, h, b, AlphaAlphaBDotN);
		tmp1.Zero();  tmp1 += NextAlpha; tmp1 += dAlpha1;  // tmp1 is NextAlpha + dAlpha1
		tmp2.Zero();  tmp2 += NextFabric; tmp2 += dFabric1;  // tmp2 is NextFabric + dFabric1
		GetStateDependent(tmp0, tmp1, alpha_in, alpha_in_p, tmp2, mFabric_in, G, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R2, mKp, alphaD, Cka, h, b, AlphaAlphaBDotN);
		// r = GetDevPart(NextStress + dSigma1) / p;
		r = GetDevPart(tmp0); r /= p;

		temp4 = mKp + 2 * G - K * D * DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma2.Zero();
			dAlpha2.Zero();
			dFabric2.Zero();
			// dPStrain2 = dDevStrain + dVolStrain * mI1;
			dPStrain2 = dPStrain1;
		}
		else {
			NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4;
			mDGamma = NextL;
			if (NextL < 0)
			{
				if (debugFlag) {
					opserr << "2 NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				// dSigma2 = 2 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1;
				tmp2 = mI1; tmp2 *= (K * dVolStrain);
				dSigma2 = ToContraviant(dDevStrain); dSigma2 *= (2.0 * G);
				dSigma2 += tmp2;
				// dAlpha2 = 2.0*(GetDevPart(NextStress + dSigma2) / GetTrace(NextStress + dSigma2) - GetDevPart(NextStress) / GetTrace(NextStress));
				dAlpha2.Zero();
				dFabric2.Zero();
				dPStrain2.Zero();
			}
			else {
				// dSigma2 = 2.0 * G * mIIcon * dDevStrain + K*dVolStrain*mI1 - Macauley(NextL)*
				// 	(2.0 * G * n + K * D * mI1);
				tmp0 = n; tmp0 *= (2.0 * G);
				tmp1 = mI1; tmp1 *= (K * D); tmp1 += tmp0; tmp1 *= (-Macauley(NextL));
				tmp2 = mI1; tmp2 *= (K * dVolStrain);
				dSigma2 = ToContraviant(dDevStrain); dSigma2 *= (2.0 * G);
				dSigma2 += tmp2; dSigma2 += tmp1;
				// update fabric
				// if (DoubleDot2_2_Contr(alphaD - (NextAlpha + dAlpha1), n) < 0.0) {
				alphaD_NextAlpha = alphaD; alphaD_NextAlpha -= NextAlpha; alphaD_NextAlpha -= dAlpha1;
				if (DoubleDot2_2_Contr(alphaD_NextAlpha, n) < 0.0) {
					// Equation 57
					//dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + NextFabric + dFabric1);
					dFabric2 = n;
					dFabric2 *= m_z_max;
					dFabric2 += NextFabric;
					dFabric2 += dFabric1;
					dFabric2 *= (-1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL) * MacauleyIndex(-D));
				}
				// dPStrain2 = NextL * mIIco * R2;
				// dAlpha2 = two3 * NextL * h * b;
				dPStrain2 = ToCovariant(R2);  dPStrain2 *= NextL;
				dAlpha2 = b; dAlpha2 *= (two3 * NextL * h);
			}
		}

		// nStress = NextStress + 0.5 * (dSigma1 + dSigma2);
		// nFabric = NextFabric + 0.5 * (dFabric1 + dFabric2);
		// // update alpha
		// dAlpha = 0.5 * (dAlpha1 + dAlpha2);
		// nAlpha = NextAlpha + dAlpha;
		nStress = dSigma1;
		nStress += dSigma2;
		nStress *= 0.5;
		nStress += NextStress;
		nFabric = dFabric1;
		nFabric += dFabric2;
		nFabric *= 0.5;
		nFabric += NextFabric;
		nAlpha = dAlpha1;
		nAlpha += dAlpha2;
		nAlpha *= 0.5;
		nAlpha += NextAlpha;

		p = 0.5 * GetTrace(nStress);
		if (p < 0)
		{
			if (dT == dT_min) {
				opserr << "Delta 2: p < 0";
				NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
				NextStress = CurStress;
				NextAlpha = CurAlpha;
				NextFabric = CurFabric;
				return;
			}
			dT = fmax(0.1 * dT, dT_min);
			continue;
		}

		stressNorm = GetNorm_Contr(NextStress);
		tmp0 = dSigma2; tmp0 -= dSigma1;
		if (stressNorm < 0.5)
			// curStepError = GetNorm_Contr(dSigma2 - dSigma1);
			curStepError = GetNorm_Contr(tmp0);
		else
			// curStepError = GetNorm_Contr(dSigma2 - dSigma1) / (2 * stressNorm);
			curStepError = GetNorm_Contr(tmp0) / (2 * stressNorm);

		if (curStepError > TolE) {
			q = fmax(0.8 * sqrt(TolE / curStepError), 0.1);
			if (dT == dT_min) {
				// opserr << "reached dT_min\n";
				// NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
				tmp0 = dPStrain1; tmp0 += dPStrain2; tmp0 *= 0.5;
				NextElasticStrain -= tmp0;
				NextStress = nStress;
				NextAlpha = nAlpha;
				Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
				//Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
				T += dT;
			}
			dT = fmax(q * dT, dT_min);
		}
		else {
			// NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
			tmp0 = dPStrain1; tmp0 += dPStrain2; tmp0 *= 0.5;
			NextElasticStrain -= tmp0;
			NextStress = nStress;
			NextAlpha = nAlpha;
			NextFabric = nFabric;
			Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
			//Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);

			T += dT;
			q = fmax(0.8 * sqrt(TolE / curStepError), 0.5);
			dT = fmax(q * dT, dT_min);
			dT = fmin(dT, 1 - T);
		}
	}
	return;

}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Runge-Kutta Integrator
/*************************************************************/
void PM4Silt::RungeKutta4(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& alpha_in_p, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	double dVolStrain, p, Cka, D, K_p, temp4, h, AlphaAlphaBDotN;
	Vector n(3), R1(3), R2(3), R3(3), R4(3), alphaD(3), dDevStrain(3), r(3), b(3);
	Vector nStress(3), nAlpha(3), nFabric(3);
	Vector dSigma1(3), dSigma2(3), dSigma3(3), dSigma4(3), dSigma(3), dAlpha1(3), dAlpha2(3),
		dAlpha3(3), dAlpha4(3), dAlpha(3), dFabric1(3), dFabric2(3), dFabric3(3), dFabric4(3),
		dFabric(3), dPStrain1(3), dPStrain2(3), dPStrain3(3), dPStrain4(3), dPStrain(3);
	double T = 0.0, dT = 0.5, dT_min = 1.0e-4, TolE = 1.0e-5;

	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextStress = CurStress;
	NextAlpha = CurAlpha;
	NextFabric = CurFabric;

	p = 0.5 * GetTrace(CurStress);
	if (p < m_Pmin / 5.0)
	{
		if (debugFlag)
			opserr << "Tag = " << this->getTag() << " : p < pmin / 5, should not happen" << endln;
		NextStress = GetDevPart(NextStress) + m_Pmin / 5.0 * mI1;
	}
	while (T < 1.0)
	{
		NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + T*(NextStrain - CurStrain));
		dVolStrain = dT * GetTrace(NextStrain - CurStrain);
		dDevStrain = dT * (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
		p = 0.5 * GetTrace(NextStress);
		// Calc Delta 1
		GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, NextFabric, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R1, K_p, alphaD, Cka, h, b, AlphaAlphaBDotN);

		r = GetDevPart(NextStress) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (temp4 < 0.0) {
			mKp = -0.5 * (2 * G - K* D *DoubleDot2_2_Contr(n, r));
			temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
			h = 1.5 * mKp / (p * AlphaAlphaBDotN);
		}
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma1.Zero();
			dAlpha1.Zero();
			dFabric1.Zero();
			dPStrain1 = dDevStrain + dVolStrain * mI1;
		}
		else {
			NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
			if (NextL < 0) {
				if (debugFlag) {
					opserr << "1 NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				dSigma1 = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
				// dAlpha1 = GetDevPart(NextStress + dSigma1) / (0.5 * GetTrace(NextStress + dSigma1))
				//	- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress));
				dAlpha1.Zero();
				dFabric1.Zero();
				dPStrain1.Zero();
			}
			else {
				dSigma1 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * mG * n + mK * D * mI1);
				// update fabric
				if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
					dFabric1 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
				}
				dPStrain1 = NextL * mIIco * R1;
				dAlpha1 = two3 * NextL * h * b;
			}
		}
		//Calc Delta 2
		p = 0.5 * GetTrace(NextStress + 0.5 * dSigma1);

		GetStateDependent(NextStress + 0.5 * dSigma1, CurAlpha + 0.5 * dAlpha1, alpha_in, alpha_in_p, NextFabric + 0.5 * dFabric1, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R2, K_p, alphaD, Cka, h, b, AlphaAlphaBDotN);
		r = GetDevPart(NextStress + 0.5 * dSigma1) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma2.Zero();
			dAlpha2.Zero();
			dFabric2.Zero();
			dPStrain2 = dDevStrain + dVolStrain * mI1;
		}
		else {
			NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
			if (NextL < 0)
			{
				if (debugFlag) {
					opserr << "2nd NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				dSigma2 = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
				// dAlpha2 = GetDevPart(NextStress + dSigma2) / (0.5 * GetTrace(NextStress + dSigma2))
				// 	- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress));
				dAlpha2.Zero();
				dFabric2.Zero();
				dPStrain2.Zero();
			}
			else {
				dSigma2 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * mG * n + mK * D * mI1);
				// update fabric
				if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
					dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric + 0.5 * dFabric1);
				}
				dPStrain2 = NextL * mIIco * R2;
				dAlpha2 = two3 * NextL * h * b;
			}
		}
		//Calc Delta 3
		p = 0.5 * GetTrace(NextStress + 0.5 * dSigma2);

		GetStateDependent(NextStress + 0.5 * dSigma2, CurAlpha + 0.5 * dAlpha2, alpha_in, alpha_in_p, NextFabric + 0.5 * dFabric2, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R3, K_p, alphaD, Cka, h, b, AlphaAlphaBDotN);
		r = GetDevPart(NextStress + 0.5 * dSigma2) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma3.Zero();
			dAlpha3.Zero();
			dFabric3.Zero();
			dPStrain3 = dDevStrain + dVolStrain * mI1;
		}
		else {
			NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
			if (NextL < 0)
			{
				if (debugFlag) {
					opserr << "3rd NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				dSigma3 = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
				// dAlpha3 = GetDevPart(NextStress + dSigma3) / (0.5 * GetTrace(NextStress + dSigma3))
				// 	- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress));
				dAlpha3.Zero();
				dFabric3.Zero();
				dPStrain3.Zero();
			}
			else {
				dSigma3 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * mG * n + mK * D * mI1);
				// update fabric
				if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
					dFabric3 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric + 0.5 * dFabric2);
				}
				dPStrain3 = NextL * mIIco * R3;
				dAlpha3 = two3 * NextL * h * b;
			}
		}
		//Calc Delta 4
		p = 0.5 * GetTrace(NextStress + dSigma3);

		GetStateDependent(NextStress + dSigma3, CurAlpha + dAlpha3, alpha_in, alpha_in_p, NextFabric + dFabric3, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R4, K_p, alphaD, Cka, h, b, AlphaAlphaBDotN);
		r = GetDevPart(NextStress + dSigma3) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma4.Zero();
			dAlpha4.Zero();
			dFabric4.Zero();
			dPStrain4 = dDevStrain + dVolStrain * mI1;
		}
		else {
			NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
			if (NextL < 0)
			{
				if (debugFlag) {
					opserr << "4th NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				dSigma4 = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
				// dAlpha4 = GetDevPart(NextStress + dSigma4) / (0.5 * GetTrace(NextStress + dSigma4))
				// 	- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress));
				dAlpha4.Zero();
				dFabric4.Zero();
				dPStrain4.Zero();
			}
			else {
				dSigma4 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * mG * n + mK * D * mI1);
				// update fabric
				if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
					dFabric4 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric + dFabric3);
				}
				dPStrain4 = NextL * mIIco * R4;
				dAlpha4 = two3 * NextL * h * b;
			}
		}

		// RK4
		dSigma = (dSigma1 + dSigma4 + 2.0 * (dSigma2 + dSigma3)) / 6.0;
		dAlpha = (dAlpha1 + dAlpha4 + 2.0 * (dAlpha2 + dAlpha3)) / 6.0;
		dFabric = (dFabric1 + dFabric4 + 2.0 * (dFabric2 + dFabric3)) / 6.0;
		dPStrain = (dPStrain1 + dPStrain4 + 2.0 * (dPStrain2 + dPStrain3)) / 6.0;

		nStress = NextStress + dSigma;
		nAlpha = NextAlpha + dAlpha;
		nFabric = NextFabric + dFabric;

		// can add error control here
		NextElasticStrain -= dPStrain;
		NextStress = nStress;
		NextAlpha = nAlpha;
		NextFabric = nFabric;
		Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
		// Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, (R1 + R4 + 2.0 * (R2 + R3)) / 6, n, r);
		T += dT;
		//q = fmax(0.8 * sqrt(TolE / curStepError), 0.5);
		//dT = fmax(q * dT, dT_min);
		//dT = fmin(dT, 1 - T);
	}
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
//            Pegasus Iterations                             //
/*************************************************************/
double
PM4Silt::IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha,
	double a0 = 0.0, double a1 = 1.0)
{
	double a = a0;
	double f, f0, f1;
	Vector dSigma(3), dSigma0(3), dSigma1(3), strainInc(3), tmp(3);

	// strainInc = NextStrain - CurStrain;
	strainInc += NextStrain;
	strainInc -= CurStrain;

	if (a0 < 0.0 || a1 > 1.0) {
		opserr << "a0 = " << a0 << "a1 = " << a1 << endln;
	}
	//GetElasticModuli(CurStress, K, G, mzcum);
	dSigma0 = a0 * DoubleDot4_2(mCe, strainInc);
	// f0 = GetF(CurStress + dSigma0, CurAlpha);
	tmp.Zero(); tmp += CurStress; tmp += dSigma0;
	f0 = GetF(tmp, CurAlpha);

	dSigma1 = a1 * DoubleDot4_2(mCe, strainInc);
	// f1 = GetF(CurStress + dSigma1, CurAlpha);
	tmp.Zero(); tmp += CurStress; tmp += dSigma1;
	f1 = GetF(tmp, CurAlpha);

	for (int i = 1; i <= 10; i++)
	{
		a = a1 - f1 * (a1 - a0) / (f1 - f0);
		dSigma = a * DoubleDot4_2(mCe, strainInc);
		// f = GetF(CurStress + dSigma, CurAlpha);
		tmp.Zero(); tmp += CurStress; tmp += dSigma;
		f = GetF(tmp, CurAlpha);
		if (fabs(f) < mTolF)
		{
			// if (debugFlag) opserr << "Found alpha in " << i << " steps" << ", alpha = " << a << endln;
			break;
		}
		if (f * f0 < 0)
		{
			a1 = a;
			f1 = f;
		}
		else {
			f1 = f1 * f0 / (f0 + f);
			a0 = a;
			f0 = f;
		}

		if (i == 10)
		{
			if (debugFlag) opserr << "Didn't find alpha!" << endln;
			a = 0;
			break;
		}
	}
	if (a > 1 - small) a = 1.0;
	if (a < small) a = 0.0;
	if (a != a) {
		if (debugFlag)
			opserr << "a is nan" << endln;
		a = 0.0;
	}
	return a;
}
/*************************************************************/
//      Pegasus Iterations  (ElastoPlastic Unloading)        //
/*************************************************************/
double
PM4Silt::IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha)
{
	double a = 0.0, a0 = 0.0, a1 = 1.0, da;
	double f, f0, f1, fs;
	int nSub = 20;
	Vector dSigma(3), dSigma0(3), dSigma1(3), strainInc(3), tmp(3);
	bool flag = false;

	// strainInc = NextStrain - CurStrain;
	strainInc += NextStrain; strainInc -= CurStrain;

	f0 = GetF(CurStress, CurAlpha);
	fs = f0;

	// GetElasticModuli(CurStress, K, G, mzcum);
	dSigma = DoubleDot4_2(mCe, strainInc);

	for (int i = 1; i < 10; i++)
	{
		da = (a1 - a0) / nSub;
		for (int k = 1; k < nSub; k++) {
			a = a0 + da;
			// f = GetF(CurStress + a * dSigma, CurAlpha);
			tmp = dSigma; tmp *= a; tmp += CurStress;
			f = GetF(tmp, CurAlpha);
			if (f > mTolF)
			{
				a1 = a;
				if (f0 < -mTolF) {
					f1 = f;
					flag = true;
					break;
				}
				else {
					a0 = 0.0;
					f0 = fs;
					break;
				}
			}
			else {
				a0 = a;
				f0 = f;
			}
			if (i == 10) {
				if (debugFlag)
					opserr << "Didn't find alpha! - Unloading" << ", a0 = " << a0 << ", a1 = " << a1 << endln;
				return 0.0;
			}
			if (flag) break;
		}
	}
	if (debugFlag)
		opserr << "Found alpha - Unloading" << ", a0 = " << a0 << ", a1 = " << a1 << endln;
	return IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1);
}
/*************************************************************/
//            Stress Correction                              //
/*************************************************************/
void
PM4Silt::Stress_Correction(Vector& NextStress, Vector& NextAlpha, const Vector& alpha_in, const Vector& alpha_in_p,
	const Vector& CurFabric, double& NextVoidRatio)
{
	Vector dSigmaP(3), dfrOverdSigma(3), dfrOverdAlpha(3), n(3), R(3), alphaD(3), b(3), aBar(3), r(3);
	Vector nAlpha(3), nStress(3), dSigma(3), tmp0(3), tmp1(3);
	double lambda, D, K_p, Cka, h, p, fr, AlphaAlphaBDotN;
	Matrix aC(3, 3);
	// Vector CurStress = NextStress;

	int maxIter = 25;
	p = 0.5 * GetTrace(NextStress);
	if (p < m_Pmin / 5.0) {
		fr = GetF(NextStress, NextAlpha);
		if (fr < mTolF) {
			// stress state inside yield surface
			NextStress += (m_Pmin / 5.0 - p)  * mI1;
		}
		else {
			// stress state ouside yield surface
			NextStress = m_Pmin / 5.0 * mI1;
			NextStress(2) = 0.8 * m_Mc * m_Pmin / 5.0;
			NextAlpha.Zero();
			NextAlpha(2) = 0.8 * m_Mc;
			return;
		}
	}
	else {
		fr = GetF(NextStress, NextAlpha);
		if (fr < mTolF) {
			// stress state inside yield surface
			return;
		}
		else {
			nStress = NextStress;
			nAlpha = NextAlpha;
			for (int i = 1; i <= maxIter; i++) {
				// r = GetDevPart(nStress) / p;
				r = GetDevPart(nStress); r /= p;
				GetStateDependent(nStress, nAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum
					, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R, K_p, alphaD, Cka, h, b, AlphaAlphaBDotN);
				aC = GetStiffness(mK, mG);
				dSigmaP = DoubleDot4_2(aC, mDGamma * ToCovariant(R));
				aBar = two3 * h * b;
				// dfrOverdSigma = n - 0.5 * DoubleDot2_2_Contr(n, r) * mI1;
				dfrOverdSigma.Zero(); dfrOverdSigma += mI1;
				dfrOverdSigma *= (-0.5 * DoubleDot2_2_Contr(n, r));	dfrOverdSigma += n;
				// dfrOverdAlpha = -p * n;
				dfrOverdAlpha = n; dfrOverdAlpha *= (-p);
				lambda = fr / (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP) - DoubleDot2_2_Contr(dfrOverdAlpha, aBar));
				// if (fabs(GetF(nStress - lambda * dSigmaP, nAlpha + lambda * aBar)) < fabs(fr))
				tmp0 = dSigmaP; tmp0 *= (-lambda); tmp0 += nStress;
				tmp1 = aBar; tmp1 *= lambda; tmp1 += nAlpha;
				if (fabs(GetF(tmp0, tmp1)) < fabs(fr))
				{
					// nStress -= lambda * dSigmaP;
					tmp0 = dSigmaP; tmp0 *= lambda; nStress -= tmp0;
					// nAlpha += lambda * aBar;
					tmp0 = aBar; tmp0 *= lambda; nAlpha += tmp0;
				}
				else {
					lambda = fr / DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma);
					// nStress -= lambda * dfrOverdSigma;
					tmp0 = dfrOverdSigma; tmp0 *= lambda; nStress -= tmp0;
				}
				fr = GetF(nStress, nAlpha);
				if (fabs(fr) < mTolF) {
					NextStress = nStress;
					NextAlpha = nAlpha;
					return;
				}

				p = fmax(0.5 * GetTrace(nStress), m_Pmin);
			}
			// if (fabs(fr) < fabs(GetF(NextStress, NextAlpha))) {
			// 	NextStress = nStress;
			// 	NextAlpha = nAlpha;
			// }
			if (debugFlag) {
				opserr << "Still outside with f =  " << fr << endln;
				opserr << "NextStress = " << NextStress;
				opserr << "nStress = " << nStress;
				opserr << "NextAlpha = " << NextAlpha;
			}

			Vector dSigma = NextStress - mSigma;
			double alpha_up = 1.0;
			double alpha_mid = 0.5;
			double alpha_down = 0.0;
			double fr_old = GetF(mSigma + alpha_mid * dSigma, NextAlpha);
			for (int jj = 0; jj < maxIter; jj++) {
				if (fr_old < 0.0) {
					alpha_down = alpha_mid;
					alpha_mid = 0.5 * (alpha_up + alpha_mid);
				}
				else {
					alpha_up = alpha_mid;
					alpha_mid = 0.5 * (alpha_down + alpha_mid);
				}

				// fr_old = GetF(mSigma + alpha_mid * dSigma, NextAlpha);
				tmp0 = dSigma; tmp0 *= alpha_mid; tmp0 += mSigma;
				fr_old = GetF(tmp0, NextAlpha);
				if (fabs(fr_old) < mTolF) {
					// NextStress = mSigma + alpha_mid * dSigma;
					NextStress = dSigma; NextStress *= alpha_mid; NextStress += mSigma;
					break;
				}
			}
		}
	}
}

/************************************************************/
/************************************************************/
void
PM4Silt::Stress_Correction(Vector& NextStress, Vector& NextAlpha, const Vector& dAlpha,
	const double m, const Vector& R, const Vector& n, const Vector& r)
{
	Vector dfrOverdSigma(3);
	double lambda;
	int maxIter = 50;
	double f = GetF(NextStress, NextAlpha);
	if (f < mTolF)
		return;
	else {
		// Method C from Potts and Gens 1983. 
		for (int i = 1; i <= maxIter; i++) {
			// GetStateDependent(CurStress, CurAlpha, alpha_in, CurFabric, mFabric_in, mG, mzcum
			// 	, mzpeak, mpzp, mMcur, ksi, CurDr, n, D, R, K_p, alphaD, Cka, h, b);
			dfrOverdSigma = n - 0.5 * DoubleDot2_2_Contr(n, r)*mI1;
			lambda = f / DoubleDot2_2_Contr(dfrOverdSigma, R);
			NextStress = NextStress - R*lambda;
			NextAlpha = NextAlpha - dAlpha*lambda;
			f = GetF(NextStress, NextAlpha);
			if (fabs(f) < mTolF) {
				break;
			}
			if (i == maxIter)
				if (debugFlag)
					opserr << "Still outside with f =  " << f << endln;
		}
	}
}
/*************************************************************/
/*************************************************************/
//            MATERIAL SPECIFIC METHODS                      //
/*************************************************************/
/*************************************************************/
// Macauley() -------------------------------------------------
double PM4Silt::Macauley(double x)
{
	// Macauley bracket
	return (x > 0 ? x : 0.0);
}
/*************************************************************/
// MacauleyIndex() --------------------------------------------
double PM4Silt::MacauleyIndex(double x)
{
	// Macauley index
	return (x > 0 ? 1.0 : 0.0);
}
/*************************************************************/
// GetF() -----------------------------------------------------
double
PM4Silt::GetF(const Vector& nStress, const Vector& nAlpha)
{
	// PM4Silt's yield function
	Vector s(3); s = GetDevPart(nStress);
	double p = 0.5 * GetTrace(nStress);
	s -= p * nAlpha;
	double f = GetNorm_Contr(s) - root12 * m_m * p;
	return f;
}
/*************************************************************/
// GetKSI() ---------------------------------------------------
double
PM4Silt::GetKsi(const double& e, const double& p)
{
	double pn = p;
	pn = (pn <= m_Pmin) ? (m_Pmin) : pn;
	double ksi = e - me0 + m_lambda * log(101.3 * pn / (m_P_atm * m_Fsu));
	return ksi;
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
PM4Silt::GetElasticModuli(const Vector& sigma, double &K, double &G, double &Mcur, const double& zcum)
// Calculates G, K, including effects of fabric and current stress ratio
{
	int msr = 4;
	double Csr0 = 0.5;
	double pn = 0.5 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;
	double qn = 2 * sqrt(pow((0.5*(sigma(0) - sigma(1))), 2) + pow(sigma(2), 2));
	//double q = sqrt(2.0 * DoubleDot2_2_Contr(GetDevPart(sigma), GetDevPart(sigma)));
	// Mcur = 2 * sqrt(2) * GetNorm_Contr(GetDevPart(sigma)) / GetTrace(sigma);
	Mcur = qn / pn;
	double Csr = 1 - Csr0 * fmin(1.0, pow((Mcur / mMb), msr));
	double temp = zcum / m_z_max;
	if (me2p == 0)
		G = m_G0 * m_P_atm;
	else {
		G = m_G0 * m_P_atm * pow(pn / m_P_atm, m_nG) * Csr * (1 + temp) / (1 + temp * m_Cgd);
		if (m_PostShake) {
			// reduce elastic shear modulus for post shaking consolidation
			double G_c_min = 8 * pn / m_lambda * (1.0 / (1 + (m_CG_consol - 1) * (mzcum / (mzcum + m_z_max))));
			double F_consol = 1 - (1 - G_c_min / G) * pow(Macauley(1 - Mcur / mMd), 0.25);
			G = G * F_consol;
		}
	}
	m_nu = (m_nu >= 0.5) ? 0.4999 : m_nu;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}
void
PM4Silt::GetElasticModuli(const Vector& sigma, double &K, double &G)
// Calculates G, K
{
	double pn = 0.5 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;

	if (me2p == 0)
		G = m_G0 * m_P_atm;
	else
		G = m_G0 * m_P_atm * sqrt(pn / m_P_atm);
	m_nu = (0.5 - m_nu < small) ? 0.4999 : m_nu;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}
/*************************************************************/
// GetStiffness() ---------------------------------------------
Matrix
PM4Silt::GetStiffness(const double& K, const double& G)
// returns the stiffness matrix in its contravarinat-contravariant form
{
	Matrix C(3, 3);
	double a = K + 4.0*one3 * G;
	double b = K - 2.0*one3 * G;
	C(0, 0) = C(1, 1) = a;
	C(2, 2) = G;
	C(0, 1) = C(1, 0) = b;
	return C;
}
/*************************************************************/
// GetCompliance() ---------------------------------------------
Matrix
PM4Silt::GetCompliance(const double& K, const double& G)
// returns the compliance matrix in its covariant-covariant form
{
	Matrix D(3, 3);
	double a = (K + 4.0 / 3.0 * G) / (4.0 * G * K + 4.0 / 3.0 * pow(G, 2));
	double b = (K - 2.0 / 3.0 * G) / (4.0 * G * K + 4.0 / 3.0 * pow(G, 2));
	double c = 1 / G;
	D(0, 0) = D(1, 1) = a;
	D(2, 2) = c;
	D(0, 1) = D(1, 0) = b;
	return D;
}
/*************************************************************/
// GetElastoPlasticTangent()---------------------------------------
Matrix
PM4Silt::GetElastoPlasticTangent(const Vector& NextStress, const Matrix& aCe, const Vector& R,
	const Vector& n, const double K_p)
{
	double p = 0.5 * GetTrace(NextStress);
	if (p < m_Pmin) p = m_Pmin;
	Vector r = GetDevPart(NextStress) / p;
	Matrix aCep(3, 3);
	aCep.Zero();
	Vector temp1 = DoubleDot4_2(aCe, R);
	Vector temp2 = DoubleDot2_4(n - 1 / 2 * DoubleDot2_2_Contr(n, r)*mI1, aCe*mIIco);
	double temp3 = DoubleDot2_2_Contr(temp2, R) + K_p;
	if (temp3 < small) {
		aCep = aCe;
	}
	else {
		aCep = aCe - 1 / temp3 * Dyadic2_2(temp1, temp2);
	}

	return aCep;
}
/*************************************************************/
// GetNormalToYield() ----------------------------------------
Vector
PM4Silt::GetNormalToYield(const Vector &stress, const Vector &alpha)
{
	// Vector devStress(3); devStress = GetDevPart(stress);
	// double p = 0.5 * GetTrace(stress);
	// Vector n(3);
	// if (fabs(p) < small) {
	// 	n.Zero();
	// }
	// else {
	// 	n = devStress - p * alpha;
	// 	double normN = GetNorm_Contr(n);
	// 	normN = (normN < small) ? 1.0 : normN;
	// 	n = n / normN;
	// }
	Vector n(3);
	double p = 0.5 * GetTrace(stress);
	if (fabs(p) < small) {
		// n.Zero();
		// change loading direction to simple shear when p is small
		n(2) = root12;
	}
	else {
		n = alpha; n *= (-p);
		n += GetDevPart(stress);
		double normN = GetNorm_Contr(n);
		normN = (normN < small) ? 1.0 : normN;
		n /= normN;
	}
	return n;
}
/*************************************************************/
// Check() ---------------------------------------------------
int
PM4Silt::Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha)
// Check if the solution of implicit integration makes sense
{
	return 0;
}
/*************************************************************/
// GetStateDependent() ----------------------------------------
void
PM4Silt::GetStateDependent(const Vector &stress, const Vector &alpha, const Vector &alpha_in, const Vector &alpha_in_p
	, const Vector &fabric, const Vector &fabric_in, const double &G, const double &zcum, const double &zpeak
	, const double &pzp, const double &Mcur, const double &CurVoidRatio, Vector &n, double &D, Vector &R, double &K_p
	, Vector &alphaD, double &Cka, double &h, Vector &b, double &AlphaAlphaBDotN)
{
	Vector alphaD_alpha(3), alphaDr_alpha(3), alpha_mAlpha_in(3), alpha_mAlpha_in_true(3), alpha_mAlpha_p(3), minusFabric(3);
	double Czpk1, Czpk2, Cpzp2, Cg1, Ckp, AlphaAlphaInDotN, AlphaAlphaInTrueDotN, Czin1, Crot1, Mdr;
	double p = 0.5 * GetTrace(stress);
	if (p <= m_Pmin) p = m_Pmin;
	double ksi = GetKsi(CurVoidRatio, p);
	n = GetNormalToYield(stress, alpha);
	mMd = fmin(1.4142136, m_Mc * exp(m_nd * ksi / m_lambda));
	if (ksi < 0) {
		// dense of critical
		mMb = m_Mc * pow((1 + mC_MB) / (p / mpcs + mC_MB), m_nbdry);
	}
	else {
		//loose of critical
		mMb = m_Mc * exp(-1.0 * m_nbwet * ksi / m_lambda);
	}
	//Vector alphaB = root12 * (mMb - m_m) * n;
	Vector alphaB(n);
	alphaB *= (root12 * (mMb - m_m));
	//alphaD = root12 * (mMd - m_m) * n;
	alphaD = n; alphaD *= (root12 * (mMd - m_m));

	Czpk1 = zpeak / (zcum + m_z_max / 5.0);
	Czpk2 = zpeak / (zcum + m_z_max / 100.0);
	Cpzp2 = Macauley((pzp - p)) / (Macauley((pzp - p)) + m_Pmin);
	Cg1 = m_h0 / 200.0;
	Ckp = 2.0;

	//b = alphaB - alpha;
	b = alphaB; b -= alpha;
	AlphaAlphaBDotN = DoubleDot2_2_Contr(b, n);
	// double AlphaAlphaInDotN = Macauley(DoubleDot2_2_Contr(alpha - mAlpha_in, n));
	alpha_mAlpha_in = alpha; alpha_mAlpha_in -= mAlpha_in;
	AlphaAlphaInDotN = Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in, n));
	// double AlphaAlphaInTrueDotN = Macauley(DoubleDot2_2_Contr(alpha - mAlpha_in_true, n));
	alpha_mAlpha_in_true = alpha; alpha_mAlpha_in_true -= mAlpha_in_true;
	AlphaAlphaInTrueDotN = Macauley(DoubleDot2_2_Contr(alpha_mAlpha_in_true, n));
	Cka = 1.0 + m_Ckaf / (1.0 + pow(2.5*AlphaAlphaInTrueDotN, 2)) * Cpzp2 * Czpk1;
	// updataed K_p formulation following PM4Silt V1. mAlpha_in is the apparent back-stress ratio.
	// if (DoubleDot2_2_Contr(alpha - alpha_in_p, n) <= 0) {
	alpha_mAlpha_p = alpha; alpha_mAlpha_p -= alpha_in_p;
	if (fabs(AlphaAlphaInDotN) < small) {
		// adding this condition to avoid division by zero error
		h = 1.0e10;
	}
	else if (DoubleDot2_2_Contr(alpha_mAlpha_p, n) <= 0) {
		h = 1.5 * G * m_h0 / p / (exp(AlphaAlphaInDotN) - 1 + Cg1) / sqrt(fabs(AlphaAlphaBDotN)) *
			Cka / (1 + Ckp * zpeak / m_z_max * Macauley(AlphaAlphaBDotN) * sqrt(1 - Czpk2));
		h = h * (AlphaAlphaInDotN + Cg1) / (AlphaAlphaInTrueDotN + Cg1);
	}
	else {
		h = 1.5 * G * m_h0 / p / (exp(AlphaAlphaInDotN) - 1 + Cg1) / sqrt(fabs(AlphaAlphaBDotN)) *
			Cka / (1 + Ckp * zpeak / m_z_max * Macauley(AlphaAlphaBDotN) * sqrt(1 - Czpk2));
	}

	K_p = two3 * h * p * AlphaAlphaBDotN;
	Czin1 = Macauley(1.0 - exp(-2.0 * fabs((DoubleDot2_2_Contr(fabric_in, n) - DoubleDot2_2_Contr(fabric, n)) / m_z_max)));
	// rotated dilatancy surface
	minusFabric = fabric; minusFabric *= (-1.0);
	double temp = Macauley(DoubleDot2_2_Contr(minusFabric, n)) * root12;
	Crot1 = fmax((1.0 + 2 * temp / m_z_max * (1 - Czin1)), 1.0);
	Mdr = mMd / Crot1;
	// Vector alphaDr = root12 * (Mdr - m_m) * n;
	alphaDr_alpha = n; alphaDr_alpha *= root12 * (Mdr - m_m); alphaDr_alpha -= alpha;
	alphaD_alpha = alphaD; alphaD_alpha -= alpha;
	// dilation
	if (DoubleDot2_2_Contr(alphaDr_alpha, n) <= 0) {
		double Cpzp, Czin2, Ad, Drot;
		Cpzp = 1.0 / (1.0 + pow((2.5* p / mpzp), 5.0));
		Czin2 = (1 + Czin1*(zcum - zpeak) / (3.0 * m_z_max)) / (1 + 3.0 * Czin1*(zcum - zpeak) / (3.0 * m_z_max));
		Ad = m_Ado * Czin2 / ((pow(zcum, 2) / m_z_max) * pow(1.0 - temp / zpeak, 3) * pow(m_ce, 2.0) * Cpzp * Czin1 + 1.0);
		D = Ad * DoubleDot2_2_Contr(alphaD_alpha, n);
		Drot = Ad * temp / m_z_max * DoubleDot2_2_Contr(alphaDr_alpha, n) / 3.0;
		if (D > Drot) {
			D = D + (Drot - D)*Macauley(mMb - Mcur) / (Macauley(mMb - Mcur) + 0.01);
		}
		if (m_Pmin <= p && p <= 2 * m_Pmin) {
			D = fmin(D, -3.5 * m_Ado * Macauley(mMb - mMd) * (2 * m_Pmin - p) / m_Pmin);
		}
	}
	else {
		//contraction
		double hp, Crot2, Cdz, Cwet, Adc, Cin, C_pmin;
		// bound K_p to non - negative, following flac practice
		K_p = fmax(0.0, K_p);
		hp = m_hpo * exp(-0.7 + 0.2 * pow(Macauley(3 - ksi / m_lambda), 2.0));
		Crot2 = 1 - Czpk2;
		Cdz = fmax((1 - Crot2 * sqrt(2.0) * zpeak / m_z_max)*(m_z_max / (m_z_max + Crot2*zcum)), 1 / (1 + m_z_max / 2.0));
		// double Cdz = (1 - Crot2 * sqrt(2.0) * zpeak / m_z_max)*(m_z_max / (m_z_max + Crot2*zcum));
		Cwet = fmin(1.0, (1.0 / (1 + pow(0.02 / AlphaAlphaBDotN, 4.0)) + 1.0 / (1 + pow(ksi / m_lambda / 0.1, 2.0))));
		Adc = m_Ado * (1 + Macauley(DoubleDot2_2_Contr(fabric, n))) / (hp * Cdz * Cwet);
		Cin = Macauley(DoubleDot2_2_Contr(fabric, n)) * sqrt(2.0) / m_z_max;
		D = fmin(Adc * pow((DoubleDot2_2_Contr(alpha_mAlpha_in, n) + Cin), 2.0), m_Ado) *
			DoubleDot2_2_Contr(alphaD_alpha, n) / (DoubleDot2_2_Contr(alphaD_alpha, n) + 0.10);
		// Apply a factor to D so it doesn't go very big when p is small
		if (p < m_Pmin * 2.0)
		{
			C_pmin = 0.0;
		}
		else if (p >= m_Pmin * 8.0)
		{
			C_pmin = 1.0;
		}
		else {
			C_pmin = (p - 2.0 * m_Pmin) / (6.0 * m_Pmin);
		}
		D *= C_pmin;
	}
	R = n + one3 * D * mI1;
	mTracker(1) = D;
}


/*************************************************************/
/*************************************************************/
//            SYMMETRIC TENSOR OPERATIONS                    //
/*************************************************************/
/*************************************************************/
// In all the functions below, by contravariant tensors, we mean stress-like tensors
// and by covariant tensors we mean strain-like tensors

//  GetTrace() ---------------------------------------------
double
PM4Silt::GetTrace(const Vector& v)
// computes the trace of the input argument
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Silt::GetTrace requires vector of size(3)!" << endln;

	return (v(0) + v(1));
}
/*************************************************************/
//  GetDevPart() ---------------------------------------------
Vector
PM4Silt::GetDevPart(const Vector& aV)
// computes the deviatoric part of the input tensor
{
	if (aV.Size() != 3)
		opserr << "\n ERROR! PM4Silt::GetDevPart requires vector of size(3)!" << endln;

	Vector result(3);
	double p = GetTrace(aV);
	result = aV;
	result(0) -= 0.5 * p;
	result(1) -= 0.5 * p;

	return result;
}
/*************************************************************/
// DoubleDot2_2_Contr() ---------------------------------------
double
PM4Silt::DoubleDot2_2_Contr(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "contravariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot2_2_Contr requires vector of size(3)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (i > 1) * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Cov() ---------------------------------------
double
PM4Silt::DoubleDot2_2_Cov(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "covariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot2_2_Cov requires vector of size(3)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) - (i > 1) * 0.5 * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Mixed() ---------------------------------------
double
PM4Silt::DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot2_2_Mixed requires vector of size(3)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// GetNorm_Contr() ---------------------------------------------
double
PM4Silt::GetNorm_Contr(const Vector& v)
// computes contravariant (stress-like) norm of input 6x1 tensor
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Silt::GetNorm_Contr requires vector of size(3)!" << endln;

	double result = 0.0;
	result = sqrt(DoubleDot2_2_Contr(v, v));

	return result;
}
/*************************************************************/
// GetNorm_Cov() ---------------------------------------------
double
PM4Silt::GetNorm_Cov(const Vector& v)
// computes covariant (strain-like) norm of input 6x1 tensor
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Silt::GetNorm_Cov requires vector of size(3)!" << endln;

	double result = 0.0;
	result = sqrt(DoubleDot2_2_Cov(v, v));

	return result;
}
/*************************************************************/
// Dyadic2_2() ---------------------------------------------
Matrix
PM4Silt::Dyadic2_2(const Vector& v1, const Vector& v2)
// computes dyadic product for two vector-storage arguments
// the coordinate form of the result depends on the coordinate form of inputs
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Silt::Dyadic2_2 requires vector of size(3)!" << endln;

	Matrix result(3, 3);

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++)
			result(i, j) = v1(i) * v2(j);
	}

	return result;
}
/*************************************************************/
// DoubleDot4_2() ---------------------------------------------
Vector
PM4Silt::DoubleDot4_2(const Matrix& m1, const Vector& v1)
// computes doubledot product for matrix-vector arguments
// caution: second coordinate of the matrix should be in opposite variant form of vector
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Silt::DoubleDot4_2 requires vector of size(3)!" << endln;
	if ((m1.noCols() != 3) || (m1.noRows() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot4_2 requires 3-by-3 matrix " << endln;

	return m1*v1;
}
/*************************************************************/
// DoubleDot2_4() ---------------------------------------------
Vector
PM4Silt::DoubleDot2_4(const Vector& v1, const Matrix& m1)
// computes doubledot product for matrix-vector arguments
// caution: first coordinate of the matrix should be in opposite 
// variant form of vector
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Silt::DoubleDot2_4 requires vector of size(3)!" << endln;
	if ((m1.noCols() != 3) || (m1.noRows() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot2_4 requires 3-by-3 matrix " << endln;

	return  m1^v1;
}
/*************************************************************/
// DoubleDot4_4() ---------------------------------------------
Matrix
PM4Silt::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
// caution: second coordinate of the first matrix should be in opposite 
// variant form of the first coordinate of second matrix
{
	if ((m1.noCols() != 3) || (m1.noRows() != 3) || (m2.noCols() != 3) || (m2.noRows() != 3))
		opserr << "\n ERROR! PM4Silt::DoubleDot4_4 requires 3-by-3 matrices " << endln;

	return m1*m2;
}
/*************************************************************/
// ToContraviant() ---------------------------------------------
Vector PM4Silt::ToContraviant(const Vector& v1)
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Silt::ToContraviant requires vector of size(3)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=12
	Vector res = v1;
	res(2) *= 0.5;

	return res;
}
/*************************************************************/
// ToCovariant() ---------------------------------------------
Vector PM4Silt::ToCovariant(const Vector& v1)
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Silt::ToCovariant requires vector of size(3)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=12
	Vector res = v1;
	res(2) *= 2.0;

	return res;
}