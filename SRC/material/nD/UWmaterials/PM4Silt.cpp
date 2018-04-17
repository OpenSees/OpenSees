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

// Written: Long Chen, Pedro Arduino
//          Jan 2018, University of Washington

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
OPS_PM4Silt(void)
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
	double oData[27];

	oData[0] = 1.0;      // Fsu
	oData[1] = 101.3;    // p_atm
	oData[2] = 0.75;     // nG
	oData[3] = 0.5;      // h0
	oData[4] = 0.9;      // e_init
	oData[5] = 0.06;     // Lambda
	oData[6] = 32.0;     // phi_cv
	oData[7] = 0.8;      // nb_wet
	oData[8] = 0.5;      // nb_dry
	oData[9] = 0.3;      // nd
	oData[10] = 0.8;     // Ado
	oData[11] = -1;      // ru_max
	oData[12] = -1;      // z_max
	oData[13] = 100.0;   // Cz
	oData[14] = -1;      // Ce
	oData[15] = 3.0;     // Cgd
	oData[16] = 4.0;     // Ckaf
	oData[17] = 0.3;     // nu
	oData[18] = 0.01;    // m
	oData[19] = 0.005;   // crhg
	oData[20] = -1;		 // chg
	oData[21] = 2.0;	 // CG_consol
	oData[22] = 0.0;     // residualP
	oData[23] = 1;       // integratoin scheme
	oData[24] = 0;       // tangent type
	oData[25] = 1.0e-7;	// TolF
	oData[26] = 1.0e-10;	// TolR

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
		oData[15], oData[16], oData[17], oData[18], oData[19], oData[20], oData[21], oData[22], (int)oData[23], (int)oData[24], oData[25], oData[26]);


	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory for nDMaterial PM4Silt material with tag: " << tag << endln;
	}

	return theMaterial;
}

// full constructor
PM4Silt::PM4Silt(int tag, int classTag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nG, double h0,
	double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd, double Ado,
	double ru_max, double z_max, double cz, double ce, double Cgd, double Ckaf, double nu, double m,
	double crhg, double chg, double CG_consol, double residualP, int integrationScheme, int tangentType,
	double TolF, double TolR) : NDMaterial(tag, classTag),
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
	mAlpha_in_true(3),
	mAlpha_in_max(3),
	mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3)
{
	m_Su = Su;
	m_Su_rate = Su_rate;
	m_G0 = G0;
	m_hpo = hpo;
	massDen = mDen;
	m_Fsu = Fsu;
	m_P_atm = P_atm;
	m_nG = nG;
	m_h0 = h0;
	m_e_init = einit;
	m_lambda = lambda;
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nbwet = nbwet;
	m_nbdry = nbdry;
	m_nd = nd;
	m_Ado = Ado;
	m_ru_max = ru_max;
	m_z_max = z_max;
	m_cz = cz;
	m_Cgd = Cgd;
	m_Ckaf = Ckaf;
	m_nu = nu;
	m_m = m;
	m_crhg = crhg;
	m_chg = chg;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = CG_consol;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;
	mresidualP = residualP;
	mIter = 0;

	initialize();
}

// full constructor
PM4Silt::PM4Silt(int tag, double Su, double Su_rate, double G0, double hpo, double mDen, double Fsu, double P_atm, double nG, double h0,
	double einit, double lambda, double phi_cv, double nbwet, double nbdry, double nd, double Ado,
	double ru_max, double z_max, double cz, double ce, double Cgd, double Ckaf, double nu, double m,
	double crhg, double chg, double CG_consol, double residualP, int integrationScheme, int tangentType,
	double TolF, double TolR) : NDMaterial(tag, ND_TAG_PM4Silt),
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
	mAlpha_in_true(3),
	mAlpha_in_max(3),
	mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3)
{
	m_Su = Su;
	m_Su_rate = Su_rate;
	m_G0 = G0;
	m_hpo = hpo;
	massDen = mDen;
	m_Fsu = Fsu;
	m_P_atm = P_atm;
	m_nG = nG;
	m_h0 = h0;
	m_e_init = einit;
	m_lambda = lambda;
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nbwet = nbwet;
	m_nbdry = nbdry;
	m_nd = nd;
	m_Ado = Ado;
	m_ru_max = ru_max;
	m_z_max = z_max;
	m_cz = cz;
	m_Cgd = Cgd;
	m_Ckaf = Ckaf;
	m_nu = nu;
	m_m = m;
	m_crhg = crhg;
	m_chg = chg;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = CG_consol;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;
	mresidualP = residualP;
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
	mAlpha_in_true(3),
	mAlpha_in_max(3),
	mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3)
{
	m_Su = 0.0;
	m_Su_rate = -1.0;
	m_G0 = 0.0;
	m_hpo = 0.0;
	massDen = 0.0;
	m_Fsu = 0.0;
	m_P_atm = 0.0;
	m_nG = 0.0;
	m_h0 = 0.0;
	m_e_init = 0.0;
	m_lambda = 0.0;
	m_Mc = 0.0;
	m_nbwet = 0.0;
	m_nbdry = 0.0;
	m_nd = 0.0;
	m_Ado = 0.0;
	m_ru_max = 0.0;
	m_z_max = 0.0;
	m_cz = 0.0;
	m_Cgd = 0.0;
	m_Ckaf = 0.0;
	m_nu = 0.0;
	m_m = 0.0;
	m_crhg = 0.0;
	m_chg = 0.0;
	m_FirstCall = 0;
	m_PostShake = 0;
	m_CG_consol = 0.0;
	mScheme = 1;
	mTangType = 0;
	mTolF = 1.0e-7;
	mTolR = 1.0e-7;
	mresidualP = 0.0;

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
		clone = new PM4Silt(this->getTag(), m_Su, m_Su_rate, m_G0, m_hpo, massDen, m_Fsu, m_P_atm, m_nG, m_h0, m_e_init,
			m_lambda, phi_cv, m_nbwet, m_nbdry, m_nd, m_Ado, m_ru_max, m_z_max, m_cz, m_ce, m_Cgd, m_Ckaf, m_nu, m_m,
			m_crhg, m_chg, m_CG_consol, mresidualP, mScheme, mTangType, mTolF, mTolR);
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
	mSigma_n = mSigma;
	mEpsilon_n = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	dFabric = mFabric - mFabric_n;
	// update cumulated fabric
	mzcum = mzcum + sqrt(DoubleDot2_2_Contr(dFabric, dFabric) / 2.0);
	mzpeak = fmax(sqrt(DoubleDot2_2_Contr(mFabric, mFabric) / 2.0), mzpeak);
	mFabric_n = mFabric;
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
	PM4Silt  *clone;
	clone = new PM4Silt();
	*clone = *this;
	return clone;
}

const char*
PM4Silt::getType(void) const
{
	return "PlaneStrain";
}

int
PM4Silt::getOrder(void) const
{
	return 3;
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
	else if (strcmp(argv[0], "dilatancy") == 0 || strcmp(argv[0], "dilatancy") == 0)
		return new MaterialResponse(this, 7, this->getTraker());
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
			*(matInfo.theVector) = getTraker();
		return 0;
	default:
		return -1;
	}
}

int
PM4Silt::sendSelf(int commitTag, Channel &theChannel)
{

	int res = 0;
	static Vector data(95);

	data(0) = this->getTag();

	// need to be implemented

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
	static Vector data(95);

	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: PM4Silt::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}
	this->setTag((int)data(0));

	// need to be implemented

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
		mAlpha_n = GetDevPart(initStress) / (p0 + mresidualP);
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
	// initialize chg
	m_chg = fmax(p0 * m_crhg, m_chg);
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
	mMd = m_Mc * exp(m_nd * ksi);
	// the maximum friction angle that can be mobilized near the origin needed to be defined
	mphi_max = 40.0;
	if (ksi < 0) {
		// dense of critical
		double Mb_max = 2 * sin(mphi_max / 180.0 * 3.14159265359);
		double C_Mb = 1 / (pow(Mb_max / m_Mc, 1 / m_nbdry) - 1.0);
		mMb = m_Mc * pow((1 + C_Mb) / (p0 / mpcs + C_Mb), m_nbdry);
	}
	else {
		//loose of critical
		mMb = m_Mc * exp(-1.0 * m_nbwet * ksi);
	}

	// check if initial stresses are inside bounding/dilatancy surface 
	double Mcut = fmax(mMb, mMd);
	double Mfin = sqrt(2) * GetNorm_Contr(GetDevPart(mSigma_n));
	Mfin = Mfin / (p0 + mresidualP);
	if (Mfin > Mcut)
	{
		Vector r = (mSigma_n - p0 * mI1) / (p0 + mresidualP) * Mcut / Mfin;
		mSigma_n = p0 * mI1 + r * (p0 + mresidualP);
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
	mAlpha_in_true = mAlpha_n;
	mAlpha_in_max = mAlpha_n;
	mAlpha_in_min = mAlpha_n;
	mFabric.Zero();
	mFabric_in.Zero();
	mFabric_n.Zero();
	mzpeak = m_z_max / 100000.0;
	mpzp = fmax(p0, m_Pmin) / 100.0;
	mzxp = 0.0;
	m_FirstCall = 1;
	m_pzpFlag = true;

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

	GetElasticModuli(mSig, mK, mG);
	mCe = mCep = mCep_Consistent = GetStiffness(mK, mG);

	return 0;
}

int
PM4Silt::setTrialStrain(const Vector &strain_from_element) {
	mEpsilon = -1.0 * strain_from_element;   // -1.0 is for geotechnical sign convention
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
	return mAlpha;
}
//send back fabric tensor
const Vector
PM4Silt::getFabric()
{
	return mFabric;
}
//send back alpha_in tensor
const Vector
PM4Silt::getAlpha_in()
{
	return mAlpha_in_n;
}
//send back internal parameter for tracking
double
PM4Silt::getTraker()
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
	return mAlpha_in_p;
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
	mSigma_r = -1.0 * (mSigma + mSigma_b);
	return  mSigma_r;  // -1.0 is for geotechnical sign convention
}
/*************************************************************/
const Vector &
PM4Silt::getStrain() {
	mEpsilon_r = -1.0 * mEpsilon;   // -1.0 is for geotechnical sign convention
	return mEpsilon_r;
}
/*************************************************************/
const Vector &
PM4Silt::getElasticStrain() {
	mEpsilonE_r = -1.0 * mEpsilonE;   // -1.0 is for geotechnical sign convention
	return mEpsilonE_r;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void PM4Silt::integrate()
{
	// update alpha_in in case of unloading
	Vector n_tr(3);
	n_tr = GetNormalToYield(mSigma_n + mCe*(mEpsilon - mEpsilon_n), mAlpha_n);

	if (DoubleDot2_2_Contr(mAlpha_n - mAlpha_in_n, n_tr) < 0.0) {
		// This is a loading reversal
		// update pzp
		double p = 0.5 * GetTrace(mSigma_n);
		p = (p <= m_Pmin) ? (m_Pmin) : p;
		double zxpTemp = GetNorm_Contr(mFabric_n) * p;
		if (zxpTemp > mzxp && p > mpzp && m_pzpFlag) {
			mzxp = zxpTemp;
			mpzp = p;
			// only update pzp once
			m_pzpFlag = false;
		}
		// track initial back-stress ratio history 
		for (int ii = 0; ii < 3; ii++) {
			if (mAlpha_in(ii) > 0.0)
				// minimum positive value
				mAlpha_in_min(ii) = fmin(mAlpha_in_min(ii), mAlpha_n(ii));
			else
				// maximum negative value
				mAlpha_in_max(ii) = fmax(mAlpha_in_max(ii), mAlpha_n(ii));
		}
		// update initial back-stress ratio
		if (DoubleDot2_2_Contr(mAlpha_n - mAlpha_in_p, n_tr) < 0.0) {
			// small unload-reload cycle, update initial back-stress ratio using apparent back-stress ratio
			mAlpha_in_p = mAlpha_in_true;
			mAlpha_in = mAlpha_n;
			// update components of apparent initial back-stress ratio
			for (int ii = 0; ii < 3; ii++) {
				if (mAlpha_in(ii) > 0.0)
					// positive loading direction
					mAlpha_in(ii) = fmax(0.0, mAlpha_in_min(ii));
				else
					// negative loading direction
					mAlpha_in(ii) = fmin(0.0, mAlpha_in_max(ii));
			}
		}
		else {
			// update initial back-stress ratio using true back-stress ratio
			mAlpha_in_p = mAlpha_in_true;
			mAlpha_in = mAlpha_n;
		}
		mAlpha_in_true = mAlpha_n;
		mFabric_in = mFabric_n;
	}
	else
		mAlpha_in = mAlpha_in_n;

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
	dStrain = NextStrain - CurStrain;
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + dStrain;
	aCep_Consistent = aCep = aC = GetStiffness(K, G);
	NextStress = CurStress + DoubleDot4_2(aC, dStrain);
	double p = 0.5 * GetTrace(NextStress) + mresidualP;
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
	Vector dSigma(3), dDevStrain(3), n(3);

	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + NextStrain - CurStrain;
	dVolStrain = GetTrace(NextStrain - CurStrain);
	dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
	aC = GetStiffness(K, G);
	dSigma = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
	NextStress = CurStress + dSigma;

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
		mTracker = 0.0;
		// Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, 
		// NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, G, K , aC, aCep, aCep_Consistent);

		return;

	}
	else if (fn < -mTolF) {
		// This is a transition from elastic to plastic
		elasticRatio = IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, 0.0, 1.0);
		dSigma = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
		(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
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
			dSigma = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
			(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
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
	double CurVoidRatio, Cka, h, p, dVolStrain, D;
	Vector n(3), R(3), alphaD(3), dPStrain(3), b(3), dDevStrain(3), r(3);
	Vector dSigma(3), dAlpha(3), dFabric(3);

	CurVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	p = 0.5 * GetTrace(CurStress) + mresidualP;
	p = p < (m_Pmin + mresidualP) ? (m_Pmin + mresidualP) : p;
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	// using NextStress instead of CurStress to get correct n
	GetStateDependent(NextStress, CurAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum
		, mzpeak, mpzp, mMcur, CurVoidRatio, n, D, R, mKp, alphaD, Cka, h, b);
	dVolStrain = GetTrace(NextStrain - CurStrain);
	dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3.0 * mI1;
	r = GetDevPart(CurStress) / p;
	double temp4 = mKp + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);

	if (fabs(temp4) < small) {
		// Neutral loading
		dSigma.Zero();
		dAlpha.Zero();
		dFabric.Zero();
		dPStrain = dDevStrain + dVolStrain * mI1;
	}
	else {
		NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
		mDGamma = NextL;
		if (NextL < 0) {
			if (debugFlag) {
				opserr << "NextL is smaller than 0\n";
				opserr << "NextL = " << NextL << endln;
			}
			dSigma = 2 * mG * ToContraviant(dDevStrain) + mK * dVolStrain * mI1;
			dAlpha = GetDevPart(NextStress + dSigma) / (0.5 * GetTrace(NextStress + dSigma) + mresidualP)
				- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress) + mresidualP);
			dFabric.Zero();
			dPStrain.Zero();
		}
		else {
			dSigma = 2.0*mG*mIIcon*dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
				(2.0 * mG * n + mK * D * mI1);
			// update fabric
			if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
				dFabric = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
				NextFabric = CurFabric + dFabric;
			}
			// update alpha
			dAlpha = two3 * NextL * h * b;
			dPStrain = NextL * mIIco * R;
		}
	}
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
	NextStress = CurStress + dSigma;
	NextAlpha = CurAlpha + dAlpha;
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
		exp_int = &PM4Silt::ForwardEuler;
		break;
	}
	Vector StrainInc(3); StrainInc = NextStrain - CurStrain;
	double maxInc = StrainInc(0);

	for (int ii = 1; ii < 3; ii++)
		if (fabs(StrainInc(ii)) > fabs(maxInc))
			maxInc = StrainInc(ii);

	if (fabs(maxInc) > maxStrainInc) {
		int numSteps = (int)floor(fabs(maxInc) / maxStrainInc) + 1;
		StrainInc = (NextStrain - CurStrain) / (double)numSteps;

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
	double dVolStrain, p, Cka, temp4, curStepError, q, stressNorm, h, D;
	Vector n(3), R1(3), R2(3), alphaD(3), dDevStrain(3), r(3), b(3);
	Vector nStress(3), nAlpha(3), nFabric(3);
	Vector dSigma1(3), dSigma2(3), dAlpha1(3), dAlpha2(3), dAlpha(3), dFabric1(3), dFabric2(3), dPStrain1(3), dPStrain2(3);
	double T = 0.0, dT = 1.0, dT_min = 1e-4, TolE = 1e-5;

	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextStress = CurStress;
	NextAlpha = CurAlpha;
	NextFabric = CurFabric;

	this->GetElasticModuli(NextStress, K, G, mMcur, mzcum);

	p = 0.5 * GetTrace(CurStress) + mresidualP;
	if (p < m_Pmin / 5.0 + mresidualP)
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
		p = 0.5 * GetTrace(NextStress) + mresidualP;
		// Calc Delta 1
		GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, NextFabric, mFabric_in, G, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R1, mKp, alphaD, Cka, h, b);

		r = GetDevPart(NextStress) / p;

		temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma1.Zero();
			dAlpha1.Zero();
			dFabric1.Zero();
			dPStrain1 = dDevStrain + dVolStrain * mI1;
		}
		else {
			NextL = (2 * G * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * K * dVolStrain) / temp4;
			if (NextL < 0) {
				if (debugFlag) {
					opserr << "1 NextL is smaller than 0\n";
					opserr << "NextL = " << NextL << endln;
				}
				dSigma1 = 2 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1;
				dAlpha1 = 2.0*(GetDevPart(NextStress + dSigma1) / GetTrace(NextStress + dSigma1) - GetDevPart(NextStress) / GetTrace(NextStress));
				dFabric1.Zero();
				dPStrain1.Zero();
			}
			else {
				dSigma1 = 2.0 * G * mIIcon * dDevStrain + K*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * G * n + K * D * mI1);
				// update fabric
				// if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
					dFabric1 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
				// }
				dPStrain1 = NextL * mIIco * R1;
				dAlpha1 = two3 * NextL * h * b;
			}
		}
		//Calc Delta 2
		p = 0.5 * GetTrace(NextStress + dSigma1) + mresidualP;
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

		GetStateDependent(NextStress + dSigma1, NextAlpha + dAlpha1, alpha_in, alpha_in_p, NextFabric + dFabric1, mFabric_in, G, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R2, mKp, alphaD, Cka, h, b);
		r = GetDevPart(NextStress + dSigma1) / p;

		temp4 = mKp + 2 * G - K* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) {
			// neutral loading
			dSigma2.Zero();
			dAlpha2.Zero();
			dFabric2.Zero();
			dPStrain2 = dDevStrain + dVolStrain * mI1;
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
				dSigma2 = 2 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1;
				dAlpha2 = 2.0*(GetDevPart(NextStress + dSigma2) / GetTrace(NextStress + dSigma2) - GetDevPart(NextStress) / GetTrace(NextStress));
				dFabric2.Zero();
				dPStrain2.Zero();
			}
			else {
				dSigma2 = 2.0 * G * mIIcon * dDevStrain + K*dVolStrain*mI1 - Macauley(NextL)*
					(2.0 * G * n + K * D * mI1);
				// update fabric
				// if (DoubleDot2_2_Contr(alphaD - (CurAlpha + dAlpha1), n) < 0.0) {
					dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric + dFabric1);
				// }
				dPStrain2 = NextL * mIIco * R2;
				dAlpha2 = two3 * NextL * h * b;
			}
		}

		nStress = NextStress + 0.5 * (dSigma1 + dSigma2);
		nFabric = NextFabric + 0.5 * (dFabric1 + dFabric2);
		// update alpha

		dAlpha = 0.5 * (dAlpha1 + dAlpha2);
		nAlpha = NextAlpha + dAlpha;

		p = 0.5 * GetTrace(nStress) + mresidualP;
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
		if (stressNorm < 0.5)
			curStepError = GetNorm_Contr(dSigma2 - dSigma1);
		else
			curStepError = GetNorm_Contr(dSigma2 - dSigma1) / (2 * stressNorm);

		if (curStepError > TolE) {
			q = fmax(0.8 * sqrt(TolE / curStepError), 0.1);
			if (dT == dT_min) {
				// opserr << "reached dT_min\n";
				NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
				NextStress = nStress;
				NextAlpha = nAlpha;
				// Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
				//Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
				T += dT;
			}
			dT = fmax(q * dT, dT_min);
		}
		else {
			NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
			NextStress = nStress;
			NextAlpha = nAlpha;
			NextFabric = nFabric;
			// Stress_Correction(NextStress, NextAlpha, alpha_in, alpha_in_p, CurFabric, NextVoidRatio);
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
	double dVolStrain, p, Cka, D, K_p, temp4, h;
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

	p = 0.5 * GetTrace(CurStress) + mresidualP;
	if (p < m_Pmin / 5.0 + mresidualP)
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
		p = 0.5 * GetTrace(NextStress) + mresidualP;
		// Calc Delta 1
		GetStateDependent(NextStress, NextAlpha, alpha_in, alpha_in_p, NextFabric, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R1, K_p, alphaD, Cka, h, b);

		r = GetDevPart(NextStress) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
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
				dAlpha1 = GetDevPart(NextStress + dSigma1) / (0.5 * GetTrace(NextStress + dSigma1) + mresidualP)
					- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress) + mresidualP);
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
		p = 0.5 * GetTrace(NextStress + 0.5 * dSigma1) + mresidualP;

		GetStateDependent(NextStress + 0.5 * dSigma1, CurAlpha + 0.5 * dAlpha1, alpha_in, alpha_in_p, NextFabric + 0.5 * dFabric1, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R2, K_p, alphaD, Cka, h, b);
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
				dAlpha2 = GetDevPart(NextStress + dSigma2) / (0.5 * GetTrace(NextStress + dSigma2) + mresidualP)
					- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress) + mresidualP);
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
		p = 0.5 * GetTrace(NextStress + 0.5 * dSigma2) + mresidualP;

		GetStateDependent(NextStress + 0.5 * dSigma2, CurAlpha + 0.5 * dAlpha2, alpha_in, alpha_in_p, NextFabric + 0.5 * dFabric2, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R3, K_p, alphaD, Cka, h, b);
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
				dAlpha3 = GetDevPart(NextStress + dSigma3) / (0.5 * GetTrace(NextStress + dSigma3) + mresidualP)
					- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress) + mresidualP);
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
		p = 0.5 * GetTrace(NextStress + dSigma3) + mresidualP;

		GetStateDependent(NextStress + dSigma3, CurAlpha + dAlpha3, alpha_in, alpha_in_p, NextFabric + dFabric3, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R4, K_p, alphaD, Cka, h, b);
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
				dAlpha4 = GetDevPart(NextStress + dSigma4) / (0.5 * GetTrace(NextStress + dSigma4) + mresidualP)
					- GetDevPart(NextStress) / (0.5 * GetTrace(NextStress) + mresidualP);
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
	Vector dSigma(3), dSigma0(3), dSigma1(3), strainInc(3);

	strainInc = NextStrain - CurStrain;

	if (a0 < 0.0 || a1 > 1.0) {
		opserr << "a0 = " << a0 << "a1 = " << a1 << endln;
	}
	//GetElasticModuli(CurStress, K, G, mzcum);
	dSigma0 = a0 * DoubleDot4_2(mCe, strainInc);
	f0 = GetF(CurStress + dSigma0, CurAlpha);

	dSigma1 = a1 * DoubleDot4_2(mCe, strainInc);
	f1 = GetF(CurStress + dSigma1, CurAlpha);

	for (int i = 1; i <= 10; i++)
	{
		a = a1 - f1 * (a1 - a0) / (f1 - f0);
		dSigma = a * DoubleDot4_2(mCe, strainInc);
		f = GetF(CurStress + dSigma, CurAlpha);
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
	Vector dSigma(3), dSigma0(3), dSigma1(3), strainInc(3);
	bool flag = false;

	strainInc = NextStrain - CurStrain;

	f0 = GetF(CurStress, CurAlpha);
	fs = f0;

	// GetElasticModuli(CurStress, K, G, mzcum);
	dSigma = DoubleDot4_2(mCe, strainInc);

	for (int i = 1; i < 10; i++)
	{
		da = (a1 - a0) / nSub;
		for (int k = 1; k < nSub; k++) {
			a = a0 + da;
			f = GetF(CurStress + a * dSigma, CurAlpha);
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
	double lambda, D, K_p, Cka, h, p, fr;
	Matrix aC(3, 3);
	// Vector CurStress = NextStress;

	int maxIter = 25;
	p = 0.5 * GetTrace(NextStress) + mresidualP;
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
			Vector nStress = NextStress;
			Vector nAlpha = NextAlpha;
			for (int i = 1; i <= maxIter; i++) {
				r = GetDevPart(nStress) / p;
				GetStateDependent(nStress, nAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum
					, mzpeak, mpzp, mMcur, NextVoidRatio, n, D, R, K_p, alphaD, Cka, h, b);
				aC = GetStiffness(mK, mG);
				dSigmaP = DoubleDot4_2(aC, mDGamma * ToCovariant(R));
				aBar = two3 * h * b;
				dfrOverdSigma = n - 0.5 * DoubleDot2_2_Contr(n, r) * mI1;
				dfrOverdAlpha = -p * n;
				lambda = fr / (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP) - DoubleDot2_2_Contr(dfrOverdAlpha, aBar));
				if (fabs(GetF(nStress - lambda * dSigmaP, nAlpha + lambda * aBar)) < fabs(fr))
				{
					nStress -= lambda * dSigmaP;
					nAlpha += lambda * aBar;
				}
				else {
					lambda = fr / DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma);
					nStress -= lambda * dfrOverdSigma;
				}
				fr = GetF(nStress, nAlpha);
				if (fabs(fr) < mTolF) {
					NextStress = nStress;
					NextAlpha = nAlpha;
					return;
				}

				p = fmax(0.5 * GetTrace(nStress) + mresidualP, m_Pmin);
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

				fr_old = GetF(mSigma + alpha_mid * dSigma, NextAlpha);
				if (fabs(fr_old) < mTolF) {
					NextStress = mSigma + alpha_mid * dSigma;
					break;
				}
			}

			// // stress state ouside yield surface
			// double CurDr = (m_emax - NextVoidRatio) / (m_emax - m_emin);
			// Vector nStress = NextStress;
			// Vector nAlpha = NextAlpha;
			// for (int i = 1; i <= maxIter; i++) {
			// 	// Sloan, Abbo, Sheng 2001, Refined explicit integration of elastoplastic models with automatic 
			// 	// error control
			// 	r = GetDevPart(nStress) / p;
			// 	GetStateDependent(nStress, nAlpha, alpha_in, alpha_in_p, CurFabric, mFabric_in, mG, mzcum
			// 		, mzpeak, mpzp, mMcur, CurDr, n, D, R, K_p, alphaD, Cka, h, b);
			// 	dSigmaP = DoubleDot4_2(mCe, ToCovariant(R));
			// 	aBar = h * b;
			// 	dfrOverdSigma = n - 0.5 * DoubleDot2_2_Contr(n, r) * mI1;
			// 	dfrOverdAlpha = -p * n;
			// 	lambda = fr / (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP) - DoubleDot2_2_Contr(dfrOverdAlpha, aBar));
			// 
			// 	if (fabs(GetF(nStress - lambda * dSigmaP, nAlpha + lambda * aBar)) < fabs(fr))
			// 	{
			// 		nStress -= lambda * dSigmaP;
			// 		nAlpha += lambda * aBar;
			// 	}
			// 	else {
			// 		lambda = fr / DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma);
			// 		if (fabs(GetF(nStress - lambda * dfrOverdSigma, nAlpha)) < fabs(fr))
			// 			nStress -= lambda * dfrOverdSigma;
			// 		else
			// 		{
			// 			if (debugFlag)
			// 				opserr << "PM4Silt::StressCorrection() Couldn't decrease the yield function." << endln;
			// 			return;
			// 		}
			// 	}
			// 
			// 	fr = GetF(nStress, nAlpha);
			// 	if (fabs(fr) < mTolF) {
			// 		NextStress = nStress;
			// 		NextAlpha = nAlpha;
			// 		break;
			// 	}
			// 
			// 	if (i == maxIter) {
			// 		if (debugFlag)
			// 			opserr << "Still outside with f =  " << fr << endln;
			// 		Vector dSigma = NextStress - CurStress;
			// 		double alpha_up = 1.0;
			// 		double alpha_mid = 0.5;
			// 		double alpha_down = 0.0;
			// 		double fr_old = GetF(CurStress + alpha_mid * dSigma, NextAlpha);
			// 		for (int jj = 0; jj < maxIter; jj++) {
			// 			if (fr_old < 0.0) {
			// 				alpha_down = alpha_mid;
			// 				alpha_mid = 0.5 * (alpha_up + alpha_mid);
			// 			}
			// 			else {
			// 				alpha_up = alpha_mid;
			// 				alpha_mid = 0.5 * (alpha_down + alpha_mid);
			// 			}
			// 
			// 			fr_old = GetF(CurStress + alpha_mid * dSigma, NextAlpha);
			// 			if (fabs(fr_old) < mTolF) {
			// 				NextStress = CurStress + alpha_mid * dSigma;
			// 				break;
			// 			}
			// 			// }
			// 			if (jj == maxIter) opserr << "stress still outside!!" << endln;
			// 		}
			// 		p = 0.5 * GetTrace(nStress) + mresidualP;
			// 	}
			// }
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
			if (abs(f) < mTolF) {
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
	double p = 0.5 * GetTrace(nStress) + mresidualP;
	s = s - p * nAlpha;
	double f = GetNorm_Contr(s) - root12 * m_m * p;
	return f;
}
/*************************************************************/
// GetPSI() ---------------------------------------------------
double
PM4Silt::GetKsi(const double& e, const double& p)
{
	double pn = p;
	pn = (pn <= m_Pmin) ? (m_Pmin) : pn;
	double ksi = e - me0 + m_lambda * log(101.3 * pn / m_P_atm / m_Fsu);
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
	double Csr = 1 - Csr0 *pow((Mcur / mMb), msr);
	Csr = (Csr < 0.4) ? 0.4 : Csr;
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
	m_nu = (m_nu == 0.5) ? 0.4999 : m_nu;
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
	m_nu = (m_nu == 0.5) ? 0.4999 : m_nu;
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
	double p = 0.5 * GetTrace(NextStress) + mresidualP;
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
	Vector devStress(3); devStress = GetDevPart(stress);
	double p = 0.5 * GetTrace(stress) + mresidualP;
	Vector n(3);
	if (fabs(p) < small) {
		n.Zero();
	}
	else {
		n = devStress - p * alpha;
		double normN = GetNorm_Contr(n);
		normN = (normN < small) ? 1.0 : normN;
		n = n / normN;
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
	, Vector &alphaD, double &Cka, double &h, Vector &b)
{
	double p = 0.5 * GetTrace(stress);
	if (p <= m_Pmin) p = m_Pmin;
	double ksi = GetKsi(CurVoidRatio, p);
	n = GetNormalToYield(stress, alpha);

	mMd = m_Mc * exp(m_nd * ksi);
	if (ksi < 0) {
		// dense of critical
		double phi_max = 40.0;   // max phi at p = 1kPa
		double Mb_max = 2 * sin(phi_max / 180.0 * 3.14159265359);
		double C_Mb = 1 / (pow(Mb_max / m_Mc, 1 / m_nbdry) - 1.0);
		mMb = m_Mc * pow((1 + C_Mb) / (p / mpcs + C_Mb), m_nbdry);
	}
	else {
		//loose of critical
		mMb = m_Mc * exp(-1.0 * m_nbwet * ksi);
	}

	Vector alphaB = root12 * (mMb - m_m) * n;
	alphaD = root12 * (mMd - m_m) * n;
	double Czpk1 = zpeak / (zcum + m_z_max / 5.0);
	double Czpk2 = zpeak / (zcum + m_z_max / 100.0);
	double Cpzp2 = Macauley((pzp - p)) / (Macauley((pzp - p)) + m_Pmin);
	double Cg1 = m_h0 / 200.0;
	double Ckp = 2.0;

	b = alphaB - alpha;
	Cka = 1.0 + m_Ckaf / (1.0 + pow(2.5*Macauley(DoubleDot2_2_Contr(alpha - mAlpha_in_true, n)), 2.0))*Cpzp2*Czpk1;
	double AlphaAlphaBDotN = fabs(DoubleDot2_2_Contr(b, n));

	// updataed K_p formulation following PM4Silt V3.1. mAlpha_in is apparent back-stress ratio. 
	if (DoubleDot2_2_Contr(alpha - alpha_in_p, n) <= 0) {
		h = 1.5 * G * m_h0 / p / (exp(DoubleDot2_2_Contr(alpha - mAlpha_in, n)) - 1 + Cg1) / sqrt(AlphaAlphaBDotN) *
			Cka / (1 + Ckp * zpeak / m_z_max * Macauley(AlphaAlphaBDotN) * sqrt(1 - Czpk2));
		h = h * (DoubleDot2_2_Contr(alpha - mAlpha_in, n) + Cg1) / (DoubleDot2_2_Contr(alpha - mAlpha_in_true, n) + Cg1);
	}
	else {
		h = 1.5 * G * m_h0 / p / (exp(DoubleDot2_2_Contr(alpha - mAlpha_in, n)) - 1 + Cg1) / sqrt(AlphaAlphaBDotN) *
			Cka / (1 + Ckp * zpeak / m_z_max * Macauley(AlphaAlphaBDotN) * sqrt(1 - Czpk2));
	}

	K_p = two3 * h * p * DoubleDot2_2_Contr(b, n);
	// bound K_p to non - negative, following flac practice
	K_p = fmax(0.0, K_p);
	double Czin1 = Macauley(1.0 - exp(-2.0*abs((DoubleDot2_2_Contr(fabric_in, n) - DoubleDot2_2_Contr(fabric, n)) / m_z_max)));
	// rotated dilatancy surface
	double Crot1 = fmax((1.0 + 2 * Macauley(DoubleDot2_2_Contr(-1.0 * fabric, n)) / (sqrt(2.0)*m_z_max)*(1 - Czin1)), 1.0);
	double Mdr = mMd / Crot1;
	Vector alphaDr = root12 * (Mdr - m_m) * n;
	// dilation
	if (DoubleDot2_2_Contr(alphaDr - alpha, n) <= 0) {
		double Cpzp = 1.0 / (1.0 + pow((2.5* p / mpzp), 5.0));
		double Czin2 = (1 + Czin1*(zcum - zpeak) / 3.0 / m_z_max) / (1 + 3.0 * Czin1*(zcum - zpeak) / 3 / m_z_max);
		double temp = pow((1.0 - Macauley(DoubleDot2_2_Contr(-1.0 * fabric, n)) * root12 / zpeak), 3);
		double Ad = m_Ado * Czin2 / ((pow(zcum, 2) / m_z_max) * temp * pow(m_ce, 2.0) * Cpzp * Czin1 + 1.0);
		D = Ad *(-1.0 * Macauley(-DoubleDot2_2_Contr(alphaD - alpha, n)));
		double Drot = Ad * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max) * DoubleDot2_2_Contr(alphaDr - alpha, n) / 3.0;
		if (D > Drot) {
			D = D + (Drot - D)*Macauley(mMb - Mcur) / (Macauley(mMb - Mcur) + 0.01);
		}
		if (p <= 2 * m_Pmin) {
			D = fmin(D, -3.5 * m_Ado * Macauley(mMb - mMd) * (2 * m_Pmin - p) / m_Pmin);
		}
	}
	else {
		//contraction
		double hp = m_hpo * exp(-0.7 + 0.2 * pow(Macauley(3 - ksi / m_lambda), 2.0));
		double Crot2 = 1 - Czpk2;
		double Cdz = fmax((1 - Crot2*sqrt(2.0)*zpeak / m_z_max)*(m_z_max / (m_z_max + Crot2*zcum)), 1 / (1 + m_z_max / 2.0));
		double Cwet = fmin(1.0, 1.0 / (1 / (1 + pow(0.02 / AlphaAlphaBDotN, 4)) + 1 / (1 + pow(ksi / m_lambda / 0.1, 2.0))));
		double Adc = m_Ado * (1 + Macauley(DoubleDot2_2_Contr(fabric, n))) / (hp * Cdz * Cwet);
		double Cin = 2.0 * Macauley(DoubleDot2_2_Contr(fabric, n)) / sqrt(2.0) / m_z_max;
		D = fmin(Adc * pow((DoubleDot2_2_Contr(alpha - mAlpha_in_true, n) + Cin), 2), 1.5 * m_Ado) *
			DoubleDot2_2_Contr(alphaD - alpha, n) / (DoubleDot2_2_Contr(alphaD - alpha, n) + 0.10);
		// Apply a factor to D so it doesn't go very big when p is small
		double C_pmin;
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