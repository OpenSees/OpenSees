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

// Written: Haoyuan Liu (TUDelft) Jose Abell (UANDES, github.com/jaabell) 
//          and Federico Pisano (TUDelft). 
//          Based on the 2016 version of ManzariDafalias by 
//          Pedro Arduino's group at UW. We acknowledge their great work 
//          and their generosity in sharing quality code in open source. 
//          We stand on the shoulders of giants!


// Description: This file contains the implementation for the SAniSandMS class.

//#include <SAniSandMS.h>
#include "SAniSandMS.h"
//#include <SAniSandMS3D.h>
#include "SAniSandMS3D.h"
//#include <SAniSandMSPlaneStrain.h>
#include "SAniSandMSPlaneStrain.h"
#include <MaterialResponse.h>
#include <string.h>

// #if defined(_WIN32) || defined(_WIN64)
#include <algorithm>
#define fmax std::max
#define fmin std::min
// #endif
#include <fenv.h>

#define INT_MAXENE_MFE    0
#define INT_ModifiedEuler 1
#define INT_BackwardEuler 2
#define INT_RungeKutta    3
#define INT_MAXENE_FE     4
#define INT_ForwardEuler  5
#define INT_MAXENE_RK     6
#define INT_MAXSTR_MFE    7
#define INT_MAXSTR_RK     8
#define INT_MAXSTR_FE     9

const double        SAniSandMS::one3 = 1.0 / 3.0;
const double        SAniSandMS::two3 = 2.0 / 3.0;
const double        SAniSandMS::root23 = sqrt(2.0 / 3.0);
const double        SAniSandMS::small = 1e-10;
const double        SAniSandMS::maxStrainInc = 1e-5;
const bool          SAniSandMS::debugFlag = false;
const char unsigned SAniSandMS::mMaxSubStep = 10;
char  unsigned      SAniSandMS::mElastFlag = 1;

Vector              SAniSandMS::mI1(6);
Matrix              SAniSandMS::mIIco(6, 6);
Matrix              SAniSandMS::mIIcon(6, 6);
Matrix              SAniSandMS::mIImix(6, 6);
Matrix              SAniSandMS::mIIvol(6, 6);
Matrix              SAniSandMS::mIIdevCon(6, 6);
Matrix              SAniSandMS::mIIdevMix(6, 6);
Matrix              SAniSandMS::mIIdevCo(6, 6);
SAniSandMS::initTensors SAniSandMS::initTensorOps;

static int numSAniSandMSMaterials = 0;

void *
OPS_SAniSandMSMaterial(void)
{
	// feenableexcept(FE_DIVBYZERO);// | FE_INVALID );//| FE_OVERFLOW);

	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numSAniSandMSMaterials == 0)
	{
		opserr << "SAniSandMS nDmaterial - \n"
			<< "          By:  Haoyuan Liu (Student, TU Delft), \n"
			<< "               Jose Abell (Prof. Universidad de los Andes, Chile) and \n"
			<< "               Federico Pisano (Prof. TU Delft) \n\n"
			<< "          From original implementation of Manzari-Dafalias by: \n"
			<< "                A.Ghofrani, P.Arduino, U.Washington\n";
	}
	numSAniSandMSMaterials++;

	NDMaterial *theMaterial = 0;


	if (numArgs < 20) {
		opserr << "Want: nDMaterial SAniSandMS tag? G0? nu? e_init? Mc? c? lambda_c? e0? ksi?" <<
			" P_atm? m? h0? Ch? nb? A0? nd? zeta? mu0? beta? Rho? < IntScheme? TanType? JacoType? TolF? TolR?>" << endln;
		return 0;
	}

	int    tag;
	double dData[19];
	int flags[5];
	double oData[2];

	flags[0] = 3;          // IntScheme
	flags[1] = 2;          // TanType
	flags[2] = 1;          // JacoType
	oData[0] = 1.0e-7;     // TolF
	oData[1] = 1.0e-7;     // TolR


	int numData = 1;
	if (OPS_GetInt(&numData, &tag) != 0) {
		opserr << "WARNING 1: invalid nDMaterial SAniSandMS material tag" << endln;
		return 0;
	}

	numData = 19;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING 2: invalid material data for nDMaterial SAniSandMS material  with tag: " << tag << endln;
		return 0;
	}

	// if(SAniSandMS::debugFlag)

	int one = 1;
	numData = numArgs - 19;
	if (numData != 0)
	{
		int pos = 0;
		while (pos < std::min(numData, 3))
		{
			OPS_GetInt(&one, &flags[pos]);
			++pos;
		}
		numData -= 5;
		pos = 0;
		while (pos < std::min(numData, 2))
		{
			OPS_GetDouble(&one, &oData[pos]);
			++pos;
		}
	}
	// if (OPS_GetDouble(&numData, oData) != 0) {
	//     opserr << "WARNING 3 invalid material data for nDMaterial SAniSandMS material  with tag: " << tag << endln;
	//     return 0;
	// }

	theMaterial = new SAniSandMS(tag, ND_TAG_SAniSandMS, dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
		dData[6], dData[7], dData[8], dData[9], dData[10], dData[11],
		dData[12], dData[13], dData[14], dData[15], dData[16],
		dData[17], dData[18], 
		flags[0], flags[1], flags[2],
		oData[0], oData[1]);


	if (theMaterial == 0) {
		opserr << "WARNING ran out of memory for nDMaterial SAniSandMS material with tag: " << tag << endln;
	}

	return theMaterial;
}

// full constructor
SAniSandMS::SAniSandMS(int tag, double G0, double nu,
	double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double zeta, double mu0,
	double beta,
	double mDen,
	int integrationScheme, int tangentType,
	int JacoType, double TolF, double TolR) : NDMaterial(tag, ND_TAG_SAniSandMS),
	mEpsilon(6),
	mEpsilon_n(6),
	mSigma(6),
	mSigma_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
	mAlpha(6),
	mAlpha_n(6),
	mAlphaM(6),
	mAlphaM_n(6),
	malpha_in(6),
	malpha_in_n(6),
	mCe(6, 6),
	mCep(6, 6),
	mCep_Consistent(6, 6)
{
	num_load_reversals = 0;
	num_load_reversals_n = 0;
	m_G0 = G0;
	m_nu = nu;
	m_e_init = e_init;
	m_Mc = Mc;
	m_c = c;
	m_lambda_c = lambda_c;
	m_e0 = e0;
	m_ksi = ksi;
	m_P_atm = P_atm;
	m_m = m;
	m_h0 = h0;
	m_ch = ch;
	m_nb = nb;
	m_A0 = A0;
	m_nd = nd;

	m_zeta = zeta;
	m_mu0 = mu0;
	m_beta = beta;

	
	mMM_plus = m;
	mMM_plus_n = m;
	mMM_minus = 0;
	mMM_minus_n = 0;

	// For debuggin
	// opserr << "SAniSandMS::SAniSandMS(1)" << endln;
	// opserr << "G0 = " << G0 << endln;
	// opserr << "nu = " << nu << endln;
	// opserr << "e_init = " << e_init << endln;
	// opserr << "Mc = " << Mc << endln;
	// opserr << "c = " << c << endln;
	// opserr << "lambda_c = " << lambda_c << endln;
	// opserr << "e0 = " << e0 << endln;
	// opserr << "ksi = " << ksi << endln;
	// opserr << "P_atm = " << P_atm << endln;
	// opserr << "m = " << m << endln;
	// opserr << "h0 = " << h0 << endln;
	// opserr << "ch = " << ch << endln;
	// opserr << "nb = " << nb << endln;
	// opserr << "A0 = " << A0 << endln;
	// opserr << "nd = " << nd << endln;
	// opserr << "zeta = " << zeta << endln;
	// opserr << "mu0 = " << mu0 << endln;
	// opserr << "beta = " << m_beta << endln;
	// opserr << "mDen = " << mDen << endln;

	// opserr << "integrationScheme = " << integrationScheme << endln;
	// opserr << "tangentType = " << tangentType << endln;
	// opserr << "JacoType = " << JacoType << endln;
	// opserr << "TolF = " << TolF << endln;
	// opserr << "TolR = " << TolR << endln;


	massDen = mDen;
	mTolF = TolF;
	mTolR = TolR;
	mJacoType = JacoType;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mUseElasticTan = false;
	mIter = 0;

	m_firstLoading = true;

	initialize();
}

// full constructor
SAniSandMS::SAniSandMS(int tag, int classTag, double G0, double nu,
	double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double zeta, double mu0,
	double beta,
	double mDen,
	int integrationScheme, int tangentType,
	int JacoType, double TolF, double TolR) : NDMaterial(tag, classTag),
	mEpsilon(6),
	mEpsilon_n(6),
	mSigma(6),
	mSigma_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
	mAlpha(6),
	mAlpha_n(6),
	mAlphaM(6),
	mAlphaM_n(6),
	malpha_in(6),
	malpha_in_n(6),
	mCe(6, 6),
	mCep(6, 6),
	mCep_Consistent(6, 6)
{
	num_load_reversals = 0;
	num_load_reversals_n = 0;
	m_G0 = G0;
	m_nu = nu;
	m_e_init = e_init;
	m_Mc = Mc;
	m_c = c;
	m_lambda_c = lambda_c;
	m_e0 = e0;
	m_ksi = ksi;
	m_P_atm = P_atm;
	m_m = m;
	m_h0 = h0;
	m_ch = ch;
	m_nb = nb;
	m_A0 = A0;
	m_nd = nd;
	m_zeta = zeta;
	m_mu0 = mu0;
	m_beta = beta;

	mMM_plus = m;
	mMM_plus_n = m;
	mMM_minus = 0;
	mMM_minus_n = 0;

//	opserr << "SAniSandMS::SAniSandMS(1)" << endln;
//	opserr << "G0 = " << G0 << endln;
//	opserr << "nu = " << nu << endln;
//	opserr << "e_init = " << e_init << endln;
//	opserr << "Mc = " << Mc << endln;
//	opserr << "c = " << c << endln;
//	opserr << "lambda_c = " << lambda_c << endln;
//	opserr << "e0 = " << e0 << endln;
//	opserr << "ksi = " << ksi << endln;
//	opserr << "P_atm = " << P_atm << endln;
//	opserr << "m = " << m << endln;
//	opserr << "h0 = " << h0 << endln;
//	opserr << "ch = " << ch << endln;
//	opserr << "nb = " << nb << endln;
//	opserr << "A0 = " << A0 << endln;
//	opserr << "nd = " << nd << endln;
//	opserr << "zeta = " << zeta << endln;
//	opserr << "mu0 = " << mu0 << endln;
//	opserr << "beta = " << m_beta << endln;
//	opserr << "mDen = " << mDen << endln;
//	opserr << "integrationScheme = " << integrationScheme << endln;
//	opserr << "tangentType = " << tangentType << endln;
//	opserr << "JacoType = " << JacoType << endln;
//	opserr << "TolF = " << TolF << endln;
//	opserr << "TolR = " << TolR << endln;
	

	massDen = mDen;
	mTolF = TolF;
	mTolR = TolR;
	mJacoType = JacoType;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mUseElasticTan = false;
	mIter = 0;

	m_firstLoading = true;


	initialize();
}

// null constructor
SAniSandMS::SAniSandMS()
	: NDMaterial(),
	mEpsilon(6),
	mEpsilon_n(6),
	mSigma(6),
	mSigma_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
	mAlpha(6),
	mAlpha_n(6),
	mAlphaM(6),
	mAlphaM_n(6),
	malpha_in(6),
	malpha_in_n(6),
	mCe(6, 6),
	mCep(6, 6),
	mCep_Consistent(6, 6)
{
	num_load_reversals = 0;
	num_load_reversals_n = 0;
	m_G0 = 0.0;
	m_nu = 0.0;
	m_e_init = 0.0;
	m_Mc = 0.0;
	m_c = 0.0;
	m_lambda_c = 0.0;
	m_e0 = 0.0;
	m_ksi = 0.0;
	m_P_atm = 0.0;
	m_m = 0.0;
	m_h0 = 0.0;
	m_ch = 0.0;
	m_nb = 0.0;
	m_A0 = 0.0;
	m_nd = 0.0;
	m_zeta = 0.0;
	m_mu0 = 0.0;
	m_beta = 0.0;



	mMM_plus = 0;
	mMM_plus_n = 0;
	mMM_minus = 0;
	mMM_minus_n = 0;
	
	massDen = 0.0;
	mTolF = 1.0e-7;
	mTolR = 1.0e-7;
	mJacoType = 1;
	mScheme = 3;
	mTangType = 2;
	mIter = 0;
	mUseElasticTan = false;

	m_firstLoading = true;

	this->initialize();
}

// destructor
SAniSandMS::~SAniSandMS()
{
}

NDMaterial*
SAniSandMS::getCopy(const char *type)
{
	if (strcmp(type, "PlaneStrain2D") == 0 || strcmp(type, "PlaneStrain") == 0) {
		SAniSandMSPlaneStrain *clone;
		clone = new SAniSandMSPlaneStrain(this->getTag(), m_G0, m_nu, m_e_init, m_Mc,
			m_c, m_lambda_c, m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0,
			m_nd, m_zeta, m_mu0,
			m_beta, massDen,
			mScheme, mTangType, mJacoType, mTolF, mTolR);
		return clone;
	}
	else if (strcmp(type, "ThreeDimensional") == 0 || strcmp(type, "3D") == 0) {
		SAniSandMS3D *clone;
		clone = new SAniSandMS3D(this->getTag(), m_G0, m_nu, m_e_init, m_Mc, m_c, m_lambda_c,
			m_e0, m_ksi, m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_zeta, m_mu0,
			m_beta, massDen,
			mScheme, mTangType, mJacoType, mTolF, mTolR);
		return clone;
	}
	else {
		opserr << "SAniSandMS::getCopy failed to get copy: " << type << endln;
		return 0;
	}
	return 0;
}

int
SAniSandMS::commitState(void)
{
	double D;
	malpha_in_n = malpha_in;

	if ((GetTrace(mSigma) / 3) > (m_P_atm / 5))
	  mUseElasticTan = false;

	// mAlpha_in_n  = mAlpha_in;
	mSigma_n = mSigma;
	mEpsilon_n = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	mAlpha_n = mAlpha;
	mAlphaM_n = mAlphaM;

	mMM_plus_n = mMM_plus;
	mMM_minus_n = mMM_minus;

	

	mDGamma_n = mDGamma;
	mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);

	num_load_reversals_n = num_load_reversals;

	this->GetElasticModuli(mSigma, mVoidRatio, mK, mG, D);

	

	return 0;
}

int SAniSandMS::revertToLastCommit(void)
{

	return 0;
}

int SAniSandMS::revertToStart(void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
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
SAniSandMS::getCopy(void)
{
	opserr << "SAniSandMS::getCopy -- subclass responsibility\n";
	exit(-1);
	return 0;
}

const char*
SAniSandMS::getType(void) const
{
	opserr << "SAniSandMS::getType -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int
SAniSandMS::getOrder(void) const
{
	opserr << "SAniSandMS::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}


Response*
SAniSandMS::setResponse(const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else if (strcmp(argv[0], "alpha") == 0 || strcmp(argv[0], "backstressratio") == 0)
		return new MaterialResponse(this, 4, this->getAlpha());
	else if (strcmp(argv[0], "alphaM") == 0)
		return new MaterialResponse(this, 5, this->getAlphaM());
	else if (strcmp(argv[0], "alpha_in") == 0 || strcmp(argv[0], "rin") == 0)
		return new MaterialResponse(this, 6, this->getalpha_in());
	else if (strcmp(argv[0], "MM") == 0 || strcmp(argv[0], "M") == 0)
		return new MaterialResponse(this, 7, this->getMM());
	else if (strcmp(argv[0], "estrain") == 0 || strcmp(argv[0], "elasticstrain") == 0)
	{
		return new MaterialResponse(this, 8, this->getEStrain());
	}
	else
	{
		opserr << "SAniSandMS::setResponse  --  Unrecognized response option \"" << argv[0] << "\"" << endln;
		return 0;
	}
}

int
SAniSandMS::getResponse(int responseID, Information &matInfo)
{
	// opserr << "responseID = " << responseID << endln;
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
			*(matInfo.theVector) = getAlphaM();
		return 0;
	case 6:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getalpha_in();
		return 0;
	case 7:
		if (matInfo.theDouble != 0)
			matInfo.theDouble = getMM();
		return 0;
	case 8:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getEStrain();
		return 0;
	default:
		return -1;
	}
}

int
SAniSandMS::sendSelf(int commitTag, Channel &theChannel)
{

	opserr << "SAniSandMS::sendSelf - not yet implemented! Contact https://github.com/jaabell" << endln;

	return 0;
}

int
SAniSandMS::recvSelf(int commitTag, Channel &theChannel,
	FEM_ObjectBroker &theBroker)
{
	opserr << "SAniSandMS::recvSelf - not yet implemented! Contact https://github.com/jaabell" << endln;
	return 0;
}

void SAniSandMS::Print(OPS_Stream &s, int flag)
{
	s << "SAniSandMS Material, tag: " << this->getTag() << endln;
	s << "Type: " << this->getType() << endln;
	s << "mSigma_n = " << mSigma_n << endln;
	s << "mEpsilon_n = " << mEpsilon_n << endln;
	s << "mEpsilonE_n = " << mEpsilonE_n << endln;
	s << "mAlpha_n = " << mAlpha_n << endln;
	s << "mAlphaM_n = " << mAlphaM_n << endln;
	s << "mMM_plus_n = " << mMM_plus_n << endln;
	s << "mMM_minus_n = " << mMM_minus_n << endln;
	s << "malpha_in_n = " << malpha_in_n << endln;
	s << "mDGamma_n = " << mDGamma_n << endln;
	s << "mVoidRatio = " << mVoidRatio << endln;


}

int
SAniSandMS::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {    // enforce elastic/elastoplastic response
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0], "materialState") == 0) {    // enforce elastic/elastoplastic response
			return param.addObject(5, this);
		}
		else if (strcmp(argv[0], "IntegrationScheme") == 0) { // change integration scheme (Explicit/Implicit)
			return param.addObject(2, this);
		}
		else if (strcmp(argv[0], "Jacobian") == 0) {         // change type of Jacobian used for newton iterations
			return param.addObject(3, this);
		}
		else if ((strcmp(argv[0], "refShearModulus") == 0) ||
			(strcmp(argv[0], "ShearModulus") == 0)) {  // change G0
			return param.addObject(6, this);
		}
		else if (strcmp(argv[0], "poissonRatio") == 0) {     // change nu
			return param.addObject(7, this);
		}
		else if (strcmp(argv[0], "voidRatio") == 0) {       // change e_init
			return param.addObject(8, this);
		}
	}
	return -1;
}

int
SAniSandMS::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		mElastFlag = info.theInt;
	}
	// called materialState in tcl file
	else if (responseID == 5) {
		mElastFlag = (int)info.theDouble;
	}
	// called update IntegrationScheme
	else if (responseID == 2) {
		mScheme = (int)info.theDouble;
	}
	// called update Jacobian
	else if (responseID == 3) {
		mJacoType = (int)info.theDouble;
	}
	// called update refShearModulus
	else if (responseID == 6) {
		m_G0 = info.theDouble;
	}
	// called update poissonRatio
	else if (responseID == 7) {
		m_nu = info.theDouble;
	}
	// called update voidRatio
	else if (responseID == 8) {
		double eps_v = GetTrace(mEpsilon);
		opserr << "(before) m_e_init = " << m_e_init << endln;
		m_e_init = (info.theDouble + eps_v) / (1 - eps_v);
		opserr << "(after) m_e_init = " << m_e_init << endln;
	}
	else {
		return -1;
	}

	return 0;
}


// Initialize SAniSandMS Material
void
SAniSandMS::initialize()
{
	// set Initial Ce with p = p_atm
	Vector mSig(6);
	mSig(0) = m_P_atm;
	mSig(1) = m_P_atm;
	mSig(2) = m_P_atm;

	// set minimum allowable p
	m_Pmin = 1.0e-4 * m_P_atm;

	// strain and stress terms
	mEpsilon.Zero();
	mEpsilon_n.Zero();
	mSigma.Zero();
	mSigma_n.Zero();

	mSigma(0) = mSigma(1) = mSigma(2) = 0.0001 * m_P_atm;
	mSigma_n(0) = mSigma_n(1) = mSigma_n(2) = 0.0001 * m_P_atm;

	mEpsilonE.Zero();
	mAlpha.Zero();
	mAlpha_n.Zero();
	malpha_in.Zero();
	malpha_in_n.Zero();
	mDGamma = 0.0;
	mVoidRatio = m_e_init;
	mAlphaM.Zero();
	mAlphaM_n.Zero();
	mMM_plus = m_m;
	mMM_plus_n = m_m;
	mMM_minus = 0;
	mMM_minus_n = 0;

	// calculate initial stiffness parameters
	GetElasticModuli(mSig, mVoidRatio, mK, mG);
	mCe = GetStiffness(mK, mG);
	mCep = mCe;
	mCep_Consistent = mCe;

	// calculate machine epsilon (used for FDM Jacobian)
	mEPS = machineEPS();

	mUseElasticTan = false;
}

//send back the state parameters to the recorders
const Vector
SAniSandMS::getState()
{
	Vector result(26);

	double p = one3 * GetTrace(mSigma_n);
	Vector n = GetNormalToYield(mSigma_n, mAlpha);
	Vector r = GetDevPart(mSigma_n) / p;

	result.Assemble(mEpsilonE, 0, 1.0);
	result.Assemble(mAlpha, 6, 1.0);
	result(21) = mMM_plus_n;
	result(22) = mMM_minus_n;
	result(23) = DoubleDot2_2_Contr(mAlpha - malpha_in_n, n);
	result(24) = mVoidRatio;
	result(25) = mDGamma;

	return result;
}

//send back alpha tensor
const Vector
SAniSandMS::getAlpha()
{
	return mAlpha;
}

//send back alphaM tensor
const Vector
SAniSandMS::getAlphaM()
{
	return mAlphaM_n;
}

//send back memory surface size
double
SAniSandMS::getMM()
{
	return mMM_plus_n - fabs(mMM_minus_n);
}

//send back alpha_in tensor
const Vector
SAniSandMS::getalpha_in()
{
	return malpha_in_n;
}


//send back elastic strain tensor
const Vector&
SAniSandMS::getEStrain()
{
	opserr << "SAniSandMS::getEStrain() - Base class being called!" << endln;
	return mEpsilonE;
}


void SAniSandMS::integrate()
{

	Vector trialDirection(6);
	//trialDirection = GetNormalToYield(mSigma_n + mCe * (mEpsilon - mEpsilon_n), mAlpha_n);
	// trialDirection = GetNormalToYield(mSigma, mAlpha_n);
	trialDirection = mCe * (mEpsilon - mEpsilon_n);
	if (DoubleDot2_2_Contr(mAlpha_n - malpha_in_n, trialDirection) < 0.0)
	//if (DoubleDot2_2_Contr(mAlpha - malpha_in_n, trialDirection) < 0.0)
		malpha_in = mAlpha_n;
	else
		malpha_in = malpha_in_n;


	if (mElastFlag == 0) {
		elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha,
			mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent, mAlphaM);
	}
	// ElastoPlastic response
	else {
		// implicit schemes
		if (mScheme == INT_BackwardEuler)
			opserr << "SAniSandMS::integrate() -- Implicit integration not available yet" << endln;
		// explicit schemes
		else
			explicit_integrator(
				mSigma_n,
				mEpsilon_n,
				mEpsilonE_n,
				mAlpha_n,
				mAlphaM_n,
				mMM_plus_n,
				mMM_minus_n,
				malpha_in,
				//malpha_in_n,
				mEpsilon,
				mEpsilonE,
				mSigma,
				mAlpha,
				mAlphaM,
				mMM_plus,
				mMM_minus,
				mDGamma,
				mVoidRatio,
				mG,
				mK,
				mCe,
				mCep,
				mCep_Consistent);

	}
}

void SAniSandMS::elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
	double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent, Vector& NextAlphaM)
{
	Vector dStrain(6);

	// calculate elastic response
	dStrain = NextStrain - CurStrain;
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + dStrain;
	GetElasticModuli(CurStress, NextVoidRatio, K, G);
	aCep_Consistent = aCep = aC = GetStiffness(K, G);
	NextStress = CurStress + DoubleDot4_2(aC, dStrain);

	//update State variables
	if (one3 * GetTrace(NextStress) > m_Pmin)
	{
		NextAlpha = 3.0 * GetDevPart(NextStress) / GetTrace(NextStress);
		NextAlphaM = NextAlpha;
	}
		
	return;
}



void SAniSandMS::explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// function pointer to the integration scheme
	void (SAniSandMS::*exp_int) (
		const Vector&,
		const Vector&,
		const Vector&,
		const Vector&,
		const Vector&,
		const double,
		const double,
		const Vector&,
		const Vector&,
		Vector&,
		Vector&,
		Vector&,
		Vector&,
		double&,
		double&,
		double&,
		double&,
		double&,
		double&,
		Matrix&,
		Matrix&,
		Matrix&);

	exp_int = &SAniSandMS::RungeKutta4;  //Default to RK4

	switch (mScheme) {
	case INT_ForwardEuler:    // Forward Euler
							  // exp_int = &SAniSandMS::ForwardEuler;
		opserr << "SAniSandMS::explicit_integrator() - Forward Euler (does not work)" << endln;
		break;

	case INT_ModifiedEuler:    // Modified Euler with error control
							   // opserr << "SAniSandMS::explicit_integrator() - Modified Euler " << endln;
							   // opserr << "SAniSandMS::ModifiedEuler() " << endln;
		exp_int = &SAniSandMS::ModifiedEuler;
		break;

	case INT_RungeKutta:    // Runge-Kutta-England 45. As in Sloan.
							// opserr << "SAniSandMS::RungeKutta4() " << endln;
		exp_int = &SAniSandMS::RungeKutta4;
		break;

	case INT_MAXSTR_FE:    // Forward Euler constraining maximum strain increment
	case INT_MAXSTR_MFE:    // Modified Euler constraining maximum strain increment
	case INT_MAXSTR_RK:    // Runge-Kutta 4-th order constraining maximum strain increment
		opserr << "SAniSandMS::explicit_integrator() - INT_MAXSTR_RK - Not yet implemented " << endln;
		exit(0);
		exp_int = &SAniSandMS::MaxStrainInc;
		break;

	case INT_MAXENE_FE:    // Forward Euler constraining maximum energy increment
	case INT_MAXENE_MFE:    // Modified Euler constraining maximum energy increment
	case INT_MAXENE_RK:    //  Runge-Kutta 4-th order constraining maximum energy increment
		opserr << "SAniSandMS::explicit_integrator() - MaxEnergyInc - Not yet implemented " << endln;
		exit(0);
		exp_int = &SAniSandMS::MaxEnergyInc;
		break;

	default:
		opserr << "SAniSandMS::explicit_integrator() - Defaulting to ModifiedEuler " << endln;
		exp_int = &SAniSandMS::RungeKutta4;
		break;
	}

	double elasticRatio, p, pn, f, fn;
	Vector dSigma(6), dStrain(6);
	bool   p_tr_pos = true;

	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	dStrain = NextStrain - CurStrain;
	NextElasticStrain = CurElasticStrain + dStrain;
	aC = GetStiffness(K, G);
	dSigma = DoubleDot4_2(aC, dStrain);
	NextStress = CurStress + dSigma;
	f = GetF(NextStress, CurAlpha);
	p = one3 * GetTrace(NextStress);

	if (GetNorm_Contr(dStrain) < small)
	{
		if (debugFlag)
			opserr << "  Non-step, returning\n";

		return;
	}


	if (p < m_Pmin)
		p_tr_pos = false;

	if (p_tr_pos && (f <= mTolF))
	{
		// This is a pure elastic loading/unloading
		if (debugFlag)
		{
			opserr << "pure elastic loading/unloading" << endln;
		}
		NextAlpha = CurAlpha;
		NextAlphaM = CurAlphaM;
		NextMM_plus = CurMM_plus;
		NextMM_minus = CurMM_minus;
		NextDGamma = 0;
		aCep_Consistent = aCep = aC;

		return;

	}
	else {
		fn = GetF(CurStress, CurAlpha);
		pn = one3 * GetTrace(CurStress);
		if (pn < 0)
		{
			if (debugFlag)
				opserr << "Manzari Dafalias (tag = " << this->getTag() << ") : p_n < 0, This should have not happened!" << endln;
			return;
		}

		if (fn > mTolF)
		{
			// This is an illegal stress state! This shouldn't happen.
			if (debugFlag) opserr << "stress state outside the yield surface!" << endln;
			if (debugFlag) opserr << "SAniSandMS : Encountered an illegal stress state! Tag: " << this->getTag() << endln;
			if (debugFlag) opserr << "                  f = " << GetF(CurStress, CurAlpha) << endln;
			(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha,
				NextAlphaM, NextMM_plus, NextMM_minus, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

		}
		else if (fn < -mTolF) {
			// This is a transition from elastic to plastic
			elasticRatio = IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, 0.0, 1.0);
			dSigma = DoubleDot4_2(aC, elasticRatio * (NextStrain - CurStrain));
			(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio * (NextStrain - CurStrain), CurElasticStrain + elasticRatio * (NextStrain - CurStrain),
				CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextAlphaM, NextMM_plus, NextMM_minus, NextDGamma, NextVoidRatio,
				G, K, aC, aCep, aCep_Consistent);

		}
		else if (fabs(fn) < mTolF) {

			if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha), dSigma) / (GetNorm_Contr(dSigma) == 0 ? 1.0 : GetNorm_Contr(dSigma)) > (-sqrt(mTolF))) {
				// This is a pure plastic step
				(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha,
					NextAlphaM, NextMM_plus, NextMM_minus, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
			}
			else {
				// This is an elastic unloding followed by plastic loading
				elasticRatio = IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha);
				dSigma = DoubleDot4_2(aC, elasticRatio * (NextStrain - CurStrain));
				(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio * (NextStrain - CurStrain), CurElasticStrain + elasticRatio * (NextStrain - CurStrain),
					CurAlpha, CurAlphaM, CurMM_plus, CurMM_minus, Curalpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextAlphaM, NextMM_plus, NextMM_minus, NextDGamma, NextVoidRatio,
					G, K, aC, aCep, aCep_Consistent);
			}
		}
	}
}


void SAniSandMS::MaxStrainInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{

}


void SAniSandMS::MaxEnergyInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	
}



void SAniSandMS::ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// double CurVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	// NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	// aC = GetStiffness(K, G);
	// Vector n(6), d(6), b(6), bM(6), R(6), dPStrain(6);
	// double Cos3Theta, h, hM, psi, rBtheta, alphaDtheta, b0, A, B, C, D, Z;
	// GetStateDependent(CurStress, CurAlpha, CurAlphaM, CurMM, Curalpha_in, CurVoidRatio,  n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, alphaDtheta, b0,
	//     A, D, B, C, R, Z);
	// double dVolStrain = GetTrace(NextStrain - CurStrain);
	// Vector dDevStrain = GetDevPart(NextStrain - CurStrain);
	// double p = one3 * GetTrace(CurStress);


	// Vector r(6);
	// if (p > small)
	//     Vector r = GetDevPart(CurStress) / p;

	// double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);

	// double temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n))))
	//     - K*D*DoubleDot2_2_Contr(n,r));

	// // TODO: if temp4 == 0, the whole step is plastic. Take correct steps here.
	// if (fabs(temp4) < small) temp4 = small;

	// NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
	// Vector dSigma   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
	//           (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
	// Vector dAlpha   = Macauley(NextDGamma) * two3 * h * b;
	// // Vector dFabric  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric);
	// dPStrain = NextDGamma * ToCovariant(R);

	// // Change!
	// Vector rM = mAlphaM + root23*n*CurMM;
	// double norm_rM_minus_r = sqrt(GetNorm_Contr(rM - r));

	// Vector nM(6);
	// if(norm_rM_minus_r < small)
	// {
	//     nM = n;
	// }
	// else
	// {
	//     nM = (rM - r) / norm_rM_minus_r;
	// }

	// Vector rB = rBtheta * n;

	// Vector r_tilde = CurAlpha - root23 * m_m * n;
	// Vector rM_tilde = CurAlphaM - root23 * CurMM * n;

	// double x1 = DoubleDot2_2_Contr(nM, rM - r);
	// double x2 = DoubleDot2_2_Contr(nM, r - r_tilde);
	// double x3 = DoubleDot2_2_Contr(nM, rM - rM_tilde);

	// double h_tilde = b0 / DoubleDot2_2_Contr(rM - Curalpha_in, n);
	// // double hM = 0.5*h_tilde + 0.5 / root23 * CurMM / m_zeta * ( 1 - (x1 + x2) / x3 )* Macauley(-D);


	// double dPStrain_v = GetTrace(dPStrain);
	// Vector dAlphaM = root23*Macauley(NextDGamma)*hM*(rB - rM);
	// double dMM = 1/root23 * DoubleDot2_2_Contr(n, dAlphaM) - CurMM / m_zeta * ( 1 - (x1 + x2) / x3 ) * Macauley(-dPStrain_v);

	// Matrix temp1 = 2.0*G*mIIdevMix + K*mIIvol;
	// Vector temp2 = 2.0*G*n - K*DoubleDot2_2_Contr(n,r)*mI1;
	// Vector temp3 = 2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1;

	// aCep = temp1 - MacauleyIndex(NextDGamma) * Dyadic2_2(temp3, temp2) / temp4;
	// aCep_Consistent = aCep;

	// NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
	// NextStress = CurStress + dSigma;
	// NextAlpha  = CurAlpha  + dAlpha;
	// // NextFabric = CurFabric + dFabric;
	// // Change!
	// NextAlphaM  = CurAlphaM  + dAlphaM;
	// NextMM  = CurMM  + dMM;

	return;
}


void SAniSandMS::ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM,
	double& NextMM_plus,
	double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{

	return;
}


void SAniSandMS::RungeKutta4(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// opserr << "rk";
	// opserr << "rk";
	double CurVoidRatio, dVolStrain;
	Vector n(6), d(6), b(6), bM(6), R(6), dDevStrain(6), r(6);
	double Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, B, C, D, p, Kp, Z;

	double bM_distance_trial, bM_distance_tilde_trial;

	double T = 0.0, dT = 1.0, dT_min = 1e-3, TolE = mTolR;

	Vector nStress(6), nAlpha(6), nAlphaM(6), ndPStrain(6), n_trial(6);
	Vector
		dSigma1(6), dSigma2(6), dSigma3(6), dSigma4(6), dSigma5(6), dSigma6(6), dSigma(6),
		dAlpha1(6), dAlpha2(6), dAlpha3(6), dAlpha4(6), dAlpha5(6), dAlpha6(6), dAlpha(6),
		dAlphaM1(6), dAlphaM2(6), dAlphaM3(6), dAlphaM4(6), dAlphaM5(6), dAlphaM6(6), dAlphaM(6),
		dPStrain1(6), dPStrain2(6), dPStrain3(6), dPStrain4(6), dPStrain5(6), dPStrain6(6), dPStrain(6);
	Matrix aCep1(6, 6), aCep2(6, 6), aCep3(6, 6), aCep4(6, 6), aCep5(6, 6), aCep6(6, 6), aCep_thisStep(6, 6), aD(6, 6);

	// double nMM_plus, dMM1_plus, dMM2_plus, dMM3_plus, dMM4_plus, dMM5_plus, dMM6_plus, dMM_plus;
	double nMM_plus, dMM1_plus, dMM2_plus, dMM3_plus, dMM4_plus, dMM5_plus,  dMM_plus;
	// double nMM_minus, dMM1_minus, dMM2_minus, dMM3_minus, dMM4_minus, dMM5_minus, dMM6_minus, dMM_minus;
	double nMM_minus, dMM1_minus, dMM2_minus, dMM3_minus, dMM4_minus, dMM5_minus,  dMM_minus;
	double temp4, q;

	Vector thisSigma(6), thisAlpha(6), thisAlphaM(6), r_alphaM_trial(6), rM_Alphatilde_trial(6), rM_trial(6);
	// double thisMM_plus, thisMM_minus, MM_trial, thisVoidRatio, nVoidRatio;
	double thisMM_plus, thisMM_minus, MM_trial, thisVoidRatio;

	CurVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);

	GetElasticModuli(CurStress, CurVoidRatio, K, G);


	aC = GetStiffness(K, G);
	aD = GetCompliance(K, G);

	NextStress = CurStress;
	NextAlpha = CurAlpha;
	NextAlphaM = CurAlphaM;
	NextMM_plus = CurMM_plus;
	NextMM_minus = CurMM_minus;

	p = one3 * GetTrace(NextStress);
	if (p < m_Pmin)
	{
		if (debugFlag)
			opserr << "Tag = " << this->getTag() << " : I have a problem (p < 0) - This should not happen!!!" << endln;
		NextStress = GetDevPart(NextStress) + m_Pmin * mI1;
		p = one3 * GetTrace(NextStress);
	}
	// Set aCep_Consistent to zero for substepping process
	aCep_Consistent.Zero();

	n_trial = GetNormalToYield(NextStress, NextAlpha);
	
	MM_trial = NextMM_plus + NextMM_minus;
	r_alphaM_trial = NextAlphaM + root23 * (MM_trial - m_m) * n_trial;
	rM_Alphatilde_trial = NextAlphaM - root23 * (MM_trial - m_m) * n_trial;
	
	
	bM_distance_trial = DoubleDot2_2_Contr(r_alphaM_trial - NextAlpha, n_trial);
	bM_distance_tilde_trial = DoubleDot2_2_Contr(NextAlpha - rM_Alphatilde_trial, n_trial);

	if (bM_distance_trial < -1.0e-3 || bM_distance_tilde_trial < -1.0e-3) 
		opserr << "bM_distance_trial = " << bM_distance_trial << "bM_distance_tilde_trial = " << bM_distance_tilde_trial << endln;
	
	if (MM_trial > 50 || MM_trial < m_m - 1.0e-6) opserr << "weired MM = " << MM_trial << endln;
	if (CurVoidRatio > 2 || CurVoidRatio < 0.1) opserr << "weired e = " << CurVoidRatio << endln;


	while (T < 1.0)
	{
		NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + T * (NextStrain - CurStrain));

		dVolStrain = dT * GetTrace(NextStrain - CurStrain);
		dDevStrain = dT * GetDevPart(NextStrain - CurStrain);

		double depsilon_pv = 0;

		// Calc Delta 1
		thisSigma = NextStress;
		thisAlpha = NextAlpha;
		thisAlphaM = NextAlphaM;
		thisMM_plus = NextMM_plus;
		thisMM_minus = NextMM_minus;
		thisVoidRatio = NextVoidRatio;
		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;


		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);


		//r = GetDevPart(NextStress) / p;
		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma1 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha1 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM1 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain1 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain1);
		dMM1_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM1, n);
		dMM1_minus = -Z * Macauley(-depsilon_pv);

		aCep1 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Calc Delta 2
		thisSigma = NextStress + 0.5 * dSigma1;
		thisAlpha = NextAlpha + 0.5 * dAlpha1;
		thisAlphaM = NextAlphaM + 0.5 * dAlphaM1;
		thisMM_plus = NextMM_plus + 0.5 * dMM1_plus;
		thisMM_minus = NextMM_minus + 0.5 * dMM1_minus;

		thisVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + 0.5 * dPStrain1);

		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;

	


		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);

		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma2 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha2 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM2 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain2 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain2);
		dMM2_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM2, n);
		dMM2_minus = -Z * Macauley(-depsilon_pv);

		aCep2 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Calc Delta 3
		thisSigma = NextStress + 0.25 * (dSigma1 + dSigma2);
		thisAlpha = NextAlpha + 0.25 * (dAlpha1 + dAlpha2);
		thisAlphaM = NextAlphaM + 0.25 * (dAlphaM1 + dAlphaM2);
		thisMM_plus = NextMM_plus + 0.25 * (dMM1_plus + dMM2_plus);
		thisMM_minus = NextMM_minus + 0.25 * (dMM1_minus + dMM2_minus);
		thisVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + 0.25 * (dPStrain1 + dPStrain2));

		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;

		

		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);

		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma3 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha3 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM3 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain3 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain3);
		dMM3_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM3, n);
		dMM3_minus = -Z * Macauley(-depsilon_pv);

		aCep3 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Calc Delta 4
		thisSigma = NextStress - dSigma2 + 2 * dSigma3;
		thisAlpha = NextAlpha - dAlpha2 + 2 * dAlpha3;
		thisAlphaM = NextAlphaM - dAlphaM2 + 2 * dAlphaM3;
		thisMM_plus = NextMM_plus - dMM2_plus + 2 * dMM3_plus;
		thisMM_minus = NextMM_minus - dMM2_minus + 2 * dMM3_minus;
		thisVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain - dPStrain2 + 2 * dPStrain3);

		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;

	

		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);

		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma4 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha4 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM4 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain4 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain4);
		dMM4_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM4, n);
		dMM4_minus = -Z * Macauley(-depsilon_pv);

		aCep4 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Calc Delta 5
		thisSigma = NextStress + (7 * dSigma1 + 10 * dSigma2 + dSigma4) / 27;
		thisAlpha = NextAlpha + (7 * dAlpha1 + 10 * dAlpha2 + dAlpha4) / 27;
		thisAlphaM = NextAlphaM + (7 * dAlphaM1 + 10 * dAlphaM2 + dAlphaM4) / 27;
		thisMM_plus = NextMM_plus + (7 * dMM1_plus + 10 * dMM2_plus + dMM4_plus) / 27;
		thisMM_minus = NextMM_minus + (7 * dMM1_minus + 10 * dMM2_minus + dMM4_minus) / 27;

		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;

		

		thisVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + (7 * dPStrain1 + 10 * dPStrain2 + dPStrain4) / 27);

		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);

		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma5 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha5 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM5 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain5 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain5);
		dMM5_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM5, n);
		dMM5_minus = -Z * Macauley(-depsilon_pv);


		aCep5 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Calc Delta 6
		thisSigma = NextStress + (28 * dSigma1 - 125 * dSigma2 + 546 * dSigma3 + 54 * dSigma4 - 378 * dSigma5) / 625;
		thisAlpha = NextAlpha + (28 * dAlpha1 - 125 * dAlpha2 + 546 * dAlpha3 + 54 * dAlpha4 - 378 * dAlpha5) / 625;
		thisAlphaM = NextAlphaM + (28 * dAlphaM1 - 125 * dAlphaM2 + 546 * dAlphaM3 + 54 * dAlphaM4 - 378 * dAlphaM5) / 625;
		thisMM_plus = NextMM_plus + (28 * dMM1_plus - 125 * dMM2_plus + 546 * dMM3_plus + 54 * dMM4_plus - 378 * dMM5_plus) / 625;
		thisMM_minus = NextMM_minus + (28 * dMM1_minus - 125 * dMM2_minus + 546 * dMM3_minus + 54 * dMM4_minus - 378 * dMM5_minus) / 625;

		p = one3 * GetTrace(thisSigma);
		if (p < m_Pmin)
		{
			thisSigma = GetDevPart(thisSigma) + m_Pmin * mI1;
			p = one3 * GetTrace(thisSigma);
		}

		r = GetDevPart(thisSigma) / p;

		

		thisVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + (28 * dPStrain1 - 125 * dPStrain2 + 546 * dPStrain3 + 54 * dPStrain4 - 378 * dPStrain5) / 625);

		GetElasticModuli(thisSigma, thisVoidRatio, K, G);
		GetStateDependent(thisSigma, thisAlpha, thisAlphaM, thisMM_plus, thisMM_minus, Curalpha_in, thisVoidRatio, n, d, b, bM, Cos3Theta, h, hM, psi, rBtheta, rDtheta, b0, A, D, B, C, R, Z);


		Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		temp4 = (Kp + 2.0 * G * (B - C * GetTrace(SingleDot(n, SingleDot(n, n))))
			- K * D * DoubleDot2_2_Contr(n, r));
		if (fabs(temp4) < small) temp4 = small;
		NextDGamma = (2.0 * G * DoubleDot2_2_Mixed(n, dDevStrain) - K * dVolStrain * DoubleDot2_2_Contr(n, r)) / temp4;
		dSigma6 = 2.0 * G * ToContraviant(dDevStrain) + K * dVolStrain * mI1 - Macauley(NextDGamma) *
			(2.0 * G * (B * n - C * (SingleDot(n, n) - 1.0 / 3.0 * mI1)) + K * D * mI1);
		dAlpha6 = Macauley(NextDGamma) * two3 * h * b;
		dAlphaM6 = Macauley(NextDGamma) * two3 * hM * bM;
		dPStrain6 = NextDGamma * ToCovariant(R);
		depsilon_pv = GetTrace(dPStrain6);
		// dMM6_plus = 1 / root23 * DoubleDot2_2_Contr(dAlphaM6, n);
		// dMM6_minus = -Z * Macauley(-depsilon_pv);

		aCep6 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


		// Update
		dSigma = (dSigma1 + 4 * dSigma3 + dSigma4) / 6;
		//dSigma = (14 * dSigma1 + 35 * dSigma4 + 162 * dSigma5 + 125 * dSigma6) / 336;
		dAlpha = (dAlpha1 + 4 * dAlpha3 + dAlpha4) / 6;
		//dAlpha = (14 * dAlpha1 + 35 * dAlpha4 + 162 * dAlpha5 + 125 * dAlpha6) / 336;
		dAlphaM = (dAlphaM1 + 4 * dAlphaM3 + dAlphaM4) / 6;
		//dAlphaM = (14 * dAlphaM1 + 35 * dAlphaM4 + 162 * dAlphaM5 + 125 * dAlphaM6) / 336;
		// dMM =      ( dMM1      + 4*dMM3      + dMM4       ) / 6;

		dMM_plus = (dMM1_plus + 4 * dMM3_plus + dMM4_plus) / 6;
		dMM_minus = (dMM1_minus + 4 * dMM3_minus + dMM4_minus) / 6;
		//dMM_plus = (14 * dMM1_plus + 35 * dMM4_plus + 162 * dMM5_plus + 125 * dMM6_plus) / 336;
		//dMM_minus = (14 * dMM1_minus + 35 * dMM4_minus + 162 * dMM5_minus + 125 * dMM6_minus) / 336;


		dPStrain = (dPStrain1 + 4 * dPStrain3 + dPStrain4) / 6;
		//dPStrain = (14 * dPStrain1 + 35 * dPStrain4 + 162 * dPStrain5 + 125 * dPStrain6) / 336;

		nStress = NextStress + dSigma;
		nAlpha = NextAlpha + dAlpha;
		nAlphaM = NextAlphaM + dAlphaM;
		nMM_plus = NextMM_plus + dMM_plus;
		nMM_minus = NextMM_minus + dMM_minus;
		// nVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + dPStrain);

		n_trial = GetNormalToYield(nStress, nAlpha);

		if (nMM_plus + nMM_minus < m_m)
		{
			nMM_minus = m_m - nMM_plus;
		}

		MM_trial = nMM_plus + nMM_minus;
		r_alphaM_trial = nAlphaM + root23 * (MM_trial - m_m) * n_trial;
		

		bM_distance_trial = DoubleDot2_2_Contr(r_alphaM_trial - nAlpha, n_trial);
		

		if (bM_distance_trial < -1.0e-4)
		{
			nAlphaM = nAlpha - root23 * (MM_trial - m_m) * n_trial;
		}

		rM_Alphatilde_trial = nAlphaM - root23 * (MM_trial - m_m) * n_trial;
		bM_distance_tilde_trial = DoubleDot2_2_Contr(nAlpha - rM_Alphatilde_trial, n_trial);
	
		if (bM_distance_tilde_trial < -1.0e-4)
		{
			nMM_plus += fabs(bM_distance_tilde_trial)/root23;
		}

		
		// Compute error

		p = one3 * GetTrace(nStress);


		if (p < 0)
		{
			if (dT == dT_min)
				return;
			dT = fmax(0.1 * dT, dT_min);
			continue;
		}

		double stressNorm = GetNorm_Contr(NextStress);
		double alphaNorm = GetNorm_Contr(NextAlpha);

		double curStepError1 = GetNorm_Contr(-42 * dSigma1 - 224 * dSigma3 - 21 * dSigma4 + 162 * dSigma5 + 125 * dSigma6) / 336;
		if (stressNorm >= 0.5) { curStepError1 /= (2 * stressNorm); }

		double curStepError2 = GetNorm_Contr(-42 * dAlpha1 - 224 * dAlpha3 - 21 * dAlpha4 + 162 * dAlpha5 + 125 * dAlpha6) / 336;
		if (alphaNorm >= 0.5) { curStepError2 /= (2 * alphaNorm); }

		double curStepError = fmax(curStepError1, curStepError2);

		if (curStepError > TolE)
		{

			//q = fmax(0.8 * sqrt(TolE / curStepError), 0.1);
			q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.1);

			if (dT == dT_min) {

				mUseElasticTan = true;

				NextElasticStrain -= dPStrain;
				NextStress = nStress;
				
				//double eta = sqrt(13.5) * GetNorm_Contr(GetDevPart(NextStress)) / GetTrace(NextStress);
				//if (eta > m_Mc)
				//	NextStress = one3 * GetTrace(NextStress) * mI1 + m_Mc / eta * GetDevPart(NextStress);
				//NextAlpha = CurAlpha + 3.0 * (GetDevPart(NextStress) / GetTrace(NextStress) - GetDevPart(CurStress) / GetTrace(CurStress));
				NextAlpha = nAlpha;
				NextAlphaM = NextAlpha;
				NextMM_plus = CurMM_plus;
				NextMM_minus = m_m - CurMM_plus;
				//NextAlphaM = nAlphaM;
				//NextAlpha = nAlpha;
				//NextMM_plus = nMM_plus;
				//NextMM_minus = nMM_minus;

				T += dT;
			}
			dT = fmax(q * dT, dT_min);
		}
		else {

			if (debugFlag)
				// if (true) 
				opserr << "+++ Successful increment: T = " << T << ", dT = " << dT << "  Error =  " << curStepError << endln;

			NextElasticStrain -= dPStrain;
			NextStress = nStress;
			NextAlpha = nAlpha;
			NextAlphaM = nAlphaM;
			NextMM_plus = nMM_plus;
			NextMM_minus = nMM_minus;

			// Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurAlphaM, CurMM, Curr_in, NextStrain, NextElasticStrain, NextStress,
			//     NextAlpha, NextAlphaM, NextMM, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

			//aCep_thisStep = (14 * aCep1 + 35 * aCep4 + 162 * aCep5 + 125 * aCep6) / 336;
			aCep_thisStep = (aCep1 + 4 * aCep3 + aCep4 ) / 6;
			aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix);

				//q = fmax(0.8 * sqrt(TolE / curStepError), 0.5);
				q = fmin(0.8 * pow(TolE / curStepError, 0.2), 2.0);
				T += dT;
				dT = fmax(q * dT, dT_min);
				dT = fmin(dT, 1 - T);
		}

	}
	return;
}


// int SAniSandMS::BackwardEuler_CPPM(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
//         const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
//         Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
//         double& NextDGamma,    double& NextVoidRatio,  double& G, double& K, Matrix& Ce, Matrix& Cep, Matrix& Cep_Consistent, int implicitLevel)
int SAniSandMS::BackwardEuler_CPPM(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	return -1;
}


double
SAniSandMS::IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha,
	double a0, double a1)
{
	double a = a0;
	double G, K, vR, f, f0, f1;
	Vector dSigma(6), dSigma0(6), dSigma1(6), strainInc(6);

	strainInc = NextStrain - CurStrain;

	vR = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + a0 * strainInc);
	GetElasticModuli(CurStress, vR, K, G);
	dSigma0 = a0 * DoubleDot4_2(GetStiffness(K, G), strainInc);
	f0 = GetF(CurStress + dSigma0, CurAlpha);

	vR = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + a1 * strainInc);
	GetElasticModuli(CurStress, vR, K, G);
	dSigma1 = a1 * DoubleDot4_2(GetStiffness(K, G), strainInc);
	f1 = GetF(CurStress + dSigma1, CurAlpha);

	for (int i = 1; i <= 10; i++)
	{
		a = a1 - f1 * (a1 - a0) / (f1 - f0);
		dSigma = a * DoubleDot4_2(GetStiffness(K, G), strainInc);
		f = GetF(CurStress + dSigma, CurAlpha);
		if (fabs(f) < mTolF)
		{
			if (debugFlag) opserr << "Found alpha in " << i << " steps" << ", alpha = " << a << endln;
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
	return a;
}


double
SAniSandMS::IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha)
{
	double a = 0.0, a0 = 0.0, a1 = 1.0, da;
	double G, K, vR, f;
	int nSub = 20;
	Vector dSigma(6), dSigma0(6), dSigma1(6), strainInc(6);

	strainInc = NextStrain - CurStrain;


	vR = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	GetElasticModuli(CurStress, vR, K, G);
	dSigma = DoubleDot4_2(GetStiffness(K, G), strainInc);

	for (int i = 1; i < nSub; i++)
	{
		da = (a1 - a0) / 2.0;
		a = a1 - da;
		f = GetF(CurStress + a * dSigma, CurAlpha);
		if (f > mTolF)
		{
			a1 = a;
		}
		else if (f < -mTolF) {
			a0 = a;
			break;
		}
		else {
			if (debugFlag)
				opserr << "Found alpha - Unloading" << ", a = " << a << endln;
			return a;
		}

		if (i == nSub) {
			if (debugFlag)
				opserr << "Didn't find alpha! - Unloading" << ", a0 = " << a0 << ", a1 = " << a1 << endln;
			return 0.0;
		}
	}
	if (debugFlag)
		opserr << "Found alpha - Unloading" << ", a0 = " << a0 << ", a1 = " << a1 << endln;
	return IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, a0, a1);
}


// void
// SAniSandMS::Stress_Correction(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
//         const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
//         Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
//         double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
void SAniSandMS::Stress_Correction(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurAlphaM, const double CurMM_plus, const double CurMM_minus, const Vector& Curalpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextAlphaM, double& NextMM_plus, double& NextMM_minus,
	double& NextDGamma, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{

	return;
}


int
SAniSandMS::NewtonIter(const Vector& xo, const Vector& inVar, Vector& x, Matrix& aCepPart)
{
	
	return -1;
}

int
SAniSandMS::NewtonIter2(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart)
{
	
	return -1;
}


Vector
SAniSandMS::NewtonRes(const Vector& x, const Vector& inVar)
{
	

	Vector error(19);//To make compilers happy...
	return error;
}


int
SAniSandMS::NewtonSol(const Vector &xo, const Vector &inVar, Vector& del, Matrix& Cep)
{
	

	return 0;
}











int
SAniSandMS::NewtonIter3(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart)
{
	
	return -1;
}


int
SAniSandMS::NewtonSol2(const Vector &xo, const Vector &inVar, Vector& res, Vector& JRes, Vector& del, Matrix& Cep)
{
	
	return -1;
}




int
SAniSandMS::NewtonSol_negP(const Vector &xo, const Vector &inVar, Vector& del, Matrix& Cep)
{
	
	return -1;
}


Vector
SAniSandMS::NewtonRes_negP(const Vector& x, const Vector& inVar)
{
	

	Vector error(20); //To make compilers happy...
	return error;
}


Vector
SAniSandMS::GetResidual(const Vector& x, const Vector& inVar)
{
	

	Vector error(19);//To make compilers happy...
	return error;
}


Matrix
SAniSandMS::GetJacobian(const Vector &x, const Vector &inVar)
{
	

	Matrix error(19, 19);//To make compilers happy...
	return error;
}



Matrix
SAniSandMS::GetFDMJacobian(const Vector& delta, const Vector& inVar)
{
	

	Matrix error(19, 19);//To make compilers happy...
	return error;
}



Vector
SAniSandMS::SetManzariComponent(const Vector& stress, const Vector& alpha,
	const Vector& fabric, const double& dGamma)
{
	// flush the all data field
	// mSize = 19;
	// Caution: Vector::Assemble() adds the number to current values
	// Vector result(19);
	// result.Assemble(stress, 0);        // Stress
	// result.Assemble(alpha, 6);        // Alpha
	// result.Assemble(fabric, 12);    // Fabric
	// result(18) = dGamma;            // DGamma
	// return result;

	Vector error(19);//To make compilers happy...
	return error;


}


Vector
SAniSandMS::SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, const Vector& cEStrain,
	const Vector& cAlpha, const Vector& cFabric, const double& cVoidRatio, const double& nVoidRatio,
	const Vector& Alpha_in)
{
	// flush the all data field
	// mSize = 44;
	// Caution: Vector::Assemble() adds the number to current values
	// Vector result(44);
	// result.Assemble(nStrain, 0);
	// result.Assemble(cStrain, 6);
	// result.Assemble(cStress, 12);
	// result.Assemble(cEStrain, 18);
	// result.Assemble(cAlpha, 24);
	// result.Assemble(cFabric, 30);
	// result(36) = cVoidRatio;
	// result(37) = nVoidRatio;
	// result.Assemble(Alpha_in, 38);
	// return result;

	Vector error(44);//To make compilers happy...
	return error;
}



double
SAniSandMS::machineEPS()
{
	double eps = 1.0;
	while (((double) 1.0 + eps) > ((double) 1.0))
		eps /= 2.0;
	return eps;
}


double SAniSandMS::Macauley(double x)
{
	// Macauley bracket
	return (x > 0 ? x : 0.0);
}


double SAniSandMS::MacauleyIndex(double x)
{
	// Macauley index
	return (x > 0 ? 1.0 : 0.0);
}


double
SAniSandMS::g(const double cos3theta, const double c)
{
	return 2 * c / ((1 + c) - (1 - c) * cos3theta);
}


double
SAniSandMS::GetF(const Vector& nStress, const Vector& nAlpha)
{
	// Manzari's yield function
	Vector s(6); s = GetDevPart(nStress);
	double p = one3 * GetTrace(nStress);
	s = s - p * nAlpha;
	return GetNorm_Contr(s) - root23 * m_m * p;
}


double
SAniSandMS::GetPSI(const double& e, const double& p)
{
	return e - (m_e0 - m_lambda_c * pow((p / m_P_atm), m_ksi));
}


double
SAniSandMS::GetLodeAngle(const Vector& n)
// Returns cos(3*theta)
{
	double Cos3Theta = sqrt(6.0) * GetTrace(SingleDot(n, SingleDot(n, n)));
	Cos3Theta = Cos3Theta > 1 ? 1 : Cos3Theta;
	Cos3Theta = Cos3Theta < -1 ? -1 : Cos3Theta;
	return Cos3Theta;
}


void
SAniSandMS::GetElasticModuli(const Vector& sigma, const double& en, const double& en1, const Vector& nEStrain,
	const Vector& cEStrain, double &K, double &G)
	// Calculates G, K
{
	double pn = one3 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;

	// this part could make problems
	/*
	//if (fabs(GetTrace(nEStrain - cEStrain)) < small)
	if (fabs(en1 - en) < small)
	{
	G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en) * sqrt(pn / m_P_atm);
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
	} else {
	double ken = pow((2.97 - en),2) / (1+en);
	double ken1= pow((2.97 - en1),2) / (1+en1);
	double pn1 = pow((sqrt(pn) + 0.5* two3 * (1 + m_nu) / (1 - 2 * m_nu) * m_G0 * sqrt(m_P_atm) * (ken1*GetTrace(nEStrain) - ken*GetTrace(cEStrain))),2);
	K = (pn1-pn) / (GetTrace(nEStrain - cEStrain));
	G = 1.5 * (1 - 2 * m_nu) / (1 + m_nu) * K;
	}
	*/
	if (mElastFlag == 0)
	{
		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en);
	}
	else
	{
		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en) * sqrt(pn / m_P_atm);
	}
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


void
SAniSandMS::GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G, const double& D)
// Calculates G, K
{
	double pn = one3 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;

	if (mElastFlag == 0)
	{
		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en);
	}
	else
	{
		// G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init) * sqrt(pn / m_P_atm);

		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en) * sqrt(pn / m_P_atm);
	}
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


void
SAniSandMS::GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G)
// Calculates G, K
{
	double pn = one3 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;

	if (mElastFlag == 0)
		// G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init);
	{
		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en);
	}
	else
	{
		G = m_G0 * m_P_atm * pow((2.97 - en), 2) / (1 + en) * sqrt(pn / m_P_atm);
	}
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


Matrix
SAniSandMS::GetStiffness(const double& K, const double& G)
// returns the stiffness matrix in its contravarinat-contravariant form
{
	Matrix C(6, 6);
	double a = K + 4.0 * one3 * G;
	double b = K - 2.0 * one3 * G;
	C(0, 0) = C(1, 1) = C(2, 2) = a;
	C(3, 3) = C(4, 4) = C(5, 5) = G;
	C(0, 1) = C(0, 2) = C(1, 2) = b;
	C(1, 0) = C(2, 0) = C(2, 1) = b;


	// if(debugFlag){
	if (false) {
		opserr << "     SAniSandMS::GetStiffness() " << endln;
		opserr << "     K = " << K << endln;
		opserr << "     G = " << G << endln;
		opserr << "     C = " << C << endln;
	}

	return C;
}


Matrix
SAniSandMS::GetCompliance(const double& K, const double& G)
// returns the compliance matrix in its covariant-covariant form
{
	Matrix D(6, 6);
	double a = 1 / (9 * K) + 1 / (3 * G);
	double b = 1 / (9 * K) - 1 / (6 * G);
	double c = 1 / G;
	D(0, 0) = D(1, 1) = D(2, 2) = a;
	D(3, 3) = D(4, 4) = D(5, 5) = c;
	D(0, 1) = D(0, 2) = D(1, 2) = b;
	D(1, 0) = D(2, 0) = D(2, 1) = b;
	return D;
}


Matrix
SAniSandMS::GetElastoPlasticTangent(const Vector& NextStress, const double& NextDGamma,
	const Vector& CurStrain, const Vector& NextStrain,
	const double& G, const double& K, const double& B,
	const double& C, const double& D, const double& h,
	const Vector& n, const Vector& d, const Vector& b)
{
	double p = one3 * GetTrace(NextStress);
	p = (p < small) ? small : p;
	Vector r = GetDevPart(NextStress) / p;
	double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);

	Matrix aC(6, 6), aCep(6, 6);
	Vector temp1(6), temp2(6), R(6);
	double temp3;

	aC = GetStiffness(K, G);
	R = ToCovariant((B * n) - (C * (SingleDot(n, n) - one3 * mI1)) + (one3 * D * mI1));
	// R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;

	temp1 = DoubleDot4_2(aC, ToCovariant(R));
	temp2 = DoubleDot2_4(ToCovariant(n - one3 * DoubleDot2_2_Contr(n, r) * mI1), aC);
	temp3 = DoubleDot2_2_Contr(temp2, R) + Kp;
	if (fabs(temp3) < small) return aC;

	aCep = (aC - (MacauleyIndex(NextDGamma) / temp3 * (Dyadic2_2(temp1, temp2))));
	return aCep;
}


Vector
SAniSandMS::GetNormalToYield(const Vector &stress, const Vector &alpha)
{
	static Vector devStress(6);
	static Vector n(6);
	devStress.Zero();
	n.Zero();
	devStress = GetDevPart(stress);

	double p = one3 * GetTrace(stress);

	if (fabs(p) < m_Pmin)
	{
		n.Zero();
	}
	else {
		n = devStress - p * alpha;
		double normN = GetNorm_Contr(n);
		normN = (normN < small) ? small : normN;
		n = n / normN;
	}
	return n;
}


int
SAniSandMS::Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha)
// Check if the solution of implicit integration makes sense
{
	int result = 1;

	if (GetTrace(stress) < 0)
	{
		if (debugFlag)
			opserr << "p < 0 !!!!!! Check this : SAniSandMS::Check()" << endln;
		//result = -2;
	}

	Vector n(6);    n = GetNormalToYield(stress, CurAlpha);
	Vector n_tr(6); n_tr = GetNormalToYield(TrialStress, CurAlpha);

	// check the direction of stress and trial stress
	if (DoubleDot2_2_Contr(n, n_tr) < 0)
	{
		if (debugFlag)
			opserr << "Direction of n and n_tr are more than 90 degrees apart!" << endln;
		result = -4;
	}

	// add any other checks here

	return result;
}


void
SAniSandMS::GetStateDependent(const Vector &stress, const Vector &alpha,
	const Vector &alphaM, const double MM_plus, const double MM_minus, const Vector &alpha_in,
	const double &e, Vector &n, Vector &d, Vector &b, Vector &bM,
	double &cos3Theta, double &h, double &hM, double &psi,
	double &rBtheta, double &rDtheta, double &b0,
	double& A, double& D, double& B, double& C, Vector& R, double &Z)
{

	double D_factor = 1.0;
	double p = one3 * GetTrace(stress);
	Vector stress_use = stress;
	if (p < m_Pmin)
	{
		stress_use = GetDevPart(stress) + m_Pmin * mI1;
		p = one3 * GetTrace(stress_use);
	}

	Vector r = GetDevPart(stress_use) / p;

	n = GetNormalToYield(stress_use, alpha);

	double AlphaAlphaInDotN;
	AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in, n);

	psi = GetPSI(e, p);

	cos3Theta = GetLodeAngle(n);

	rBtheta = g(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
	double rBthetaPLUSpi = g(-cos3Theta, m_c)* m_Mc * exp(-1.0 * m_nb * psi) - m_m;

	Vector alphaBtheta = root23 * rBtheta * n;

	Vector alphaBthetaPLUSpi = -root23 * rBthetaPLUSpi * n;

	rDtheta = g(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;

	Vector alphaDtheta = root23 * rDtheta * n;

	b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);

	d = alphaDtheta - alpha;

	b = alphaBtheta - alpha;

	// the memory surface and image points on it
	double MM = fmax(0.01, MM_plus + MM_minus);
	
	Vector r_alphaM = alphaM + root23 * (MM - m_m) * n;
	Vector rM_tilde = alphaM - root23 * MM * n;
	Vector rM_Alphatilde = alphaM - root23 * (MM - m_m) * n;
	Vector rM = alphaM + root23 * MM * n;

	double bM_distance = DoubleDot2_2_Contr(r_alphaM - alpha, n);
	double bM_distance_tilde = DoubleDot2_2_Contr(alpha - rM_Alphatilde, n);



	if (bM_distance < 0) 
		//|| bM_distance_tilde < -1.0e-4)
	{
		//if (bM_distance < -1.0e-3)
		//{
		//	opserr << "happened...... " << bM_distance << endln;
		//}

		r_alphaM = alpha;
		rM = r;
		bM_distance = 0;
	}


	Vector r_tilde = alpha - root23 * m_m * n;
	//Vector rM_tilde = alphaM - root23 * MM * n;
	double f_shr;
	double x2 = 0;
	double x3;
	// double x1 = DoubleDot2_2_Contr(n, r_alphaM - alpha);
	if (bM_distance_tilde < 0)
	{
		rM_tilde = r_tilde;
		rM_Alphatilde = alpha;
		bM_distance = 2 * root23 * (MM-m_m);
		
	}

	x3 = DoubleDot2_2_Contr(n, rM - rM_tilde);

	x2 = DoubleDot2_2_Contr(n, r_tilde - rM_tilde);
	if (x2 < 0)
	{
		if (x2 < -1.0e-4) opserr << "x2 = " << x2 << endln;
		f_shr = 0;
	}
	else
	{
		f_shr = x2 / x3;
		//f_shr = (1 - (x1 + x2 + 2*root23 *m_m) / (2*root23*MM));
	}

	if (f_shr < 0) {
		//if (f_shr < -1.0e-7)
		//{
		//	opserr << "f_shr = " << f_shr << endln;
		//	opserr << " MM = " << MM << endln;
		//	opserr << " MM - m_m = " << DoubleDot2_2_Contr(n, r_alphaM - rM_Alphatilde) << endln;
		//}
		f_shr = 0;
	}
	if (f_shr > 1) {
		//if (f_shr > 1 +1.0e-7)
		//{
		//	opserr << "f_shr = " << f_shr << endln;
		//	opserr << " MM = " << MM << endln;
		//	opserr << " MM - m_m = " << DoubleDot2_2_Contr(n, r_alphaM - rM_Alphatilde) << endln;
		//}
		f_shr = 1;
	}

	//  Change!
	double gthetaPlusPi = 2 * m_c / ((1 + m_c) + (1 - m_c) * cos3Theta);
	Vector alphad_tilde = root23 * (gthetaPlusPi * m_Mc * exp(m_nd * psi) - m_m) * ((-1) * n);

	double b_dM_tilde = DoubleDot2_2_Contr(alphad_tilde - alpha_in, n);
	double bref = fmax(2 * root23 * m_m, DoubleDot2_2_Contr(alphaBtheta - alphaBthetaPLUSpi, n));
	double bref_D = fmax(2 * root23 * m_m, GetNorm_Contr(alpha_in));

	double b_d_r = DoubleDot2_2_Contr(d, n);

	//if (b_d_r > 0)
	//{
		A = m_A0 * exp(m_beta * Macauley(b_dM_tilde) / bref_D);
	//}
	//else {
	//	A = m_A0;
	//}

	D = A * b_d_r;

	// Apply a factor to D so it doesn't go very big when p is small
	if (p < 0.001 * m_P_atm)
	{
		D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
	}
	else {
		D_factor = 1.0;
	}

	D *= D_factor;

	B = 1.0 + 1.5 * (1 - m_c) / m_c * g(cos3Theta, m_c) * cos3Theta;

	C = 3.0 * sqrt(1.5) * (1 - m_c) / m_c * g(cos3Theta, m_c);

	R = B * n - C * (SingleDot(n, n) - one3 * mI1) + one3 * D * mI1;
	

	//f_shr = (1 - (x1 + x2) / x3);
	
	
	Z = MM / m_zeta * f_shr;

	bM = alphaBtheta - r_alphaM;
	double b_bM = DoubleDot2_2_Contr(bM, n);
	//if (b_bM < 0.00001) b_bM = 0;

	// compute hM
	double b_rM_rin = DoubleDot2_2_Contr(r_alphaM - alpha_in, n);
	if (b_rM_rin < 1.0e-7)
		b_rM_rin = 1.0e-7;
	hM = fmin(1.0e10, 0.5 * b0 / b_rM_rin + 0.5 / root23 * Z * Macauley(-D) / b_bM);

	// compute h
	if (AlphaAlphaInDotN < small)
	{
		AlphaAlphaInDotN = small;
	}
	
	h = fmin(1.0e7, b0 / AlphaAlphaInDotN * exp(m_mu0 * pow(p / m_P_atm, 0.5) * pow(bM_distance / bref, 2)));

	//if (hM > h) hM = h;
	

}





/*************************************************************/
/*************************************************************/
//            SYMMETRIC TENSOR OPERATIONS                    //
/*************************************************************/
/*************************************************************/
// In all the functions below, contravariant means a stress-like tensor
// and covariant means a strain-like tensor

double
SAniSandMS::GetTrace(const Vector& v)
// computes the trace of the input argument
{
	if (v.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::GetTrace requires vector of size(6)!" << endln;

	return (v(0) + v(1) + v(2));
}

Vector
SAniSandMS::GetDevPart(const Vector& aV)
// computes the deviatoric part of the input tensor
{
	if (aV.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::GetDevPart requires vector of size(6)!" << endln;

	static Vector result(6); result.Zero();

	double p = GetTrace(aV);
	result = aV;
	result(0) -= one3 * p;
	result(1) -= one3 * p;
	result(2) -= one3 * p;

	return result;
}

Vector
SAniSandMS::SingleDot(const Vector& v1, const Vector& v2)
// computes v1.v2, v1 and v2 should be both in their "contravariant" form
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! SAniSandMS::SingleDot requires vector of size(6)!" << endln;

	Vector result(6);
	result(0) = v1(0) * v2(0) + v1(3) * v2(3) + v1(5) * v2(5);
	result(1) = v1(3) * v2(3) + v1(1) * v2(1) + v1(4) * v2(4);
	result(2) = v1(5) * v2(5) + v1(4) * v2(4) + v1(2) * v2(2);
	result(3) = 0.5 * (v1(0) * v2(3) + v1(3) * v2(0) + v1(3) * v2(1) + v1(1) * v2(3) + v1(5) * v2(4) + v1(4) * v2(5));
	result(4) = 0.5 * (v1(3) * v2(5) + v1(5) * v2(3) + v1(1) * v2(4) + v1(4) * v2(1) + v1(4) * v2(2) + v1(2) * v2(4));
	result(5) = 0.5 * (v1(0) * v2(5) + v1(5) * v2(0) + v1(3) * v2(4) + v1(4) * v2(3) + v1(5) * v2(2) + v1(2) * v2(5));
	return result;
}

double
SAniSandMS::DoubleDot2_2_Contr(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "contravariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot2_2_Contr requires vector of size(6)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (i > 2) * v1(i) * v2(i);
	}

	return result;
}

double
SAniSandMS::DoubleDot2_2_Cov(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "covariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot2_2_Cov requires vector of size(6)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) - (i > 2) * 0.5 * v1(i) * v2(i);
	}

	return result;
}

double
SAniSandMS::DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot2_2_Mixed requires vector of size(6)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}

double
SAniSandMS::GetNorm_Contr(const Vector& v)
// computes contravariant (stress-like) norm of input 6x1 tensor
{
	if (v.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::GetNorm_Contr requires vector of size(6)!" << endln;

	double result = 0.0;
	result = sqrt(DoubleDot2_2_Contr(v, v));

	return result;
}

double
SAniSandMS::GetNorm_Cov(const Vector& v)
// computes covariant (strain-like) norm of input 6x1 tensor
{
	if (v.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::GetNorm_Cov requires vector of size(6)!" << endln;

	double result = 0.0;
	result = sqrt(DoubleDot2_2_Cov(v, v));

	return result;
}

Matrix
SAniSandMS::Dyadic2_2(const Vector& v1, const Vector& v2)
// computes dyadic product for two vector-storage arguments
// the coordinate form of the result depends on the coordinate form of inputs
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! SAniSandMS::Dyadic2_2 requires vector of size(6)!" << endln;

	Matrix result(6, 6);

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++)
			result(i, j) = v1(i) * v2(j);
	}

	return result;
}

Vector
SAniSandMS::DoubleDot4_2(const Matrix& m1, const Vector& v1)
// computes doubledot product for matrix-vector arguments
// caution: second coordinate of the matrix should be in opposite variant form of vector
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::DoubleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot4_2 requires 6-by-6 matrix " << endln;

	return m1 * v1;
}

Vector
SAniSandMS::DoubleDot2_4(const Vector& v1, const Matrix& m1)
// computes doubledot product for matrix-vector arguments
// caution: first coordinate of the matrix should be in opposite
// variant form of vector
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::DoubleDot2_4 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot2_4 requires 6-by-6 matrix " << endln;

	return  m1 ^ v1;
}

Matrix
SAniSandMS::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
// caution: second coordinate of the first matrix should be in opposite
// variant form of the first coordinate of second matrix
{
	if ((m1.noCols() != 6) || (m1.noRows() != 6) || (m2.noCols() != 6) || (m2.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::DoubleDot4_4 requires 6-by-6 matrices " << endln;

	return m1 * m2;
}

Matrix
SAniSandMS::SingleDot4_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::SingleDot4_2 requires 6-by-6 matrix " << endln;

	Matrix result(6, 6);
	for (int i = 0; i < 6; i++) {
		result(i, 0) = m1(i, 0) * v1(0) + m1(i, 3) * v1(3) + m1(i, 5) * v1(5);
		result(i, 1) = m1(i, 3) * v1(3) + m1(i, 1) * v1(1) + m1(i, 4) * v1(4);
		result(i, 2) = m1(i, 5) * v1(5) + m1(i, 4) * v1(4) + m1(i, 2) * v1(2);
		result(i, 3) = 0.5 * (m1(i, 0) * v1(3) + m1(i, 3) * v1(1) + m1(i, 5) * v1(4)
			+ m1(i, 3) * v1(0) + m1(i, 1) * v1(3) + m1(i, 4) * v1(5));
		result(i, 4) = 0.5 * (m1(i, 3) * v1(5) + m1(i, 1) * v1(4) + m1(i, 4) * v1(2)
			+ m1(i, 5) * v1(3) + m1(i, 4) * v1(1) + m1(i, 2) * v1(4));
		result(i, 5) = 0.5 * (m1(i, 0) * v1(5) + m1(i, 3) * v1(4) + m1(i, 5) * v1(2)
			+ m1(i, 5) * v1(0) + m1(i, 4) * v1(3) + m1(i, 2) * v1(5));
	}

	return result;
}

Matrix
SAniSandMS::SingleDot2_4(const Vector& v1, const Matrix& m1)
// computes singledot product for vector-matrix arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::SingleDot2_4 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::SingleDot2_4 requires 6-by-6 matrix " << endln;

	Matrix result(6, 6);
	for (int i = 0; i < 6; i++) {
		result(0, i) = m1(0, i) * v1(0) + m1(3, i) * v1(3) + m1(5, i) * v1(5);
		result(1, i) = m1(3, i) * v1(3) + m1(1, i) * v1(1) + m1(4, i) * v1(4);
		result(2, i) = m1(5, i) * v1(5) + m1(4, i) * v1(4) + m1(2, i) * v1(2);
		result(3, i) = 0.5 * (m1(0, i) * v1(3) + m1(3, i) * v1(1) + m1(5, i) * v1(4)
			+ m1(3, i) * v1(0) + m1(1, i) * v1(3) + m1(4, i) * v1(5));
		result(4, i) = 0.5 * (m1(3, i) * v1(5) + m1(1, i) * v1(4) + m1(4, i) * v1(2)
			+ m1(5, i) * v1(3) + m1(4, i) * v1(1) + m1(2, i) * v1(4));
		result(5, i) = 0.5 * (m1(0, i) * v1(5) + m1(3, i) * v1(4) + m1(5, i) * v1(2)
			+ m1(5, i) * v1(0) + m1(4, i) * v1(3) + m1(2, i) * v1(5));
	}
	return result;
}

Matrix
SAniSandMS::Trans_SingleDot4T_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::SingleDot4_2 requires 6-by-6 matrix " << endln;
	Matrix result(6, 6);
	for (int i = 0; i < 6; i++) {
		result(0, i) = m1(0, i) * v1(0) + m1(3, i) * v1(3) + m1(5, i) * v1(5);
		result(1, i) = m1(3, i) * v1(3) + m1(1, i) * v1(1) + m1(4, i) * v1(4);
		result(2, i) = m1(5, i) * v1(5) + m1(4, i) * v1(4) + m1(2, i) * v1(2);
		result(3, i) = 0.5 * (m1(0, i) * v1(3) + m1(3, i) * v1(1) + m1(5, i) * v1(4)
			+ m1(3, i) * v1(0) + m1(1, i) * v1(3) + m1(4, i) * v1(5));
		result(4, i) = 0.5 * (m1(3, i) * v1(5) + m1(1, i) * v1(4) + m1(4, i) * v1(2)
			+ m1(5, i) * v1(3) + m1(4, i) * v1(1) + m1(2, i) * v1(4));
		result(5, i) = 0.5 * (m1(0, i) * v1(5) + m1(3, i) * v1(4) + m1(5, i) * v1(2)
			+ m1(5, i) * v1(0) + m1(4, i) * v1(3) + m1(2, i) * v1(5));
	}
	return result;
}
double SAniSandMS::Det(const Vector& aV)
{
	if (aV.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::Det requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	return (aV[0] * aV[1] * aV[2]
		+ 2 * aV[3] * aV[4] * aV[5]
		- aV[0] * aV[5] * aV[5]
		- aV[2] * aV[3] * aV[3]
		- aV[1] * aV[4] * aV[4]);
}

Vector SAniSandMS::Inv(const Vector& aV)
{
	if (aV.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::Inv requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	double det = Det(aV);
	if (det == 0)
	{
		opserr << "\n Error! SAniSandMS::Inv - Singular tensor - return 0 tensor" << endln;
		return aV;
	}
	Vector res(6);
	res(0) = aV(1) * aV(2) - aV(4) * aV(4);
	res(1) = aV(0) * aV(2) - aV(5) * aV(5);
	res(2) = aV(0) * aV(1) - aV(3) * aV(3);
	res(3) = aV(4) * aV(5) - aV(2) * aV(3);
	res(4) = aV(3) * aV(5) - aV(0) * aV(4);
	res(5) = aV(3) * aV(4) - aV(1) * aV(5);
	res = res / det;

	return res;
}

Vector SAniSandMS::ToContraviant(const Vector& v1)
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::ToContraviant requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	Vector res = v1;
	res(3) *= 0.5;
	res(4) *= 0.5;
	res(5) *= 0.5;

	return res;
}

Vector SAniSandMS::ToCovariant(const Vector& v1)
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! SAniSandMS::ToCovariant requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	Vector res = v1;
	res(3) *= 2.0;
	res(4) *= 2.0;
	res(5) *= 2.0;

	return res;
}

Matrix SAniSandMS::ToContraviant(const Matrix& m1)
{
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::ToContraviant requires 6-by-6 matrix " << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	Matrix res = m1;
	for (int ii = 0; ii < 6; ii++)
	{
		res(3, ii) *= 0.5;
		res(4, ii) *= 0.5;
		res(5, ii) *= 0.5;
	}

	return res;
}

Matrix SAniSandMS::ToCovariant(const Matrix& m1)
{
	if ((m1.noCols() != 6) || (m1.noRows() != 6))
		opserr << "\n ERROR! SAniSandMS::ToCovariant requires 6-by-6 matrix " << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	Matrix res = m1;
	for (int ii = 0; ii < 6; ii++)
	{
		res(3, ii) *= 2.0;
		res(4, ii) *= 2.0;
		res(5, ii) *= 2.0;
	}

	return res;
}


//Jose A. Abell. 2021.  "Deo omnis gloria"
