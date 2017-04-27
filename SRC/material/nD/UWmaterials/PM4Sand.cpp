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
//          Nov 2016, University of Washington
                                                                      
// Description: This file contains the implementation for the PM4Sand class.

#include <PM4Sand.h>
#include <MaterialResponse.h>

// #include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <algorithm>
#define fmax std::max
#define fmin std::min
#endif

#define INT_ForwardEuler  2
#define INT_ModifiedEuler 1

const double		PM4Sand::root12         = sqrt(1.0 / 2.0);
const double		PM4Sand::one3			= 1.0/3.0;
const double		PM4Sand::two3			= 2.0/3.0;
const double		PM4Sand::root23			= sqrt(2.0/3.0);
const double		PM4Sand::small			= 1e-10;
const double		PM4Sand::maxStrainInc   = 1e-5;
const bool  		PM4Sand::debugFlag		= true;
const char unsigned	PM4Sand::mMaxSubStep	= 10;
char  unsigned		PM4Sand::mElastFlag		= 1;

Vector 			PM4Sand::mI1(3);
Matrix  		PM4Sand::mIIco(3,3);
Matrix 			PM4Sand::mIIcon(3,3);
Matrix 			PM4Sand::mIImix(3,3);
Matrix 			PM4Sand::mIIvol(3,3);
Matrix 			PM4Sand::mIIdevCon(3,3);
Matrix 			PM4Sand::mIIdevMix(3,3);
Matrix 			PM4Sand::mIIdevCo(3,3);
PM4Sand::initTensors PM4Sand::initTensorOps;

static int numPM4SandMaterials = 0;

void *
OPS_PM4Sand(void)
{
  if (numPM4SandMaterials == 0) {
    numPM4SandMaterials++;
    opserr << "PM4Sand nDmaterial - Written: L.Chen, P.Arduino, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: nDMaterial PM4Sand tag? DR? G0? hpo? rho?" << endln;
    return 0;	
  }
  
  int tag;
  double dData[4];
  double oData[24];
	
  oData[0] = 101.3;    // P_atm
  oData[1] = -1;       // h0
  oData[2] = 0.8;      // emax
  oData[3] = 0.5;      // emin
  oData[4] = 0.5;      // nb
  oData[5] = 0.1;      // nd
  oData[6] = -1;       // Ado
  oData[7] = -1;       // z_max
  oData[8] = 250.0;    // cz
  oData[9] = -1;       // ce
  oData[10] = 33.0;    // phi_cv
  oData[11] = 0.3;     // nu
  oData[12] = 2.0;     // Cgd
  oData[13] = -1;      // Ckaf
  oData[14] = 10.0;    // Q
  oData[15] = 1.5;     // R
  oData[16] = 0.01;    // m
  oData[17] = 0.005;   // crhg
  oData[18] = -1;      // chg
  oData[19] = 0;       // FirstCall

  oData[20] = 1;		// IntScheme
  oData[21] = 0;		// TanType
  oData[22] = 1.0e-9;	// TolF
  oData[23] = 1.0e-10;	// TolR

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial PM4Sand material tag" << endln;
    return 0;
  }

  numData = 4;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial PM4Sand material  with tag: " << tag << endln;
    return 0;
  }

  numData = numArgs - 5;
  if (numData != 0)
	if (OPS_GetDouble(&numData, oData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial PM4Sand material  with tag: " << tag << endln;
		return 0;
	}

	theMaterial = new PM4Sand(tag, ND_TAG_PM4Sand, dData[0] , dData[1] , dData[2] , dData[3], oData[0] , oData[1] , oData[2] ,
					     oData[3] , oData[4] , oData[5] , oData[6] , oData[7], oData[8],
					     oData[9], oData[10], oData[11], oData[12], oData[13], oData[14], 
					     oData[15], oData[16], oData[17], oData[18], (int)oData[19], (int)oData[20], (int)oData[21], oData[22], oData[23]);

  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial PM4Sand material with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
PM4Sand::PM4Sand(int tag, int classTag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
	double emin, double nb, double nd, double Ado, double z_max, double cz,
	double ce, double phi_cv, double nu, double Cgd, double Ckaf, double Q,
	double R, double m, double crhg, double chg, int FirstCall, int integrationScheme, int tangentType,
	double TolF, double TolR) : NDMaterial(tag, classTag),
	mEpsilon(3), 
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_n(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
    mAlpha_in_max(3),
    mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3,3),
	mCep(3,3),
	mCep_Consistent(3,3)
{
	m_Dr = Dr;
	m_G0 = G0;
	m_hpo = hp0;
	massDen = mDen;
	m_P_atm = P_atm;
	m_h0 = (h0 < 0) ? fmax(0.3,(0.25+ m_Dr)/2) : h0;
	m_emax = emax;
	m_emin = emin;
	m_nb = nb;
	m_nd = nd;
	m_Ado = Ado;
	m_z_max = z_max;
	m_cz = cz;
	if (ce > 0)
		m_ce = ce;
	else {
		if (m_Dr > 0.75)
			m_ce = 1.0;
		else if (m_Dr < 0.35)
			m_ce = 5.0;
		else
			m_ce = 5.0 - (m_Dr - 0.35) / 0.4 * 4;
	}
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nu = nu;
	m_Cgd = Cgd;
	m_Ckaf = (Ckaf < 0) ? (5.0 + 220.0 *pow((m_Dr - 0.26), 3)) : Ckaf;
	m_Ckaf = m_Ckaf > 35 ? 35 : m_Ckaf;
	m_Ckaf = m_Ckaf < 4 ? 4 : m_Ckaf;
	m_Q = Q;
	m_R = R;
	m_m = m;
	m_crhg = crhg;
	m_chg = chg;
	m_FirstCall = 0;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;

	m_e_init	= emax - (emax - emin) * m_Dr;
		
	m_isSmallp	= false;

	mOrgTangType	= tangentType;
	mIter			= 0;
	
	initialize();
}

// full constructor
PM4Sand::PM4Sand(int tag, double Dr, double G0, double hp0, double mDen, double P_atm, double h0, double emax,
	double emin, double nb, double nd, double Ado, double z_max, double cz,
	double ce, double phi_cv, double nu, double Cgd, double Ckaf, double Q,
	double R, double m, double crhg, double chg, int FirstCall, int integrationScheme, int tangentType,
	double TolF, double TolR) : NDMaterial(tag, ND_TAG_PM4Sand),
	mEpsilon(3),
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_n(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
	mAlpha_in_max(3),
	mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3)
{
	m_Dr = Dr;
	m_G0 = G0;
	m_hpo = hp0;
	massDen = mDen;
	m_P_atm = P_atm;
	m_h0 = (h0 < 0) ? fmax(0.3, (0.25 + m_Dr) / 2) : h0;
	m_emax = emax;
	m_emin = emin;
	m_nb = nb;
	m_nd = nd;
	m_Ado = Ado;
	m_z_max = z_max;
	m_cz = cz;
	if (ce > 0)
		m_ce = ce;
	else {
		if (m_Dr > 0.75)
			m_ce = 1.0;
		else if (m_Dr < 0.35)
			m_ce = 5.0;
		else
			m_ce = 5.0 - (m_Dr - 0.35) / 0.4 * 4;
	}
	m_Mc = 2 * sin(phi_cv / 180.0 * 3.14159265359);
	m_nu = nu;
	m_Cgd = Cgd;
	m_Ckaf = (Ckaf < 0) ? (5.0 + 220.0 *pow((m_Dr - 0.26), 3)) : Ckaf;
	m_Ckaf = m_Ckaf > 35 ? 35 : m_Ckaf;
	m_Ckaf = m_Ckaf < 4 ? 4 : m_Ckaf;
	m_Q = Q;
	m_R = R;
	m_m = m;
	m_crhg = crhg;
	m_chg = chg;
	m_FirstCall = 0;
	mScheme = integrationScheme;
	mTangType = tangentType;
	mTolF = TolF;
	mTolR = TolR;

	m_e_init = emax - (emax - emin) * m_Dr;

	m_isSmallp = false;

	mOrgTangType	= tangentType;
	mIter = 0;

	initialize();
}

// null constructor
PM4Sand ::PM4Sand() 
    : NDMaterial(),
	mEpsilon(3),
	mEpsilon_n(3),
	mEpsilonE(3),
	mEpsilonE_n(3),
	mSigma(3),
	mSigma_n(3),
	mAlpha(3),
	mAlpha_n(3),
	mAlpha_in(3),
	mAlpha_in_n(3),
	mAlpha_in_p(3),
	mAlpha_in_max(3),
	mAlpha_in_min(3),
	mFabric(3),
	mFabric_n(3),
	mFabric_in(3),
	mCe(3, 3),
	mCep(3, 3),
	mCep_Consistent(3, 3)
{
	m_Dr = 0.0;
	m_G0 = 0.0;
	m_hpo = 0.0;
	massDen = 0.0;
	m_P_atm = 0.0;
	m_h0 = 0.0;
	m_emax = 0.0;
	m_emin = 0.0;
	m_nb = 0.0;
	m_nd = 0.0;
	m_Ado = 0.0;
	m_z_max = 0.0;
	m_cz = 0.0;
	m_ce = 0.0;
	m_Mc = 0.0;
	m_nu = 0.0;
	m_Cgd = 0.0;
	m_Ckaf = 0.0;
	m_Q = 0.0;
	m_R = 0.0;
	m_m = 0.0;
	m_crhg = 0.0;
	m_chg = 0.0;
	m_FirstCall = 0;
	mScheme = 2;
	mTangType = 0;
	mTolF = 1.0e-9;
	mTolR = 1.0e-10;

	m_e_init = 0.0;

	m_isSmallp = false;

	mIter = 0;

	this->initialize();
}

// destructor
PM4Sand::~PM4Sand()
{
}

NDMaterial*
PM4Sand::getCopy(const char *type)
{
	if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		PM4Sand *clone;
		double phi_cv = asin(m_Mc / 2.0) * 180.0 / 3.14159265359;
		clone = new PM4Sand(this->getTag(), m_Dr, m_G0, m_hpo, massDen, m_P_atm, m_h0, m_emax,
			m_emin, m_nb, m_nd, m_Ado, m_z_max, m_cz, m_ce, phi_cv, m_nu, m_Cgd, m_Ckaf, m_Q,
			m_R, m_m, m_crhg, m_chg, m_FirstCall, mScheme, mTangType, mTolF, mTolR);
		return clone;
	}
	else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		opserr << "This is a 2D model and it's not compatible with " << type << endln;
	  	return 0;
	}
	else {
	  	opserr << "PM4Sand::getCopy failed to get copy: " << type << endln;
		return 0;
	}
}

int 
PM4Sand::commitState(void)
{
	Vector n(3), R(3), alphaD(3), dFabric(3);
	double Cka, D, K_p, dr;

	mAlpha_in_n = mAlpha_in;

	// update commited state variables
	// if ((GetTrace(mSigma) / 2) > (m_P_atm / 5))
	// 	mTangType   = mOrgTangType;
	mSigma_n    = mSigma;
	mEpsilon_n  = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	mAlpha_n	= mAlpha;
	mAlpha_in_max = Vector2Max(mAlpha_in_max, mAlpha_n);
	mAlpha_in_min = Vector2Min(mAlpha_in_min, mAlpha_n);
	dFabric = mFabric - mFabric_n;
	// update cumulated fabric
	mzcum = mzcum + sqrt(DoubleDot2_2_Contr(dFabric, dFabric) / 2.0);
	mzpeak = fmax(sqrt(DoubleDot2_2_Contr(mFabric, mFabric) / 2.0), mzpeak);
	mFabric_n	= mFabric;
	mDGamma_n   = mDGamma;
	mVoidRatio  = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);
	dr = (m_emax - mVoidRatio) / (m_emax - m_emin);;

	this->GetElasticModuli(mSigma, mK, mG, mMcur, mzcum);
	GetStateDependent(mSigma, mAlpha, mAlpha_in,  mFabric, mFabric_in, mG, mzcum
		, mzpeak, mpzp, mMcur, mksi, dr, n, D, R, K_p, alphaD, Cka);

	mCe = GetStiffness(mK,mG);
	mCep = GetElastoPlasticTangent(mSigma_n, mCe, R, n, K_p);
	mCep_Consistent = mCe;
	return 0;
}

int PM4Sand::revertToLastCommit (void)
{
	// need to be added
    return 0;
}

int PM4Sand::revertToStart(void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->setStress(mSigma);
	}

    return 0;
}

NDMaterial*
PM4Sand::getCopy (void)
{
	opserr << "PM4Sand::getCopy -- subclass responsibility\n"; 
  	exit(-1);
  	return 0;
}

const char*
PM4Sand::getType (void) const
{
	opserr << "PM4Sand::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}

int
PM4Sand::getOrder (void) const
{
	opserr << "PM4Sand::getOrder -- subclass responsibility\n";
	exit(-1);
	return 0;
}


Response*
PM4Sand::setResponse (const char **argv, int argc, OPS_Stream &output)
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
	else
		return 0;
}

int
PM4Sand::getResponse(int responseID, Information &matInfo)
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
	default:
		return -1;
	}
}

int
PM4Sand::sendSelf(int commitTag, Channel &theChannel)
{
	
	int res = 0;
	static Vector data(86);

	data(0) = this->getTag();

	data(1) = m_Dr;
	data(2) = m_G0;
	data(3) = m_hpo;
	data(4) = massDen;
	data(5) = m_P_atm;
	data(6) = m_Pmin;
	data(7) = m_h0;
	data(8) = m_emax;
	data(9) = m_emin;
	data(10) = m_nb;
	data(11) = m_nd;
	data(12) = m_Ado;
	data(13) = m_z_max;
	data(14) = m_cz;
	data(15) = m_ce;
	data(16) = m_Mc;
	data(17) = m_nu;
	data(18) = m_Cgd;
	data(19) = m_Ckaf;
	data(20) = m_Q;
	data(21) = m_R;
	data(22) = m_m;
	data(23) = m_crhg;
	data(24) = m_chg;
	data(25) = m_FirstCall;

	data(26) = mScheme;
	data(27) = mTangType;
	data(28) = mTolF;
	data(29) = mTolR;
	data(30) = mElastFlag;
	data(31) = mOrgTangType;

	data(32) = mEpsilon(0);		data(35) = mEpsilon_n(0);	data(38) = mSigma(0);	data(41) = mSigma_n(0);
	data(33) = mEpsilon(1);		data(36) = mEpsilon_n(1);	data(39) = mSigma(1);	data(42) = mSigma_n(1);
	data(34) = mEpsilon(2);		data(37) = mEpsilon_n(2);	data(40) = mSigma(2);	data(43) = mSigma_n(2);

	data(44) = mEpsilonE(0);	data(47) = mEpsilonE_n(0);	data(50) = mAlpha(0);	data(53) = mAlpha_n(0);
	data(45) = mEpsilonE(1);	data(48) = mEpsilonE_n(1);	data(51) = mAlpha(1);	data(54) = mAlpha_n(1);
	data(46) = mEpsilonE(2);	data(49) = mEpsilonE_n(2);	data(52) = mAlpha(2);	data(55) = mAlpha_n(2);

	data(56) = mFabric(0);		data(59) = mFabric_n(0);	data(62) = mAlpha_in_n(0);	   data(65) = mAlpha_in_p(0);
	data(57) = mFabric(1);		data(60) = mFabric_n(1);	data(63) = mAlpha_in_n(1);	   data(66) = mAlpha_in_p(1);
	data(58) = mFabric(2);		data(61) = mFabric_n(2);	data(64) = mAlpha_in_n(2);	   data(67) = mAlpha_in_p(2);

	data(68) = mAlpha_in_max(0);      data(71) = mAlpha_in_min(0);		data(74) = mFabric_in(0);
	data(69) = mAlpha_in_max(1);      data(72) = mAlpha_in_min(1);		data(75) = mFabric_in(1);
	data(70) = mAlpha_in_max(2);      data(73) = mAlpha_in_min(2);		data(76) = mFabric_in(2);

	data(77) = mDGamma_n;
	data(78) = mDGamma;
	data(79) = mK;
	data(80) = mG;
	data(81) = m_Pmin;

	data(82) = mVoidRatio;
	data(83) = mzcum;
	data(84) = mzpeak;
	data(85) = mpzp;
	data(86) = mMcur;

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: PM4Sand::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}

	return 0;
}

int 
PM4Sand::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
	int res = 0;
	static Vector data(86);

	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: PM4Sand::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}
	this->setTag((int)data(0));

	m_Dr = data(1);
	m_G0 = data(2);
	m_hpo = data(3);
	massDen = data(4);
	m_P_atm = data(5);
	m_Pmin = data(6);
	m_h0 = data(7);
	m_emax = data(8);
	m_emin = data(9);
	m_nb = data(10);
	m_nd = data(11);
	m_Ado= data(12);
	m_z_max= data(13);
	m_cz= data(14);
	m_ce= data(15);
	m_Mc= data(16);
	m_nu= data(17);
	m_Cgd= data(18);
	m_Ckaf= data(19);
	m_Q= data(20);
	m_R= data(21);
	m_m= data(22);
	m_crhg= data(23);
	m_chg= data(24);
	m_FirstCall = data(25); 

	mScheme = data(26);
	mTangType = data(27);
	mTolF = data(28);
	mTolR = data(29);
	mElastFlag = data(30);
	mOrgTangType = data(31);

	mEpsilon(0) = data(32);		mEpsilon_n(0)= data(35);	mSigma(0) = data(38);	mSigma_n(0) = data(41);
	mEpsilon(1) = data(33);		mEpsilon_n(1) = data(36);	mSigma(1) = data(39);	mSigma_n(1) = data(42);
	mEpsilon(2) = data(34); 	mEpsilon_n(2) = data(37);	mSigma(2) = data(40);	mSigma_n(2) = data(43);

	mEpsilonE(0) = data(44);	mEpsilonE_n(0) = data(47);	mAlpha(0) = data(50);	mAlpha_n(0) = data(53);
	mEpsilonE(1) = data(45);	mEpsilonE_n(1) = data(48);	mAlpha(1) = data(51);	mAlpha_n(1) = data(54);
	mEpsilonE(2) = data(46);	mEpsilonE_n(2) = data(49);	mAlpha(2) = data(52);	mAlpha_n(2) = data(55);

	mFabric(0) = data(56);		mFabric_n(0) = data(59);	mAlpha_in_n(0) = data(62);	   mAlpha_in_p(0) = data(65);
	mFabric(1) = data(57);		mFabric_n(1) = data(60);	mAlpha_in_n(1) = data(63);	   mAlpha_in_p(1) = data(66);
	mFabric(2) = data(58);		mFabric_n(2) = data(61);	mAlpha_in_n(2) = data(64);	   mAlpha_in_p(2) = data(67);

	mAlpha_in_max(0) = data(68);      mAlpha_in_min(0) = data(71);		mFabric_in(0) = data(74);
	mAlpha_in_max(1) = data(69);      mAlpha_in_min(1) = data(72);		mFabric_in(1) = data(75);
	mAlpha_in_max(2) = data(70);      mAlpha_in_min(2) = data(73);		mFabric_in(2) = data(76);

	mDGamma_n = data(77);
	mDGamma = data(78);
	mK = data(79);
	mG = data(80);
	m_Pmin = data(81);

	mVoidRatio = data(82);
	mzcum = data(83);
	mzpeak = data(84);
	mpzp = data(85);
	mMcur = data(86);

	return 0;
}

void PM4Sand::Print(OPS_Stream &s, int flag )
{
	s << "PM4Sand Material, tag: " << this->getTag() << endln;
	s << "Type: " << this->getType() << endln;
}

int
PM4Sand::setParameter(const char **argv, int argc, Parameter &param)
{
	if (argc < 2)
		return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0], "updateMaterialStage") == 0) {     // enforce elastic/elastoplastic response
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
		else if (strcmp(argv[0], "firstCall") == 0) {       // update first call, remove fabric
			return param.addObject(8, this);
		}
		//else if (strcmp(argv[0], "voidRatio") == 0) {        // change e_init
		//	return param.addObject(8, this);
		//}
	}
	return -1;
}

int
PM4Sand::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		mElastFlag = info.theInt;
		if (mElastFlag == 1) {
			// set zmax and Ado when changing from elastic to plastic stage
			setStress(mSigma);
		}
	}
	// called materialState in tcl file
	else if (responseID == 5) {
		mElastFlag = (int)info.theDouble;
		if (mElastFlag == 1) {
			// set zmax and Ado when changing from elastic to plastic stage
			setStress(mSigma);
		}
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
		setStress(mSigma_n);
	}
	// called update voidRatio
	//else if (responseID == 8) {
	//	double eps_v = GetTrace(mEpsilon);
	//	m_e_init = (info.theDouble + eps_v) / (1 - eps_v);
	//}
	else {
		return -1;
	}

	return 0;
}

int
PM4Sand::setStress(const Vector&initStress)
{	
	// set minimum allowable p
    m_Pmin = m_P_atm / 1000;
	if (m_FirstCall == 0) {
		double p = fmax(0.5 * GetTrace(mSigma), m_Pmin);
		if (mSigma(0) < 0.0) {
			//stress tensile
			p = m_P_atm / 20.0;
			mSigma = p * mI1;
		}
		mSigma_n = mSigma;
		m_chg = fmax(p*m_crhg, m_chg);
		mksi = GetKsi(m_Dr, p);
		if (m_z_max < 0) {
			m_z_max = fmin(0.7 * exp(-6.1 * mksi), 20.0);
		}
		if (m_Ado < 0) {
			if (mksi < 0) {
				// dense of critical
				mMb = m_Mc * exp(-1.0 * m_nb * mksi);
				mMd = m_Mc * exp(m_nd * mksi);
				m_Ado = 2.5 * (asin(mMb / 2.0) - asin(m_Mc / 2.0)) / (mMb - mMd);
			}
			else {
				mMb = m_Mc * exp(-1.0 * m_nb / 4.0 * mksi);
				mMd = m_Mc * exp(m_nd * 4.0 * mksi);
				m_Ado = 1.24;
			}
		}
		double Mcut = fmax(mMb, mMd);
		double Mfin = -2.0 / p * sqrt(0.5) * GetNorm_Contr(mSigma - p * mI1);
		if (Mfin > Mcut) {
			Vector r = (mSigma - p * mI1) / p * Mcut / Mfin;
			mSigma = p * mI1 + r * p;
			mAlpha = r * (Mcut - m_m) / Mcut;
		}
		else {
			mAlpha = (mSigma - p * mI1) / p;
		}
		GetElasticModuli(mSigma, mK, mG);
		mCe = mCep = mCep_Consistent = GetStiffness(mK, mG);
		mKp = 100 * mG;
		mAlpha_in = mAlpha;
		mAlpha_in_max = mAlpha;
		mAlpha_in_min = mAlpha;
		mFabric.Zero();
		mFabric_in.Zero();
		mFabric_n.Zero();
		mzcum = 0.0;
		mzpeak = m_z_max / 100000.0;
		mpzp = p / 100.0;
		m_FirstCall = 1;
	}
	return 0;
}

int
PM4Sand::setTrialStrain(const Vector &strain_from_element){
	mEpsilon = -1.0 * strain_from_element;   // -1.0 is for geotechnical sign convention
	integrate();
	return 0;
}

// unused trial strain functions
int
PM4Sand::setTrialStrain(const Vector &v, const Vector &r)
{
	return this->setTrialStrain(v);
}

// Initialize PM4Sand Material
void 
PM4Sand::initialize()
{
	// set Initial Ce with p = p_atm
	Vector mSig(3);
	mSig(0) = m_P_atm;
	mSig(1) = m_P_atm;
	mSig(2) = m_P_atm;
	
	setStress(mSig);
}

//send back the state parameters to the recorders
const Vector 
PM4Sand::getState()
 {
	Vector result(14);
	result.Assemble(mEpsilonE, 0, 1.0);
	result.Assemble(mAlpha, 3, 1.0);
	result.Assemble(mFabric, 6, 1.0);
	result.Assemble(mAlpha_in, 9, 1.0);
	result(12) = mVoidRatio;
	result(13) = mDGamma;

	return result;
 }
//send back alpha tensor
const Vector
PM4Sand::getAlpha()
{
	return mAlpha;
}
//send back fabric tensor
const Vector
PM4Sand::getFabric()
{
	return mFabric;
}
//send back alpha_in tensor
const Vector
PM4Sand::getAlpha_in()
{
return mAlpha_in_n;
}
//send back previous alpha_in tensor
const Vector
PM4Sand::getAlpha_in_p()
{
	return mAlpha_in_p;
}
//send back previous L
const double
PM4Sand::getDGamma()
{
	return mDGamma;
}
/*************************************************************/
const Matrix&
PM4Sand::getTangent() {
	if (mTangType == 0)
		return mCe;
	else if (mTangType == 1)
		return mCep;
	else return mCep_Consistent;
}
/*************************************************************/
const Matrix &
PM4Sand::getInitialTangent() {
	return mCe;
}
/*************************************************************/
const Vector &
PM4Sand::getStress() {
	mSigma_r = -1.0 *mSigma;
	return  mSigma_r;  // -1.0 is for geotechnical sign convention
}
/*************************************************************/
const Vector &
PM4Sand::getStrain() {
	mEpsilon_r =  -1.0 * mEpsilon;   // -1.0 is for geotechnical sign convention
	return mEpsilon_r;
}
/*************************************************************/
const Vector &
PM4Sand::getElasticStrain() {
	mEpsilonE_r = -1.0 * mEpsilonE;   // -1.0 is for geotechnical sign convention
	return mEpsilonE_r;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void PM4Sand::integrate()
{
	mAlpha_in = mAlpha_in_n;

	// Force elastic response
	if (mElastFlag == 0) {
		elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha,
			mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent);
	}
	// ElastoPlastic response
	else {
		// explicit schemes
		explicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in,
			mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG,
			mK, mCe, mCep, mCep_Consistent);
	}

}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Elastic Integrator
/*************************************************************/
void PM4Sand::elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
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

	//update State variables
	if (0.5 * GetTrace(NextStress) > m_Pmin)
		NextAlpha = 2.0 * GetDevPart(NextStress) / GetTrace(NextStress);
	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Explicit Integrator
/*************************************************************/
void PM4Sand::explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
	const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
	Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
	double& NextL, double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
	// function pointer to the integration scheme
	void (PM4Sand::*exp_int) (const Vector&, const Vector&, const Vector&, const Vector&, const Vector&, const Vector&,
		const Vector&, Vector&, Vector&, Vector&, Vector&, double&, double&, double&, double&,
		Matrix&, Matrix&, Matrix&);

	switch (mScheme) {
	case INT_ForwardEuler:	// Forward Euler
		exp_int = &PM4Sand::ForwardEuler;
		break;

	case INT_ModifiedEuler:	// Modified Euler with error control
		exp_int = &PM4Sand::ModifiedEuler;
		break;

	default:
		exp_int = &PM4Sand::ForwardEuler;
		break;
	}

	double elasticRatio, p, pn, f, fn;
	Vector dSigma(3), dStrain(3), n(3);

	NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	dStrain = NextStrain - CurStrain;
	NextElasticStrain = CurElasticStrain + dStrain;
	aC = GetStiffness(K, G);
	dSigma = DoubleDot4_2(aC, dStrain);
	NextStress = CurStress + dSigma;

	f = GetF(NextStress, CurAlpha);
	p = 0.5 * GetTrace(NextStress);

	fn = GetF(CurStress, CurAlpha);
	pn = 0.5 * GetTrace(CurStress);

	if (fn > mTolF)
		if (debugFlag) {
			opserr << "stress state outside the yield surface!" << endln;
		}

	n = GetNormalToYield(NextStress, CurAlpha);
	if (DoubleDot2_2_Contr(CurAlpha - alpha_in, n) < 0.0) {
		// This is stress reverse
		if (DoubleDot2_2_Contr(CurAlpha - mAlpha_in_p, n) < 0.0) {
		// update alpha initial using true back-stress ratio
			mAlpha_in_p = alpha_in;
			mAlpha_in = CurAlpha;
			mFabric_in = CurFabric;
		}
		else {
		// update alpha initial using apparent back-stress ratio
			mAlpha_in_p = alpha_in;
			for (int kk = 0; kk < 3; kk++) {
				if (CurAlpha(kk) > 0.0) {
					mAlpha_in(kk) = fmax(mAlpha_in_min(kk), 0.0);
				}
				else {
					mAlpha_in(kk) = fmin(mAlpha_in_max(kk), 0.0);
				}
			}
			mFabric_in = CurFabric;
		}
}

	if (f <= mTolF)
	{
		// This is a pure elastic loading/unloading
		NextAlpha		= CurAlpha;
		NextFabric		= CurFabric;
		NextL		= 0;
		aCep_Consistent = aCep = aC;
		// Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, 
			// NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, G, K , aC, aCep, aCep_Consistent);

		return;

	} else if (fn < -mTolF) {
		// This is a transition from elastic to plastic
		elasticRatio = IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, 0.0, 1.0);
		dSigma		 = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
		(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
			CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, 
			G, K, aC, aCep, aCep_Consistent);

		return;
	} else if (fabs(fn) < mTolF){
		if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha),dSigma)/(GetNorm_Contr(dSigma) == 0 ? 1.0 : GetNorm_Contr(dSigma)) > (- sqrt(mTolF))) {
			// This is a pure plastic step
			(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, 
				NextFabric, NextL, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

			return;
		} else {
			// This is an elastic unloding followed by plastic loading
			elasticRatio = IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha);
			dSigma		 = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
			(this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
				CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextL, NextVoidRatio, 
				G, K, aC, aCep, aCep_Consistent);

			return;
		}
	} else {
		// This is an illegal stress state! This shouldn't happen.
		if (debugFlag) opserr << "PM4Sand : Encountered an illegal stress state! Tag: " << this->getTag() << endln;
		if (debugFlag) opserr << "                  f = " << GetF(CurStress, CurAlpha) << endln;
		(this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, 
				NextFabric, NextL, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
		return;
	}
}

// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Forward-Euler Integrator
/*************************************************************/
void PM4Sand::ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextL, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{	
	double CurVoidRatio, CurDr, Cka, D, K_p, ksi;
	
	CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	CurDr = (m_emax - CurVoidRatio) / (m_emax - m_emin);
	double p = 0.5 * GetTrace(CurStress);
	p = p < m_Pmin ? m_Pmin : p;
	NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);

	// GetElasticModuli(CurStress,K,G,mzcum);
	// aC = GetStiffness(mK, mG);
	Vector n(3), R(3), alphaD(3), dPStrain(3);
	GetStateDependent(CurStress,  CurAlpha, alpha_in, CurFabric, mFabric_in, mG, mzcum
		, mzpeak, mpzp, mMcur, ksi, CurDr, n, D, R, K_p, alphaD, Cka);
	double dVolStrain = GetTrace(NextStrain - CurStrain);
	Vector dDevStrain = (NextStrain - CurStrain) - dVolStrain / 3 * mI1;
	Vector r = GetDevPart(CurStress) / p;
	
	double temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
	if (fabs(temp4) < small) temp4 = small;

	NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain)/temp4;
	Vector dSigma   = 2.0*mG*mIIcon*dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
			  (2.0 * mG * n + mK * D * mI1);
	// update fabric
	if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
		Vector dFabric = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
		// Update pzp, the mean stress at which max z were formed
		// 	if DoubleDotT2_T2_Contr((z + dz), GetTrace(NextStress)*Identity(2)) > DoubleDotT2_T2_Contr(z, GetTrace(CurStress)*Identity(2))
		// 		pzp = 0.5 * GetTrace(NextStress);
		// end
		NextFabric = CurFabric + dFabric;
	}
	// update alpha
	Vector ds = dSigma - 0.5 * GetTrace(dSigma) * mI1;
	Vector s = CurStress - 0.5 * GetTrace(CurStress) * mI1;
	double dp = 0.5*GetTrace(dSigma);
	Vector dAlpha = (s + ds) / (p + dp) - s / p;

	dPStrain = NextL * mIIco * R;

	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
	NextStress = CurStress + dSigma;
	NextAlpha  = CurAlpha  + dAlpha;
	Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, R, n, r);
	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Modified-Euler Integrator
/*************************************************************/
void PM4Sand::ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
		double& NextL, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{	
	double NextDr, dVolStrain, p, dp, Cka, D, K_p, ksi, temp4, curStepError, q;
	Vector n(3), R1(3), R2(3), alphaD(3), dDevStrain(3), r(3);
	Vector nStress(3), nAlpha(3), nFabric(3);
	Vector ds(3), s(3), dSigma1(3), dSigma2(3), dAlpha(3), dFabric1(3), dFabric2(3), dPStrain1(3), dPStrain2(3);
	double T = 0.0, dT = 1.0, dT_min = 1e-3, TolE = 1e-4;

	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextStress = CurStress;
	NextAlpha = CurAlpha;
	NextFabric = CurFabric;

	while (T < 1.0)
	{
		NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + T*(NextStrain - CurStrain));
		NextDr = (m_emax - NextVoidRatio) / (m_emax - m_emin);
		dVolStrain = dT * GetTrace(NextStrain - CurStrain);
		dDevStrain = dT * (NextStrain - CurStrain) - dVolStrain / 3 * mI1;

		// Calc Delta 1
		GetStateDependent(NextStress, NextAlpha, alpha_in, NextFabric, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, ksi, NextDr, n, D, R1, K_p, alphaD, Cka);
		p = 0.5 * GetTrace(NextStress);
		p = p < m_Pmin ? m_Pmin : p;
		r = GetDevPart(NextStress) / p;

		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) temp4 = small;

		NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
		dSigma1 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
			(2.0 * mG * n + mK * D * mI1);
		// update fabric
		if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
			dFabric1 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric);
		}
		dPStrain1 = NextL * mIIco * R1;
		dp = 0.5 * GetTrace(dSigma1);
		ds = dSigma1 - dp * mI1;
		s = NextStress - p * mI1;
		dAlpha = (s + ds) / (p + dp) - (s / p);

		//Calc Delta 2
		GetStateDependent(NextStress + dSigma1, NextAlpha + dAlpha, alpha_in, NextFabric + dFabric1, mFabric_in, mG, mzcum
			, mzpeak, mpzp, mMcur, ksi, NextDr, n, D, R2, K_p, alphaD, Cka);
		p = 0.5 * GetTrace(NextStress + dSigma1);
		p = p < m_Pmin ? m_Pmin : p;
		r = GetDevPart(NextStress + dSigma1) / p;
		
		temp4 = K_p + 2 * mG - mK* D *DoubleDot2_2_Contr(n, r);
		if (fabs(temp4) < small) temp4 = small;
		
		NextL = (2 * mG * DoubleDot2_2_Mixed(n, dDevStrain) - DoubleDot2_2_Contr(n, r) * mK * dVolStrain) / temp4;
		dSigma2 = 2.0 * mG * mIIcon * dDevStrain + mK*dVolStrain*mI1 - Macauley(NextL)*
			(2.0 * mG * n + mK * D * mI1);
		// update fabric
		if (DoubleDot2_2_Contr(alphaD - CurAlpha, n) < 0.0) {
			dFabric2 = -1.0 * m_cz / (1 + Macauley(mzcum / 2.0 / m_z_max - 1.0)) * Macauley(NextL)*MacauleyIndex(-D)*(m_z_max * n + CurFabric + dFabric1);
		}
		dPStrain2 = NextL * mIIco * R2;

		nStress = NextStress + 0.5 * (dSigma1 + dSigma2);
		nFabric = NextFabric + 0.5 * (dFabric1 + dFabric2);
		// update alpha
		p = 0.5 * GetTrace(NextStress);
		s = NextStress - p * mI1;
		dp = 0.5 * GetTrace(0.5 * (dSigma1 + dSigma2));
		ds = 0.5 * (dSigma1 + dSigma2) - dp * mI1;
		if (p == 0 ){
			dAlpha = (s + ds) / (p + dp);
		} else{
			dAlpha = (s + ds) / (p + dp) - s / p;}
		nAlpha = NextAlpha + dAlpha;
		
		curStepError = GetNorm_Contr(dSigma2 - dSigma1) / (2 * GetNorm_Contr(nStress));
		if (curStepError > TolE) {
			q = fmax(0.9*sqrt(TolE / curStepError), 0.1);
			if (dT == dT_min) {

				NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
				NextStress = nStress;
				NextAlpha = nAlpha;
				NextFabric = nFabric;

				Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
				T += dT;
			}
			dT = fmax(q * dT, dT_min);
			//dT = fmin(dT, 1 - T);
		}
		else {
			NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
			NextStress = nStress;
			NextAlpha = nAlpha;
			NextFabric = nFabric;
			Stress_Correction(NextStress, NextAlpha, dAlpha, m_m, 0.5 * (R1 + R2), n, r);
			q = fmin(0.9*sqrt(TolE / curStepError), 1.1);
			dT = q*dT;
			T += dT;
			dT = fmax(dT, dT_min);
			dT = fmin(dT, 1 - T);
		}
	}
	return;
	
}

/*************************************************************/
//            Pegasus Iterations                             //
/*************************************************************/
double
PM4Sand::IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha,
	double a0 = 0.0, double a1 = 1.0)
{
	double a = a0;
	double f, f0, f1;
	Vector dSigma(3), dSigma0(3), dSigma1(3), strainInc(3); 

	strainInc = NextStrain - CurStrain;

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
/*************************************************************/
//      Pegasus Iterations  (ElastoPlastic Unloading)        //
/*************************************************************/
double
PM4Sand::IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha)
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
PM4Sand::Stress_Correction(Vector& NextStress, Vector& NextAlpha, const Vector& dAlpha, 
	const double m, const Vector& R, const Vector& n, const Vector& r)
{
	double f = GetF(NextStress,NextAlpha);
	if (f < 1e-9)
		return;
	else {
		Vector dfrOverdSigma = n - 0.5 * DoubleDot2_2_Contr(n, r)*mI1;
		double lamda = f / DoubleDot2_2_Contr(dfrOverdSigma, R);
		NextStress = NextStress - R*lamda;
		NextAlpha = NextAlpha - dAlpha*lamda;
	}
}

/************************************************************/
/************************************************************/


/*************************************************************/
/*************************************************************/
//            MATERIAL SPECIFIC METHODS                      //
/*************************************************************/
/*************************************************************/
// Macauley() -------------------------------------------------
double PM4Sand::Macauley(double x)
{
	// Macauley bracket
	return (x > 0 ? x : 0.0);
}
/*************************************************************/
// MacauleyIndex() --------------------------------------------
double PM4Sand::MacauleyIndex(double x)
{
	// Macauley index
	return (x > 0 ? 1.0 : 0.0);
}
/*************************************************************/
// GetF() -----------------------------------------------------
double 
PM4Sand::GetF(const Vector& nStress, const Vector& nAlpha)
{
	// PM4Sand's yield function
	Vector s(3); s = GetDevPart(nStress);
	double p = 0.5 * GetTrace(nStress);
	s = s - p * nAlpha;
	return GetNorm_Contr(s) - root12 * m_m * p;
}
/*************************************************************/
// GetPSI() ---------------------------------------------------
double 
PM4Sand::GetKsi(const double& dr, const double& p)
{
	//Bolton
	return m_R / (m_Q - log(100.0 * p / m_P_atm)) - dr;
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
PM4Sand::GetElasticModuli(const Vector& sigma, double &K, double &G, double &Mcur, const double& zcum)
// Calculates G, K, including effects of fabric and current stress ratio
{
	int msr = 4;
	double Csr0 = 0.5;
	double pn = 0.5 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;
	double qn = 2 * sqrt(pow((0.5*(sigma(0) - sigma(1))), 2) + pow(sigma(2), 2));
	Mcur = qn / pn;
	double Csr = 1 - Csr0 *pow((Mcur / mMb), msr);
	double temp = zcum / m_z_max;
	if (mElastFlag == 0)
		G = m_G0 * m_P_atm;
	else
		G = m_G0 * m_P_atm * sqrt(pn / m_P_atm) * Csr * (1+ temp) / (1 + temp * m_Cgd);
	m_nu = (m_nu == 0.5) ? 0.4999 : m_nu;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}
void
PM4Sand::GetElasticModuli(const Vector& sigma, double &K, double &G)
// Calculates G, K
{
	double pn = 0.5 * GetTrace(sigma);
	pn = (pn <= m_Pmin) ? m_Pmin : pn;

	if (mElastFlag == 0) 
		G = m_G0 * m_P_atm;
	else
		G = m_G0 * m_P_atm * sqrt(pn / m_P_atm);
	m_nu = (m_nu == 0.5) ? 0.4999 : m_nu;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}
/*************************************************************/
// GetStiffness() ---------------------------------------------
Matrix
PM4Sand::GetStiffness(const double& K, const double& G)
// returns the stiffness matrix in its contravarinat-contravariant form
{
	Matrix C(3,3);
	double a = K + 4.0*one3 * G;
	double b = K - 2.0*one3 * G;
	C(0,0) = C(1,1) = a;
	C(2,2) = G;
	C(0,1) = C(1,0) = b;
	return C;
}
/*************************************************************/
// GetCompliance() ---------------------------------------------
Matrix
PM4Sand::GetCompliance(const double& K, const double& G)
// returns the compliance matrix in its covariant-covariant form
{
	Matrix D(3,3);
	double a = (K + 4.0 / 3.0 * G) / (4.0 * G * K + 4.0 / 3.0 * pow(G,2));
	double b = (K - 2.0 / 3.0 * G) / (4.0 * G * K + 4.0 / 3.0 * pow(G, 2));
	double c = 1 / G;
	D(0,0) = D(1,1) = a;
	D(2,2) = c;
	D(0,1) = D(1,0) = b;
	return D;
}
/*************************************************************/
// GetElastoPlasticTangent()---------------------------------------
Matrix
PM4Sand::GetElastoPlasticTangent(const Vector& NextStress, const Matrix& aCe, const Vector& R, 
					const Vector& n, const double K_p) 
{	
	double p = 0.5 * GetTrace(NextStress);
	if (p < m_Pmin) p = m_Pmin;
	Vector r = GetDevPart(NextStress) / p;
	Matrix aCep(3,3);
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
PM4Sand::GetNormalToYield(const Vector &stress, const Vector &alpha)
{
	Vector devStress(3); devStress = GetDevPart(stress);
	double p = 0.5 * GetTrace(stress);
	Vector n(3); n = devStress - p * alpha;
	if (GetNorm_Contr(n) != 0)
		n = n / GetNorm_Contr(n);

	return n;
}
/*************************************************************/
// Check() ---------------------------------------------------
int
PM4Sand::Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha)
// Check if the solution of implicit integration makes sense
{
	return 0;
}
/*************************************************************/
// GetStateDependent() ----------------------------------------
void 
PM4Sand::GetStateDependent(const Vector &stress, const Vector &alpha, const Vector &alpha_in
	,  const Vector &fabric, const Vector &fabric_in, const double &G, const double &zcum, const double &zpeak
	, const double &pzp, const double &Mcur, double &ksi, const double &CurDr, Vector &n, double &D, Vector &R, double &K_p
	, Vector &alphaD, double &Cka)
{
	double p = 0.5 * GetTrace(stress);
	if (p <= m_Pmin) p = m_Pmin;
	ksi = GetKsi(CurDr, p);
	n = GetNormalToYield(stress, alpha);
	if (ksi <= 0.0) {
		// dense of critical
		mMb = m_Mc * exp(-1.0 * m_nb * ksi);
		mMd = m_Mc * exp(m_nd * ksi);
	}
	else {
		mMb = m_Mc * exp(-1.0 * m_nb / 4.0 * ksi);
		mMd = m_Mc * exp(m_nd * 4.0 * ksi);
	}
	Vector alphaB = root12 * (mMb - m_m) * n;
	alphaD = root12 * (mMd - m_m) * n;
	double Czpk1 = zpeak / (zcum + m_z_max / 5.0);
	double Czpk2 = zpeak / (zcum + m_z_max / 100.0);
	double Cpzp2 = 1.0;
	//double Cpzp2 = Macauley((pzp-p))/(Macauley((pzp-p))+m_Pmin);
	double Cg1 = m_h0 / 200.0;
	double Ckp = 2.0;

	Cka = 1.0 + m_Ckaf / (1.0 + pow(2.5*Macauley(DoubleDot2_2_Contr(alpha - alpha_in, n)), 2.0))*Cpzp2*Czpk1;
	K_p = G * m_h0 * sqrt(abs(DoubleDot2_2_Contr((alphaB - alpha), n))) /
		(exp(DoubleDot2_2_Contr(alpha - alpha_in, n)) - 1 + Cg1) *
		Cka / (1 + Ckp * zpeak / m_z_max * Macauley(DoubleDot2_2_Contr(alphaB - alpha, n)) * sqrt(1 - Czpk2));
	// loose of critical state
	if (DoubleDot2_2_Contr(alphaB - alpha, n) < 0.0) {
		K_p = -K_p;
	}
	// bound K_p to non - negative, following flac practice
	K_p = fmax(0.0, K_p);

	// dilation
	if (DoubleDot2_2_Contr(alphaD - alpha, n) <= 0) {
		double C_epsilon = m_ce;
		//double Cpzp = (pzp == 0.0) ? 1.0 : 1.0 / (1.0 + pow((2.5*p / pzp), 5.0));
		double Cpzp = 1.0;
		double Cpmin = 1.0 / (1.0 + pow((5.0 * m_Pmin / p), 3));
		double Czin1 = 1.0 - exp(-2.0*abs((DoubleDot2_2_Contr(fabric_in, n) - DoubleDot2_2_Contr(fabric, n)) / m_z_max));
		double Czin2 = (1 + Czin1*(zcum - zpeak) / 3.0 / m_z_max) / (1 + 3.0 * Czin1*(zcum - zpeak) / 3 / m_z_max);
		double temp = pow((1.0 - Macauley(DoubleDot2_2_Contr(-1.0 * fabric, n)) * root12 / zpeak), 3);
		double Ad = m_Ado * Czin2 / ((pow(zcum, 2) / m_z_max)*temp*C_epsilon*Cpzp*Cpmin*Czin1 + 1.0);
		D = Ad * DoubleDot2_2_Contr(alphaD - alpha, n);
		// rotated dilatancy surface
		double Crot1 = fmax(1.0 + 2 * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max)*(1 - Czin1), 1.0);
		double Mdr = mMd / Crot1;
		Vector alphaDr = root12 * (Mdr - m_m) * n;
		double Cdr = fmin(5 + 25 * Macauley(m_Dr - 0.35), 10.0);
		double Drot = Ad * Macauley(DoubleDot2_2_Contr(-1.0*fabric, n)) / (sqrt(2.0)*m_z_max) * DoubleDot2_2_Contr(alphaDr - alpha, n) / Cdr;
		if (D > Drot) {
			D = D + (Drot - D)*Macauley(mMb - Mcur) / Macauley(mMb - Mcur + 0.01);
		}
		if (ksi > 0) {
			D = 0;
		}
	}
	else {
		//contraction
		double hp;
		if (ksi <= 0.5) {
			// change exponent from 2.5 to 2 to improve contraction rate
			hp = m_hpo * exp(-0.7 + 7.0 * pow((0.5 - ksi), 2.5));
		}
		else {
			hp = m_hpo * exp(-0.7);
		}
		double Crot2 = 1 - Czpk2;
		double Cdz = fmax((1 - Crot2*sqrt(2.0)*zpeak / m_z_max)*(m_z_max / (m_z_max + Crot2*zcum)), 1 / (1 + m_z_max / 2.0));
		double Adc = m_Ado * (1 + Macauley(DoubleDot2_2_Contr(fabric, n))) / hp / Cdz;
		double Cin = 2.0 * DoubleDot2_2_Contr(fabric, n) / sqrt(2.0) / m_z_max;
		D = fmin(Adc * pow((DoubleDot2_2_Contr(alpha - alpha_in, n) + Cin),2), 1.5 * m_Ado) *
			DoubleDot2_2_Contr(alphaD - alpha, n) / (DoubleDot2_2_Contr(alphaD - alpha, n) + 0.1);
		//    D = D / 1.66;
		if (p <= m_Pmin) D = D * 1e-4;
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
PM4Sand::GetTrace(const Vector& v) 
// computes the trace of the input argument
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Sand::GetTrace requires vector of size(3)!" << endln;

	return (v(0) + v(1));
}
/*************************************************************/
//  GetDevPart() ---------------------------------------------
Vector 
PM4Sand::GetDevPart(const Vector& aV)
// computes the deviatoric part of the input tensor
{
	if (aV.Size() != 3)
		opserr << "\n ERROR! PM4Sand::GetDevPart requires vector of size(3)!" << endln;

	Vector result(3);
	double p = GetTrace(aV);
	result = aV;
	result(0) -= 0.5 * p;
	result(1) -= 0.5 * p;

	return result;
}
/*************************************************************/
//  SingleDot() ---------------------------------------------
Vector 
PM4Sand::SingleDot(const Vector& v1, const Vector& v2)
// computes v1.v2, v1 and v2 should be both in their "contravariant" form
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! PM4Sand::SingleDot requires vector of size(6)!" << endln;

    Vector result(6);
	result(0) = v1(0)*v2(0) + v1(3)*v2(3) + v1(5)*v2(5);
	result(1) = v1(3)*v2(3) + v1(1)*v2(1) + v1(4)*v2(4);
	result(2) = v1(5)*v2(5) + v1(4)*v2(4) + v1(2)*v2(2);
	result(3) = 0.5*(v1(0)*v2(3) + v1(3)*v2(0) + v1(3)*v2(1) + v1(1)*v2(3) + v1(5)*v2(4) + v1(4)*v2(5));
	result(4) = 0.5*(v1(3)*v2(5) + v1(5)*v2(3) + v1(1)*v2(4) + v1(4)*v2(1) + v1(4)*v2(2) + v1(2)*v2(4));
	result(5) = 0.5*(v1(0)*v2(5) + v1(5)*v2(0) + v1(3)*v2(4) + v1(4)*v2(3) + v1(5)*v2(2) + v1(2)*v2(5));
	return result;
}
/*************************************************************/
// DoubleDot2_2_Contr() ---------------------------------------
double
PM4Sand::DoubleDot2_2_Contr(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "contravariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Sand::DoubleDot2_2_Contr requires vector of size(3)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (i>1) * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Cov() ---------------------------------------
double
PM4Sand::DoubleDot2_2_Cov(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "covariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Sand::DoubleDot2_2_Cov requires vector of size(3)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) - (i>1) * 0.5 * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Mixed() ---------------------------------------
double
PM4Sand::DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Sand::DoubleDot2_2_Mixed requires vector of size(3)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// GetNorm_Contr() ---------------------------------------------
double
PM4Sand::GetNorm_Contr(const Vector& v)
// computes contravariant (stress-like) norm of input 6x1 tensor
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Sand::GetNorm_Contr requires vector of size(3)!" << endln;

	double result=0.0;	
	result = sqrt(DoubleDot2_2_Contr(v,v));

	return result;
}
/*************************************************************/
// GetNorm_Cov() ---------------------------------------------
double
PM4Sand::GetNorm_Cov(const Vector& v)
// computes covariant (strain-like) norm of input 6x1 tensor
{
	if (v.Size() != 3)
		opserr << "\n ERROR! PM4Sand::GetNorm_Cov requires vector of size(3)!" << endln;
	
	double result=0.0;	
	result = sqrt(DoubleDot2_2_Cov(v,v));

	return result;
}
/*************************************************************/
// Dyadic2_2() ---------------------------------------------
Matrix 
PM4Sand::Dyadic2_2(const Vector& v1, const Vector& v2)
// computes dyadic product for two vector-storage arguments
// the coordinate form of the result depends on the coordinate form of inputs
{
	if ((v1.Size() != 3) || (v2.Size() != 3))
		opserr << "\n ERROR! PM4Sand::Dyadic2_2 requires vector of size(3)!" << endln;

	Matrix result(3,3);

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}
/*************************************************************/
// DoubleDot4_2() ---------------------------------------------
Vector
PM4Sand::DoubleDot4_2(const Matrix& m1, const Vector& v1)
// computes doubledot product for matrix-vector arguments
// caution: second coordinate of the matrix should be in opposite variant form of vector
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Sand::DoubleDot4_2 requires vector of size(3)!" << endln;
	if ((m1.noCols() != 3) || (m1.noRows() != 3)) 
		opserr << "\n ERROR! PM4Sand::DoubleDot4_2 requires 3-by-3 matrix " << endln;

	return m1*v1;
}
/*************************************************************/
// DoubleDot2_4() ---------------------------------------------
Vector
PM4Sand::DoubleDot2_4(const Vector& v1, const Matrix& m1)
// computes doubledot product for matrix-vector arguments
// caution: first coordinate of the matrix should be in opposite 
// variant form of vector
{
	if (v1.Size() != 3)
		opserr << "\n ERROR! PM4Sand::DoubleDot2_4 requires vector of size(3)!" << endln;
	if ((m1.noCols() != 3) || (m1.noRows() != 3) )
		opserr << "\n ERROR! PM4Sand::DoubleDot2_4 requires 3-by-3 matrix " << endln;

	return  m1^v1;
}
/*************************************************************/
// DoubleDot4_4() ---------------------------------------------
Matrix
PM4Sand::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
// caution: second coordinate of the first matrix should be in opposite 
// variant form of the first coordinate of second matrix
{
	if ((m1.noCols() != 3) || (m1.noRows() != 3) || (m2.noCols() != 3) || (m2.noRows() != 3)) 
		opserr << "\n ERROR! PM4Sand::DoubleDot4_4 requires 3-by-3 matrices " << endln;

	return m1*m2;
}
/*************************************************************/
// SingleDot4_2() ---------------------------------------------
Matrix
PM4Sand::SingleDot4_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! PM4Sand::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
		opserr << "\n ERROR! PM4Sand::SingleDot4_2 requires 6-by-6 matrix " << endln;

	Matrix result(6,6);
	for (int i = 0; i < 6; i++){
		result(i,0) = m1(i,0) * v1(0) + m1(i,3) * v1(3) + m1(i,5) * v1(5);
		result(i,1) = m1(i,3) * v1(3) + m1(i,1) * v1(1) + m1(i,4) * v1(4);
		result(i,2) = m1(i,5) * v1(5) + m1(i,4) * v1(4) + m1(i,2) * v1(2);
		result(i,3) = 0.5 * (m1(i,0) * v1(3) + m1(i,3) * v1(1) + m1(i,5) * v1(4)
							+m1(i,3) * v1(0) + m1(i,1) * v1(3) + m1(i,4) * v1(5));
		result(i,4) = 0.5 * (m1(i,3) * v1(5) + m1(i,1) * v1(4) + m1(i,4) * v1(2)
							+m1(i,5) * v1(3) + m1(i,4) * v1(1) + m1(i,2) * v1(4));
		result(i,5) = 0.5 * (m1(i,0) * v1(5) + m1(i,3) * v1(4) + m1(i,5) * v1(2)
							+m1(i,5) * v1(0) + m1(i,4) * v1(3) + m1(i,2) * v1(5));
	}
	
	return result;
}
/*************************************************************/
// SingleDot2_4() ---------------------------------------------
Matrix
PM4Sand::SingleDot2_4(const Vector& v1, const Matrix& m1)
// computes singledot product for vector-matrix arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
		opserr << "\n ERROR! PM4Sand::SingleDot2_4 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
		opserr << "\n ERROR! PM4Sand::SingleDot2_4 requires 6-by-6 matrix " << endln;

	Matrix result(6,6);
	for (int i = 0; i < 6; i++){
		result(0,i) = m1(0,i) * v1(0) + m1(3,i) * v1(3) + m1(5,i) * v1(5);
		result(1,i) = m1(3,i) * v1(3) + m1(1,i) * v1(1) + m1(4,i) * v1(4);
		result(2,i) = m1(5,i) * v1(5) + m1(4,i) * v1(4) + m1(2,i) * v1(2);
		result(3,i) = 0.5 * (m1(0,i) * v1(3) + m1(3,i) * v1(1) + m1(5,i) * v1(4)
							+m1(3,i) * v1(0) + m1(1,i) * v1(3) + m1(4,i) * v1(5));
		result(4,i) = 0.5 * (m1(3,i) * v1(5) + m1(1,i) * v1(4) + m1(4,i) * v1(2)
							+m1(5,i) * v1(3) + m1(4,i) * v1(1) + m1(2,i) * v1(4));
		result(5,i) = 0.5 * (m1(0,i) * v1(5) + m1(3,i) * v1(4) + m1(5,i) * v1(2)
							+m1(5,i) * v1(0) + m1(4,i) * v1(3) + m1(2,i) * v1(5));
	}
	return result;
}
/*************************************************************/
// Trans_SingleDot4T_2() ---------------------------------------------
Matrix
PM4Sand::Trans_SingleDot4T_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! PM4Sand::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! PM4Sand::SingleDot4_2 requires 6-by-6 matrix " << endln;
	Matrix result(6,6);
	for (int i = 0; i < 6; i++){
		result(0,i) = m1(0,i) * v1(0) + m1(3,i) * v1(3) + m1(5,i) * v1(5);
		result(1,i) = m1(3,i) * v1(3) + m1(1,i) * v1(1) + m1(4,i) * v1(4);
		result(2,i) = m1(5,i) * v1(5) + m1(4,i) * v1(4) + m1(2,i) * v1(2);
		result(3,i) = 0.5 * (m1(0,i) * v1(3) + m1(3,i) * v1(1) + m1(5,i) * v1(4)
							+m1(3,i) * v1(0) + m1(1,i) * v1(3) + m1(4,i) * v1(5));
		result(4,i) = 0.5 * (m1(3,i) * v1(5) + m1(1,i) * v1(4) + m1(4,i) * v1(2)
							+m1(5,i) * v1(3) + m1(4,i) * v1(1) + m1(2,i) * v1(4));
		result(5,i) = 0.5 * (m1(0,i) * v1(5) + m1(3,i) * v1(4) + m1(5,i) * v1(2)
							+m1(5,i) * v1(0) + m1(4,i) * v1(3) + m1(2,i) * v1(5));
	}
	return result;
}
/*************************************************************/
// Det() ------------------------------------------------------
double PM4Sand::Det(const Vector& aV)
{
	if (aV.Size() != 6)
		opserr << "\n ERROR! PM4Sand::Det requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	return ( +   aV[0] * aV[1] * aV[2] 
		     + 2*aV[3] * aV[4] * aV[5] 
			 -   aV[0] * aV[5] * aV[5] 
			 -   aV[2] * aV[3] * aV[3] 
			 - 	 aV[1] * aV[4] * aV[4]);
}
/*************************************************************/
// Inv() ------------------------------------------------------
Vector PM4Sand::Inv(const Vector& aV)
{
	if (aV.Size() != 6)
		opserr << "\n ERROR! PM4Sand::Inv requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	double det = Det(aV);
	if (det == 0) 
	{
		opserr << "\n Error! PM4Sand::Inv - Singular tensor - return 0 tensor" << endln;
		return aV;
	}
	Vector res(6);
	res(0) = aV(1)*aV(2)-aV(4)*aV(4);
	res(1) = aV(0)*aV(2)-aV(5)*aV(5);
	res(2) = aV(0)*aV(1)-aV(3)*aV(3);
	res(3) = aV(4)*aV(5)-aV(2)*aV(3);
	res(4) = aV(3)*aV(5)-aV(0)*aV(4);
	res(5) = aV(3)*aV(4)-aV(1)*aV(5);
	res = res / det;
	
	return res;
}






/*************************************************************/
/*************************************************************/
//            OTHER AUXILLARY MATH FUNCTIONS                 //
/*************************************************************/
/*************************************************************/
//------------------------------------------------------------
double
MatrixMax_Rows(const Matrix& mat, int rowNo)
{
	int n = mat.noCols();
	double Rmax = fabs(mat(rowNo, 0));
	for (int i = 1; i < n; i++){
		if (Rmax < fabs(mat(rowNo, i)))
			Rmax = fabs(mat(rowNo, i));
	}
	return Rmax;
}
//------------------------------------------------------------
double
MatrixMax_Cols(const Matrix& mat, int colNo)
{
	int n = mat.noRows();
	double Cmax = fabs(mat(0, colNo));
	for (int i = 1; i < n; i++){
		if (Cmax < fabs(mat(i, colNo)))
			Cmax = fabs(mat(i, colNo));
	}
	return Cmax;
}
//------------------------------------------------------------
double
MatrixMin_Rows(const Matrix& mat, int rowNo)
{
	int n = mat.noCols();
	double Rmin = fabs(mat(rowNo, 0));
	for (int i = 1; i < n; i++){
		if (Rmin > fabs(mat(rowNo, i)))
			Rmin = fabs(mat(rowNo, i));
	}
	return Rmin;
}
//------------------------------------------------------------
double
MatrixMin_Cols(const Matrix& mat, int colNo)
{
	int n = mat.noRows();
	double Cmin = fabs(mat(0, colNo));
	for (int i = 1; i < n; i++){
		if (Cmin > fabs(mat(i, colNo)))
			Cmin = fabs(mat(i, colNo));
	}
	return Cmin;
}
//------------------------------------------------------------
double
Sgn(const double& x)
{
	return x == 0 ? 1.0 : (fabs(x)/x);
}
//------------------------------------------------------------
double
VectorMax(const Vector& v)
{
	int n = v.Size();
	double Vmax = fabs(v(0));
	for (int i = 1; i < n; i++){
		if (Vmax < fabs(v(i)))
			Vmax = fabs(v(i));
	}
	return Vmax;
}
//------------------------------------------------------------
// max of two vectors
Vector 
Vector2Max(const Vector& v1, const Vector& v2) {
	if (v1.Size() != v2.Size())
		opserr << "\n ERROR! PM4Sand::Vector2Max requires vector of size(3)!" << endln;
	Vector result(3);
	for (int i = 0; i < v1.Size(); i++) {
		result(i) = (v1(i) > v2(i)) ? v1(i) : v2(i);
	}
	return result;
}
//------------------------------------------------------------
// min of two vectors
Vector 
Vector2Min(const Vector& v1, const Vector& v2) {
	if (v1.Size() != v2.Size())
		opserr << "\n ERROR! PM4Sand::Vector2Min requires vector of size(3)!" << endln;
	Vector result(3);
	for (int i = 0; i < v1.Size(); i++) {
		result(i) = (v1(i) > v2(i)) ? v2(i) : v1(i);
	}
	return result;
}