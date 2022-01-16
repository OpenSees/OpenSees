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
                                                                        
// Written: Alborz Ghofrani, Pedro Arduino
//			Nov 2014, University of Washington
                                                                      
// Description: This file contains the implementation for the ManzariDafaliasRO class.

#include <ManzariDafaliasRO.h>
#include <ManzariDafalias3DRO.h>
#include <ManzariDafaliasPlaneStrainRO.h>
#include <MaterialResponse.h>

#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <algorithm>
#define fmax std::max
#define fmin std::min
#endif

static int numManzariDafaliasMaterials = 0;

void * OPS_ADD_RUNTIME_VPV(OPS_ManzariDafaliasMaterialRO)
{
  if (numManzariDafaliasMaterials == 0) {
    numManzariDafaliasMaterials++;
    opserr << "ManzariDafaliasRO nDmaterial - Written: A.Ghofrani, P.Arduino, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 22) {
    opserr << "Want: nDMaterial ManzariDafaliasRO tag? G0? nu? B? a1? gamma1? e_init? Mc? c? lambda_c? e0? ksi?" <<
		" P_atm? m? h0? Ch? nb? A0? nd? z_max? cz? Rho? <kappa? IntScheme? TanType? JacoType? TolF? TolR?>" << endln;
    return 0;	
  }
  
  int tag;
  double dData[21];
  double oData[6];

  oData[0] = 2.0;		// kappa (Ramberg Osgood)
  oData[1] = 2;			// IntScheme
  oData[2] = 2;			// TanType
  oData[3] = 1;			// JacoType
  oData[4] = 1.0e-7;	// TolF
  oData[5] = 1.0e-7;	// TolR

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial ManzariDafaliasRO material tag" << endln;
    return 0;
  }

  numData = 21;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ManzariDafaliasRO material  with tag: " << tag << endln;
    return 0;
  }

  numData = numArgs - 22;
  if (numData != 0)
	if (OPS_GetDouble(&numData, oData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial ManzariDafaliasRO material  with tag: " << tag << endln;
		return 0;
	}

	theMaterial = new ManzariDafaliasRO(tag, ND_TAG_ManzariDafaliasRO , dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
					     dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], 
						 dData[14], dData[15], dData[16], dData[17], dData[18], dData[19], dData[20], 
					     oData[0],(int)oData[1], (int)oData[2], (int)oData[3], oData[4], oData[5]);

  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafaliasRO material with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafaliasRO::ManzariDafaliasRO(int tag, int classTag, double G0, double nu, 
	double B, double a1, double gamma1, double e_init, double Mc, double c, 
	double lambda_c, double e0, double ksi,	double P_atm, double m, double h0, 
	double ch, double nb, double A0, double nd,	double z_max, double cz, double mDen,
	double kappa, int integrationScheme, int tangentType, int JacoType, double TolF, 
	double TolR): ManzariDafalias(tag, classTag, G0, nu, e_init, Mc, c, lambda_c, e0, 
								ksi, P_atm, m, h0, ch, nb, A0, nd, z_max, cz, mDen, 
								integrationScheme, tangentType, JacoType, TolF, TolR),
	mDevEpsSR(6),
	mSigmaSR(6)
{
	m_B			= B;
	m_a1		= a1;
	m_gamma1	= gamma1;
	m_kappa		= kappa;
	
	initialize();
}

// null constructor
ManzariDafaliasRO::ManzariDafaliasRO() : ManzariDafalias(),
	mDevEpsSR(6),
	mSigmaSR(6)
{
	initialize();
}

// destructor
ManzariDafaliasRO::~ManzariDafaliasRO()
{
}

NDMaterial*
ManzariDafaliasRO::getCopy(const char *type)
{
	if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		ManzariDafaliasPlaneStrainRO *clone;
		clone = new ManzariDafaliasPlaneStrainRO(this->getTag(), m_G0,  m_nu, m_B, m_a1, m_gamma1, m_e_init,  m_Mc,  
				m_c, m_lambda_c,  m_e0,  m_ksi,  m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, 
				m_nd, m_z_max, m_cz, massDen, m_kappa, mScheme, mTangType, mJacoType, mTolF, mTolR);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		ManzariDafalias3DRO *clone;
     		clone = new ManzariDafalias3DRO(this->getTag(), m_G0,  m_nu, m_B, m_a1, m_gamma1, m_e_init,  m_Mc,  m_c, m_lambda_c,
				m_e0,  m_ksi,  m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz, massDen, 
				m_kappa, mScheme, mTangType, mJacoType, mTolF, mTolR);
	 	return clone;
  	} else {
	  	opserr << "ManzariDafaliasRO::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

void 
ManzariDafaliasRO::integrate()
{
	double chi_e, chi_en;
	Vector devEps(6), devEps_n(6);
	//
	devEps		= GetDevPart(mEpsilon);
	devEps_n	= GetDevPart(mEpsilon_n);
	chi_e		= sqrt(0.5 * DoubleDot2_2_Cov(devEps   - mDevEpsSR, devEps   - mDevEpsSR));
	chi_en		= sqrt(0.5 * DoubleDot2_2_Cov(devEps_n - mDevEpsSR, devEps_n - mDevEpsSR));
	if (mIsFirstShear && fabs(chi_e - chi_en) < 1.0e-10) { // This is required in case of consolidation (mEta1 should be updated)
		// how small should 1.0e-10 be?
		double p = one3 * GetTrace(mSigma_n);
		double Gmax  = m_B * m_P_atm / (0.3 + 0.7 * mVoidRatio*mVoidRatio) * sqrt(p / m_P_atm);
		mEta1 = m_a1 * Gmax * m_gamma1 / p;
		chi_e = chi_en = 0.0;
	}
	if ((chi_e - chi_en) * mDChi_e < -1.0e-14) { // how small should be -1.0e-14?
		mSigmaSR  = mSigma_n;

		//mSigmaSR(3) = mSigmaSR(4) = mSigmaSR(5) -= 0.1;
		mDevEpsSR = GetDevPart(mEpsilon_n);
		double pSR = one3 * GetTrace(mSigmaSR);
		double GmaxSR  = m_B * m_P_atm / (0.3 + 0.7 * mVoidRatio*mVoidRatio) * sqrt(pSR / m_P_atm);
		mEta1 = m_a1 * GmaxSR * m_gamma1 / pSR;
		mIsFirstShear = false;
		GetElasticModuli(mSigma_n, mVoidRatio, mK, mG);
	}
	ManzariDafalias::integrate();
}

int 
ManzariDafaliasRO::commitState(void)
{
	double chi_e, chi_en;
    Vector devEps(6), devEps_n(6);
    
    devEps          = GetDevPart(mEpsilon);
    devEps_n        = GetDevPart(mEpsilon_n);
    chi_e           = sqrt(0.5 * DoubleDot2_2_Cov(devEps   - mDevEpsSR, devEps   - mDevEpsSR));
    chi_en          = sqrt(0.5 * DoubleDot2_2_Cov(devEps_n - mDevEpsSR, devEps_n - mDevEpsSR));

    mDChi_e			= chi_e - chi_en;
	
	int res = ManzariDafalias::commitState();
	GetElasticModuli(mSigma, mVoidRatio, mK, mG);

	return res;
}

NDMaterial*
ManzariDafaliasRO::getCopy (void)
{
	opserr << "ManzariDafaliasRO::getCopy -- subclass responsibility\n"; 
  	exit(-1);
  	return 0;
}

const char*
ManzariDafaliasRO::getType (void) const
{
    opserr << "ManzariDafaliasRO::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
ManzariDafaliasRO::getOrder (void) const
{
    opserr << "ManzariDafaliasRO::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
ManzariDafaliasRO::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int 
ManzariDafaliasRO::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
    return 0;
}

void ManzariDafaliasRO::Print(OPS_Stream &s, int flag )
{
	s << "ManzariDafaliasRO Material, tag: " << this->getTag() << endln;
	s << "Type: " << this->getType() << endln;
}

// Initialize Manzari Dafalias Material
void 
ManzariDafaliasRO::initialize()
{
	mSigma_n = mSigma = mSigmaSR = m_Pmin * mI1;
	
	mDChi_e = 0.0;
	double GmaxSR  = m_B * m_P_atm / (0.3 + 0.7 * mVoidRatio*mVoidRatio) * sqrt(m_Pmin / m_P_atm);
	mEta1 = m_a1 * GmaxSR * m_gamma1 / m_Pmin;
	mIsFirstShear = true;

	// calculate initial stiffness parameters
	GetElasticModuli(mSigma_n,mVoidRatio,mK,mG);
	mCe = GetStiffness(mK,mG);
	mCep = mCe;
	mCep_Consistent = mCe;
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
ManzariDafaliasRO::GetElasticModuli(const Vector& sigma, const double& en, const double& en1, const Vector& nEStrain, 
				const Vector& cEStrain, double &K, double &G)
// Calculates G, K
{
	double p, pSR, Gmax, T, temp;
	Vector r(6), rSR(6);

	p = one3 * GetTrace(sigma);
	p = (p <= m_Pmin) ? m_Pmin : p;
	r = GetDevPart(sigma) / p;

	pSR = one3 * GetTrace(mSigmaSR);
	pSR = (pSR <= m_Pmin) ? m_Pmin : pSR;
	rSR = GetDevPart(mSigmaSR) / pSR;

	Gmax = m_B * m_P_atm / (0.3 + 0.7 * en * en) * sqrt(p / m_P_atm);
	if (mElastFlag == 0) {
		mIsFirstShear = true;
		T = 1.0;
	} else {
		mChi_r = sqrt(0.5 * DoubleDot2_2_Contr(r-rSR, r-rSR));
		temp = m_kappa * (1.0 / m_a1 - 1);
		if (mIsFirstShear)
			T = 1 + temp * pow(mChi_r / mEta1, m_kappa - 1);
		else
			T = 1 + temp * pow(mChi_r / mEta1 / 2.0, m_kappa - 1);
		T = (T < (1.0+temp)) ? T : (1.0+temp);
		T = (T < 1.0) ? 1.0 : T;
	}

	G = Gmax / T;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;	
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
ManzariDafaliasRO::GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G)
// Calculates G, K
{
	double p, pSR, Gmax, T, temp;
	Vector r(6), rSR(6);

	p = one3 * GetTrace(sigma);
	p = (p <= m_Pmin) ? m_Pmin : p;
	r = GetDevPart(sigma) / p;

	pSR = one3 * GetTrace(mSigmaSR);
	pSR = (pSR <= m_Pmin) ? m_Pmin : pSR;
	rSR = GetDevPart(mSigmaSR) / pSR;

	Gmax = m_B * m_P_atm / (0.3 + 0.7 * en * en) * sqrt(p / m_P_atm);
	if (mElastFlag == 0) {
		mIsFirstShear = true;
		T = 1.0;
	} else {
		mChi_r = sqrt(0.5 * DoubleDot2_2_Contr(r-rSR, r-rSR));
		temp = m_kappa * (1.0 / m_a1 - 1);
		if (mIsFirstShear)
			T = 1 + temp * pow(mChi_r / mEta1, m_kappa - 1);
		else
			T = 1 + temp * pow(mChi_r / mEta1 / 2.0, m_kappa - 1);
		T = (T < (1.0+temp)) ? T : (1.0+temp);
		T = (T < 1.0) ? 1.0 : T;
	}

	G = Gmax / T;
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;	
}
