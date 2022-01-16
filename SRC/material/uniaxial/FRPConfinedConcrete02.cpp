// Written/implemented by: J.Y. LU,							Southeast University, Nanjing, China
//						   G. LIN (guan.lin@polyu.edu.hk),	The Hong Kong Polytechnic University, Hong Kong, China
//
// Description: This file contains the class definition for Uniaxial material FRPConfinedConcrete02. 
//		FRPConfinedConcrete02.cpp

/* *********************************************************************************************
** Compression Part:
** Cyclic behavior based on Ref: 
**		L. Lam and J.G. Teng (2009), 
**		"Stress-strain model for FRP-confined concrete under cyclic axial compression",
**		Engineering Structures, 31, 308-321.										
** Envelop curve based on Ref:
**		J.G. Teng, T. Jiang, L. Lam, Y.Z. Luo (2009),
**		"Refinement of a design-oriented stress-strain model for FRP-confined concrete",
**		Journal of Composites for Construction, ASCE, 13(4), 269-278
** Column simulation based on Ref:
**		J.G. Teng, L. Lam, G. Lin, J.Y. Lu and Q.G. Xiao (2015),
**		Numerical Simulation of FRP-Jacketed RC Columns Subjected to Cyclic and Seismic Loading,
**		Journal of Composites for Construction, ASCE, 20(1), 04015021.
** Tension Part:
** Based on Ref:
**		M.H.M. Yassin (1994),
**		Nonlinear analysis of prestressed concrete structures under monotonic and cyclic loads. 
**		Dissertation.University of California. Berkeley, California.
**		and made some improvements
** *********************************************************************************************/

#include <elementAPI.h>
#include "FRPConfinedConcrete02.h"

#include <Vector.h>
#include <Channel.h>
#include <math.h>
#include <float.h>

#define TAG_LOADING 1
#define TAG_UNLOADING -1

void * OPS_ADD_RUNTIME_VPV(OPS_FRPConfinedConcrete02)
{
  // Pointer to a uniaxial material that will be returned
  UniaxialMaterial *theMaterial = 0;
  
  int    iData[1];
  int numData = 1;
  
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid uniaxialMaterial FRPConfinedConcrete02 tag" << endln;
    return 0;
  }
  
  numData = OPS_GetNumRemainingInputArgs(); // get number of input args excluding the material tag
  
  if (numData != 6 && numData != 9 && numData != 11) {
    opserr << "Incorrect # args, want: uniaxialMaterial FRPConfinedConcrete02 tag? fc0? ec0? Ec? ft? Ets? Unit?" << endln;
    opserr << "Or: uniaxialMaterial FRPConfinedConcrete02 tag? fc0? ec0? Ec? -Ultimate fcc? ecu? ft? Ets? Unit?" << endln;
    opserr << "Or: uniaxialMaterial FRPConfinedConcrete02 tag? fc0? ec0? Ec? -JacketC t? Efrp? eps_h_rup? R? ft? Ets? Unit?" << endln;
    return 0;
  }
  
  if (numData == 6){ // Input for unconfined concrete
    int numData1 = 6;
    double dData[6];
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << "fc0? ec0? Ec? ft? Ets? Unit?" << endln;
      return 0;	
    }
    theMaterial = new FRPConfinedConcrete02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], (int)dData[5]);
    
  } else if (numData == 9) { // Ultimate stress/strain input by users
    double dData[8];
    int numData1 = 3;
    int numData2 = 5;
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << "fc0? ec0? Ec? -Ultimate fcc? ecu? ft? Ets? Unit?" << endln;
      return 0;	
    }
    
    const char *str = OPS_GetString();
    // OPS_GetStringCopy(&str);
    if (strcmp(str, "-Ultimate") == 0) {
      if (OPS_GetDoubleInput(&numData2, dData+3) != 0) {
	opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << "fc0? ec0? Ec? -Ultimate fcc? ecu? ft? Ets? Unit?" << endln;
	return 0;
      }
    } else {
      opserr << "Invalid input parameter for uniaxialMaterial FRPConfinedConcrete02 with tag  " << iData[0] << ", want: -Ultimate"<< endln;
      return 0;
    }
    theMaterial = new FRPConfinedConcrete02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], (int)dData[7]);     
  } else { // FRP-confined concrete in circular columns, FRP jacket properties input by users
    double dData[10];
    int numData1 = 3;
    int numData2 = 7;
    if (OPS_GetDoubleInput(&numData1, dData) != 0) {
      opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << "fc0? ec0? Ec? -JacketC tfrp? Efrp? erup? R? ft? Ets? Unit?" << endln;
      return 0;	
    }
    
    const char *str = OPS_GetString();
    //OPS_GetStringCopy(&str);
    if (strcmp(str, "-JacketC") == 0) {
      if (OPS_GetDoubleInput(&numData2, dData+3) != 0) {
	opserr << "Invalid #args, want: uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << "fc0? ec0? Ec? -JacketC tfrp? Efrp? erup? R? ft? Ets? Unit?" << endln;
	return 0;
      }
    } else {
      opserr << "Invalid input parameter for uniaxialMaterial FRPConfinedConcrete02 with tag " << iData[0] << ", want: -JacketC"<< endln;
      return 0;
    }
    theMaterial = new FRPConfinedConcrete02(iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5], dData[6], dData[7], dData[8], (int)dData[9]);    
  }
  
  if (theMaterial == 0) {
    opserr << "WARNING could not create uniaxialMaterial FRPConfinedConcrete02 " << iData[0] << endln;
    return 0;
  }
  
  return theMaterial;
}

// Constructor for unconfined concrete
FRPConfinedConcrete02::FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double ft, double Ets, int Unit)
:UniaxialMaterial(tag,MAT_TAG_FRPConfinedConcrete02),
m_fc0(-fc0), m_Ec(Ec), m_epsc0(-ec0), m_Ets(Ets), m_ft(-ft), m_Unit(Unit)
{
	m_Tstrain = 0.0;
	m_Tstress = 0.0;
	m_trialTangent = Ec;
	m_Unitscale = 1.0;
	if (m_Unit == 0)
		m_Unitscale = 6.895;
	//////////////////////////////////////////////////////////////////////////
	m_fcc = m_fc0*0.85;
	m_epscu = m_epsc0*1.75;  // Unconfined concrete, Teng et al. (2009), Journal of Composites for Construction, 13(4), 269-278

	m_E2 = (m_fcc-m_fc0)/m_epscu;
	m_epst = 2.0*m_fc0/(m_Ec-m_E2);

	m_Eun0 = m_Ec;

	m_Etr1 = m_Ec;
	m_Etr2 = m_Ec;
	m_epstn = (1.0)*m_ft/m_Etr1;
	m_epstu = (1.0)*(m_epstn + m_ft/m_Ets);

	m_fi = 1.0;
	m_fiful = 1.0;
	m_gamare = 0;
	m_betaun = 0;

	m_loadingflag = TAG_LOADING;
	m_n=0;
	m_ne=1;
	m_epsunenv = 0.0;
	m_sigunenv = 0.0;
	m_trialStrainlast = 0.0;
	m_trialStresslast = 0.0;
	m_epsretenv = 0.0;
	m_epspl=0.0;
	m_bSmallStress = false;
	m_bUltimate = false;
	m_trialTangentlast = m_trialTangent;

	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////
}

// Constructor for the case of <-Ultimate> (ultimate stress/strain input by users)
FRPConfinedConcrete02::FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double fcc, double ecu, double ft, double Ets, int Unit)
:UniaxialMaterial(tag,MAT_TAG_FRPConfinedConcrete02)
{
	m_fc0 = -fc0;
	m_Ec = Ec;
	m_epsc0 = -ec0;
	m_Ets = Ets;
	m_ft = -ft; 
	m_Unit = Unit;

	m_Tstrain = 0.0;
	m_Tstress = 0.0;
	m_trialTangent = Ec;
	m_Unitscale = 1.0;
	if (m_Unit == 0)
		m_Unitscale = 6.895;
	//////////////////////////////////////////////////////////////////////////
	m_fcc = -fcc;
	m_epscu = -ecu;
	m_E2 = (m_fcc-m_fc0)/m_epscu;
	m_epst = 2.0*m_fc0/(m_Ec-m_E2);

	m_Eun0 = m_Ec;

	m_Etr1 = m_Ec;
	m_Etr2 = m_Ec;
	m_epstn = (1.0)*m_ft/m_Etr1;
	m_epstu = (1.0)*(m_epstn + m_ft/m_Ets);

	m_fi = 1.0;
	m_fiful = 1.0;
	m_gamare = 0;
	m_betaun = 0;

	m_loadingflag = TAG_LOADING;
	m_n=0;
	m_ne=1;
	m_epsunenv = 0.0;
	m_sigunenv = 0.0;
	m_trialStrainlast = 0.0;
	m_trialStresslast = 0.0;
	m_epsretenv = 0.0;
	m_epspl=0.0;
	m_bSmallStress = false;
	m_bUltimate = false;
	m_trialTangentlast = m_trialTangent;

	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////
}

// Constructor for the case of <-JacketC> (FRP-confined concrete in circular columns, FRP jacket properties input by users)
FRPConfinedConcrete02::FRPConfinedConcrete02(int tag, double fc0, double Ec, double ec0, double t, double Efrp, double eps_h_rup, double R, double ft, double Ets, int Unit)
:UniaxialMaterial(tag,MAT_TAG_FRPConfinedConcrete02)
{
	m_fc0 = -fc0;
	m_Ec = Ec;
	m_epsc0 = -ec0;
	m_t = t;
	m_Efrp = Efrp;
	m_eps_h_rup = eps_h_rup;
	m_R = R;
	m_Ets = Ets;
	m_ft = -ft; 
	m_Unit = Unit;

	m_Tstrain = 0.0;
	m_Tstress = 0.0;
	m_trialTangent = Ec;
	m_Unitscale = 1.0;
	if (m_Unit == 0)
		m_Unitscale = 6.895;
	//////////////////////////////////////////////////////////////////////////
	// Jiang & Teng (2009) envelop
	m_fl = m_Efrp*m_t*m_eps_h_rup/m_R;
	m_fcc = m_fc0*(1.0 + 3.5*m_fl/m_fc0 - 0.035*m_eps_h_rup/m_epsc0);
	m_epscu = m_epsc0*(1.75 + 6.5*pow(m_fl/m_fc0,0.8)*pow(m_eps_h_rup/m_epsc0,0.65));
	m_E2 = (m_fcc-m_fc0)/m_epscu;
	m_epst = 2.0*m_fc0/(m_Ec-m_E2);

	m_Eun0 = m_Ec;

	m_Etr1 = m_Ec;
	m_Etr2 = m_Ec;
	m_epstn = (1.0)*m_ft/m_Etr1;
	m_epstu = (1.0)*(m_epstn + m_ft/m_Ets);

	m_fi = 1.0;
	m_fiful = 1.0;
	m_gamare = 0;
	m_betaun = 0;

	m_loadingflag = TAG_LOADING;
	m_n=0;
	m_ne=1;
	m_epsunenv = 0.0;
	m_sigunenv = 0.0;
	m_trialStrainlast = 0.0;
	m_trialStresslast = 0.0;
	m_epsretenv = 0.0;
	m_epspl=0.0;
	m_bSmallStress = false;
	m_bUltimate = false;
	m_trialTangentlast = m_trialTangent;

	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////
}

FRPConfinedConcrete02::FRPConfinedConcrete02()
:UniaxialMaterial(0,MAT_TAG_FRPConfinedConcrete02),
m_fc0(0.0), m_Ec(0.0), m_fcc(0.0), m_epscu(0.0), m_Ets(0.0), m_ft(0.0), m_Unit(1),
m_Tstrain(0.0), m_Tstress(0.0), m_trialTangent(0.0)
{

	// AddingSensitivity:BEGIN /////////////////////////////////////
	parameterID = 0;
	SHVs = 0;
	// AddingSensitivity:END //////////////////////////////////////
}

FRPConfinedConcrete02::~FRPConfinedConcrete02()
{
}

int 
FRPConfinedConcrete02::setTrialStrain(double strain, double strainRate)
{

	m_Tstrain = -strain;

	if (m_bUltimate)
	{
		m_Tstress = 0.0;
		m_trialStrainlast = m_Tstrain;
		m_trialStresslast = m_Tstress;
		m_trialTangent = 1e-15;
		return 0;
	}

	if (m_Tstrain == m_trialStrainlast && m_trialStrainlast == 0.0)
	{
		m_Tstress = 0.0;
		m_trialStrainlast = m_Tstrain;
		m_trialStresslast = m_Tstress;
		return 0;
	}
	if (fabs(m_Tstrain - m_trialStrainlast) <= 1e-15) //Lu JY add  2009.06.06
	{
		m_Tstrain = m_trialStrainlast;
		m_Tstress = m_trialStresslast;
		m_trialTangent = m_trialTangentlast;

		return 0;
	}

	if (m_Tstrain >= m_epspl) // Compression part
	{
		if (m_Tstrain > m_epscu)
		{
			m_Tstress = 0.0;
			m_trialStrainlast = m_Tstrain;
			m_trialStresslast = m_Tstress;
			m_trialTangent = 1e-15;
			m_bUltimate = true;
			return 0;
		}
		// Define loading types & parameters
		if ((m_Tstrain - m_trialStrainlast) > 0) // loading & reloading
		{
			if ((m_loadingflag == TAG_UNLOADING) || (m_trialStrainlast < m_epspl)) // Begin loading point: Show define Reloading curve parameters here
			{
				// define reloading point
				m_epsre = m_trialStrainlast;
				m_sigre = m_trialStresslast;
				if (m_trialStrainlast < m_epspl)
				{
					m_epsre = m_epspl;
					m_sigre = 0.0;
				}
				if (m_n == 1 && m_sigunenv != 0) // Eq.28 (Lam and Teng 2009)
					m_betaun = (m_sigunenv - m_sigre)/m_sigunenv;
				else if (m_n >= 2 && m_signew != 0)
					m_betaun = (m_sigun - m_sigre)/m_signew; // m_signew: last step
				GetDeterioratedStress();

			}
			m_loadingflag = TAG_LOADING; // Set current loading type
		}
		else
		{
			if (m_loadingflag == TAG_LOADING) // Begin unloading point: Show define Unloading curve parameters here
			{
				if (m_trialStresslast > m_sigunenv)
					m_n = 1;
				else if (m_trialStresslast <= m_sigunenv)
					m_n+=1;
				// define unloading point
				m_epsun = m_trialStrainlast;
				m_sigun = m_trialStresslast;

				if (m_n == 1)
				{
					m_epsunenv = m_epsun;
					m_sigunn1 = m_sigunenv = m_sigun;
				}
				// Calculate reference point
				GetRefPoint();
				if (m_n >= 2 && (m_epsreflast - m_epspl)!=0)
				{
					m_gamare = (m_epsun - m_epspl)/(m_epsreflast - m_epspl); // Eq.29 (Lam and Teng 2009)
					// Calculate effective cyclic loading num
					if (m_gamare > 0.7 && m_betaun > 0.7) // Eq.30 (Lam and Teng 2009)
						m_ne +=1;
				}
				Compr_GetPlasticStrain(); // calculate new m_epspl
			}
			m_loadingflag = TAG_UNLOADING; // Set current loading type
		}

		if (m_loadingflag == TAG_LOADING)
		{
			if (m_n == 0) // Loading on the Envelop curve
				Compr_Envlp(m_Tstrain,m_Tstress,m_trialTangent);
			else if (m_n >= 1) // Internal Reloading
				Compr_ReloadingPath(m_Tstrain,m_Tstress,m_trialTangent);
		} 
		else if (m_loadingflag == TAG_UNLOADING)
		{
			if (m_n == 1) // Unloading from envelop
			{
				m_ne = 1;
				Compr_UnloadingPath(m_Tstrain,m_Tstress,m_trialTangent);
			}
			else if (m_n >= 2) // Internal unloading 
				Compr_UnloadingPath(m_Tstrain,m_Tstress,m_trialTangent);
		}
	}
	else // Tension part
	{
		if (m_Tstrain <= (m_epstu + m_epspl) || fabs(m_Etr2)<1)
		{
			m_Etr1 = m_Etr2 = 0.0;
			m_Tstress = 0;
			m_trialTangent = 1e-15;		
			return 0;
		}
		if ((m_Tstrain - m_trialStrainlast) > 0)
		{
			if (m_loadingflag == TAG_UNLOADING) // Begin loading point: Show define Reloading curve parameters here
			{
				if (m_trialStrainlast <= (m_epstn+m_epspl) && (m_epspl != m_trialStrainlast))
					m_Etr2 = (0.0 - m_trialStresslast)/(m_epspl - m_trialStrainlast);
				m_Etr1 = m_Etr2;
			}
			m_loadingflag = TAG_LOADING; // Set current loading type
		}
		else
			m_loadingflag = TAG_UNLOADING; // Set current loading type
		double epstemp;
		epstemp = m_Tstrain - m_epspl;

		if (m_loadingflag == TAG_LOADING)
		{
			m_Tstress = epstemp*m_Etr2;
			m_trialTangent = m_Etr2;
		}
		else if(m_loadingflag == TAG_UNLOADING)
		{
			m_Etr1 = m_Etr1<m_Eun0 ? m_Etr1:m_Eun0;
			m_epstn = m_epstu/(m_Etr1/m_Ets + 1.0);
			Tens_Envlp(epstemp,m_Tstress,m_trialTangent);
		}
	}
	return 0;
}

double 
FRPConfinedConcrete02::getStrain(void)
{
	return -m_Tstrain;
}

double 
FRPConfinedConcrete02::getStress(void)
{
	return -m_Tstress;
}


double 
FRPConfinedConcrete02::getTangent(void)
{
	return m_trialTangent;
}

int 
FRPConfinedConcrete02::commitState(void)
{
	m_epstnlast = m_epstn;
	m_epstulast = m_epstu;
	m_Etr1last = m_Etr1;
	m_Etr2last = m_Etr2;

	// History
	m_nlast = m_n;
	m_nelast = m_ne;
	m_loadingflaglast = m_loadingflag;

	m_Erelast = m_Ere;
	m_epsunenvlast = m_epsunenv;
	m_sigunenvlast = m_sigunenv;
	m_sigunn1last = m_sigunn1;
	m_epsretenvlast = m_epsretenv;
	m_epsunlast = m_epsun;
	m_sigunlast = m_sigun;
	m_epsrelast = m_epsre;
	m_sigrelast = m_sigre;
	m_epsreflast = m_epsref;
	m_sigreflast = m_sigref;
	m_Eun0last = m_Eun0; 


	m_betaunlast = m_betaun;
	m_filast = m_fi;
	m_fifullast = m_fiful;
	m_signewlast = m_signew;

	m_gamarelast = m_gamare;
	m_omglast = m_omg;
	m_omgfullast = m_omgful;
	m_epspllast = m_epspl;

	m_bSmallStresslast = m_bSmallStress;
	m_bUltimatelast = m_bUltimate;

	m_trialStrainlast = m_Tstrain;
	m_trialStresslast = m_Tstress;

	m_trialTangentlast = m_trialTangent;
	return 0;
}	


int 
FRPConfinedConcrete02::revertToLastCommit(void)
{
	m_epstn = m_epstnlast;
	m_epstu = m_epstulast;
	m_Etr1 = m_Etr1last;
	m_Etr2 = m_Etr2last;

	// History
	m_n = m_nlast;
	m_ne = m_nelast;
	m_loadingflag = m_loadingflaglast;

	m_Ere = m_Erelast;
	m_epsunenv = m_epsunenvlast;
	m_sigunenv = m_sigunenvlast;
	m_sigunn1 = m_sigunn1last;
	m_epsretenv = m_epsretenvlast;
	m_epsun = m_epsunlast;
	m_sigun = m_sigunlast;
	m_epsre = m_epsrelast;
	m_sigre = m_sigrelast;
	m_epsref = m_epsreflast;
	m_sigref = m_sigreflast;
	m_Eun0 = m_Eun0last; 

	m_betaun = m_betaunlast;
	m_fi = m_filast;
	m_fiful = m_fifullast;
	m_signew = m_signewlast;

	m_gamare = m_gamarelast;
	m_omg = m_omglast;
	m_omgful = m_omgfullast;
	m_epspl = m_epspllast;

	m_bSmallStress = m_bSmallStresslast;
	m_bUltimate = m_bUltimatelast;

	m_Tstrain = m_trialStrainlast;
	m_Tstress = m_trialStresslast;

	m_trialTangent = m_trialTangentlast;
	return 0;
}


int 
FRPConfinedConcrete02::revertToStart(void)
{
	m_Eun0 = m_Ec;
	m_Etr1 = m_Ec;
	m_Etr2 = m_Ec;
	m_epstn = m_ft/m_Etr1;
	m_epstu = m_epstn + m_ft/m_Ets;

	m_n=0;
	m_ne=1;
	m_loadingflag = TAG_LOADING;
	m_Ere = 0.0;
	m_epsunenv = 0.0;
	m_sigunenv = 0.0;
	m_sigunn1 = 0.0;
	m_epsretenv = 0.0;

	m_epsun = 0.0;
	m_sigun = 0.0;
	m_epsre = 0.0;
	m_sigre = 0.0;
	m_epsref = 0.0;
	m_sigref = 0.0;

	m_betaun = 0.0;
	m_fi = 1.0;
	m_fiful = 0.0;
	m_signew = 0.0;

	m_gamare = 0.0;
	m_omg = 0.0;
	m_omgful = 0.0;

	m_epspl = 0.0;
	m_bSmallStress = false;
	m_bUltimate = false;

	m_Tstrain = 0.0;
	m_Tstress = 0.0;

	m_trialTangent = m_Ec;

	//////////////////////////////////////////////////////////////////////////
	m_Etr1last = m_Ec;
	m_Etr2last = m_Ec;
	m_epstnlast = m_ft/m_Etr1;
	m_epstulast = m_epstn + m_ft/m_Ets;

	m_nlast = 0;
	m_nelast = 1;
	m_loadingflaglast = TAG_LOADING;
	m_Erelast = 0.0;
	m_epsunenvlast = 0.0;
	m_sigunenvlast = 0.0;
	m_sigunn1last = 0.0;
	m_epsretenvlast = 0.0;

	m_epsunlast = 0.0;
	m_sigunlast = 0.0;
	m_epsrelast = 0.0;
	m_sigrelast = 0.0;
	m_epsreflast = 0.0;;
	m_sigreflast = 0.0;;

	m_betaunlast = 0.0;
	m_filast = 0.0;
	m_fifullast = 0.0;
	m_signewlast = 0.0;

	m_gamarelast = 0.0;
	m_omglast = 0.0;
	m_omgfullast = 0.0;

	m_epspllast = 0.0;
	m_bSmallStresslast = false;
	m_bUltimatelast = false;

	m_trialStrainlast = 0.0;
	m_trialStresslast = 0.0;

	m_trialTangentlast = m_Ec;
	return 0;
}


UniaxialMaterial *
FRPConfinedConcrete02::getCopy(void)
{
	FRPConfinedConcrete02* theCopy =
		new FRPConfinedConcrete02(this->getTag(),-m_fc0,m_Ec,-m_epsc0,-m_fcc,-m_epscu,-m_ft,m_Ets,m_Unit);

	//////////////////////////////////////////////////////////////////////////
	theCopy->m_epstn = this->m_epstn;
	theCopy->m_epstu = this->m_epstu;
	theCopy->m_Etr1 = this->m_Etr1;
	theCopy->m_Etr2 = this->m_Etr2;

	theCopy->m_n = this->m_n;
	theCopy->m_ne = this->m_ne;
	theCopy->m_loadingflag = this->m_loadingflag;

	theCopy->m_Ere = this->m_Ere;
	theCopy->m_epsunenv = this->m_epsunenv;
	theCopy->m_sigunenv = this->m_sigunenv;
	theCopy->m_sigunn1 = this->m_sigunn1;
	theCopy->m_epsretenv = this->m_epsretenv;
	theCopy->m_epsun = this->m_epsun;
	theCopy->m_sigun = this->m_sigun;
	theCopy->m_epsre = this->m_epsre;
	theCopy->m_sigre = this->m_sigre;
	theCopy->m_epsref = this->m_epsref;
	theCopy->m_sigref = this->m_sigref;
	theCopy->m_Eun0 = this->m_Eun0;

	theCopy->m_betaun = this->m_betaun;
	theCopy->m_fi = this->m_fi;
	theCopy->m_fiful = this->m_fiful;
	theCopy->m_signew = this->m_signew;

	theCopy->m_gamare = this->m_gamare;
	theCopy->m_omg = this->m_omg;
	theCopy->m_omgful = this->m_omgful;
	theCopy->m_epspl = this->m_epspl;

	theCopy->m_bSmallStress = this->m_bSmallStress;
	theCopy->m_bUltimate = this->m_bUltimate;

	theCopy->m_Tstrain = this->m_Tstrain;
	theCopy->m_Tstress = this->m_Tstress;
	theCopy->m_trialTangent = this->m_trialTangent;
	//////////////////////////////////////////////////////////////////////////
	theCopy->m_epstnlast = this->m_epstnlast;
	theCopy->m_epstulast = this->m_epstulast;
	theCopy->m_Etr1last = this->m_Etr1last;
	theCopy->m_Etr2last = this->m_Etr2last;

	theCopy->m_nlast = this->m_nlast;
	theCopy->m_nelast = this->m_nelast;
	theCopy->m_loadingflaglast = this->m_loadingflaglast;

	theCopy->m_Erelast = this->m_Erelast;
	theCopy->m_epsunenvlast = this->m_epsunenvlast;
	theCopy->m_sigunenvlast = this->m_sigunenvlast;
	theCopy->m_sigunn1last = this->m_sigunn1last;
	theCopy->m_epsretenvlast = this->m_epsretenvlast;
	theCopy->m_epsunlast = this->m_epsunlast;
	theCopy->m_sigunlast = this->m_sigunlast;
	theCopy->m_epsrelast = this->m_epsrelast;
	theCopy->m_sigrelast = this->m_sigrelast;
	theCopy->m_epsreflast = this->m_epsreflast;
	theCopy->m_sigreflast = this->m_sigreflast;
	theCopy->m_Eun0last = this->m_Eun0last;

	theCopy->m_betaunlast = this->m_betaunlast;
	theCopy->m_filast = this->m_filast;
	theCopy->m_fifullast = this->m_fifullast;
	theCopy->m_signewlast = this->m_signewlast;

	theCopy->m_gamarelast = this->m_gamarelast;
	theCopy->m_omglast = this->m_omglast;
	theCopy->m_omgfullast = this->m_omgfullast;
	theCopy->m_epspllast = this->m_epspllast;

	theCopy->m_bSmallStresslast = this->m_bSmallStresslast;
	theCopy->m_bUltimatelast = this->m_bUltimatelast;
	theCopy->m_trialStrainlast = this->m_trialStrainlast;
	theCopy->m_trialStresslast = this->m_trialStresslast;
	theCopy->m_trialTangentlast = this->m_trialTangentlast;
	return theCopy;
}


int 
FRPConfinedConcrete02::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;
	static Vector data(49);
	data(0) = this->getTag();
	data(1) = m_fc0;
	data(2) = m_Ec;
	data(3) = m_t;
	data(4) = m_Efrp;
	data(5) = m_eps_h_rup;
	data(6) = m_R;
	data(7) = m_Ets;
	data(8) = m_ft;
	data(9) = m_fl;
	data(10) = m_epsc0;
	data(11) = m_fcc;
	data(12) = m_epscu;
	data(13) = m_E2;
	data(14) = m_epst;

	data(15) = m_Unit;
	data(16) = m_Unitscale;
	//////////////////////////////////////////////////////////////////////////
	data(17) = m_epstnlast;
	data(18) = m_epstulast;
	data(19) = m_Etr1last;
	data(20) = m_Etr2last;

	// History
	data(21) = m_nlast;
	data(22) = m_nelast;
	data(23) = m_loadingflaglast;

	data(24) = m_Erelast;
	data(25) = m_epsunenvlast;
	data(26) = m_sigunenvlast;
	data(27) = m_sigunn1last;
	data(28) = m_epsretenvlast;
	data(29) = m_epsunlast;
	data(30) = m_sigunlast;
	data(31) = m_epsrelast;
	data(32) = m_sigrelast;
	data(33) = m_epsreflast;
	data(34) = m_sigreflast;

	data(35) = m_betaunlast;
	data(36) = m_filast;
	data(37) = m_fifullast;
	data(38) = m_signewlast;

	data(39) = m_gamarelast;
	data(40) = m_omglast;
	data(41) = m_omgfullast;
	data(42) = m_epspllast;

	data(43) = m_bSmallStresslast;
	data(44) = m_bUltimatelast;

	data(45) = m_trialStrainlast;
	data(46) = m_trialStresslast;

	data(47) = m_trialTangentlast;
	data(48) = m_Eun0last;

	res = theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "FRPConfinedConcrete02::sendSelf() - failed to send data\n";

	opserr << "sendSelf\n";
	return res;
}

int 
FRPConfinedConcrete02::recvSelf(int cTag, Channel &theChannel, 
								FEM_ObjectBroker &theBroker)
{
	int res = 0;
	static Vector data(49);
	res = theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) 
		opserr << "FRPConfinedConcrete02::recvSelf() - failed to recv data\n";
	else {
		this->setTag(data(0));

		m_fc0 = data(1);
		m_Ec = data(2);
		m_t = data(3);
		m_Efrp = data(4);
		m_eps_h_rup = data(5);
		m_R = data(6);
		m_Ets = data(7);
		m_ft = data(8);
		m_fl = data(9);
		m_epsc0 = data(10);
		m_fcc = data(11);
		m_epscu = data(12);
		m_E2 = data(13);
		m_epst = data(14);

		m_Unit = data(15);
		m_Unitscale = data(16);
		//////////////////////////////////////////////////////////////////////////
		m_epstnlast = data(17);
		m_epstulast = data(18);
		m_Etr1last = data(19);
		m_Etr2last = data(20);

		// History
		m_nlast = data(21);
		m_nelast = data(22);
		m_loadingflaglast = data(23);

		m_Erelast = data(24);
		m_epsunenvlast = data(25);
		m_sigunenvlast = data(26);
		m_sigunn1last = data(27);
		m_epsretenvlast = data(28);
		m_epsunlast = data(29);
		m_sigunlast = data(30);
		m_epsrelast = data(31);
		m_sigrelast = data(32);
		m_epsreflast = data(33);
		m_sigreflast = data(34);

		m_betaunlast = data(35);
		m_filast = data(36);
		m_fifullast = data(37);
		m_signewlast = data(38);

		m_gamarelast = data(39);
		m_omglast = data(40);
		m_omgfullast = data(41);
		m_epspllast = data(42);

		m_bSmallStresslast = data(43);
		m_bUltimatelast = data(44);

		m_trialStrainlast = data(45);
		m_trialStresslast = data(46);

		m_trialTangentlast = data(47);

		m_Eun0last = data(48);

		m_Tstrain = m_trialStrainlast;
		m_Tstress = m_trialStresslast;

		m_trialTangent = m_trialTangentlast;

	}
	opserr << "recvSelf\n";

	return res;
}

void 
FRPConfinedConcrete02::Print(OPS_Stream &s, int flag)
{
	s << "FRPConfinedConcrete02 tag: " << this->getTag() << endln;
	s << "  fc0: " << m_fc0 << endln;
	s << "  Ec: " << m_Ec << endln;
	s << "  ec0: " << m_epsc0 << endln;
	s << "  Ets: " << m_Ets << endln;
	s << "  ft: " << m_ft << endln;

	s << "  stress: " << m_Tstress << " tangent: " << m_trialTangent << endln;
}


void
FRPConfinedConcrete02::Tens_Envlp(double epsc, double &sigc, double &Ect)
{
	if (epsc<=0.0 && epsc>=m_epstn)		//   linear ascending branch between 0 and m_epstn
	{
		sigc = m_Etr1*epsc;
		Ect  = m_Etr1;
	} 
	else if (epsc<=0.0 && epsc>m_epstu)		//   linear descending branch between m_epstn and epstu
	{
		sigc = m_Etr1*m_epstn - m_Ets*(epsc-m_epstn);
		Ect  = -m_Ets;
	}
	else if (epsc<=m_epstu)
	{
		sigc = 0.0;
		Ect  = 1e-15;
	}
	return;
}

void
FRPConfinedConcrete02::Compr_Envlp(double epsc, double &sigc, double &Ect)
{
	m_bSmallStress = false;
	if (epsc>=0.0 && epsc<=m_epst) // Eq.1 (Lam and Teng 2009)
	{
		sigc = m_Ec*epsc - pow((m_Ec-m_E2)*epsc,2)/4.0/m_fc0;
		Ect  = m_Ec - pow(m_Ec-m_E2,2)*epsc/2.0/m_fc0;
	} 
	else if (epsc>m_epst && epsc<=m_epscu) // Eq.2 (Lam and Teng 2009)
	{
		//   linear ascending branch between epst and epscu
		sigc = m_fc0 + m_E2*epsc;
		Ect  = m_E2;
	}
	else if (epsc>m_epscu)
	{
		sigc = 0.0;
		Ect  = 1e-15;
	}

	return;
}

void
FRPConfinedConcrete02::Compr_UnloadingPath(double epsc, double &sigc, double &Ect)
{
	double eta;
	double Eun0,Eun01,Eun02;
	double a,b,c;

	eta = 350*m_epsun + 3;

	Eun01 = Eun02 = m_Ec;

	if (m_epsun !=0)
	{
		Eun01 = 0.5*m_fc0/m_epsun;
	}
	if (m_epsun != m_epspl)
	{
		Eun02 = m_sigun/(m_epsun - m_epspl);
	}

	Eun0 = Eun01<Eun02 ? Eun01:Eun02;  // Eq.13 (Lam and Teng 2009)

	a = (m_sigun-Eun0*(m_epsun-m_epspl))/(pow(m_epsun,eta)-pow(m_epspl,eta)-eta*pow(m_epspl,eta-1)*(m_epsun-m_epspl));
	b = Eun0 - eta*pow(m_epspl,eta-1)*a;
	c = -a*pow(m_epspl,eta) - b*m_epspl;

	sigc = a*pow(epsc,eta) + b*epsc + c; // Eq.8-12 (Lam and Teng 2009)

	Ect = a*eta*pow(epsc,eta-1) + b;
	m_Eun0 = Eun0;
}

void
FRPConfinedConcrete02::Compr_GetPlasticStrain() 
{

	if (m_n == 1) // epspl1 Eq.25 (Lam and Teng 2009)
	{
		if (m_epsunenv > 0 && m_epsunenv <= 0.001)
			m_epspl = 0.0;
		else if (m_epsunenv >= 0.001 && m_epsunenv < 0.0035) //Lu JY Add m_epsunenv >= 0.001 && 2009.06.04
			m_epspl = (1.4*(0.87 - 0.004*m_fc0*m_Unitscale) - 0.64)*(m_epsunenv - 0.001);
		else if (m_epsunenv >= 0.0035 && m_epsunenv <= m_epscu) //Lu JY Add m_epsunenv >= 0.0035 && 2009.06.04
			m_epspl = (0.87 - 0.004*m_fc0*m_Unitscale)*m_epsunenv - 0.0016;
	}
	else if (m_n >= 2) // Eq.34 (Lam and Teng 2009)
	{
		Compr_GetStrainRecoveryRatio();
		m_epspl = (1.0 - m_omg)*m_epsun + m_omg*m_epspl;
	}
}

void
FRPConfinedConcrete02::Compr_GetStrainRecoveryRatio()
{
	if (m_n >= 2)
	{
		if (m_ne == 1)
			m_omgful = 1.0;
		if (m_ne >= 2 && m_ne <=5)
		{
			if (m_epsunenv >= 0.0 && m_epsunenv <= 0.001) // Get Omg_n,ful Eq.33 (Lam and Teng 2009)
				m_omgful = 1.0;
			else if (m_epsunenv > 0.001 && m_epsunenv < 0.0035)
				m_omgful = 1.0 + 400.0*(0.0212*m_ne - 0.12)*(m_epsunenv - 0.001);
			else if (m_epsunenv > 0.0035 && m_epsunenv <= m_epscu) //Lu JY Add if (m_epsunenv > 0.0035) && m_epsunenv <= m_epscu 2009.06.04

				m_omgful = 0.0212*m_ne + 0.88;
		} 
		else if (m_ne >= 6)
			m_omgful = 1.0;

		double omgtemp;
		omgtemp = m_omgful - 0.25*(m_gamare-1);
		m_omg = 1.0<omgtemp ? 1.0:omgtemp; // Get Omgn Eq.32 (Lam and Teng 2009)
	}
}

void
FRPConfinedConcrete02::Compr_ReloadingPath(double epsc, double &sigc, double &Ect)
{
	if (epsc >= m_epsre && epsc <= m_epsref)	// Linear part
	{
		if ((m_epsunenv != m_epsre) && ((m_epsunenv <= 0.001) || ((m_n == 1) && (m_epsunenv > 0.001) && (m_sigre > 0.85*m_sigunenv)) || ((m_n > 1) && (m_epsunenv > 0.001) && (m_sigre > 0.85*m_sigunenv) && (m_epsretenv == m_epsunenv))))	// Eq.17
		{
			m_epsretenv = m_epsunenv; // return to the unloading point on the envelop
			m_Ere = (m_sigunenv - m_sigre)/(m_epsunenv - m_epsre);
			m_bSmallStress = true;
		}
		else if (m_epsref != m_epsre)	// Eq.15-Eq.16 (Lam and Teng 2009)
		{
			m_Ere = (m_signew - m_sigre)/(m_epsref - m_epsre);
			m_bSmallStress = false;
		}
		sigc = m_sigre + m_Ere*(epsc - m_epsre);
		Ect = m_Ere;
	}
	else if (epsc > m_epsref && m_bSmallStress == false)	// Parabolic part Eq.18-Eq.24 (Lam and Teng 2009)
	{	
		double A,B,C;
		if (m_epsunenv < m_epst)
		{

			A = (pow(m_Ec-m_E2,2)*(m_Ere*m_epsref-m_signew) + pow(m_Ec-m_Ere,2)*m_fc0)/(4.0*(m_signew - m_Ec*m_epsref)*m_fc0 + pow((m_Ec-m_E2)*m_epsref,2));
			B = m_Ere - 2.0*A*m_epsref;
			C = m_signew - A*pow(m_epsref,2) - B*m_epsref;

			m_epsretenv = (m_Ec - B)/(2.0*A + pow(m_Ec-m_E2,2)/m_fc0/2.0);
			if (m_epsretenv >= m_epst && (m_signew - m_fc0 - m_E2*m_epsref) != 0)
			{

				A = pow(m_Ere-m_E2,2)/(4.0*(m_signew - m_fc0 - m_E2*m_epsref));
				B = m_Ere - 2.0*A*m_epsref;
				C = m_signew - A*pow(m_epsref,2) - B*m_epsref;

				m_epsretenv = (m_E2 - B)/(2.0*A);			
			} 
		} 
		else if ((m_signew - m_fc0 - m_E2*m_epsref) != 0)
		{
			A = pow(m_Ere-m_E2,2)/(4.0*(m_signew - m_fc0 - m_E2*m_epsref));
			B = m_Ere - 2.0*A*m_epsref;
			C = m_signew - A*pow(m_epsref,2) - B*m_epsref;

			m_epsretenv = (m_E2 - B)/(2.0*A);
		}
		if (epsc > m_epsref && epsc <= m_epsretenv)
		{
			sigc = A*pow(epsc,2) + B*epsc + C;
			Ect = 2*A*epsc + B;
		}
		else  // return to envelope
		{
			m_n = 0;
			Compr_Envlp(epsc,sigc,Ect);
		}
	}
	else if (epsc > m_epsref && m_bSmallStress == true)  // return to envelope
	{
		m_n = 0;
		Compr_Envlp(epsc,sigc,Ect);
	}
}

void
FRPConfinedConcrete02::GetRefPoint()
{
	if (m_n == 1) // Eq.14 (Lam and Teng 2009)
	{
		m_epsref = m_epsunenv;
		m_sigref = m_sigunenv;

		m_epsreflast = m_epsref;
		m_sigreflast = m_sigref;
	}
	else if (m_n >= 2)
	{
		if (m_epsun <= m_epsreflast)
		{
			m_epsref = m_epsreflast;
			m_sigref = m_signew;
		}
		else
		{
			m_epsref = m_epsun;	
			m_sigref = m_sigun;	
		}
		m_epsreflast = m_epsref;
		m_sigreflast = m_sigref;
	}
}

void
FRPConfinedConcrete02::GetDeterioratedStress()
{
	GetStressDeteriorationRatio();
	m_signew = m_fi*m_sigref;  // m_signew Eq.26 (Lam and Teng 2009)
}

void
FRPConfinedConcrete02::GetStressDeteriorationRatio()
{
	if (m_n == 1)  // Get Fi_1 Eq.27 (Lam and Teng 2009)
	{
		if (m_epsunenv>=0 && m_epsunenv<=0.001)
			m_fi = 1.0;
		else if (m_epsunenv>0.001 && m_epsunenv<0.002)
			m_fi = 1.0 - 80*(m_epsunenv - 0.001);
		else
			m_fi = 0.92;
	}
	else if (m_n >= 2) // Get Fi_n Eq.35 (Lam and Teng 2009)
	{
		if (m_ne == 1)
			m_fiful = 1.0;
		if (m_ne >= 2 && m_ne <=5)
		{
			if (m_epsunenv <= 0.001) // Get Fi_n,ful Eq.36 (Lam and Teng 2009)
				m_fiful = 1.0;
			else if (m_epsunenv > 0.001 && m_epsunenv < 0.002)
				m_fiful = 1.0 + 1000.0*(0.013*m_ne - 0.075)*(m_epsunenv - 0.001);
			else
				m_fiful = 0.013*m_ne + 0.925;
		} 
		else if (m_ne >= 6)
			m_fiful = 1.0;

		double fitemp;
		fitemp = m_fiful - 0.2*(m_betaun-1);

		m_fi = 1.0<fitemp ? 1.0:fitemp;
	}
}

// AddingSensitivity:BEGIN ///////////////////////////////////
int
FRPConfinedConcrete02::setParameter(const char **argv, int argc, Information &info)
{
	return -1;
}


int
FRPConfinedConcrete02::updateParameter(int parameterID, Information &info)
{
	return 0;
}


int
FRPConfinedConcrete02::activateParameter(int passedParameterID)
{
	return 0;
}

double
FRPConfinedConcrete02::getStressSensitivity(int gradNumber, bool conditional)
{
	return 0;
}


int
FRPConfinedConcrete02::commitSensitivity(double TstrainSensitivity, int gradNumber, int numGrads)
{
	return 0;
}
// AddingSensitivity:END /////////////////////////////////////////////
