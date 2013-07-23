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
//			May 2013, University of Washington
                                                                      
// Description: This file contains the implementation for the ManzariDafalias class.

#include <ManzariDafalias.h>
#include <ManzariDafalias3D.h>
#include <ManzariDafaliasPlaneStrain.h>

#include <string.h>

const double ManzariDafalias::one3   = 1.0/3.0 ;
const double ManzariDafalias::two3   = 2.0/3.0;
const double ManzariDafalias::root23 = sqrt(2.0/3.0);
const double ManzariDafalias::small  = 1e-10;
double       ManzariDafalias::mElastFlag = 1;

static int numManzariDafaliasMaterials = 0;

void *
OPS_NewManzariDafaliasMaterial(void)
{
  if (numManzariDafaliasMaterials == 0) {
    numManzariDafaliasMaterials++;
    opserr << "ManzariDafalias nDmaterial - Written: P.Arduino, A.Ghofrani, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 19) {
    opserr << "Want: nDMaterial ManzariDafalias tag? G0? nu? e_init? Mc? c? lambda_c? e0? ksi? P_atm? m? h0? Ch? nb? A0? nd? z_max? cz? Rho?" << endln;
    return 0;	
  }
  
  int tag;
  double dData[22];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial ManzariDafalias material  tag" << endln;
    return 0;
  }

  if (numArgs == 19) {
      numData = 18;
  } else if (numArgs == 20) {
      numData = 19;
  } else if (numArgs == 21) {
      numData = 20;
  } else if (numArgs == 22) {
      numData = 21;
  } else {
	  numData = 22;
  }

  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
    return 0;
  }

 if (numArgs == 19) {
	 theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
		 dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], 
		 dData[16], dData[17]);
  } else if (numArgs == 20) {
	 theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
		 dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], 
		 dData[16], dData[17], dData[18]);
  } else if (numArgs == 21) {
	 theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
		 dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], 
		 dData[16], dData[17], dData[18], dData[19]);
  } else if (numArgs == 22) {
	 theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
		 dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], 
		 dData[16], dData[17], dData[18], dData[19], (int)dData[20]);
  } else {
	 theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
		 dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11], dData[12], dData[13], dData[14], dData[15], 
		 dData[16], dData[17], dData[18], dData[19], (int)dData[20], (int)dData[21]);
  }
  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, int classTag, double G0, double nu, 
	double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double z_max, double cz, double mDen, double TolF, double TolR, int JacoType,
	int integrationScheme)
  : NDMaterial(tag,ND_TAG_ManzariDafalias),
    mEpsilon(6), 
	mEpsilon_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
    mSigma(6),
    mSigma_n(6),
	mAlpha(6),
	mAlpha_n(6),
	mAlpha_in(6),
	mFabric(6),
	mFabric_n(6),
    mCe(6,6),
    mCep(6,6),
    mI1(6),
    mIIco(6,6),
	mIIcon(6,6),
	mIImix(6,6),
    mIIvol(6,6),
	mIIdevCon(6,6),
	mIIdevMix(6,6)
{
	m_G0		= G0;
	m_nu		= nu;
	m_e_init	= e_init;
	m_Mc		= Mc;
	m_c			= c;
	m_lambda_c  = lambda_c;
	m_e0		= e0;
	m_ksi		= ksi;
	m_P_atm		= P_atm;
	m_m			= m;
	m_h0		= h0;
	m_ch		= ch;
	m_nb		= nb;
	m_A0		= A0;
	m_nd		= nd;
	m_z_max		= z_max;
	m_cz		= cz;

	massDen     = mDen;
	mTolF		= TolF;
	mTolR		= TolR;
	mJacoType	= JacoType;
	mScheme		= integrationScheme;
	mIter		= 0;
	
	this->initialize();
}

// null constructor
ManzariDafalias ::ManzariDafalias() 
    : NDMaterial(),
    mEpsilon(6), 
	mEpsilon_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
    mSigma(6),
    mSigma_n(6),
	mAlpha(6),
	mAlpha_in(6),
	mFabric(6),
    mCe(6,6),
    mCep(6,6),
    mI1(6),
    mIIco(6,6),
	mIIcon(6,6),
	mIImix(6,6),
    mIIvol(6,6),
	mIIdevCon(6,6),
	mIIdevMix(6,6)
{
	this->initialize();
}

// destructor
ManzariDafalias::~ManzariDafalias()
{
}

NDMaterial*
ManzariDafalias::getCopy(const char *type)
{
	if (strcmp(type,"PlanStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		ManzariDafaliasPlaneStrain *clone;
		clone = new ManzariDafaliasPlaneStrain(this->getTag(), m_G0,  m_nu,  m_e_init,  m_Mc,  m_c, m_lambda_c,  m_e0,  m_ksi,  m_P_atm,  
															   m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz, massDen);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		ManzariDafalias3D *clone;
     	clone = new ManzariDafalias3D(this->getTag(), m_G0,  m_nu,  m_e_init,  m_Mc,  m_c, m_lambda_c,  m_e0,  m_ksi,  m_P_atm,  
															   m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz, massDen);
	 	return clone;
  	} else {
	  	opserr << "ManzariDafalias::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int 
ManzariDafalias::commitState(void)
{
	Vector n(6);
	n = GetDevPart(mSigma) -  one3 * GetTrace(mSigma) * mAlpha;
    if (GetNorm_Contr(n)!=0) n = n / GetNorm_Contr(n);
    if (DoubleDot2_2_Contr(mAlpha - mAlpha_in,n) < 0)
        mAlpha_in = mAlpha;

	mSigma_n    = mSigma;
	mEpsilon_n  = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	mAlpha_n	= mAlpha;
	mFabric_n	= mFabric;
	mVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);
	return 0;
}

int ManzariDafalias::revertToLastCommit (void)
{
    return 0;
}

int ManzariDafalias::revertToStart(void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
    	this->initialize();
	}

    return 0;
}


NDMaterial*
ManzariDafalias::getCopy (void)
{
	opserr << "ManzariDafalias::getCopy -- subclass responsibility\n"; 
  	exit(-1);
  	return 0;
}

const char*
ManzariDafalias::getType (void) const
{
    opserr << "ManzariDafalias::getType -- subclass responsibility\n";
    exit(-1);
    return 0;
}

int
ManzariDafalias::getOrder (void) const
{
    opserr << "ManzariDafalias::getOrder -- subclass responsibility\n";
    exit(-1);
    return 0;
}


Response*
ManzariDafalias::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->getState());
	else
		return 0;
}

int
ManzariDafalias::getResponse(int responseID, Information &matInfo)
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

int
ManzariDafalias::sendSelf(int commitTag, Channel &theChannel)
{
  return 0;
}

int 
ManzariDafalias::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
  return 0;
}

void ManzariDafalias::Print(OPS_Stream &s, int flag )
{
  s << "ManzariDafalias" << endln;
}

int
ManzariDafalias::setParameter(const char **argv, int argc, Parameter &param)
{
  	if (argc < 2)
    	return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0],"updateMaterialStage") == 0) {
			return param.addObject(1, this);
		}
		else if (strcmp(argv[0],"IntegrationScheme") == 0) {
			return param.addObject(2, this);
		}
		else if (strcmp(argv[0],"Jacobian") == 0) {
			return param.addObject(3, this);
		}
	}
    return -1;
}

int
ManzariDafalias::updateParameter(int responseID, Information &info)
{
	// called updateMaterialStage in tcl file
	if (responseID == 1) {
		mElastFlag = info.theInt;
	}
	// called materialState in tcl file
	if (responseID == 5) {
		mElastFlag = info.theDouble;
	}
	if (responseID == 2) {
		mScheme = info.theInt;
	}
	if (responseID == 3) {
		mJacoType = info.theInt;
	}
	return 0;
}


// Initialize Manzari Dafalias Material
void 
ManzariDafalias::initialize()
{
	// strain and stress terms
	mEpsilon.Zero();
	mEpsilon_n.Zero();
	mSigma.Zero();
	mSigma_n.Zero();

	mEpsilonE.Zero();
	mAlpha.Zero();
	mAlpha_n.Zero();
	mAlpha_in.Zero();
	mDGamma = 0.0;
	mFabric.Zero();
	mFabric_n.Zero();
	mVoidRatio = m_e_init;

	mSigma_n(0) = m_P_atm / 200; // P_atm / 200
	mSigma_n(1) = m_P_atm / 200;
	mSigma_n(2) = m_P_atm / 200;
	GetElasticModuli(mSigma_n,mVoidRatio,mVoidRatio,mEpsilon,mEpsilonE,mK,mG);
	mCe = GetStiffness(mK,mG);
	mCep = mCe;

	mEPS = machineEPS();

	initializeState = true;
    
	// 2nd order identity tensor
	mI1.Zero();
	mI1(0) = 1.0;
	mI1(1) = 1.0;
	mI1(2) = 1.0;

	// 4th order mixed variant identity tensor
	mIImix.Zero();
	for (int i = 0; i<6; i++) {
		mIImix(i,i) = 1.0;
	}

	// 4th order covariant identity tensor
	mIIco = mIImix;
	mIIco(3,3) = 2.0;
	mIIco(4,4) = 2.0;
	mIIco(5,5) = 2.0;

	// 4th order contravariant identity tensor
	mIIcon = mIImix;
	mIIcon(3,3) = 0.5;
	mIIcon(4,4) = 0.5;
	mIIcon(5,5) = 0.5;

	// 4th order volumetric tensor, IIvol = I1 tensor I1
	mIIvol.Zero();
	for (int i = 0; i<3; i++) {
		mIIvol(i,0) = 1.0;
		mIIvol(i,1) = 1.0;
		mIIvol(i,2) = 1.0;
	}
	
	// 4th order contravariant deviatoric tensor
	mIIdevCon = mIIcon - one3*mIIvol;

	// 4th order covariant deviatoric tensor
	mIIdevCo  = mIIco  - one3*mIIvol;

	// 4th order mixed variant deviatoric tensor
	mIIdevMix = mIImix - one3*mIIvol;

}

//send back the state parameters
const Vector 
ManzariDafalias::getState()
 {
	 Vector result(26);
	 result.Assemble(mEpsilonE,0,1.0);
	 result.Assemble(mAlpha,6,1.0);
	 result.Assemble(mFabric,12,1.0);
	 result.Assemble(mAlpha_in,18,1.0);
	 result(24) = mVoidRatio;
	 result(25) = mDGamma;

	 return result;
 }

// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void ManzariDafalias::plastic_integrator() 
{	
	Vector CurStress(6);
	Vector CurStrain(6);
	Vector CurElasticStrain(6);
	Vector alpha_in(6);
	double CurVoidRatio;

	Vector NextStress(6);
	Vector NextStrain(6);
	Vector NextElasticStrain(6);
	Vector NextAlpha(6);
	Vector NextFabric(6);
	Matrix aC(6,6);
	Matrix aCep(6,6);
	double NextDGamma, NextVoidRatio;
	
	CurStress	      = mSigma_n;
	CurStrain	      = mEpsilon_n;
	CurElasticStrain  = mEpsilonE_n;
	alpha_in		  = mAlpha_in;
	CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);

	NextStrain	      = mEpsilon;
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextAlpha	      = mAlpha_n;
	NextFabric	      = mFabric_n;
	NextDGamma	      = 0.0;
	NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	//Trial Elastic Stress
	double G, K;
	GetElasticModuli(CurStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
	aC = GetStiffness(K, G);

	if (mElastFlag == 0) {    // Force elastic response
		NextStress = CurStress + DoubleDot4_2(aC,(NextElasticStrain - CurElasticStrain));
		aCep = aC;
	} else {                  // ElastoPlastic response
		if (mScheme == 0){        //Explicit Integration
			Vector n(6), d(6), b(6), dPStrain(6); 
			double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
			GetStateDependent(CurStress, NextAlpha, CurVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
			double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(NextFabric, n)));
			double D = A * DoubleDot2_2_Contr(d, n);
			double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
			double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
			Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
			double dVolStrain = GetTrace(NextStrain - CurStrain);
			Vector dDevStrain = GetDevPart(NextStrain - CurStrain);
			double p = one3 * GetTrace(CurStress);
			Vector r = GetDevPart(CurStress) / p;
			double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
		
			NextDGamma = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/(Kp + 2.0*G*(B-C*
				GetTrace(SingleDot(n,SingleDot(n,n)))) - K*D*DoubleDot2_2_Contr(n,r));
			Vector dSigma = 2.0*G*mIIcon*dDevStrain + K*dVolStrain*mI1 - Macauley(NextDGamma)*(2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
			Vector dAlpha = Macauley(NextDGamma) * two3 * h * b;
			Vector dFabric = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + NextFabric);
			dPStrain = NextDGamma * mIIco * R;

			Matrix temp1 = 2.0*G*mIIdevCo + K*mIIvol;
			Vector temp2 = 2.0*G*n - DoubleDot2_2_Contr(n,r)*mI1;
			Vector temp3 = 2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1;
			double temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
				- K*D*DoubleDot2_2_Contr(n,r));
			aCep = temp1 + MacauleyIndex(NextDGamma) * Dyadic2_2(temp3, temp2) / temp4;
		
			NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
			NextStress = CurStress + dSigma;
			NextAlpha = NextAlpha + dAlpha;
			NextFabric = NextFabric + dFabric;

		} 
		else {                //Implicit Integration
			NextStress = CurStress + DoubleDot4_2(aC,(NextElasticStrain - CurElasticStrain));
			//In case of pure elastic response
			aCep = aC;
			//Trial yield function
			double NextF = GetF(NextStress, NextAlpha);
	
			if (NextF > mTolF)
			{
				Vector Delta0(19), InVariants(44), Delta(19);
				Delta0 = SetManzariComponent(NextStress, NextAlpha, NextFabric, NextDGamma);
				InVariants = SetManzariStateInVar(NextStrain, CurStrain, CurStress, CurElasticStrain, NextAlpha, NextFabric,
					CurVoidRatio, NextVoidRatio, alpha_in);
				Delta = NewtonSolve(Delta0, InVariants, aCep);
		
				NextStress.Extract(Delta, 0, 1.0);
				NextAlpha.Extract(Delta, 6, 1.0);
				NextFabric.Extract(Delta, 12, 1.0);
				NextDGamma = Delta(18);
		
				Vector n(6), d(6), b(6), dPStrain(6); 
				double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
				GetStateDependent(NextStress, NextAlpha, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
				double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(NextFabric, n)));
				double D = A * DoubleDot2_2_Contr(d, n);
				double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
				double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
				Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	
				dPStrain          = NextDGamma * mIIco * R;
				NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
				GetElasticModuli(CurStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
				aC                = GetStiffness(K, G);
			}
		}
	}
	//update State variables
	
	
	mSigma     = NextStress;
	mEpsilonE  = NextElasticStrain;
	mAlpha     = NextAlpha;
	mFabric    = NextFabric;
	mDGamma    = NextDGamma;
	mCe        = aC;
	mCep       = aCep;
	mVoidRatio = NextVoidRatio;
	return;
}
/*************************************************************/
//            NewtonSolve                                    //
/*************************************************************/
Vector
ManzariDafalias::NewtonSolve(const Vector& xo, const Vector& inVar, Matrix& aCepPart)
{
	mIter = 0;
	int MaxIter = 200;
	// Declare variables to be used
	Vector sol(19), R(19), dX(19);
	Matrix  jaco(19,19), jInv(19,19);
	sol = xo;
	R = GetResidual(sol, inVar);
	while ((R.Norm() > mTolR) && (mIter <= MaxIter))
	{
		mIter++;
		if (mJacoType == 1)
			jaco = GetJacobian(sol, inVar);
		else
			jaco = GetFDMJacobian(sol, inVar);
		if (jaco.Invert(jInv) != 0)
			opserr << "ManzariDafalias - Jacobian Matrix is Singular" << endln;
		dX   = -1.0*jInv * R;
		sol += dX;
		R    = GetResidual(sol, inVar);
	}
	aCepPart.Extract(jInv, 0, 0, 1.0);
	return sol;
}
/*************************************************************/
//            GetResidual                                    //
/*************************************************************/
Vector
ManzariDafalias::GetResidual(const Vector& x, const Vector& inVar)
{
	// Inside this function, if head of variable doesn't have 'Cur', it indicate 'Next' State
	// Initialization of Variables to be used
	// Material Parameters

	Vector Res(19);   // Residual Vector
	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
	Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
	Vector fabric(6), curFabric(6);
	double dGamma, curVoidRatio, voidRatio;

	// Next Iteration variable setting
	stress.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	fabric.Extract(x, 12, 1.0);
	dGamma = x(18);

	// Current Iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curFabric.Extract(inVar, 30, 1.0);
	curVoidRatio = inVar[36];
	voidRatio = inVar[37];
	alpha_in.Extract(inVar,38,1.0);
	

	// Hypo elastic Stress calculation
	TrialElasticStrain = curEStrain + (strain - curStrain);
	// State Dependent variables
	Vector n(6), d(6), b(6); 
	double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
	GetStateDependent(stress, alpha, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
	Vector devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	if (p < small)
		p = small;
	double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(fabric, n)));
	double D = A * DoubleDot2_2_Contr(d, n);
	double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
	double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
	Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	Vector aBar(6); aBar = two3 * h * b;
	Vector zBar(6); zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
	double G, K;
	GetElasticModuli(curStress, curVoidRatio, voidRatio, TrialElasticStrain, curEStrain, K, G);
	Matrix De = GetCompliance(K, G);
	Vector dEstrain(6);
	dEstrain = De * (stress - curStress);
	eStrain = curEStrain + dEstrain;

	// Get Residual
	Vector g1(6); Vector g2(6); Vector g3(6); double g4;

	g1 = eStrain - TrialElasticStrain + dGamma * mIIco * R;
	g2 = alpha   - curAlpha     - dGamma * aBar;
	g3 = fabric  - curFabric    - dGamma * zBar;
	g4 = GetF(stress, alpha);

	// Arrange Residual as Row vector
	Res.Assemble(g1,  0, 1.0);
	Res.Assemble(g2,  6, 1.0);
	Res.Assemble(g3, 12, 1.0);
	Res(18) = g4;
	return Res;
}
/*************************************************************/
//            GetJacobian                                     //
/*************************************************************/
Matrix 
ManzariDafalias::GetJacobian(const Vector &x, const Vector &inVar)
{
	// Inside this function, if head of variable doesn't have 'Cur', it indicate 'Next' State
	// Initialization of Variables to be used
	// Material Parameters

	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
	Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
	Vector fabric(6), curFabric(6);
	double dGamma, curVoidRatio, voidRatio;

	// Next Iteration variable setting
	stress.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	fabric.Extract(x, 12, 1.0);
	dGamma = x(18);

	// Current Iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curFabric.Extract(inVar, 30, 1.0);
	curVoidRatio = inVar[36];
	voidRatio = inVar[37];
	alpha_in.Extract(inVar,38,1.0);

	// Hypo elastic Stress calculation
	TrialElasticStrain = curEStrain + (strain - curStrain);
	// State Dependent variables
	Vector n(6), n2(6), d(6), b(6); 
	double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
	GetStateDependent(stress, alpha, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
	n2 = SingleDot(n,n);
	Vector devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	if (p < small)
		p = small;
	Vector r(6); r = devStress - p * alpha;
	double normR = GetNorm_Contr(r);
	double gc = g(Cos3Theta, m_c);
	double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(fabric, n)));
	double D = A * DoubleDot2_2_Contr(d, n);
	double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
	double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
	Vector R(6); R = B * n - C * (n2 - one3 * mI1) + one3 * D * mI1;
	Vector aBar(6); aBar = two3 * h * b;
	Vector zBar(6); zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);

	double G, K;
	GetElasticModuli(curStress, curVoidRatio, voidRatio, TrialElasticStrain, curEStrain, K, G);
	Matrix aD(6,6);	aD = GetCompliance(K, G);

	double AlphaAlphaInDotN;
	if (abs(DoubleDot2_2_Contr(alpha - alpha_in,n)) <= small)
		AlphaAlphaInDotN = small;
	else
		AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in,n);

// Start calculation of Jacobian
	// Differentials of quantities with respect to Sigma
	Matrix dnOverdSigma(6,6), dAbarOverdSigma(6,6), dROverdSigma(6,6), dZbarOverdSigma(6,6);
	Vector dPsiOverdSigma(6), db0OverdSigma(6), dCos3ThetaOverdSigma(6), dAdOverdSigma(6), dhOverdSigma(6),
		dgOverdSigma(6), dAlphaDOverdSigma(6), dCOverdSigma(6), dBOverdSigma(6), dAlphaBOverdSigma(6), dDOverdSigma(6);
	// Differentials of quantities with respect to Alpha
	Matrix dnOverdAlpha(6,6), dAbarOverdAlpha(6,6), dROverdAlpha(6,6), dZbarOverdAlpha(6,6);
	Vector dCos3ThetaOverdAlpha(6), dAdOverdAlpha(6), dhOverdAlpha(6), dgOverdAlpha(6), dAlphaDOverdAlpha(6), 
		dCOverdAlpha(6), dBOverdAlpha(6), dAlphaBOverdAlpha(6), dDOverdAlpha(6);
	// Differentials of quantities with respect to Alpha
	Matrix dZbarOverdFabric(6,6), dROverdFabric(6,6);
	Vector dAdOverdFabric(6), dDOverdFabric(6);

	// d...OverdSigma : Arranged by order of dependence
	dnOverdSigma          = 1.0 / normR * (mIIdevCon - one3*Dyadic2_2(alpha,mI1) - 
		Dyadic2_2(n,n) + one3*DoubleDot2_2_Contr(alpha,n)*Dyadic2_2(n,mI1));
	dPsiOverdSigma        = one3 * m_ksi * m_lambda_c / m_P_atm * pow(p/m_P_atm, m_ksi-1) * mI1;
	db0OverdSigma         = -b0 / (6.0*p) * mI1;

	dCos3ThetaOverdSigma  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, mIIco*dnOverdSigma);
	dAdOverdSigma         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
		DoubleDot2_4(fabric, mIIco*dnOverdSigma);
	dhOverdSigma          = 1.0 / AlphaAlphaInDotN * (db0OverdSigma - 
		h*DoubleDot2_4(alpha-alpha_in, mIIco*dnOverdSigma));

	dgOverdSigma          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdSigma;

	dAlphaDOverdSigma     = m_Mc * exp(m_nd * psi) * (dgOverdSigma + m_nd * gc * dPsiOverdSigma);
	dCOverdSigma          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdSigma;
	dBOverdSigma          = 1.5 * (1.0 - m_c)/m_c * (dgOverdSigma * Cos3Theta + gc * dCos3ThetaOverdSigma);
	dAlphaBOverdSigma     = m_Mc * exp(-1.0*m_nb*psi) * (dgOverdSigma - m_nb * gc * dPsiOverdSigma);

	dDOverdSigma          = dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
		A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, mIIco*dnOverdSigma));
	dAbarOverdSigma       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdSigma) + 
		root23 * h * (Dyadic2_2(n, dAlphaBOverdSigma)+alphaBtheta * dnOverdSigma));

	dROverdSigma          = B * dnOverdSigma + Dyadic2_2(n, dBOverdSigma) - C * 
		(Trans_SingleDot4T_2(dnOverdSigma,n) + SingleDot2_4(n, dnOverdSigma)) -
		Dyadic2_2((n2 - one3 * mI1),dCOverdSigma) + one3 * Dyadic2_2(mI1, dDOverdSigma);
	dZbarOverdSigma       = -1.0 * m_cz * MacauleyIndex(-1.0*D) * 
		(-1.0*Dyadic2_2(m_z_max*n + fabric, dDOverdSigma) - m_z_max * D * dnOverdSigma);

	// d...OverdAlpha : Arranged by order of dependence
	dnOverdAlpha          = p / normR * (Dyadic2_2(n,n) - mIIcon);

	dCos3ThetaOverdAlpha  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, mIIco*dnOverdAlpha);
	dAdOverdAlpha         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
		DoubleDot2_4(fabric, mIIco*dnOverdAlpha);
	dhOverdAlpha          = -1.0*h / AlphaAlphaInDotN * (n + 
		DoubleDot2_4(alpha-alpha_in,mIIco*dnOverdAlpha));

	dgOverdAlpha          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdAlpha;

	dAlphaDOverdAlpha     = m_Mc * exp(m_nd * psi) * dgOverdAlpha;
	dCOverdAlpha          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdAlpha;
	dBOverdAlpha          = 1.5 * (1.0 - m_c)/m_c * (dgOverdAlpha * Cos3Theta + gc * dCos3ThetaOverdAlpha);
	dAlphaBOverdAlpha     = m_Mc * exp(-1.0*m_nb*psi) * dgOverdAlpha;

	dDOverdAlpha          = dAdOverdAlpha * (root23 * alphaDtheta - 
		DoubleDot2_2_Contr(alpha, n)) + A * (root23 * dAlphaDOverdAlpha -
		n - DoubleDot2_4(alpha, mIIco*dnOverdAlpha));
	dAbarOverdAlpha       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdAlpha) +
		root23 * h * (Dyadic2_2(n, dAlphaBOverdAlpha)+alphaBtheta * dnOverdAlpha) - h * mIIcon);

	dROverdAlpha          = B * dnOverdAlpha + Dyadic2_2(n, dBOverdAlpha) - C * 
		(Trans_SingleDot4T_2(dnOverdAlpha,n) + SingleDot2_4(n, dnOverdAlpha)) -
		Dyadic2_2((n2 - one3 * mI1),dCOverdAlpha) + one3 * Dyadic2_2(mI1, dDOverdAlpha);
	dZbarOverdAlpha       = -1.0*m_cz*MacauleyIndex(-1.0*D) * (-1.0*
		Dyadic2_2(m_z_max*n + fabric, dDOverdAlpha) - m_z_max * D * dnOverdAlpha);

	// d...OverdFabric : Arranged by order of dependence
	dAdOverdFabric        = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * n;

	dDOverdFabric         = dAdOverdFabric * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n));
	
	dROverdFabric         = one3 * mIIcon * Dyadic2_2(mI1, dDOverdFabric);

	dZbarOverdFabric      = -1.0*m_cz* MacauleyIndex(-1.0*D) * 
		(-1.0*Dyadic2_2(m_z_max * n + fabric, dDOverdFabric) - D * mIIcon);

	// Derivatives of residuals
	Matrix dR1OverdSigma(6,6), dR2OverdSigma(6,6), dR3OverdSigma(6,6); Vector dR1OverdDGamma(6);
	Matrix dR1OverdAlpha(6,6), dR2OverdAlpha(6,6), dR3OverdAlpha(6,6); Vector dR2OverdDGamma(6);
	Matrix dR1OverdFabric(6,6), dR2OverdFabric(6,6), dR3OverdFabric(6,6); Vector dR3OverdDGamma(6);
	Vector dR4OverdSigma(6), dR4OverdAlpha(6), dR4OverdFabric(6);
	double dR4OverdDGamma;
	
	dR1OverdSigma		= aD + dGamma * mIIco * dROverdSigma * mIIco;
	dR1OverdAlpha		= dGamma * mIIco * dROverdAlpha  * mIIco;
	dR1OverdFabric		= dGamma * mIIco * dROverdFabric * mIIco;
	dR1OverdDGamma		= mIIco * R;

	dR2OverdSigma		=  -1.0*dGamma * dAbarOverdSigma * mIIco;
	dR2OverdAlpha		= mIImix - dGamma * dAbarOverdAlpha * mIIco;
	dR2OverdFabric.Zero();
	dR2OverdDGamma		=  -1.0 *  aBar;

	dR3OverdSigma		=  -1.0*dGamma * dZbarOverdSigma * mIIco;
	dR3OverdAlpha		=  -1.0*dGamma * dZbarOverdAlpha * mIIco;
	dR3OverdFabric		=  mIImix - dGamma * dZbarOverdFabric * mIIco;
	dR3OverdDGamma		=  -1.0 * zBar;

	dR4OverdSigma		= mIIco * (n - one3 * DoubleDot2_2_Contr(devStress/p,n) * mI1);
	dR4OverdAlpha		= -1.0 * mIIco * p * n;
	dR4OverdFabric.Zero();  
	dR4OverdDGamma		= 0;


// End calculation of Jacobian

	// Arrange Jacobian
	Matrix j(19, 19);
	j.Zero();
    j.Assemble(dR1OverdSigma , 0,  0, 1.0);
	j.Assemble(dR1OverdAlpha , 0,  6, 1.0);
	j.Assemble(dR1OverdFabric, 0, 12, 1.0);
	j.Assemble(dR1OverdDGamma, 0, 18, 1.0);
	
    j.Assemble(dR2OverdSigma , 6,  0, 1.0);
	j.Assemble(dR2OverdAlpha , 6,  6, 1.0);
	j.Assemble(dR2OverdFabric, 6, 12, 1.0);
	j.Assemble(dR2OverdDGamma, 6, 18, 1.0);

    j.Assemble(dR3OverdSigma , 12,  0, 1.0);
	j.Assemble(dR3OverdAlpha , 12,  6, 1.0);
	j.Assemble(dR3OverdFabric, 12, 12, 1.0);
	j.Assemble(dR3OverdDGamma, 12, 18, 1.0);

	j.AssembleTranspose(dR4OverdSigma , 18,  0, 1.0);
	j.AssembleTranspose(dR4OverdAlpha , 18,  6, 1.0);
	j.AssembleTranspose(dR4OverdFabric, 18, 12, 1.0);
	j(18,18)          = dR4OverdDGamma;

	return j;
}

/*************************************************************/
//            GetFDMJacobian                                 //
/*************************************************************/
Matrix
ManzariDafalias::GetFDMJacobian(const Vector delta, const Vector inVar)
{
	int aSize = delta.Size();
	Matrix j(aSize, aSize);
	Vector x(aSize), fn(aSize), f(aSize);
	x = delta;
	fn = GetResidual(delta, inVar);
	double temp, h;
	for (int i=0; i < aSize; i++)
	{
		temp = x(i);
		h = sqrt(mEPS);
		if (h == 0.0) h =mEPS;
		x(i) = temp + h;
		h = x(i) - temp;
		f = GetResidual(x, inVar);
		x(i) = temp;
		j.Assemble(((f - fn) / h),0,i,1.0);
	}
	return j;
}

/*************************************************************/
/*************************************************************/
// ------------------------------------------------------------
// SetComponent for Manzari Model
Vector 
ManzariDafalias::SetManzariComponent(const Vector& stress, const Vector& alpha,
							 const Vector& fabric, const double& dGamma)
{
	// flush the all data field
	//mSize = 19;
	Vector result(19);
	result.Assemble(stress, 0);		// Stress
	result.Assemble(alpha, 6);		// Alpha
	result.Assemble(fabric, 12);	// Fabric
	result(18) = dGamma;			// DGamma
	return result;
}
// -------------------------------------------------------------
// Set Manzari State Invariants Component for Manzari Model
Vector 
ManzariDafalias::SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, 
									  const Vector& cEStrain, const Vector& cAlpha, const Vector& cFabric,
									  const double& cVoidRatio, const double& nVoidRatio, const Vector& Alpha_in)
{
	// flush the all data field
	//mSize = 44;
	Vector result(44);
	result.Assemble(nStrain, 0);
	result.Assemble(cStrain, 6);
	result.Assemble(cStress, 12);
	result.Assemble(cEStrain, 18);
	result.Assemble(cAlpha, 24);
	result.Assemble(cFabric, 30);
	result(36) = cVoidRatio;
	result(37) = nVoidRatio;
	result.Assemble(Alpha_in, 38);
	return result;
}
/************************************************************/
/************************************************************/

/*************************************************************/
//             SetMachineEPS                                 //
/*************************************************************/
double
ManzariDafalias::machineEPS()
{
	double eps = 1.0;
    int counter = 1;
    while ( ((double) 1.0 + eps) > ((double) 1.0) )
	{
		eps = eps/2.0;
		counter++;
      }
    eps = eps*2.0;
    counter--;
    return eps;
}
/************************************************************/
/************************************************************/


/*************************************************************/
/*************************************************************/
//            MATERIAL SPECIFIC METHODS                      //
/*************************************************************/
/*************************************************************/
// Macauley() -------------------------------------------------
double ManzariDafalias::Macauley(double x)
{
	// Macauley bracket
	//return ((double)(x > 0.0)) * x;
	if (x >= 0 ) return x;
	else return 0.0;
}
/*************************************************************/
// MacauleyIndex() --------------------------------------------
double ManzariDafalias::MacauleyIndex(double x)
{
	//return ((double)(x > 0.0)) * 1.0;
	if (x >= 0 ) return 1.0;
	else return 0.0;
}
/*************************************************************/
// g() (Manzari) ----------------------------------------------
double 
ManzariDafalias::g(const double cos3theta, const double c)
{
	return 2 * c / ((1 + c) - (1 - c) * cos3theta);
}
/*************************************************************/
// GetF() -----------------------------------------------------
double 
ManzariDafalias::GetF(const Vector& nStress, const Vector& nAlpha)
{
	// Manzari's yield function
	Vector s(6); s = GetDevPart(nStress);
	double p = one3 * GetTrace(nStress);
	s = s - p * nAlpha;
	double f = GetNorm_Contr(s) - root23 * m_m * p;
	return f;
}
/*************************************************************/
// GetPSI() ---------------------------------------------------
double 
ManzariDafalias::GetPSI(const double& e, const double& p)
{
	double ec = m_e0 - m_lambda_c * pow((p / m_P_atm),m_ksi);
	return e - ec;
}
/*************************************************************/
// GetLodeAngle() ---------------------------------------------------
double 
ManzariDafalias::GetLodeAngle(const Vector& n)
// Returns cos(3*theta)
{
	Vector n3 = SingleDot(n,SingleDot(n,n));
	double Cos3Theta = sqrt(6.0) * GetTrace(n3);
	Cos3Theta = Cos3Theta > 1 ? 1 : Cos3Theta;
	Cos3Theta = Cos3Theta < -1 ? -1 : Cos3Theta;
	return Cos3Theta;
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
ManzariDafalias::GetElasticModuli(const Vector& sigma, const double& en, const double& en1,
								const Vector& nEStrain, const Vector& cEStrain, double &K, double &G)
// Calculates G, K
{
	double pn = one3 * GetTrace(sigma), pn1, ken, ken1;
	if (pn <= small)
		pn = small;
	if (GetTrace(nEStrain - cEStrain) < small)
	{
		G = m_G0 * m_P_atm * pow((2.97 - en),2) / (1 + en) * sqrt(pn / m_P_atm);
		K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
	} else {
		ken = pow((2.97 - en),2) / (1+en);
		ken1= pow((2.97 - en1),2) / (1+en1);
		pn1 = pow((sqrt(pn) + 0.5* two3 * (1 + m_nu) / (1 - 2 * m_nu) * m_G0 * sqrt(m_P_atm) * (ken1*GetTrace(nEStrain) - ken*GetTrace(cEStrain))),2);
		K = (pn1-pn) / (GetTrace(nEStrain - cEStrain));
		G = 1.5 * (1 - 2 * m_nu) / (1 + m_nu) * K;
	}
}
/*************************************************************/
// GetStiffness() ---------------------------------------------
Matrix
ManzariDafalias::GetStiffness(const double& K, const double& G)
// returns the stiffness matrix in its covariant-covarint form
{
	Matrix C(6,6);
	double a = K + 4*one3 * G;
	double b = K - 2*one3 * G;
	C(0,0) = C(1,1) = C(2,2) = a;
	C(3,3) = C(4,4) = C(5,5) = G;
	C(0,1) = C(0,2) = C(1,2) = b;
	C(1,0) = C(2,0) = C(2,1) = b;
	return C;
}
/*************************************************************/
// GetCompliance() ---------------------------------------------
Matrix
ManzariDafalias::GetCompliance(const double& K, const double& G)
// returns the compliance matrix in its contravarinat-contravariant form
{
	Matrix D(6,6);
	double a = 1 / (9*K) + 1 / (3*G);
	double b = 1 / (9*K) - 1 / (6*G);
	double c = 1 / G;
	D(0,0) = D(1,1) = D(2,2) = a;
	D(3,3) = D(4,4) = D(5,5) = c;
	D(0,1) = D(0,2) = D(1,2) = b;
	D(1,0) = D(2,0) = D(2,1) = b;
	return D;
}
/*************************************************************/
// GetStateDependent() ----------------------------------------
void 
ManzariDafalias::GetStateDependent( const Vector &stress, const Vector &alpha, const double &e
							, const Vector &alpha_in, Vector &n, Vector &d, Vector &b
							, double &cos3Theta, double &h, double &psi, double &alphaBtheta
							, double &alphaDtheta, double &b0)
{
	Vector devStress(6); devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	if (p <= small)
		p = small;
	Vector r(6); r = devStress - p * alpha;
	if (GetNorm_Contr(r) == 0)
		n = r;
	else
		n = r / GetNorm_Contr(r);
	psi = GetPSI(e, p);
	cos3Theta = GetLodeAngle(n);
	alphaBtheta = g(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
	alphaDtheta = g(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;
	b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);
	d    = root23 * alphaDtheta * n - alpha;
	b    = root23 * alphaBtheta * n - alpha;
	double alphaAlphaInDotN;
	if (abs(DoubleDot2_2_Contr(alpha-alpha_in,n)) <= small)
		alphaAlphaInDotN = small;
	else
		alphaAlphaInDotN = DoubleDot2_2_Contr(alpha-alpha_in,n);
	h = b0 / alphaAlphaInDotN;
}






/*************************************************************/
/*************************************************************/
//            SYMMETRIC TENSOR OPERATIONS                    //
/*************************************************************/
/*************************************************************/
//  GetTrace() ---------------------------------------------
double
ManzariDafalias::GetTrace(const Vector& v) 
// computes the trace of the input argument
{
	if (v.Size() != 6)
		opserr << "\n ERROR! ManzariDafalias::GetTrace requires vector of size(6)!" << endln;

	return (v(0) + v(1) + v(2));
}
/*************************************************************/
//  GetDevPart() ---------------------------------------------
Vector 
ManzariDafalias::GetDevPart(const Vector& aV)
// computes the deviatoric part of the input tensor
{
	if (aV.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::GetDevPart requires vector of size(6)!" << endln;
	Vector result(6);
	double p = GetTrace(aV);
	result = aV;
	result(0) -= one3 * p;
	result(1) -= one3 * p;
	result(2) -= one3 * p;

	return result;
}
/*************************************************************/
//  SingleDot() ---------------------------------------------
Vector 
ManzariDafalias::SingleDot(const Vector& v1, const Vector& v2)
// computes v1.v2, v1 and v2 should be both in their "contravariant" form
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
	opserr << "\n ERROR! ManzariDafalias::SingleDot requires vector of size(6)!" << endln;
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
ManzariDafalias::DoubleDot2_2_Contr(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "contravariant"
{
	double result = 0.0;
	
	if ((v1.Size() != 6) || (v2.Size() != 6))
	opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Contr requires vector of size(6)!" << endln;

	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (double)(i>2) * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Cov() ---------------------------------------
double
ManzariDafalias::DoubleDot2_2_Cov(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "covariant"
{
	double result = 0.0;
	
	if ((v1.Size() != 6) || (v2.Size() != 6))
	opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Cov requires vector of size(6)!" << endln;

	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) - (double)(i>2) * 0.5 * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Mixed() ---------------------------------------
double
ManzariDafalias::DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
{
	double result = 0.0;
	
	if ((v1.Size() != 6) || (v2.Size() != 6))
	opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Mixed requires vector of size(6)!" << endln;

	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// GetNorm_Contr() ---------------------------------------------
double
ManzariDafalias::GetNorm_Contr(const Vector& v)
// computes contravariant (stress-type) norm of input 6x1 tensor
{
	double result=0.0;	
	if (v.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::GetNorm_Contr requires vector of size(6)!" << endln;
	result = DoubleDot2_2_Contr(v,v);

	return sqrt(result);
}
/*************************************************************/
// GetNorm_Cov() ---------------------------------------------
double
ManzariDafalias::GetNorm_Cov(const Vector& v)
// computes covariant (strain-type) norm of input 6x1 tensor
{
	double result=0.0;	
	if (v.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::GetNorm_Cov requires vector of size(6)!" << endln;
	result = DoubleDot2_2_Cov(v,v);

	return sqrt(result);
}
/*************************************************************/
// Dyadic2_2() ---------------------------------------------
Matrix 
ManzariDafalias::Dyadic2_2(const Vector& v1, const Vector& v2)
// computes dyadic product for two vector-storage arguments
// the coordinate form of the result depends on the coordinate form of inputs
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
	opserr << "\n ERROR! ManzariDafalias::Dyadic2_2 requires vector of size(6)!" << endln;
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}
/*************************************************************/
// DoubleDot4_2() ---------------------------------------------
Vector
ManzariDafalias::DoubleDot4_2(const Matrix& m1, const Vector& v1)
// computes doubledot product for matrix-vector arguments
// caution: second coordinate of the matrix should be in opposite variant form of vector
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::DoubleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::DoubleDot4_2 function requires 6-by-6 matrix " << endln;
	Vector result(6);
	result = m1*v1;

	return result;
}
/*************************************************************/
// DoubleDot2_4() ---------------------------------------------
Vector
ManzariDafalias::DoubleDot2_4(const Vector& v1, const Matrix& m1)
// computes doubledot product for matrix-vector arguments
// caution: first coordinate of the matrix should be in opposite 
// variant form of vector
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::DoubleDot2_4 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::DoubleDot2_4 function requires 6-by-6 matrix " << endln;
	Vector result(6);
	result = m1^v1;

	return result;
}
/*************************************************************/
// DoubleDot4_4() ---------------------------------------------
Matrix
ManzariDafalias::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
// caution: second coordinate of the matrix should be in opposite 
// variant form of the first coordinate of second matrix
{
	if ((m1.noCols() != 6) || (m1.noRows() != 6) || (m2.noCols() != 6) || (m2.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::DoubleDot4_4 function requires 6-by-6 matrix " << endln;
	Matrix result(6,6);
	result.Zero();
	result = m1*m2;
	return result;
}
/*************************************************************/
// SingleDot4_2() ---------------------------------------------
Matrix
ManzariDafalias::SingleDot4_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 function requires 6-by-6 matrix " << endln;
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
ManzariDafalias::SingleDot2_4(const Vector& v1, const Matrix& m1)
// computes singledot product for vector-matrix arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::SingleDot2_4 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::SingleDot2_4 function requires 6-by-6 matrix " << endln;
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
ManzariDafalias::Trans_SingleDot4T_2(const Matrix& m1, const Vector& v1)
// computes singledot product for matrix-vector arguments
// caution: this implementation is specific for contravariant forms
{
	if (v1.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 requires vector of size(6)!" << endln;
	if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
	opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 function requires 6-by-6 matrix " << endln;
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
double ManzariDafalias::Det(const Vector& aV)
{
	if (aV.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::Det requires vector of size(6)!" << endln;
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	return( +   aV[0] * aV[1] * aV[2] 
		    + 2*aV[3] * aV[4] * aV[5] 
			-   aV[0] * aV[5] * aV[5] 
			-   aV[2] * aV[3] * aV[3] 
			- 	aV[1] * aV[4] * aV[4]);
}
