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
#include <MaterialResponse.h>

#include <string.h>

const double		ManzariDafalias::one3		= 1.0/3.0 ;
const double		ManzariDafalias::two3		= 2.0/3.0;
const double		ManzariDafalias::root23		= sqrt(2.0/3.0);
const double		ManzariDafalias::small		= 1e-10;
const double		ManzariDafalias::maxStrainInc	= 1e-7;
const bool		ManzariDafalias::debugFlag	= false;
const char unsigned	ManzariDafalias::mMaxSubStep	= 10;
char  unsigned		ManzariDafalias::mElastFlag	= 1;

Vector 			ManzariDafalias::mI1(6);
Matrix  		ManzariDafalias::mIIco(6,6);
Matrix 			ManzariDafalias::mIIcon(6,6);
Matrix 			ManzariDafalias::mIImix(6,6);
Matrix 			ManzariDafalias::mIIvol(6,6);
Matrix 			ManzariDafalias::mIIdevCon(6,6);
Matrix 			ManzariDafalias::mIIdevMix(6,6);
Matrix 			ManzariDafalias::mIIdevCo(6,6);
ManzariDafalias::initTensors ManzariDafalias::initTensorOps;

static int numManzariDafaliasMaterials = 0;

void *
OPS_NewManzariDafaliasMaterial(void)
{
  if (numManzariDafaliasMaterials == 0) {
    numManzariDafaliasMaterials++;
    opserr << "ManzariDafalias nDmaterial - Written: A.Ghofrani, P.Arduino, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 19) {
    opserr << "Want: nDMaterial ManzariDafalias tag? G0? nu? e_init? Mc? c? lambda_c? e0? ksi? P_atm? m? h0? Ch? nb? A0? nd? z_max? cz? Rho? <IntScheme? TanType? JacoType? TolF? TolR?>" << endln;
    return 0;	
  }
  
  int tag;
  double dData[18];
  double oData[5];

  oData[0] = 2;
  oData[1] = 2;
  oData[2] = 1;
  oData[3] = 1.0e-7;
  oData[4] = 1.0e-7;

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial ManzariDafalias material tag" << endln;
    return 0;
  }

  numData = 18;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
    return 0;
  }

  numData = numArgs - 19;
  if (numData != 0)
	if (OPS_GetDouble(&numData, oData) != 0) {
		opserr << "WARNING invalid material data for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
		return 0;
	}

  theMaterial = new ManzariDafalias(tag, 0 , dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
					     dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11],
					     dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], 
					     (int)oData[0],(int)oData[1], (int)oData[2], oData[3], oData[4]);

  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafalias material with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, int classTag, double G0, double nu, 
	double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
	double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
	double z_max, double cz, double mDen, int integrationScheme, int tangentType, 
	int JacoType, double TolF, double TolR): NDMaterial(tag,ND_TAG_ManzariDafalias),
	mEpsilon(6), 
	mEpsilon_n(6),
	mEpsilonE(6),
	mEpsilonE_n(6),
	mSigma(6),
	mSigma_n(6),
	mAlpha(6),
	mAlpha_n(6),
	mAlpha_in(6),
	mAlpha_in_n(6),
	mFabric(6),
	mFabric_n(6),
	mCe(6,6),
	mCep(6,6),
	mCep_Consistent(6,6)
{
	m_G0		= G0;
	m_nu		= nu;
	m_e_init	= e_init;
	m_Mc		= Mc;
	m_c		= c;
	m_lambda_c	= lambda_c;
	m_e0		= e0;
	m_ksi		= ksi;
	m_P_atm		= P_atm;
	m_m		= m;
	m_h0		= h0;
	m_ch		= ch;
	m_nb		= nb;
	m_A0		= A0;
	m_nd		= nd;
	m_z_max		= z_max;
	m_cz		= cz;

	massDen		= mDen;
	mTolF		= TolF;
	mTolR		= TolR;
	mJacoType	= JacoType;
	mScheme		= integrationScheme;
	mTangType	= tangentType;
	mOrgTangType	= tangentType;
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
	mAlpha_n(6),
	mAlpha_in(6),
	mAlpha_in_n(6),
	mFabric(6),
	mFabric_n(6),
	mCe(6,6),
	mCep(6,6),
	mCep_Consistent(6,6)
{
	m_G0		= 0.0;
	m_nu		= 0.0;
	m_e_init	= 0.0;
	m_Mc		= 0.0;
	m_c		= 0.0;
	m_lambda_c	= 0.0;
	m_e0		= 0.0;
	m_ksi		= 0.0;
	m_P_atm		= 0.0;
	m_m		= 0.0;
	m_h0		= 0.0;
	m_ch		= 0.0;
	m_nb		= 0.0;
	m_A0		= 0.0;
	m_nd		= 0.0;
	m_z_max		= 0.0;
	m_cz		= 0.0;

	massDen		= 0.0;
	mTolF		= 1.0e-7;
	mTolR		= 1.0e-7;
	mJacoType	= 1;
	mScheme		= 2;
	mTangType	= 2;
	mOrgTangType	= 2;
	mIter		= 0;
	
	this->initialize();
}

// destructor
ManzariDafalias::~ManzariDafalias()
{
}

NDMaterial*
ManzariDafalias::getCopy(const char *type)
{
	if (strcmp(type,"PlaneStrain2D") == 0 || strcmp(type,"PlaneStrain") == 0) {
		ManzariDafaliasPlaneStrain *clone;
		clone = new ManzariDafaliasPlaneStrain(this->getTag(), m_G0,  m_nu,  m_e_init,  m_Mc,  
				m_c, m_lambda_c,  m_e0,  m_ksi,  m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, 
				m_nd, m_z_max, m_cz, massDen, mScheme, mTangType, mJacoType, mTolF, mTolR);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		ManzariDafalias3D *clone;
     		clone = new ManzariDafalias3D(this->getTag(), m_G0,  m_nu,  m_e_init,  m_Mc,  m_c, m_lambda_c,
				m_e0,  m_ksi,  m_P_atm, m_m, m_h0, m_ch, m_nb, m_A0, m_nd, m_z_max, m_cz, massDen, 
				mScheme, mTangType, mJacoType, mTolF, mTolR);
	 	return clone;
  	} else {
	  	opserr << "ManzariDafalias::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int 
ManzariDafalias::commitState(void)
{
	// update alpha_in if needed
	Vector n(6);
	n = GetNormalToYield(mSigma, mAlpha);
	if (DoubleDot2_2_Contr(mAlpha - mAlpha_in_n,n) < 0) 
		mAlpha_in_n     = mAlpha;

	// update commited state variables
	if ((GetTrace(mSigma) / 3) > (m_P_atm / 5))
		mTangType   = mOrgTangType;
	mSigma_n    = mSigma;
	mEpsilon_n  = mEpsilon;
	mEpsilonE_n = mEpsilonE;
	mAlpha_n	= mAlpha;
	mFabric_n	= mFabric;
	mDGamma_n   = mDGamma;
	mVoidRatio  = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);
	return 0;
}

int ManzariDafalias::revertToLastCommit (void)
{
	// need to be added
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
	else if (strcmp(argv[0], "alpha") == 0 || strcmp(argv[0],"backstressratio") == 0)
		return new MaterialResponse(this, 4, this->getAlpha());
	else if (strcmp(argv[0], "fabric") == 0)
		return new MaterialResponse(this, 5, this->getFabric());
	else if (strcmp(argv[0], "alpha_in") == 0 || strcmp(argv[0],"alphain") == 0)
		return new MaterialResponse(this, 6, this->getAlpha_in());
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
ManzariDafalias::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

    // place data in a vector
	static Vector data(93);

	data(0) = this->getTag();

	data(1)  = m_G0;
	data(2)  = m_nu;
	data(3)  = m_e_init;
	data(4)  = m_Mc;
	data(5)  = m_c;
	data(6)  = m_lambda_c;
	data(7)  = m_e0;
	data(8)  = m_ksi;
	data(9)  = m_P_atm;
	data(10) = m_m;
	data(11) = m_h0;
	data(12) = m_ch;
	data(13) = m_nb;
	data(14) = m_A0;
	data(15) = m_nd;
	data(16) = m_z_max;
	data(17) = m_cz;	
	data(18) = massDen;
	
	data(19) = mTolF;
	data(20) = mTolR;
	data(21) = mJacoType;
	data(22) = mScheme;
	data(23) = mTangType;
	data(24) = mOrgTangType;
	data(25) = mElastFlag;

	data(26) = mEpsilon(0);		data(32) = mEpsilon_n(0);	data(38) = mSigma(0);	data(44) = mSigma_n(0);
	data(27) = mEpsilon(1);		data(33) = mEpsilon_n(1);	data(39) = mSigma(1);	data(45) = mSigma_n(1);
	data(28) = mEpsilon(2);		data(34) = mEpsilon_n(2);	data(40) = mSigma(2);	data(46) = mSigma_n(2);
	data(29) = mEpsilon(3);		data(35) = mEpsilon_n(3);	data(41) = mSigma(3);	data(47) = mSigma_n(3);
	data(30) = mEpsilon(4);		data(36) = mEpsilon_n(4);	data(42) = mSigma(4);	data(48) = mSigma_n(4);
	data(31) = mEpsilon(5);		data(37) = mEpsilon_n(5);	data(43) = mSigma(5);	data(49) = mSigma_n(5);

	data(50) = mEpsilonE(0);	data(56) = mEpsilonE_n(0);	data(62) = mAlpha(0);	data(68) = mAlpha_n(0);
	data(51) = mEpsilonE(1);	data(57) = mEpsilonE_n(1);	data(63) = mAlpha(1);	data(69) = mAlpha_n(1);
	data(52) = mEpsilonE(2);	data(58) = mEpsilonE_n(2);	data(64) = mAlpha(2);	data(70) = mAlpha_n(2);
	data(53) = mEpsilonE(3);	data(59) = mEpsilonE_n(3);	data(65) = mAlpha(3);	data(71) = mAlpha_n(3);
	data(54) = mEpsilonE(4);	data(60) = mEpsilonE_n(4);	data(66) = mAlpha(4);	data(72) = mAlpha_n(4);
	data(55) = mEpsilonE(5);	data(61) = mEpsilonE_n(5);	data(67) = mAlpha(5);	data(73) = mAlpha_n(5);

	data(74) = mFabric(0);		data(80) = mFabric_n(0);	data(86) = mAlpha_in_n(0);
	data(75) = mFabric(1);		data(81) = mFabric_n(1);	data(87) = mAlpha_in_n(1);
	data(76) = mFabric(2);		data(82) = mFabric_n(2);	data(88) = mAlpha_in_n(2);
	data(77) = mFabric(3);		data(83) = mFabric_n(3);	data(89) = mAlpha_in_n(3);
	data(78) = mFabric(4);		data(84) = mFabric_n(4);	data(90) = mAlpha_in_n(4);
	data(79) = mFabric(5);		data(85) = mFabric_n(5);	data(91) = mAlpha_in_n(5);

	data(92) = mDGamma_n;

	
	res = theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: ManzariDafalias::sendSelf - failed to send vector to channel" << endln;
		return -1;
	}
    
    return 0;
}

int 
ManzariDafalias::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
	int res = 0;

	// receive data
	static Vector data(93);
	res = theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		opserr << "WARNING: ManzariDafalias::recvSelf - failed to receive vector from channel" << endln;
		return -1;
	}

    // set member variables
	this->setTag((int)data(0));

	m_G0		= data(1);	
	m_nu		= data(2);
	m_e_init	= data(3);
	m_Mc		= data(4);
	m_c		= data(5);
	m_lambda_c	= data(6);
	m_e0		= data(7);
	m_ksi		= data(8);
	m_P_atm		= data(9);
	m_m		= data(10);
	m_h0		= data(11);
	m_ch		= data(12);
	m_nb		= data(13);
	m_A0		= data(14);
	m_nd		= data(15);
	m_z_max		= data(16);
	m_cz		= data(17);
	massDen		= data(18);

	mTolF		= data(19); 
	mTolR		= data(20); 
	mJacoType	= (int)data(21); 
	mScheme		= (int)data(22); 
	mTangType	= (int)data(23); 
	mOrgTangType	= (int)data(24); 
	mElastFlag	= (int)data(25); 

	mEpsilon(0)	 = data(26);	mEpsilon_n(0)  = data(32); 	mSigma(0) = data(38); 	 mSigma_n(0) = data(44); 
	mEpsilon(1)	 = data(27);	mEpsilon_n(1)  = data(33); 	mSigma(1) = data(39); 	 mSigma_n(1) = data(45); 
	mEpsilon(2)	 = data(28);	mEpsilon_n(2)  = data(34); 	mSigma(2) = data(40); 	 mSigma_n(2) = data(46); 
	mEpsilon(3)	 = data(29);	mEpsilon_n(3)  = data(35); 	mSigma(3) = data(41); 	 mSigma_n(3) = data(47); 
	mEpsilon(4)	 = data(30);	mEpsilon_n(4)  = data(36); 	mSigma(4) = data(42); 	 mSigma_n(4) = data(48); 
	mEpsilon(5)	 = data(31);	mEpsilon_n(5)  = data(37); 	mSigma(5) = data(43); 	 mSigma_n(5) = data(49); 
																	  
	mEpsilonE(0) = data(50);	mEpsilonE_n(0) = data(56); 	mAlpha(0) = data(62); 	 mAlpha_n(0) = data(68); 
	mEpsilonE(1) = data(51);	mEpsilonE_n(1) = data(57); 	mAlpha(1) = data(63); 	 mAlpha_n(1) = data(69); 
	mEpsilonE(2) = data(52);	mEpsilonE_n(2) = data(58); 	mAlpha(2) = data(64); 	 mAlpha_n(2) = data(70); 
	mEpsilonE(3) = data(53);	mEpsilonE_n(3) = data(59); 	mAlpha(3) = data(65); 	 mAlpha_n(3) = data(71); 
	mEpsilonE(4) = data(54);	mEpsilonE_n(4) = data(60); 	mAlpha(4) = data(66); 	 mAlpha_n(4) = data(72); 
	mEpsilonE(5) = data(55);	mEpsilonE_n(5) = data(61); 	mAlpha(5) = data(67); 	 mAlpha_n(5) = data(73); 

	mFabric(0)	 = data(74);	mFabric_n(0)   = data(80); 	mAlpha_in_n(0) = data(86);  
	mFabric(1)	 = data(75);	mFabric_n(1)   = data(81); 	mAlpha_in_n(1) = data(87);  
	mFabric(2)	 = data(76);	mFabric_n(2)   = data(82); 	mAlpha_in_n(2) = data(88);  
	mFabric(3)	 = data(77);	mFabric_n(3)   = data(83); 	mAlpha_in_n(3) = data(89);  
	mFabric(4)	 = data(78);	mFabric_n(4)   = data(84); 	mAlpha_in_n(4) = data(90);  
	mFabric(5)	 = data(79);	mFabric_n(5)   = data(85); 	mAlpha_in_n(5) = data(91);  

	mDGamma_n = data(92);
	mVoidRatio  = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);

	GetElasticModuli(mSigma, mVoidRatio, mVoidRatio, mEpsilonE, mEpsilonE_n, mK, mG);
	mCe  = GetStiffness(mK, mG);
	mCep = mCe;
	mCep_Consistent = mCe;

    return 0;
}

void ManzariDafalias::Print(OPS_Stream &s, int flag )
{
	s << "ManzariDafalias Material, tag: " << this->getTag() << endln;
	s << "Type: " << this->getType() << endln;
}

int
ManzariDafalias::setParameter(const char **argv, int argc, Parameter &param)
{
  	if (argc < 2)
    	return -1;

	int theMaterialTag;
	theMaterialTag = atoi(argv[1]);

	if (theMaterialTag == this->getTag()) {

		if (strcmp(argv[0],"updateMaterialStage") == 0) {     // enforce elastic/elastoplastic response
			return param.addObject(1, this);
		}
        else if (strcmp(argv[0],"materialState") == 0) {     // enforce elastic/elastoplastic response
            return param.addObject(5, this);
        }
		else if (strcmp(argv[0],"IntegrationScheme") == 0) { // change integration scheme (Explicit/Implicit)
			return param.addObject(2, this);
		}
		else if (strcmp(argv[0],"Jacobian") == 0) {          // change type of Jacobian used for newton iterations
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
	else {
        	return -1;
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
	mAlpha_in_n.Zero();
	mDGamma = 0.0;
	mFabric.Zero();
	mFabric_n.Zero();
	mVoidRatio = m_e_init;

	// calculate initial stiffness parameters
	GetElasticModuli(mSigma_n,mVoidRatio,mVoidRatio,mEpsilon,mEpsilonE,mK,mG);
	mCe = GetStiffness(mK,mG);
	mCep = mCe;
	mCep_Consistent = mCe;

	// calculate machine epsilon (used for FDM Jacobian)
	mEPS = machineEPS();
}

//send back the state parameters to the recorders
const Vector 
ManzariDafalias::getState()
 {
	 Vector result(26);
	 result.Assemble(mEpsilonE,0,1.0);
	 result.Assemble(mAlpha,6,1.0);
	 result.Assemble(mFabric,12,1.0);
	 result.Assemble(mAlpha_in_n,18,1.0);
	 result(24) = mVoidRatio;
	 result(25) = mDGamma;

	 return result;
 }
//send back alpha tensor
const Vector
ManzariDafalias::getAlpha()
{
	return mAlpha;
}
//send back fabric tensor
const Vector
ManzariDafalias::getFabric()
{
	return mFabric;
}
//send back alpha_in tensor
const Vector
ManzariDafalias::getAlpha_in()
{
	return mAlpha_in_n;
}

// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void ManzariDafalias::integrate() 
{
	mAlpha_in = mAlpha_in_n;

	// Force elastic response
	if (mElastFlag == 0) {
		elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha, 
				mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent);
	} 
	// ElastoPlastic response
	else {                  
		//Explicit Integration
		if (mScheme == 0){
			explicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in_n,
				mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mAlpha_in, mDGamma, mVoidRatio, mG, 
				mK, mCe, mCep, mCep_Consistent);
		} 
		else if (mScheme == 3) { // Constrain strain increment to 1.0e-7
			explicit_aux(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in_n,
				mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mAlpha_in, mDGamma, mVoidRatio, mG, 
				mK, mCe, mCep, mCep_Consistent);
		}
		//Implicit Integration
		else {
			implicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in_n,
				mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mAlpha_in, mDGamma, mVoidRatio, mG, 
				mK, mCe, mCep, mCep_Consistent);
		}
	}
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Elastic Integrator
/*************************************************************/
void ManzariDafalias::elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
		double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{	
	double CurVoidRatio;
	
	// calculate elastic response
	CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	GetElasticModuli(CurStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
	aC = GetStiffness(K, G);
	NextStress = CurStress + DoubleDot4_2(aC,(NextElasticStrain - CurElasticStrain));
	aCep = aC;
	aCep_Consistent = aC;

	//update State variables
	NextStress         = NextStress;
	NextElasticStrain  = NextElasticStrain;
	if (GetTrace(NextStress) > small)
		NextAlpha	   = 3.0 * GetDevPart(NextStress)/GetTrace(NextStress);
	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Explicit Integrator (Constrained strain increment)
/*************************************************************/
void ManzariDafalias::explicit_aux(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, Vector& NextAlpha_in, 
		double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{	
	NextDGamma = 0;

	Vector StrainInc(6); StrainInc = NextStrain - CurStrain;
	double maxInc = StrainInc(0);
	for(int ii=1; ii<6; ii++)
		if(fabs(StrainInc(ii)) > fabs(maxInc)) 
			maxInc = StrainInc(ii);
	if (fabs(maxInc) > maxStrainInc){
		int numSteps = (int)floor(fabs(maxInc) / maxStrainInc) + 1;
		StrainInc = (NextStrain - CurStrain) / numSteps;	
	
		Vector cStress(6), cStrain(6), cAlpha(6), cFabric(6), cAlpha_in(6), cEStrain(6);
		Vector nStrain(6) ,nEStrain(6), nStress(6), nAlpha(6), nFabric(6), nAlpha_in(6);
		Matrix nCe(6,6), nCep(6,6), nCepC(6,6);
		double nDGamma, nVoidRatio, nG, nK;
				
		// create temporary variables
		cStress = CurStress; cStrain = CurStrain; cAlpha = CurAlpha; cFabric = CurFabric;
		cAlpha_in = alpha_in; cEStrain = CurElasticStrain;

		
		for(int ii=1; ii <= numSteps; ii++)
		{
			nStrain = cStrain + StrainInc;

			explicit_integrator(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
			nEStrain, nStress, nAlpha, nFabric, nAlpha_in, nDGamma, nVoidRatio, nG, nK, nCe, nCep, nCepC);

			cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; cAlpha_in = nAlpha_in;
		}

		NextElasticStrain	= nEStrain;
		NextStress		= nStress;
		NextAlpha		= nAlpha;
		NextFabric		= nFabric;
		NextDGamma		= nDGamma;
		NextAlpha_in		= nAlpha_in;
		NextVoidRatio		= nVoidRatio;
		G = nG; K = nK;

		Vector n(6), d(6), b(6), dPStrain(6); 
		double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
		GetStateDependent(NextStress, NextAlpha, NextVoidRatio, NextAlpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, 
				alphaDtheta, b0);
		double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(NextFabric, n)));
		double D = A * DoubleDot2_2_Contr(d, n);
		double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
		double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
		Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	
		dPStrain     = CurElasticStrain + (NextStrain - CurStrain) - NextElasticStrain;
		NextDGamma   = dPStrain.Norm() / R.Norm();

		aC    = nCe;
		aCep  = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
		aCep_Consistent = aCep;

	} else {
		explicit_integrator(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain,
			NextElasticStrain, NextStress, NextAlpha, NextFabric, NextAlpha_in, NextDGamma, NextVoidRatio, 
			G, K, aC, aCep, aCep_Consistent);
	}
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Explicit Integrator
/*************************************************************/
void ManzariDafalias::explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, Vector& NextAlpha_in, 
		double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{	
	double CurVoidRatio, dn;
	
	CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextAlpha_in      = alpha_in;

	GetElasticModuli(CurStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
	aC = GetStiffness(K, G);
	Vector n(6), d(6), b(6), dPStrain(6); 
	double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
	GetStateDependent(CurStress, CurAlpha, CurVoidRatio, NextAlpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
	double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(CurFabric, n)));
	dn = DoubleDot2_2_Contr(d, n);
	if ((psi>0)&&(dn<0)) dn = 0; // check dialtion if the material is loose (this should not happen)
	double D = A * dn;
	double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
	double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
	Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	double dVolStrain = GetTrace(NextStrain - CurStrain);
	Vector dDevStrain = GetDevPart(NextStrain - CurStrain);
	double p = one3 * GetTrace(CurStress);
	p = p < small ? small : p;
	Vector r = GetDevPart(CurStress) / p;
	double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
	
	double temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
		- K*D*DoubleDot2_2_Contr(n,r));
	if (fabs(temp4) < small) temp4 = Sgn(temp4) * small;

	NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
	Vector dSigma   = 2.0*G*mIIcon*dDevStrain + K*dVolStrain*mI1 - Macauley(NextDGamma)*
			  (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
	Vector dAlpha   = Macauley(NextDGamma) * two3 * h * b;
	Vector dFabric  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric);
	       dPStrain = NextDGamma * mIIco * R;

	Matrix temp1 = 2.0*G*mIIdevMix + K*mIIvol;
	Vector temp2 = 2.0*G*n - DoubleDot2_2_Contr(n,r)*mI1;
	Vector temp3 = 2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1;

	aCep = temp1 - MacauleyIndex(NextDGamma) * Dyadic2_2(temp3, temp2) / temp4;
	aCep_Consistent = aCep;

	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
	NextStress = CurStress + dSigma;
	NextAlpha  = CurAlpha  + dAlpha;
	NextFabric = CurFabric + dFabric;
	
	//if (DoubleDot2_2_Contr(NextAlpha-alpha_in,n)<0)
		//NextAlpha_in = CurAlpha;

	return;
}
// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Implicit Integrator
/*************************************************************/
int ManzariDafalias::implicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
		const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
		Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, Vector& NextAlpha_in,
		double& NextDGamma,	double& NextVoidRatio,  double& G, double& K, Matrix& Ce, Matrix& Cep, Matrix& Cep_Consistent, int implicitLevel) 
{
	// check if max number of substepping is reached
	if (implicitLevel > mMaxSubStep){
		if (debugFlag) opserr << " !!! MORE THAN 10 LEVELS OF SUBSTEPPING !!!" << endln;
		return -3;
	}

	int errFlag = 0, SchemeControl = 1;
	Vector TrialStress(6);
	Matrix aC(6,6), aCep(6,6), aCepConsistent(6,6);
	double CurVoidRatio;

	CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
	NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
	// elastic trial strain
	NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	NextAlpha         = CurAlpha;
	NextFabric        = CurFabric;
	NextDGamma        = 0.0;

	// elastic trial stress
	GetElasticModuli(CurStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
	aC = GetStiffness(K, G);
	TrialStress = CurStress + DoubleDot4_2(aC,(NextElasticStrain - CurElasticStrain));

	//In case of pure elastic response
	NextStress = TrialStress;
	aCep = aC;
	aCepConsistent = aC;

	//Trial yield function
	double NextF = GetF(NextStress, NextAlpha);
	
	if (NextF > mTolF) // elastoplastic response
	{
		Vector Delta0(19), InVariants(44), Delta(19);
		Delta0 = SetManzariComponent(NextStress, NextAlpha, NextFabric, NextDGamma);
		InVariants = SetManzariStateInVar(NextStrain, CurStrain, CurStress, CurElasticStrain, CurAlpha, CurFabric,
						CurVoidRatio, NextVoidRatio, NextAlpha_in);

		// do newton iterations
		errFlag = NewtonSolve(Delta0, InVariants, Delta, aCepConsistent);
		
		// check if newton converged
		if (errFlag == 0)
		{
			NextStress.Extract(Delta, 0, 1.0);
			NextAlpha.Extract(Delta, 6, 1.0);
			NextFabric.Extract(Delta, 12, 1.0);
			NextDGamma = Delta(18);
			// check if the results are acceptible
			errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
		}

		// try the considerations to get an acceptible solution
		if (mScheme == 2)
		{
			// try different approaches and continue until a solution is found
			while(errFlag != 0)
			{
				if (errFlag == -4) SchemeControl = 4;
				Vector StrainInc(6), cStress(6), cStrain(6), cAlpha(6), cFabric(6), cAlpha_in(6), cEStrain(6);
				Vector nStrain(6) ,nEStrain(6), nStress(6), nAlpha(6), nFabric(6), nAlpha_in(6);
				Matrix nCe(6,6), nCep(6,6), nCepC(6,6);
				double nDGamma, nVoidRatio, nG, nK;
				int numSteps;

				// original strain increment
				StrainInc = NextStrain - CurStrain;
				// create temporary variables
				cStress = CurStress; cStrain = CurStrain; cAlpha = CurAlpha; cFabric = CurFabric;
				cAlpha_in = NextAlpha_in; cEStrain = CurElasticStrain;

				if (SchemeControl < 3)
				{
					// use numSteps steps of explicit integration as an initial guess to the newton iteration
					if (SchemeControl == 1)
					{
						numSteps = 10;
						if (debugFlag) opserr << "******* Explicit 10" << endln;

					} else if (SchemeControl == 2)
					{
						numSteps = 50;
						if (debugFlag) opserr << "******* Explicit 50" << endln;
					}
					for(int ii=1; ii <= numSteps; ii++)
					{
						nStrain = cStrain + StrainInc / numSteps;
						explicit_integrator(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain,
								nEStrain, nStress, nAlpha, nFabric, nAlpha_in, nDGamma, nVoidRatio, 
								nG, nK, nCe, nCep, nCepC);
						cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; 
						cAlpha_in = nAlpha_in;
					}
					// do newton iterations
					Delta0 = SetManzariComponent(nStress, nAlpha,nFabric, nDGamma);
					errFlag = NewtonSolve(Delta0, InVariants, Delta, aCepConsistent);
					// check if newton converged
					if (errFlag == 0) 
					{
						NextStress.Extract(Delta, 0, 1.0);
						NextAlpha.Extract(Delta, 6, 1.0);
						NextFabric.Extract(Delta, 12, 1.0);
						NextDGamma = Delta(18);
						//NextAlpha_in = nAlpha_in;

						// check validity of the solution
						errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
						// if not valid do implicit substepping
						if (errFlag != 0)
						{
							SchemeControl = 3;
							continue;
						}
					} else {
						SchemeControl += 1;
						continue;
					}
				// explicit guess did not work, now do implicit substepping
				} else if (SchemeControl == 3)
				{
					if (debugFlag) opserr << "******* Implicit Substepping" << endln;
					implicitLevel++;
					nStrain = cStrain + StrainInc / 2;
					// do a recursive implicit_integrator on the first half of the strain increment
					errFlag = implicit_integrator(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain,
								nEStrain, nStress, nAlpha, nFabric, nAlpha_in, nDGamma, nVoidRatio, 
								nG, nK,nCe, nCep, nCepC, implicitLevel);
					// check if maximum number of substepping levels is reached
					if (errFlag == -3){
						// do explicit integration over this strain increment
						SchemeControl += 1;
						continue;
					}

					cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; 
					cAlpha_in = nAlpha_in;
						
					nStrain = cStrain + StrainInc / 2;
					// do a recursive implicit_integrator on the second half of the strain increment
					errFlag = implicit_integrator(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
						nEStrain, nStress, nAlpha, nFabric, nAlpha_in, nDGamma, nVoidRatio, nG, nK,nCe, nCep, 
						nCepC, implicitLevel);
					// check if maximum number of substepping levels is reached
					if (errFlag == -3){
						// do explicit integration over this strain increment
						SchemeControl += 1;
						continue;
					}

					// Update results from substepping
					if (errFlag == 0)
					{
						NextStress = nStress;
						NextAlpha  = nAlpha;
						NextFabric = nFabric;
						NextDGamma = nDGamma;
						//NextAlpha_in = nAlpha_in;
						aC = nCe;
						aCep = nCep;
						aCepConsistent = nCepC;
						
						/*
						Delta0 = SetManzariComponent(nStress, nAlpha,nFabric, nDGamma);
						errFlag = NewtonSolve(Delta0, InVariants, Delta, aCep);
						if (errFlag == 0) 
						{
							NextStress.Extract(Delta, 0, 1.0);
							NextAlpha.Extract(Delta, 6, 1.0);
							NextFabric.Extract(Delta, 12, 1.0);
							NextDGamma = Delta(18);
							NextAlpha_in = nAlpha_in;

							//errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
							if (errFlag != 0)
							{
								SchemeControl = 3;
								continue;
							}
						} else {
							SchemeControl += 1;
							continue;
						}
						*/

					} 
					else {
						SchemeControl += 1;
						continue;
					}

				} 
				// do explicit integration on this strain increment
				else {
					if (debugFlag) opserr << "******* 1000 Exlpicit steps" << endln;
					numSteps = 1000;
					for(int ii=1; ii <= numSteps; ii++)
					{
						nStrain = cStrain + StrainInc / numSteps;
						explicit_integrator(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
								nEStrain, nStress, nAlpha, nFabric, nAlpha_in,nDGamma, nVoidRatio, 
								nG, nK, nCe, nCep, nCepC);

						cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; 
						cAlpha_in = nAlpha_in;
					}

					// update results from explicit integration
					NextStress = nStress;
					NextAlpha  = nAlpha;
					NextFabric = nFabric;
					NextDGamma = nDGamma;
					//NextAlpha_in = nAlpha_in;
					aC = nCe;
					aCep = nCe;
					aCepConsistent = nCe;

					errFlag = 0;
				}
			}
		}


		Vector n(6), d(6), b(6), dPStrain(6); 
		double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
		GetStateDependent(NextStress, NextAlpha, NextVoidRatio, NextAlpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, 
				alphaDtheta, b0);
		double A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(NextFabric, n)));
		double D = A * DoubleDot2_2_Contr(d, n);
		double B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(Cos3Theta, m_c) * Cos3Theta;
		double C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(Cos3Theta, m_c);
		Vector R(6); R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	
		dPStrain          = NextDGamma * mIIco * R;
		NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
		GetElasticModuli(NextStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
		aC                = GetStiffness(K, G);
		aCep              = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
		//if ((DoubleDot2_2_Contr(NextAlpha-alpha_in,n)<0)&&(implicitLevel>1))
			//NextAlpha_in = CurAlpha;
	}

	Ce = aC;
	Cep = aCep;
	Cep_Consistent = aCepConsistent;

	return errFlag;
}

/*************************************************************/
//            NewtonSolve                                    //
/*************************************************************/
int
ManzariDafalias::NewtonSolve(const Vector& xo, const Vector& inVar, Vector& x, Matrix& aCepPart)
{
	int MaxIter = 20;
	int errFlag = -1;
	bool jacoFlag = true;
	// Declare variables to be used
	static Vector sol(19);
	static Vector R(19);
	static Vector dX(19);
	static Vector norms(20);
	static Matrix jaco(19,19);
	static Matrix jInv(19,19);

	sol = xo;
	R = GetResidual(sol, inVar);
	for(mIter = 1; mIter <= MaxIter; mIter++)
	{
		if (jacoFlag)
			if (mJacoType == 1)
				jaco = GetJacobian(sol, inVar);
			else
				jaco = GetFDMJacobian(sol, inVar);
		else
		{	
			static Vector aux(19);
			aux = SetManzariComponent(mSigma_n, mAlpha_n, mFabric_n, mDGamma_n);
			if (mJacoType == 1)
				jaco = GetJacobian(aux, inVar);
			else
				jaco = GetFDMJacobian(aux, inVar);
		}

		/*	
		norms = NormalizeJacobian(jaco);
		if (jaco.Invert(jInv) != 0) 
		{
			if (jacoFlag)
			{
				if (debugFlag) opserr << "Try last commited Jacobian" << endln;
				//if ((GetTrace(mSigma)/3.0) < (m_P_atm / 20))
				//	mTangType = 0;
				jacoFlag = false;
				break;
			} else {
				errFlag = -4;
				if (debugFlag) opserr << "Last comitted Jacobian didn't work!" << endln;
				break;
			}
		}
		DenormalizeJacobian(jInv, norms);
		dX   = -1.0*jInv * R;
		sol += dX;
		*/
		
		errFlag = jaco.Solve(R, dX);
		if (errFlag != 0) 
		{
			errFlag = -4;
			break;
		}
		sol -= dX;
		
		R    = GetResidual(sol, inVar);

		if (R.Norm() < mTolR) 
		{
			errFlag = 0;
			break;
		}
	}
	if (errFlag == 0)
	{	
		
		norms = NormalizeJacobian(jaco);
		if (jaco.Invert(jInv) != 0) 
		{
			errFlag = -4;
			if (debugFlag) opserr << "Singular Matrix!!! - Jacobian" << endln;
		}
		DenormalizeJacobian(jInv, norms);
		
		aCepPart.Extract(jInv, 0, 0, 1.0);
		x = sol;
	}

	return errFlag;
}
/*************************************************************/
//            GetResidual                                    //
/*************************************************************/
Vector
ManzariDafalias::GetResidual(const Vector& x, const Vector& inVar)
{
	// stress: NextStress, also for other variables

	Vector Res(19);   // Residual Vector
	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
	Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6); // Stress and Hardening
	Vector fabric(6), curFabric(6); // Fabric
	double dGamma, curVoidRatio, voidRatio;

	// read the trial variables from newton iterations
	stress.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	fabric.Extract(x, 12, 1.0);
	dGamma = x(18);

	// current iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curFabric.Extract(inVar, 30, 1.0);
	curVoidRatio = inVar[36];
	voidRatio = inVar[37];
	alpha_in.Extract(inVar,38,1.0);
	

	// elastic trial strain
	TrialElasticStrain = curEStrain + (strain - curStrain);
	
	// state dependent variables
	Vector n(6), d(6), b(6); 
	double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
	GetStateDependent(stress, alpha, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
	Vector devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	p = p < small ? small : p;
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

	// residuals
	Vector g1(6); Vector g2(6); Vector g3(6); double g4;

	g1 = eStrain - TrialElasticStrain + dGamma * mIIco * R;
	g2 = alpha   - curAlpha			  - dGamma * aBar;
	g3 = fabric  - curFabric		  - dGamma * zBar;
	g4 = GetF(stress, alpha);

	// put residuals in a one vector
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
	// stress: NextStress, also for other variables

	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
	Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
	Vector fabric(6), curFabric(6);
	double dGamma, curVoidRatio, voidRatio;
	
	// read the trial variables from newton iterations
	stress.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	fabric.Extract(x, 12, 1.0);
	dGamma = x(18);
	
	// current iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curFabric.Extract(inVar, 30, 1.0);
	curVoidRatio = inVar[36];
	voidRatio = inVar[37];
	alpha_in.Extract(inVar,38,1.0);

	// elastic trial strain
	TrialElasticStrain = curEStrain + (strain - curStrain);
	
	// state dependent variables
	Vector n(6), n2(6), d(6), b(6); 
	double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
	GetStateDependent(stress, alpha, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0);
	n2 = SingleDot(n,n);
	Vector devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	p = p < small ? small : p;
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
	if (fabs(DoubleDot2_2_Contr(alpha - alpha_in,n)) <= small)
		AlphaAlphaInDotN = Sgn(AlphaAlphaInDotN) * small;
	else
		AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in,n);

// analytical Jacobian
	// Differentials of quantities with respect to Sigma
	Matrix dnOverdSigma(6,6), dAbarOverdSigma(6,6), dROverdSigma(6,6), dZbarOverdSigma(6,6);
	Vector dPsiOverdSigma(6), db0OverdSigma(6), dCos3ThetaOverdSigma(6), dAdOverdSigma(6), dhOverdSigma(6),
		dgOverdSigma(6), dAlphaDOverdSigma(6), dCOverdSigma(6), dBOverdSigma(6), dAlphaBOverdSigma(6), dDOverdSigma(6);
	// Differentials of quantities with respect to Alpha
	Matrix dnOverdAlpha(6,6), dAbarOverdAlpha(6,6), dROverdAlpha(6,6), dZbarOverdAlpha(6,6);
	Vector dCos3ThetaOverdAlpha(6), dAdOverdAlpha(6), dhOverdAlpha(6), dgOverdAlpha(6), dAlphaDOverdAlpha(6), 
		dCOverdAlpha(6), dBOverdAlpha(6), dAlphaBOverdAlpha(6), dDOverdAlpha(6);
	// Differentials of quantities with respect to Fabric
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
ManzariDafalias::SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, const Vector& cEStrain, 
				const Vector& cAlpha, const Vector& cFabric, const double& cVoidRatio, const double& nVoidRatio, 
				const Vector& Alpha_in)
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
//------------------------------------------------------------
Vector
ManzariDafalias::NormalizeJacobian(Matrix& Jaco)
{
	Vector norms(20);
	for (int i = 0; i < 19; i++){
		norms(i) = 0.5 * (MatrixMax_Rows(Jaco, i)+MatrixMin_Rows(Jaco, i));
		if (fabs(norms(i)) < small)
			norms(i) = 1.0;
		for (int j = 0; j < 19; j++)
			Jaco(i,j) /= norms(i);
	}
	norms(19) = 0.5 * (MatrixMax_Cols(Jaco, 18)+MatrixMin_Cols(Jaco, 18));
	if (fabs(norms(19)) < small)
		norms(19) = 1.0;
	for (int j = 0; j < 19; j++)
		Jaco(j,18) /= norms(19);
	return norms;
}
//------------------------------------------------------------
void
ManzariDafalias::DenormalizeJacobian(Matrix& JInv, const Vector& norms){
	for(int i = 0; i < 19; i++){
		
		for (int j = 0; j < 19; j++)
			JInv(i,j) /= norms(j);

		JInv(18,i) /= norms(19);
	}
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
	return (x > 0 ? x : 0.0);
}
/*************************************************************/
// MacauleyIndex() --------------------------------------------
double ManzariDafalias::MacauleyIndex(double x)
{
	// Macauley index
	return (x > 0 ? 1.0 : 0.0);
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
	return GetNorm_Contr(s) - root23 * m_m * p;
}
/*************************************************************/
// GetPSI() ---------------------------------------------------
double 
ManzariDafalias::GetPSI(const double& e, const double& p)
{
	return e - (m_e0 - m_lambda_c * pow((p / m_P_atm),m_ksi));
}
/*************************************************************/
// GetLodeAngle() ---------------------------------------------------
double 
ManzariDafalias::GetLodeAngle(const Vector& n)
// Returns cos(3*theta)
{
	double Cos3Theta = sqrt(6.0) * GetTrace(SingleDot(n,SingleDot(n,n)));
	Cos3Theta = Cos3Theta > 1 ? 1 : Cos3Theta;
	Cos3Theta = Cos3Theta < -1 ? -1 : Cos3Theta;
	return Cos3Theta;
}
/*************************************************************/
// GetElasticModuli() ---------------------------------------------
void
ManzariDafalias::GetElasticModuli(const Vector& sigma, const double& en, const double& en1, const Vector& nEStrain, 
				const Vector& cEStrain, double &K, double &G)
// Calculates G, K
{
	double pn = one3 * GetTrace(sigma);
	pn = (pn <= m_P_atm / 20.0) ? m_P_atm / 20.0 : pn;

	// this part could make problems
	/*
	//if (fabs(GetTrace(nEStrain - cEStrain)) < small)
	if (fabs(en1 - en) < small)
	{
		G = m_G0 * m_P_atm * pow((2.97 - en),2) / (1 + en) * sqrt(pn / m_P_atm);
		K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
	} else {
		double ken = pow((2.97 - en),2) / (1+en);
		double ken1= pow((2.97 - en1),2) / (1+en1);
		double pn1 = pow((sqrt(pn) + 0.5* two3 * (1 + m_nu) / (1 - 2 * m_nu) * m_G0 * sqrt(m_P_atm) * (ken1*GetTrace(nEStrain) - ken*GetTrace(cEStrain))),2);
		K = (pn1-pn) / (GetTrace(nEStrain - cEStrain));
		G = 1.5 * (1 - 2 * m_nu) / (1 + m_nu) * K;
	}
	*/
	G = m_G0 * m_P_atm * pow((2.97 - en),2) / (1 + en) * sqrt(pn / m_P_atm);
	K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
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
// GetElastoPlasticTangent()---------------------------------------
Matrix
ManzariDafalias::GetElastoPlasticTangent(const Vector& NextStress, const double& NextDGamma, 
					const Vector& CurStrain, const Vector& NextStrain,
					const double& G, const double& K, const double& B, 
					const double& C,const double& D, const double& h, 
					const Vector& n, const Vector& d, const Vector& b) 
{	
	// double dVolStrain = GetTrace(NextStrain - CurStrain);
	Vector dDevStrain = GetDevPart(NextStrain - CurStrain);
	double p = one3 * GetTrace(NextStress);
	p = (p < small) ? small : p;
	Vector r = GetDevPart(NextStress) / p;
	double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);

	Matrix temp1 = 2.0*G*mIIdevMix + K*mIIvol;
	Vector temp2 = 2.0*G*n - DoubleDot2_2_Contr(n,r)*mI1;
	Vector temp3 = 2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1;
	double temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
		- K*D*DoubleDot2_2_Contr(n,r));
	if (fabs(temp4) < small) temp4 = Sgn(temp4) * small;

	return temp1 - MacauleyIndex(NextDGamma) * Dyadic2_2(temp3, temp2) / temp4;
}
/*************************************************************/
// GetNormalToYield() ----------------------------------------
Vector
ManzariDafalias::GetNormalToYield(const Vector &stress, const Vector &alpha)
{
	Vector devStress(6); devStress = GetDevPart(stress);
	double p = one3 * GetTrace(stress);
	Vector n(6); n = devStress - p * alpha;
	if (GetNorm_Contr(n) != 0)
		n = n / GetNorm_Contr(n);

	return n;
}
/*************************************************************/
// Check() ---------------------------------------------------
int
ManzariDafalias::Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha)
// Check if the solution of implicit integration makes sense
{
	int result = 0;
	Vector n(6);    n    = GetNormalToYield(stress, CurAlpha);
	Vector n_tr(6); n_tr = GetNormalToYield(TrialStress, CurAlpha);

	// check the direction of stress and trial stress
	if (DoubleDot2_2_Contr(n, n_tr) < 0) result = -1;

	// add any other checks here

	return result;
}
/*************************************************************/
// GetStateDependent() ----------------------------------------
void 
ManzariDafalias::GetStateDependent( const Vector &stress, const Vector &alpha, const double &e
				, const Vector &alpha_in, Vector &n, Vector &d, Vector &b
				, double &cos3Theta, double &h, double &psi, double &alphaBtheta
				, double &alphaDtheta, double &b0)
{
	double p = one3 * GetTrace(stress);
	p = (p < small) ? small : p;

	n = GetNormalToYield(stress, alpha);
	psi = GetPSI(e, p);
	cos3Theta = GetLodeAngle(n);
	alphaBtheta = g(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
	alphaDtheta = g(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;
	b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);
	d    = root23 * alphaDtheta * n - alpha;
	b    = root23 * alphaBtheta * n - alpha;
	double alphaAlphaInDotN;
	if (fabs(DoubleDot2_2_Contr(alpha-alpha_in,n)) <= small)
		alphaAlphaInDotN = Sgn(alphaAlphaInDotN) * small;
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
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Contr requires vector of size(6)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (i>2) * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Cov() ---------------------------------------
double
ManzariDafalias::DoubleDot2_2_Cov(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "covariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Cov requires vector of size(6)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) - (i>2) * 0.5 * v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// DoubleDot2_2_Mixed() ---------------------------------------
double
ManzariDafalias::DoubleDot2_2_Mixed(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, one "covariant" and the other "contravariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! ManzariDafalias::DoubleDot2_2_Mixed requires vector of size(6)!" << endln;
	
	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}
/*************************************************************/
// GetNorm_Contr() ---------------------------------------------
double
ManzariDafalias::GetNorm_Contr(const Vector& v)
// computes contravariant (stress-like) norm of input 6x1 tensor
{
	if (v.Size() != 6)
		opserr << "\n ERROR! ManzariDafalias::GetNorm_Contr requires vector of size(6)!" << endln;

	double result=0.0;	
	result = sqrt(DoubleDot2_2_Contr(v,v));

	return result;
}
/*************************************************************/
// GetNorm_Cov() ---------------------------------------------
double
ManzariDafalias::GetNorm_Cov(const Vector& v)
// computes covariant (strain-like) norm of input 6x1 tensor
{
	if (v.Size() != 6)
		opserr << "\n ERROR! ManzariDafalias::GetNorm_Cov requires vector of size(6)!" << endln;
	
	double result=0.0;	
	result = sqrt(DoubleDot2_2_Cov(v,v));

	return result;
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
		opserr << "\n ERROR! ManzariDafalias::DoubleDot4_2 requires 6-by-6 matrix " << endln;

	return m1*v1;
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
		opserr << "\n ERROR! ManzariDafalias::DoubleDot2_4 requires 6-by-6 matrix " << endln;

	return  m1^v1;
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
		opserr << "\n ERROR! ManzariDafalias::DoubleDot4_4 requires 6-by-6 matrices " << endln;

	return m1*m2;
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
		opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 requires 6-by-6 matrix " << endln;

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
		opserr << "\n ERROR! ManzariDafalias::SingleDot2_4 requires 6-by-6 matrix " << endln;

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
	opserr << "\n ERROR! ManzariDafalias::SingleDot4_2 requires 6-by-6 matrix " << endln;
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
	return ( +   aV[0] * aV[1] * aV[2] 
		     + 2*aV[3] * aV[4] * aV[5] 
			 -   aV[0] * aV[5] * aV[5] 
			 -   aV[2] * aV[3] * aV[3] 
			 - 	 aV[1] * aV[4] * aV[4]);
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
	double max = fabs(mat(rowNo, 0));
	for (int i = 1; i < n; i++){
		if (max < fabs(mat(rowNo, i)))
			max = fabs(mat(rowNo, i));
	}
	return max;
}
//------------------------------------------------------------
double
MatrixMax_Cols(const Matrix& mat, int colNo)
{
	int n = mat.noRows();
	double max = fabs(mat(0, colNo));
	for (int i = 1; i < n; i++){
		if (max < fabs(mat(i, colNo)))
			max = fabs(mat(i, colNo));
	}
	return max;
}
//------------------------------------------------------------
double
MatrixMin_Rows(const Matrix& mat, int rowNo)
{
	int n = mat.noCols();
	double min = fabs(mat(rowNo, 0));
	for (int i = 1; i < n; i++){
		if (min > fabs(mat(rowNo, i)))
			min = fabs(mat(rowNo, i));
	}
	return min;
}
//------------------------------------------------------------
double
MatrixMin_Cols(const Matrix& mat, int colNo)
{
	int n = mat.noRows();
	double min = fabs(mat(0, colNo));
	for (int i = 1; i < n; i++){
		if (min > fabs(mat(i, colNo)))
			min = fabs(mat(i, colNo));
	}
	return min;
}
//------------------------------------------------------------
double
Sgn(const double& x)
{
	return x == 0 ? 1.0 : (fabs(x)/x);
}
