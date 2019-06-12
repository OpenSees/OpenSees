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
//            May 2013, University of Washington
// Last Modified: May 2019, Long Chen
// Description: This file contains the implementation for the ManzariDafalias class.

#include <ManzariDafalias.h>
#include <ManzariDafalias3D.h>
#include <ManzariDafaliasPlaneStrain.h>
#include <MaterialResponse.h>

#include <string.h>

#if defined(_WIN32) || defined(_WIN64)
#include <algorithm>
#define fmax std::max
#define fmin std::min
#endif

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
#define INT_RungeKutta45  45 // By Jose Abell @ UANDES

const double        ManzariDafalias::one3            = 1.0/3.0 ;
const double        ManzariDafalias::two3            = 2.0/3.0;
const double        ManzariDafalias::root23          = sqrt(2.0/3.0);
const double        ManzariDafalias::small           = 1e-10;
const double        ManzariDafalias::maxStrainInc    = 1e-5;
const bool          ManzariDafalias::debugFlag       = false;
const char unsigned ManzariDafalias::mMaxSubStep     = 10;
char  unsigned      ManzariDafalias::mElastFlag      = 1;

Vector              ManzariDafalias::mI1(6);
Matrix              ManzariDafalias::mIIco(6,6);
Matrix              ManzariDafalias::mIIcon(6,6);
Matrix              ManzariDafalias::mIImix(6,6);
Matrix              ManzariDafalias::mIIvol(6,6);
Matrix              ManzariDafalias::mIIdevCon(6,6);
Matrix              ManzariDafalias::mIIdevMix(6,6);
Matrix              ManzariDafalias::mIIdevCo(6,6);
ManzariDafalias::initTensors ManzariDafalias::initTensorOps;

static int numManzariDafaliasMaterials = 0;

void *
OPS_ManzariDafaliasMaterial(void)
{
  if (numManzariDafaliasMaterials == 0) 
    opserr << "ManzariDafalias nDmaterial - Written: A.Ghofrani, P.Arduino, U.Washington\n";					
	numManzariDafaliasMaterials++;		

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 19) {
    opserr << "Want: nDMaterial ManzariDafalias tag? G0? nu? e_init? Mc? c? lambda_c? e0? ksi?" <<
        " P_atm? m? h0? Ch? nb? A0? nd? z_max? cz? Rho? <IntScheme? TanType? JacoType? TolF? TolR?>" << endln;
    return 0;    
  }
  
  int    tag;
  double dData[18];
  double oData[5];

  oData[0] = 1;          // IntScheme
  oData[1] = 0;          // TanType
  oData[2] = 1;          // JacoType
  oData[3] = 1.0e-7;     // TolF
  oData[4] = 1.0e-7;     // TolR

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

    theMaterial = new ManzariDafalias(tag, ND_TAG_ManzariDafalias , dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] ,
                         dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11],
                         dData[12], dData[13], dData[14], dData[15], dData[16], dData[17], 
                         (int)oData[0],(int)oData[1], (int)oData[2], oData[3], oData[4]);

  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafalias material with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, double G0, double nu, 
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
    m_G0        = G0;
    m_nu        = nu;
    m_e_init    = e_init;
    m_Mc        = Mc;
    m_c         = c;
    m_lambda_c  = lambda_c;
    m_e0        = e0;
    m_ksi       = ksi;
    m_P_atm     = P_atm;
    m_m         = m;
    m_h0        = h0;
    m_ch        = ch;
    m_nb        = nb;
    m_A0        = A0;
    m_nd        = nd;
    m_z_max     = z_max;
    m_cz        = cz;

    massDen                = mDen;
    mTolF                  = TolF;
    mTolR                  = TolR;
    mJacoType              = JacoType;
    mScheme                = integrationScheme;
    mTangType              = tangentType;
    mIter                  = 0;
    mUseElasticTan         = false;
    mStressCorrectionInUse = true;
    
    initialize();
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, int classTag, double G0, double nu, 
    double e_init, double Mc, double c, double lambda_c, double e0, double ksi,
    double P_atm, double m, double h0, double ch, double nb, double A0, double nd,
    double z_max, double cz, double mDen, int integrationScheme, int tangentType, 
    int JacoType, double TolF, double TolR): NDMaterial(tag, classTag),
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
    m_G0        = G0;
    m_nu        = nu;
    m_e_init    = e_init;
    m_Mc        = Mc;
    m_c         = c;
    m_lambda_c  = lambda_c;
    m_e0        = e0;
    m_ksi       = ksi;
    m_P_atm     = P_atm;
    m_m         = m;
    m_h0        = h0;
    m_ch        = ch;
    m_nb        = nb;
    m_A0        = A0;
    m_nd        = nd;
    m_z_max     = z_max;
    m_cz        = cz;

    massDen                = mDen;
    mTolF                  = TolF;
    mTolR                  = TolR;
    mJacoType              = JacoType;
    mScheme                = integrationScheme;
    mTangType              = tangentType;
    mIter                  = 0;
    mUseElasticTan         = false;
    mStressCorrectionInUse = true;
    
    initialize();
}


// specific-type full constructor... needed in parallel mode [J.Abell jaabell@miuandes.cl]
ManzariDafalias ::ManzariDafalias(int classTag) 
    : NDMaterial(0, classTag),
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
    m_G0        = 0.0;
    m_nu        = 0.0;
    m_e_init    = 0.0;
    m_Mc        = 0.0;
    m_c         = 0.0;
    m_lambda_c  = 0.0;
    m_e0        = 0.0;
    m_ksi       = 0.0;
    m_P_atm     = 0.0;
    m_m         = 0.0;
    m_h0        = 0.0;
    m_ch        = 0.0;
    m_nb        = 0.0;
    m_A0        = 0.0;
    m_nd        = 0.0;
    m_z_max     = 0.0;
    m_cz        = 0.0;

    massDen      = 0.0;
    mTolF        = 1.0e-7;
    mTolR        = 1.0e-7;
    mJacoType    = 1;
    mScheme      = 2;
    mTangType    = 2;
    mIter        = 0;
    mUseElasticTan= false;
    
    mElastFlag = 0;

    this->initialize();
}



// null constructor
ManzariDafalias ::ManzariDafalias() 
    : NDMaterial(0, ND_TAG_ManzariDafalias),
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
    m_G0        = 0.0;
    m_nu        = 0.0;
    m_e_init    = 0.0;
    m_Mc        = 0.0;
    m_c         = 0.0;
    m_lambda_c  = 0.0;
    m_e0        = 0.0;
    m_ksi       = 0.0;
    m_P_atm     = 0.0;
    m_m         = 0.0;
    m_h0        = 0.0;
    m_ch        = 0.0;
    m_nb        = 0.0;
    m_A0        = 0.0;
    m_nd        = 0.0;
    m_z_max     = 0.0;
    m_cz        = 0.0;

    massDen                = 0.0;
    mTolF                  = 1.0e-7;
    mTolR                  = 1.0e-7;
    mJacoType              = 1;
    mScheme                = 2;
    mTangType              = 2;
    mIter                  = 0;
    mUseElasticTan         = false;
    mStressCorrectionInUse = true;
    
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
    Vector n(6), d(6), b(6), R(6);
    double cos3Theta, h, psi, aB, aD, b0, A, D, B, C;

    mAlpha_in_n = mAlpha_in;

    mSigma_n    = mSigma;
    mEpsilon_n  = mEpsilon;
    mEpsilonE_n = mEpsilonE;
    mAlpha_n    = mAlpha;
    mFabric_n   = mFabric;
    mDGamma_n   = mDGamma;
    mVoidRatio  = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);

	// This is needed for the ManzariDafaliasRO subclass
    GetStateDependent(mSigma, mAlpha, mFabric, mVoidRatio, mAlpha_in, 
              n, d, b, cos3Theta, h, psi, aB, aD, b0, A, D, B, C, R);
    this->GetElasticModuli(mSigma, mVoidRatio, mK, mG, D);

//      double q = sqrt(1.5 * DoubleDot2_2_Contr(GetDevPart(mSigma), GetDevPart(mSigma)));
//      double p = one3*GetTrace(mSigma);
//      opserr << "Committed stress (tag = " << this->getTag() << ") = " << mSigma << "Yield = " << GetF(mSigma, mAlpha) << endln << "p = " << p << ", q = " << q << ", eta = " << q/p << endln;
 
// opserr << "psi = " << psi << endln;
// opserr << "alpha_b = " << aB << endln;
// opserr << "alpha_d = " << aD << endln;
// opserr << "b0 = " << b0 << endln;
// opserr << "d = " << d;
// opserr << "b = " << b;
// opserr << "h = " << h << endln;
// opserr << "A = " << A << endln;
// opserr << "D = " << D << endln;
// opserr << "B = " << B << endln;
// opserr << "C = " << C << endln;
// opserr << "R = " << R;
// opserr << "n = " << n;
//  opserr << endln;

    if (GetTrace(mSigma) > 0.01 * m_P_atm)
        mUseElasticTan = false;

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
    else if (strcmp(argv[0], "elasticstrains") == 0 || strcmp(argv[0],"estrains") == 0)
    {    
        return new MaterialResponse(this, 7, this->getEStrain());
    }
    else if (strcmp(argv[0], "plasticstrains") == 0 || strcmp(argv[0],"pstrains") == 0)
    {
        return new MaterialResponse(this, 8, this->getPStrain());
    }
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
        case 7://elasticstrains
            if (matInfo.theVector != 0)
                *(matInfo.theVector) = getEStrain();
            return 0;
        case 8://plasticstrains
            if (matInfo.theVector != 0)
                *(matInfo.theVector) = getPStrain();
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
    static Vector data(97);

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
    data(24) = 0; // used to be filled with mOrgTanType
    data(25) = mElastFlag;

    data(26) = mEpsilon(0);        data(32) = mEpsilon_n(0);    data(38) = mSigma(0);    data(44) = mSigma_n(0);
    data(27) = mEpsilon(1);        data(33) = mEpsilon_n(1);    data(39) = mSigma(1);    data(45) = mSigma_n(1);
    data(28) = mEpsilon(2);        data(34) = mEpsilon_n(2);    data(40) = mSigma(2);    data(46) = mSigma_n(2);
    data(29) = mEpsilon(3);        data(35) = mEpsilon_n(3);    data(41) = mSigma(3);    data(47) = mSigma_n(3);
    data(30) = mEpsilon(4);        data(36) = mEpsilon_n(4);    data(42) = mSigma(4);    data(48) = mSigma_n(4);
    data(31) = mEpsilon(5);        data(37) = mEpsilon_n(5);    data(43) = mSigma(5);    data(49) = mSigma_n(5);

    data(50) = mEpsilonE(0);       data(56) = mEpsilonE_n(0);   data(62) = mAlpha(0);    data(68) = mAlpha_n(0);
    data(51) = mEpsilonE(1);       data(57) = mEpsilonE_n(1);   data(63) = mAlpha(1);    data(69) = mAlpha_n(1);
    data(52) = mEpsilonE(2);       data(58) = mEpsilonE_n(2);   data(64) = mAlpha(2);    data(70) = mAlpha_n(2);
    data(53) = mEpsilonE(3);       data(59) = mEpsilonE_n(3);   data(65) = mAlpha(3);    data(71) = mAlpha_n(3);
    data(54) = mEpsilonE(4);       data(60) = mEpsilonE_n(4);   data(66) = mAlpha(4);    data(72) = mAlpha_n(4);
    data(55) = mEpsilonE(5);       data(61) = mEpsilonE_n(5);   data(67) = mAlpha(5);    data(73) = mAlpha_n(5);

    data(74) = mFabric(0);         data(80) = mFabric_n(0);     data(86) = mAlpha_in_n(0);
    data(75) = mFabric(1);         data(81) = mFabric_n(1);     data(87) = mAlpha_in_n(1);
    data(76) = mFabric(2);         data(82) = mFabric_n(2);     data(88) = mAlpha_in_n(2);
    data(77) = mFabric(3);         data(83) = mFabric_n(3);     data(89) = mAlpha_in_n(3);
    data(78) = mFabric(4);         data(84) = mFabric_n(4);     data(90) = mAlpha_in_n(4);
    data(79) = mFabric(5);         data(85) = mFabric_n(5);     data(91) = mAlpha_in_n(5);

    data(92) = mDGamma_n;
    data(93) = mDGamma;
    data(94) = mK;
    data(95) = mG;
    data(96) = m_Pmin;


    
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
    static Vector data(97);
    res = theChannel.recvVector(this->getDbTag(), commitTag, data);
    if (res < 0) {
        opserr << "WARNING: ManzariDafalias::recvSelf - failed to receive vector from channel" << endln;
        return -1;
    }

    // set member variables
    this->setTag((int)data(0));

    m_G0        = data(1);    
    m_nu        = data(2);
    m_e_init    = data(3);
    m_Mc        = data(4);
    m_c         = data(5);
    m_lambda_c  = data(6);
    m_e0        = data(7);
    m_ksi       = data(8);
    m_P_atm     = data(9);
    m_m         = data(10);
    m_h0        = data(11);
    m_ch        = data(12);
    m_nb        = data(13);
    m_A0        = data(14);
    m_nd        = data(15);
    m_z_max     = data(16);
    m_cz        = data(17);
    massDen     = data(18);

    mTolF        = data(19); 
    mTolR        = data(20); 
    mJacoType    = (int)data(21); 
    mScheme      = (int)data(22); 
    mTangType    = (int)data(23); 
    // mOrgTangType = (int)data(24); 
    mElastFlag   = (int)data(25); 

    mEpsilon(0)  = data(26);    mEpsilon_n(0)  = data(32);     mSigma(0) = data(38);      mSigma_n(0) = data(44); 
    mEpsilon(1)  = data(27);    mEpsilon_n(1)  = data(33);     mSigma(1) = data(39);      mSigma_n(1) = data(45); 
    mEpsilon(2)  = data(28);    mEpsilon_n(2)  = data(34);     mSigma(2) = data(40);      mSigma_n(2) = data(46); 
    mEpsilon(3)  = data(29);    mEpsilon_n(3)  = data(35);     mSigma(3) = data(41);      mSigma_n(3) = data(47); 
    mEpsilon(4)  = data(30);    mEpsilon_n(4)  = data(36);     mSigma(4) = data(42);      mSigma_n(4) = data(48); 
    mEpsilon(5)  = data(31);    mEpsilon_n(5)  = data(37);     mSigma(5) = data(43);      mSigma_n(5) = data(49); 
                                                                      
    mEpsilonE(0) = data(50);    mEpsilonE_n(0) = data(56);     mAlpha(0) = data(62);      mAlpha_n(0) = data(68); 
    mEpsilonE(1) = data(51);    mEpsilonE_n(1) = data(57);     mAlpha(1) = data(63);      mAlpha_n(1) = data(69); 
    mEpsilonE(2) = data(52);    mEpsilonE_n(2) = data(58);     mAlpha(2) = data(64);      mAlpha_n(2) = data(70); 
    mEpsilonE(3) = data(53);    mEpsilonE_n(3) = data(59);     mAlpha(3) = data(65);      mAlpha_n(3) = data(71); 
    mEpsilonE(4) = data(54);    mEpsilonE_n(4) = data(60);     mAlpha(4) = data(66);      mAlpha_n(4) = data(72); 
    mEpsilonE(5) = data(55);    mEpsilonE_n(5) = data(61);     mAlpha(5) = data(67);      mAlpha_n(5) = data(73); 

    mFabric(0)   = data(74);    mFabric_n(0)   = data(80);     mAlpha_in_n(0) = data(86);  
    mFabric(1)   = data(75);    mFabric_n(1)   = data(81);     mAlpha_in_n(1) = data(87);  
    mFabric(2)   = data(76);    mFabric_n(2)   = data(82);     mAlpha_in_n(2) = data(88);  
    mFabric(3)   = data(77);    mFabric_n(3)   = data(83);     mAlpha_in_n(3) = data(89);  
    mFabric(4)   = data(78);    mFabric_n(4)   = data(84);     mAlpha_in_n(4) = data(90);  
    mFabric(5)   = data(79);    mFabric_n(5)   = data(85);     mAlpha_in_n(5) = data(91);  

    mDGamma_n = data(92);
    mDGamma   = data(93); 
    mK        = data(94); 
    mG        = data(95); 
    m_Pmin    = data(96); 

    mVoidRatio  = m_e_init - (1 + m_e_init) * GetTrace(mEpsilon);

    //GetElasticModuli(mSigma, mVoidRatio, mK, mG);
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
        else if ((strcmp(argv[0],"refShearModulus") == 0) ||
                 (strcmp(argv[0],"ShearModulus") == 0))    {   // change G0
            return param.addObject(6, this);
        }
        else if (strcmp(argv[0],"poissonRatio") == 0) {      // change nu
            return param.addObject(7, this);
        }
        else if (strcmp(argv[0],"voidRatio") == 0) {        // change e_init
            return param.addObject(8, this);
        }
        else if (strcmp(argv[0],"stressCorrection") == 0) {        // change stress correction state
            return param.addObject(9, this);
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
        if (mElastFlag == 1)
            Elastic2Plastic();
    }
    // called materialState in tcl file
    else if (responseID == 5) {
        mElastFlag = (int)info.theDouble;
        if (mElastFlag == 1)
            Elastic2Plastic();
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
        m_e_init = (info.theDouble + eps_v) / (1 - eps_v);
    }
	// flag to apply stress correction
    else if (responseID == 9) {
        mStressCorrectionInUse = (info.theInt == 0) ? false : true;
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
    // set Initial Ce with p = p_atm
    Vector mSig(6);
    mSig(0) = m_P_atm;
    mSig(1) = m_P_atm;
    mSig(2) = m_P_atm;

    // set minimum allowable p
    m_Pmin      = 1.0e-4 * m_P_atm;
    m_Presidual = 1.0e-2 * m_P_atm;

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
    GetElasticModuli(mSig,mVoidRatio,mK,mG);
    mCe = GetStiffness(mK,mG);
    mCep = mCe;
    mCep_Consistent = mCe;

    // calculate machine epsilon (used for FDM Jacobian)
    mEPS = machineEPS();

    mUseElasticTan = false;
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


void ManzariDafalias::integrate() 
{
    // update alpha_in in case of unloading
	// I assume full elastic step and check if the new stress direction is "dramatically" 
	// different from the stress path (in reference to the center of the yield surface). 
	// Another method is to use the change in the stress direction.
    Vector trialDirection(6), tmp(6);
	// trialDirection = GetNormalToYield(mSigma_n + mCe*(mEpsilon - mEpsilon_n), mAlpha_n);
	// trialDirection = mCe * (mEpsilon - mEpsilon_n);
	tmp = mEpsilon; tmp -= mEpsilon_n;
	trialDirection = mCe * tmp;

    // if (DoubleDot2_2_Contr(mAlpha_n - mAlpha_in_n, trialDirection) < 0.0)
	tmp = mAlpha_n; tmp -= mAlpha_in_n;
	if (DoubleDot2_2_Contr(tmp, trialDirection) < 0.0)
        mAlpha_in = mAlpha_n;
    else
        mAlpha_in = mAlpha_in_n;

    // Force elastic response
    if (mElastFlag == 0) {
        elastic_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mEpsilon, mEpsilonE, mSigma, mAlpha, 
                mVoidRatio, mG, mK, mCe, mCep, mCep_Consistent);
    } 
    // ElastoPlastic response
    else {  
        // implicit schemes
        if ((mScheme == INT_BackwardEuler))
            BackwardEuler_CPPM(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in,
                mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, 
                mK, mCe, mCep, mCep_Consistent);
        // explicit schemes
        else
            explicit_integrator(mSigma_n, mEpsilon_n, mEpsilonE_n, mAlpha_n, mFabric_n, mAlpha_in,
                mEpsilon, mEpsilonE, mSigma, mAlpha, mFabric, mDGamma, mVoidRatio, mG, 
                mK, mCe, mCep, mCep_Consistent);
    }
}

void ManzariDafalias::elastic_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& NextStrain, Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha,
        double& NextVoidRatio, double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    Vector dStrain(6);
    
    // calculate elastic response
    // dStrain               = NextStrain - CurStrain;
	dStrain = NextStrain; dStrain -= CurStrain;
    NextVoidRatio         = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
    // NextElasticStrain     = CurElasticStrain + dStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;
    GetElasticModuli(CurStress, NextVoidRatio, K, G); 
    aCep_Consistent       = aCep = aC = GetStiffness(K, G);
    // NextStress            = CurStress + DoubleDot4_2(aC,dStrain);
	NextStress = CurStress; NextStress += DoubleDot4_2(aC, dStrain);

    //update State variables
    double p = one3 * GetTrace(NextStress)+ m_Presidual;
    if (p > small)
        NextAlpha = GetDevPart(NextStress) / p;
    return;
}


void ManzariDafalias::explicit_integrator(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    // function pointer to the integration scheme
    void (ManzariDafalias::*exp_int) (const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , 
        const Vector& ,    Vector& , Vector& , Vector& , Vector& , double& , double& ,  double& , double& , 
        Matrix& , Matrix& , Matrix& ) ;
    
    switch (mScheme) {
        case INT_ForwardEuler     :    // Forward Euler
            exp_int = &ManzariDafalias::ForwardEuler;
            break;

        case INT_ModifiedEuler    :    // Modified Euler with error control
            exp_int = &ManzariDafalias::ModifiedEuler;
            break;

        case INT_RungeKutta       :    // Runge Kutta 4th order
            exp_int = &ManzariDafalias::RungeKutta4;
            break;

        case INT_RungeKutta45        :    // Runge Kutta 45 with error control after Sloan (J. Abell @ UANDES added)
            exp_int = &ManzariDafalias::RungeKutta45;
            break;        

        case INT_MAXSTR_FE        :    // Forward Euler constraining maximum strain increment
        case INT_MAXSTR_MFE       :    // Modified Euler constraining maximum strain increment
        case INT_MAXSTR_RK        :    // Runge-Kutta 4-th order constraining maximum strain increment
            exp_int = &ManzariDafalias::MaxStrainInc;
            break;

        case INT_MAXENE_FE        :    // Forward Euler constraining maximum energy increment
        case INT_MAXENE_MFE       :    // Modified Euler constraining maximum energy increment
        case INT_MAXENE_RK        :    //  Runge-Kutta 4-th order constraining maximum energy increment
            exp_int = &ManzariDafalias::MaxEnergyInc;
            break;
            
        default :
            exp_int = &ManzariDafalias::ModifiedEuler;
            break;
    }
    double elasticRatio, p, pn, f, fn;
    Vector dSigma(6), dStrain(6), dElasStrain(6);
    bool   p_tr_pos = true;

    NextVoidRatio          = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
    // dStrain                = NextStrain - CurStrain;
	dStrain = NextStrain; dStrain -= CurStrain;
	// NextElasticStrain     = CurElasticStrain + dStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;
    aC                     = GetStiffness(K, G);
    dSigma                 = DoubleDot4_2(aC, dStrain);
    // NextStress             = CurStress + dSigma;
	NextStress = CurStress; NextStress += dSigma;
    f                      = GetF(NextStress, CurAlpha);
    p                      = one3 * GetTrace(NextStress) + m_Presidual;

    if (p < m_Presidual)
        p_tr_pos = false;

    if (p_tr_pos && (f <= mTolF))
    {
        // This is a pure elastic loading/unloading
        NextAlpha         = CurAlpha;
        NextFabric        = CurFabric;
        NextDGamma        = 0;
        aCep_Consistent = aCep = aC;

        return;

    } else {
        fn = GetF(CurStress, CurAlpha);
        pn = one3 * GetTrace(CurStress) + m_Presidual;
        if (pn < m_Presidual)
        {
            if (debugFlag) 
                opserr << "Manzari Dafalias (tag = " << this->getTag() << ") : p_n < 0, This should have not happened!" << endln;
            NextStress = m_Pmin * mI1;
            NextAlpha.Zero();
            return;
        }
        
        if (fn > mTolF)
        {
            // This is an illegal stress state! This shouldn't happen.
            if (debugFlag) opserr << "stress state outside the yield surface!" << endln;
            if (debugFlag) opserr << "ManzariDafalias : Encountered an illegal stress state! Tag: " << this->getTag() << endln;
            if (debugFlag) opserr << "                  f = " << GetF(CurStress, CurAlpha) << endln;
            (this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, 
                    NextFabric, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

        } else if (fn < -mTolF) {
            // This is a transition from elastic to plastic
            elasticRatio = IntersectionFactor(CurStress, CurStrain, NextStrain, CurAlpha, 0.0, 1.0);
            // dSigma         = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
			dElasStrain = dStrain; dElasStrain *= elasticRatio;
			dSigma = DoubleDot4_2(aC, dElasStrain);
            // (this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
            //     CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio,
            //     G, K, aC, aCep, aCep_Consistent);
			(this->*exp_int)(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain,
				CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio,
				G, K, aC, aCep, aCep_Consistent);

        } else if (fabs(fn) < mTolF) {

            if (DoubleDot2_2_Contr(GetNormalToYield(CurStress, CurAlpha),dSigma)/(GetNorm_Contr(dSigma) == 0 ? 1.0 : GetNorm_Contr(dSigma)) > (- sqrt(mTolF))) {
                // This is a pure plastic step
                (this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, 
                    NextFabric, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
            } else {
                // This is an elastic unloding followed by plastic loading
                elasticRatio = IntersectionFactor_Unloading(CurStress, CurStrain, NextStrain, CurAlpha);
                // dSigma         = DoubleDot4_2(aC, elasticRatio*(NextStrain - CurStrain));
				dElasStrain = dStrain; dElasStrain *= elasticRatio;
				dSigma = DoubleDot4_2(aC, dElasStrain);
                // (this->*exp_int)(CurStress + dSigma, CurStrain + elasticRatio*(NextStrain - CurStrain), CurElasticStrain + elasticRatio*(NextStrain - CurStrain),
                //     CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio,
                //     G, K, aC, aCep, aCep_Consistent);
				(this->*exp_int)(CurStress + dSigma, CurStrain + dElasStrain, CurElasticStrain + dElasStrain,
					CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio,
					G, K, aC, aCep, aCep_Consistent);
            }
        }
    }
}


void ManzariDafalias::MaxStrainInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric, 
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{        
    // function pointer to the integration scheme
    void (ManzariDafalias::*exp_int) (const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , 
        const Vector& ,    Vector& , Vector& , Vector& , Vector& , double& , double& ,  double& , double& , 
        Matrix& , Matrix& , Matrix& ) ;
    
    switch (mScheme) {
        case INT_MAXSTR_FE    : // Forward Euler constraining maximum strain increment
            exp_int = &ManzariDafalias::ForwardEuler;
            break;

        default :
            exp_int = &ManzariDafalias::ForwardEuler;
            break;
    }
    
    NextDGamma = 0;

    Vector StrainInc(6); StrainInc = NextStrain - CurStrain;
    double maxInc = StrainInc(0);
    for(int ii=1; ii < 6; ii++)
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

        
        for(int ii = 1; ii <= numSteps; ii++)
        {
            nStrain = cStrain + StrainInc;

            (this->*exp_int)(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
            nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, nG, nK, nCe, nCep, nCepC);

            cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric;
        }

        NextElasticStrain    = nEStrain;
        NextStress            = nStress;
        NextAlpha            = nAlpha;
        NextFabric            = nFabric;

        Vector n(6), d(6), b(6), R(6), dPStrain(6); 
        double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
        GetStateDependent(NextStress, NextAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, 
                alphaDtheta, b0,A, D, B, C, R);
    
        dPStrain     = CurElasticStrain + (NextStrain - CurStrain) - NextElasticStrain;
        NextDGamma   = dPStrain.Norm() / R.Norm();

        aC    = nCe;
        aCep  = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
        aCep_Consistent = aCep;

    } else {
        (this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain,
            NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio, 
            G, K, aC, aCep, aCep_Consistent);
    }
}


void ManzariDafalias::MaxEnergyInc(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    // function pointer to the integration scheme
    void (ManzariDafalias::*exp_int) (const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , const Vector& , 
        const Vector& ,    Vector& , Vector& , Vector& , Vector& , double& , double& ,  double& , double& , 
        Matrix& , Matrix& , Matrix& ) ;
    
    switch (mScheme) {
        case INT_MAXENE_FE    : // Forward Euler constraining maximum energy increment
            exp_int = &ManzariDafalias::ForwardEuler;
            break;

        case INT_MAXENE_RK    : // Runge-Kutta 4-th order scheme
            exp_int = &ManzariDafalias::RungeKutta4;
            break;

        case INT_MAXENE_MFE    : // Modified Euler constraining maximum energy increment
            exp_int = &ManzariDafalias::ModifiedEuler;
            break;
            
        default :
            exp_int = &ManzariDafalias::ModifiedEuler;
            break;
    }

    double TolE = 1.0e-4;

    (this->*exp_int)(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain,
            NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio, 
            G, K, aC, aCep, aCep_Consistent);
    
    
    

    if ((DoubleDot2_2_Mixed(NextStrain - CurStrain, NextStress - CurStress) > TolE))     // || (DoubleDot2_2_Mixed(NextStress - CurStress, NextStress - CurStress) > TolE))
    {
        if (debugFlag) opserr << "******* Energy Inc > tol --> use sub-stepping" << endln;
        Vector StrainInc(6); StrainInc = NextStrain - CurStrain;
        StrainInc = (NextStrain - CurStrain) / 2;    
    
        Vector cStress(6), cStrain(6), cAlpha(6), cFabric(6), cAlpha_in(6), cEStrain(6);
        Vector nStrain(6) ,nEStrain(6), nStress(6), nAlpha(6), nFabric(6), nAlpha_in(6);
        Matrix nCe(6,6), nCep(6,6), nCepC(6,6);
        double nDGamma, nVoidRatio, nG, nK;
        Vector n(6), d(6), b(6), R(6), dPStrain(6); 
        //double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
                
        // create temporary variables
        cStress = CurStress; cStrain = CurStrain; cAlpha = CurAlpha; cFabric = CurFabric;
        cAlpha_in = alpha_in; cEStrain = CurElasticStrain;

        
        for(int ii=1; ii <= 2; ii++)
        {
            nStrain = cStrain + StrainInc;

            (this->*exp_int)(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
            nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, nG, nK, nCe, nCep, nCepC);

            cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric;
        }

        NextElasticStrain    = nEStrain;
        NextStress            = nStress;
        NextAlpha            = nAlpha;
        NextFabric            = nFabric;
        
        //GetStateDependent(NextStress, NextAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, 
        //        alphaDtheta, b0, A, D, B, C, R);
        //
        //dPStrain     = CurElasticStrain + (NextStrain - CurStrain) - NextElasticStrain;
        //NextDGamma   = dPStrain.Norm() / (R.Norm() == 0 ? 1.0 : R.Norm());
        //GetElasticModuli(NextStress, NextVoidRatio, K, G);
        //aC = GetStiffness(K, G);
        //aCep = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
        //aCep_Consistent = aCep;
        aC = nCe;
        aCep = nCep;
        aCep_Consistent = nCepC;
    }
}


void ManzariDafalias::ForwardEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    double CurVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
    NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
    aC = GetStiffness(K, G);
    Vector n(6), d(6), b(6), R(6), dPStrain(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
    GetStateDependent(CurStress, CurAlpha, CurFabric, CurVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0,
        A, D, B, C, R);
    double dVolStrain = GetTrace(NextStrain - CurStrain);
    Vector dDevStrain = GetDevPart(NextStrain - CurStrain);
    double p = one3 * GetTrace(CurStress) + m_Presidual;

    Vector r(6);
    if (p > small)
        Vector r = GetDevPart(CurStress) / p;

    double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
    
    double temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
        - K*D*DoubleDot2_2_Contr(n,r));

    // TODO: if temp4 == 0, the whole step is plastic. Take correct steps here.
    if (fabs(temp4) < small) temp4 = small;

    NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
    Vector dSigma   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
              (2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1);
    Vector dAlpha   = Macauley(NextDGamma) * two3 * h * b;
    Vector dFabric  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric);
           dPStrain = NextDGamma * ToCovariant(R);

    Matrix temp1 = 2.0*G*mIIdevMix + K*mIIvol;
    Vector temp2 = 2.0*G*n - DoubleDot2_2_Contr(n,r)*mI1;
    Vector temp3 = 2.0*G*(B*n-C*(SingleDot(n,n)-one3*mI1)) + K*D*mI1;

    aCep = temp1 - MacauleyIndex(NextDGamma) * Dyadic2_2(temp3, temp2) / temp4;
    aCep_Consistent = aCep;

    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
    NextStress = CurStress + dSigma;
    NextAlpha  = CurAlpha  + dAlpha;
    NextFabric = CurFabric + dFabric;

    return;
}


void ManzariDafalias::ModifiedEuler(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    double dVolStrain;
    Vector n(6), d(6), b(6), R(6), dDevStrain(6), r(6), dStrain(6), tmp0(6), tmp1(6), tmp2(6), tmp3(6);
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0,A, B, C, D, p, Kp;

    double T = 0.0, dT = 1.0, dT_min = 1e-6 , TolE = 1e-4;
    
    Vector nStress(6), nAlpha(6), nFabric(6), ndPStrain(6);
    Vector dSigma1(6), dSigma2(6), dAlpha1(6), dAlpha2(6), dFabric1(6), dFabric2(6),
           dPStrain1(6), dPStrain2(6);
    Matrix aCep1(6,6), aCep2(6,6), aCep_thisStep(6,6), aD(6,6);
    double temp4, curStepError, q = 1.0;

    // NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);
	dStrain = NextStrain; dStrain -= CurStrain;
	NextElasticStrain = CurElasticStrain; NextElasticStrain += dStrain;

    aC = GetStiffness(K, G);
    aD = GetCompliance(K, G);

    NextStress = CurStress;
    NextAlpha = CurAlpha;
    NextFabric = CurFabric;

    p = one3 * GetTrace(NextStress) + m_Presidual;
    if (p < m_Pmin + m_Presidual)
    {
        if (debugFlag)
            opserr << "Tag = " << this->getTag() << " : I have a problem (p < 0) - This should not happen!!!" << endln;        
        NextStress = GetDevPart(NextStress) + m_Pmin * mI1;
		p = m_Pmin;
    }
    // Set aCep_Consistent to zero for substepping process
    aCep_Consistent.Zero();

    while (T < 1.0)
    {
        // NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + T * (NextStrain - CurStrain));
		tmp0 = dStrain; tmp0 *= T; tmp0 += NextStrain;
		NextVoidRatio = m_e_init - (1 + m_e_init) * GetTrace(tmp0);
        
        // dVolStrain = dT * GetTrace(NextStrain - CurStrain);
        // dDevStrain = dT * GetDevPart(NextStrain - CurStrain);
		dVolStrain = dT * GetTrace(dStrain);
		dDevStrain = dT * GetDevPart(dStrain);

        // Calc Delta 1
        p = one3 * GetTrace(NextStress) + m_Presidual;
        GetStateDependent(NextStress, NextAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                b0, A, D, B, C, R);

        // r = GetDevPart(NextStress) / p;
		r = GetDevPart(NextStress); r /= p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);

        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));

        if (fabs(temp4) < small) 
        {
            // Neutral loading
            dSigma1.Zero();
            dAlpha1.Zero();
            dFabric1.Zero();
            dPStrain1 = dDevStrain + dVolStrain*mI1;
            
        } else {
            NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
             
            if (NextDGamma < -small)
            {
               if (debugFlag)
                    opserr << "dGamma cannot be negative! This should not happen. Setting dGamma = 0." << endln;
                NextDGamma = 0.0;
                dSigma1   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1;
                dAlpha1   = 3.0*(GetDevPart(NextStress + dSigma1) / GetTrace(NextStress + dSigma1) - GetDevPart(NextStress) / GetTrace(NextStress)) ;
                dFabric1.Zero();
                dPStrain1.Zero();
                mUseElasticTan = true;
            } else {
                // dSigma1   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
                //   (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
				tmp0 = mI1; tmp0 *= (K * dVolStrain);
				tmp1 = n; tmp1 *= B;
				tmp2 = mI1; tmp2 *= (-1.0 / 3.0); tmp2 += SingleDot(n, n); tmp2 *= C;
				tmp1 -= tmp2; tmp1 *= (2.0 *G);
				tmp3 = mI1; tmp3 *= (K * D); tmp1 += tmp3; tmp1 *= (-Macauley(NextDGamma));
				dSigma1 = ToContraviant(dDevStrain); dSigma1 *= (2.0 * G);
				dSigma1 += tmp0; dSigma1 += tmp1;

                // dAlpha1   = Macauley(NextDGamma) * two3 * h * b;
				dAlpha1 = b; dAlpha1 *= (Macauley(NextDGamma) * two3 * h);
                // dFabric1  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + NextFabric);
				dFabric1 = n;
				dFabric1 *= m_z_max;
				dFabric1 += NextFabric;
				dFabric1 *= -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0 * D);
                // dPStrain1 = NextDGamma * ToCovariant(R);
				dPStrain1 = ToCovariant(R); dPStrain1 *= NextDGamma;
            }
            aCep1 = GetElastoPlasticTangent(NextStress + dSigma1, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
        }

        // Calc Delta 2
        // p = one3 * GetTrace(NextStress + dSigma1) + m_Presidual;
		tmp0 = NextStress; tmp0 += dSigma1;    // tmp0 is NextStress + dSigma1 until calculating dSigma2
		p = one3 * GetTrace(tmp0) + m_Presidual;

        if (p < m_Presidual)
        {
            if (dT == dT_min)
                return;
            dT = fmax(0.1 * dT, dT_min);
            continue;
        }
            
        // GetStateDependent(NextStress + dSigma1, NextAlpha + dAlpha1, NextFabric + dFabric1, NextVoidRatio, alpha_in, n, d, b,
        //         Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
		tmp1.Zero();  tmp1 += NextAlpha; tmp1 += dAlpha1;  // tmp1 is NextAlpha + dAlpha1 until calculating dSigma2
		tmp2.Zero();  tmp2 += NextFabric; tmp2 += dFabric1;  // tmp2 is NextFabric + dFabric1 until calculating dSigma2
		GetStateDependent(tmp0, tmp1, tmp2, NextVoidRatio, alpha_in, n, d, b,
			Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        // r = GetDevPart(NextStress + dSigma1) / p;
		r = GetDevPart(tmp0); r /= p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));

        if (fabs(temp4) < small) 
        {
            // Neutral loading
            dSigma2.Zero();
            dAlpha2.Zero();
            dFabric2.Zero();
            dPStrain2 = dDevStrain + dVolStrain*mI1;
        
        } else {

            NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;

            if (NextDGamma < 0.0)
            {
                NextDGamma = 0.0;
                dSigma2   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1;
                dAlpha2   = 3.0*(GetDevPart(NextStress + dSigma2) / GetTrace(NextStress + dSigma2) - GetDevPart(NextStress) / GetTrace(NextStress)) ;
                dFabric2.Zero();
                dPStrain2.Zero();
                mUseElasticTan = true;
            } else {
                // dSigma2   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
                //   (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
				tmp0 = mI1; tmp0 *= (K * dVolStrain);
				tmp1 = n; tmp1 *= B;
				tmp2 = mI1; tmp2 *= (-1.0 / 3.0); tmp2 += SingleDot(n, n); tmp2 *= C;
				tmp1 -= tmp2; tmp1 *= (2.0 *G);
				tmp3 = mI1; tmp3 *= (K * D); tmp1 += tmp3; tmp1 *= (-Macauley(NextDGamma));
				dSigma2 = ToContraviant(dDevStrain); dSigma2 *= (2.0 * G);
				dSigma2 += tmp0; dSigma2 += tmp1;

                // dAlpha2   = Macauley(NextDGamma) * two3 * h * b;
				dAlpha2 = b; dAlpha2 *= (Macauley(NextDGamma) * two3 * h);
                // dFabric2  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + NextFabric + dFabric1);
				dFabric2 = n;
				dFabric2 *= m_z_max;
				dFabric2 += NextFabric;
				dFabric2 += dFabric1;
				dFabric2 *= (-1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0 * D));
                // dPStrain2 = NextDGamma * ToCovariant(R);
				dPStrain2 = ToCovariant(R);  dPStrain2 *= NextDGamma;
            }
        }
        
        aCep2 = GetElastoPlasticTangent(NextStress + dSigma1, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);

        // nStress = NextStress + 0.5 * (dSigma1 + dSigma2);
        // nAlpha  = NextAlpha  + 0.5 * (dAlpha1 + dAlpha2);
        // nFabric = NextFabric + 0.5 * (dFabric1 + dFabric2);
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

        p = one3 * GetTrace(nStress) + m_Presidual;
        
        if (p < m_Presidual)
        {
            if (dT == dT_min)
                return;
            dT = fmax(0.1 * dT, dT_min);
            continue;
        }

            double stressNorm = GetNorm_Contr(NextStress);
			tmp0 = dSigma2; tmp0 -= dSigma1;
            if (stressNorm < 0.5)
                // curStepError = GetNorm_Contr(dSigma2 - dSigma1);
				curStepError = GetNorm_Contr(tmp0);
            else 
                // curStepError = GetNorm_Contr(dSigma2 - dSigma1) / (2 * stressNorm);
				curStepError = GetNorm_Contr(tmp0) / (2 * stressNorm);
        
        if (curStepError > TolE)
        {
            if (debugFlag) 
                opserr << "--- Unsuccessful increment (tag = " << this->getTag() << "): Error =  " << curStepError << endln;
            if (debugFlag) 
                opserr << "                           T = " << T << ", dT = " << dT << endln;
            q = fmax(0.8 * sqrt(TolE / curStepError), 0.1);

            if (dT == dT_min) {
                mUseElasticTan = true;

				// NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
				tmp0 = dPStrain1; tmp0 += dPStrain2; tmp0 *= 0.5;
				NextElasticStrain -= tmp0;
                NextStress = nStress;
                double eta = sqrt(13.5) * GetNorm_Contr(GetDevPart(NextStress)) / GetTrace(NextStress);
                if (eta > m_Mc)
                    NextStress = one3 * GetTrace(NextStress) * mI1 + m_Mc / eta * GetDevPart(NextStress);
                NextAlpha  = CurAlpha + 3.0 * (GetDevPart(NextStress)/GetTrace(NextStress) - GetDevPart(CurStress)/GetTrace(CurStress));
                
                T += dT;
            }
            dT = fmax(q * dT, dT_min);
        } else {
            
            if (debugFlag) 
                opserr << "+++ Successful increment: T = " << T << ", dT = " << dT << endln;

			// NextElasticStrain -= 0.5* (dPStrain1 + dPStrain2);
			tmp0 = dPStrain1; tmp0 += dPStrain2; tmp0 *= 0.5;
			NextElasticStrain -= tmp0;
            NextStress = nStress;
            NextAlpha  = nAlpha;
            NextFabric = nFabric;

            Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress,
                NextAlpha, NextFabric, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
            
            T += dT;

            // aCep_thisStep = 0.5 * (aCep1 + aCep2);
			aCep_thisStep = aCep1; aCep_thisStep += aCep2;
			aCep_thisStep *= 0.5;
            aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix);
        
            q = fmax(0.8 * sqrt(TolE / curStepError), 0.5);
            dT = fmax(q * dT, dT_min);
            dT = fmin(dT, 1 - T);
        }
    }
    return;
}


void ManzariDafalias::RungeKutta4(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    double CurVoidRatio, dVolStrain;
    Vector n(6), d(6), b(6), R(6), dDevStrain(6), r(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0,A, B, C, D, p, Kp;

    double T = 0.0, dT = 1.0;
    Vector nStress(6), nAlpha(6), nFabric(6), ndPStrain(6);
    Vector dSigma1(6), dSigma2(6), dSigma3(6), dSigma4(6), dSigma(6), 
        dAlpha1(6), dAlpha2(6), dAlpha3(6), dAlpha4(6), dAlpha(6), 
        dFabric1(6), dFabric2(6), dFabric3(6), dFabric4(6), dFabric(6),
        dPStrain1(6), dPStrain2(6), dPStrain3(6), dPStrain4(6), dPStrain(6);
    double temp4, q;
    
    CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
    NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);

    GetElasticModuli(CurStress, CurVoidRatio, K, G);
    aC = GetStiffness(K, G);

    NextStress = CurStress;
    NextAlpha = CurAlpha;
    NextFabric = CurFabric;

    while (T < 1.0)
    {
        NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + T * (NextStrain - CurStrain));
        
        dVolStrain = dT * GetTrace(NextStrain - CurStrain);
        dDevStrain = dT * GetDevPart(NextStrain - CurStrain);

        // Calc Delta 1
        GetStateDependent(CurStress, CurAlpha, CurFabric , CurVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
            b0, A, D, B, C, R);
        dVolStrain = GetTrace(NextStrain - CurStrain);
        dDevStrain = GetDevPart(NextStrain - CurStrain);
        p = one3 * GetTrace(CurStress) + m_Presidual;
        p = p < small ? small : p;
        r = GetDevPart(CurStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;

        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma1   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha1   = Macauley(NextDGamma) * two3 * h * b;
        dFabric1  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric);
        dPStrain1 = NextDGamma * ToCovariant(R);

        // Calc Delta 2
        GetElasticModuli(CurStress + 0.5 * dSigma1, CurVoidRatio, K, G);
        aC = GetStiffness(K, G);
        GetStateDependent(CurStress + 0.5 * dSigma1, CurAlpha + 0.5 * dAlpha1, CurFabric + 0.5 * dFabric1, CurVoidRatio, alpha_in, 
            n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        p = one3 * GetTrace(CurStress + 0.5 * dSigma1) + m_Presidual;
        p = p < small ? small : p;
        r = GetDevPart(CurStress + 0.5 * dSigma1) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;

        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,0.5*dDevStrain) - K*0.5*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma2   = 2.0*G*0.5* ToContraviant(dDevStrain) + K*0.5*dVolStrain*mI1 - Macauley(NextDGamma)*
              (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha2   = Macauley(NextDGamma) * two3 * h * b;
        dFabric2  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric + 0.5 * dFabric1);
        dPStrain2 = NextDGamma * ToCovariant(R);

        // Calc Delta 3
        GetElasticModuli(CurStress + 0.5 * dSigma2, CurVoidRatio, K, G);
        aC = GetStiffness(K, G);
        GetStateDependent(CurStress + 0.5 * dSigma2, CurAlpha + 0.5 * dAlpha2, CurFabric + 0.5 * dFabric2, CurVoidRatio, alpha_in, 
            n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        p = one3 * GetTrace(CurStress + 0.5 * dSigma2) + m_Presidual;
        p = p < small ? small : p;
        r = GetDevPart(CurStress + 0.5 * dSigma2) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;

        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,0.5*dDevStrain) - K*0.5*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma3   = 2.0*G*0.5* ToContraviant(dDevStrain) + K*0.5*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha3   = Macauley(NextDGamma) * two3 * h * b;
        dFabric3  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric + 0.5 * dFabric2);
        dPStrain3 = NextDGamma * ToCovariant(R);

        // Calc Delta 4
        GetElasticModuli(CurStress + dSigma3, CurVoidRatio, K, G);
        aC = GetStiffness(K, G);
        GetStateDependent(CurStress + dSigma3, CurAlpha + dAlpha3, CurFabric + dFabric3, CurVoidRatio, alpha_in, 
            n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        p = one3 * GetTrace(CurStress + dSigma3) + m_Presidual;
        p = p < small ? small : p;
        r = GetDevPart(CurStress + dSigma3) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;

        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma4   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha4   = Macauley(NextDGamma) * two3 * h * b;
        dFabric4  = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + CurFabric + dFabric3);
        dPStrain4 = NextDGamma * ToCovariant(R);
        
        // RK
        dSigma = (dSigma1 + dSigma4 + 2.0 * (dSigma2 + dSigma3)) / 6.0;
        dAlpha = (dAlpha1 + dAlpha4 + 2.0 * (dAlpha2 + dAlpha3)) / 6.0;
        dFabric = (dFabric1 + dFabric4 + 2.0 * (dFabric2 + dFabric3)) / 6.0;
        dPStrain = (dPStrain1 + dPStrain4 + 2.0 * (dPStrain2 + dPStrain3)) / 6.0;

        nStress = NextStress + dSigma;
        nAlpha  = NextAlpha  + dAlpha;
        nFabric = NextFabric + dFabric;


        if (false){ // Add a condition to make an adaptive integration increment
            //if (debugFlag) opserr << "---Unsuccessful increment: Error =  " << curStepError << endln;
            //if (debugFlag) opserr << "                           T = " << T << ", dT = " << dT << endln;
            q = 0.5;
            if (dT == 1e-4) {
                NextElasticStrain -= dPStrain;
                NextStress = nStress;
                NextAlpha  = nAlpha;
                NextFabric = nFabric;

                
                //Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress,
                //    NextAlpha, NextFabric, NextAlpha_in, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

                T += dT;
            }
            dT = fmax(q * dT, 1e-4);
        } else {
            
            //if (debugFlag) opserr << "+++Successful increment: T = " << T << ", dT = " << dT << endln;
            NextElasticStrain -= dPStrain;
            NextStress = nStress;
            NextAlpha  = nAlpha;
            NextFabric = nFabric;

            
            //Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress,
            //    NextAlpha, NextFabric, NextAlpha_in, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);
        
            q = 1.1;
            T += dT;
            dT = fmax(q * dT, 1e-4);
            dT = fmin(dT, 1 - T);
        }
        
    }
    return;
}




// RK45 integrator added by Jose Abell @ UANDES. (www.joseabell.com)
// Theory by Scott Sloan. 
// Sloan, S. W., Abbo, A. J., & Sheng, D. (2001). 
//  " Refined explicit integration of elastoplastic models with automatic error control. ""
//  Engineering Computations, 18(1/2), 121–194. https://doi.org/10.1108/02644400110365842
void ManzariDafalias::RungeKutta45(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent) 
{    
    static bool do_once = true;

    if (do_once) {
        opserr << 
        "ManzariDafalias::RungeKutta45 - RK45 integrator added by Jose Abell @ UANDES. (www.joseabell.com). \n\
                                         Theory by Scott Sloan. " << endln;
        do_once = false;
    }

    double CurVoidRatio, dVolStrain;
    double T = 0.0, dT = 1.0, dT_min = 1e-3 , TolE = mTolR;
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0,A, B, C, D, p, Kp;
    double temp4, q;

    //Static allocation so we can avoid mallocs and get maximum speed. This is not thread-safe
    static Vector n(6), d(6), b(6), R(6), dDevStrain(6), r(6); 
    static Vector nStress(6), nAlpha(6), nFabric(6), ndPStrain(6);
    static Vector dSigma1(6), dSigma2(6), dSigma3(6), dSigma4(6), dSigma5(6), dSigma6(6), dSigma(6), 
        dAlpha1(6), dAlpha2(6), dAlpha3(6), dAlpha4(6), dAlpha5(6), dAlpha6(6), dAlpha(6), 
        dFabric1(6), dFabric2(6), dFabric3(6), dFabric4(6), dFabric5(6), dFabric6(6), dFabric(6),
        dPStrain1(6), dPStrain2(6), dPStrain3(6), dPStrain4(6), dPStrain5(6), dPStrain6(6), dPStrain(6);
    static Matrix aCep1(6,6), aCep2(6,6), aCep3(6,6), aCep4(6,6), aCep5(6,6), aCep6(6,6), aCep_thisStep(6,6), aD(6,6);
    static Vector thisSigma(6), thisAlpha(6), thisFabric(6);    
    
    // Zero everything out for good measure. 
    n.Zero(); d.Zero(); b.Zero(); R.Zero(); dDevStrain.Zero(); r.Zero();
    nStress.Zero(); nAlpha.Zero(); nFabric.Zero(); ndPStrain.Zero();
    dSigma1.Zero(); dSigma2.Zero(); dSigma3.Zero(); dSigma4.Zero(); dSigma5.Zero(); dSigma6.Zero(); dSigma.Zero(); 
    dAlpha1.Zero(); dAlpha2.Zero(); dAlpha3.Zero(); dAlpha4.Zero(); dAlpha5.Zero(); dAlpha6.Zero(); dAlpha.Zero(); 
    dFabric1.Zero(); dFabric2.Zero(); dFabric3.Zero(); dFabric4.Zero(); dFabric5.Zero(); dFabric6.Zero(); dFabric.Zero();
    dPStrain1.Zero(); dPStrain2.Zero(); dPStrain3.Zero(); dPStrain4.Zero(); dPStrain5.Zero(); dPStrain6.Zero(); dPStrain.Zero();
    aCep1.Zero(); aCep2.Zero(); aCep3.Zero(); aCep4.Zero(); aCep5.Zero(); aCep6.Zero(); aCep_thisStep.Zero(); aD.Zero();
    thisSigma.Zero(); thisAlpha.Zero(); thisFabric.Zero();    
    
    CurVoidRatio      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain);
    NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain);
    NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain);

    GetElasticModuli(CurStress, CurVoidRatio, K, G);
    aC = GetStiffness(K, G);
    aD = GetCompliance(K, G);

    NextStress = CurStress;
    NextAlpha = CurAlpha;
    NextFabric = CurFabric;

    p = one3 * GetTrace(NextStress);
    if (p < m_Pmin)
    {
        if (debugFlag)
            opserr << "ManzariDafalias::RungeKutta45() - Tag = " << this->getTag() << " : I have a problem (p < 0) - This should not happen!!!" << endln;        
        NextStress = GetDevPart(NextStress) + m_Pmin * mI1;
        p = one3 * GetTrace(NextStress);
    }

    // Set aCep_Consistent to zero for substepping process
    aCep_Consistent.Zero();

    while (T < 1.0)
    {
        NextVoidRatio     = m_e_init - (1 + m_e_init) * GetTrace(NextStrain + T * (NextStrain - CurStrain));
        
        dVolStrain = dT * GetTrace(NextStrain - CurStrain);
        dDevStrain = dT * GetDevPart(NextStrain - CurStrain);


        // Calc Delta 1
        thisSigma = NextStress;
        thisAlpha = NextAlpha;
        thisFabric = NextFabric;

        GetStateDependent(thisSigma, thisAlpha, thisFabric,  NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        GetElasticModuli(thisSigma , NextVoidRatio, K, G);

        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma1   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha1   = Macauley(NextDGamma) * two3 * h * b;
        dFabric1   = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain1 = NextDGamma * ToCovariant(R);
        
        aCep1 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);



        // Calc Delta 2
        thisSigma = NextStress  + 0.5*dSigma1;
        thisAlpha = NextAlpha   + 0.5*dAlpha1;
        thisFabric = NextFabric + 0.5*dFabric1;

        GetStateDependent(thisSigma, thisAlpha, thisFabric, NextVoidRatio,  alpha_in,  n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma2   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha2   = Macauley(NextDGamma) * two3 * h * b;
        dFabric2   = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain2 = NextDGamma * ToCovariant(R);

        aCep2 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


        // Calc Delta 3
        thisSigma = NextStress  + 0.25*(dSigma1  + dSigma2);
        thisAlpha = NextAlpha   + 0.25*(dAlpha1  + dAlpha2);
        thisFabric = NextFabric + 0.25*(dFabric1 + dFabric2);
        GetElasticModuli(thisSigma , NextVoidRatio, K, G);
        GetStateDependent(thisSigma, thisAlpha, thisFabric, NextVoidRatio,  alpha_in,  n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma3   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dFabric3   = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain3 = NextDGamma * ToCovariant(R);

        aCep3 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


        // Calc Delta 4
        thisSigma = NextStress  - dSigma2  + 2*dSigma3;
        thisAlpha = NextAlpha   - dAlpha2  + 2*dAlpha3;
        thisFabric = NextFabric - dFabric2 + 2*dFabric3;
        GetElasticModuli(thisSigma , NextVoidRatio, K, G);
        GetStateDependent(thisSigma, thisAlpha, thisFabric, NextVoidRatio,  alpha_in,  n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma4   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dFabric4   = -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain4 = NextDGamma * ToCovariant(R);

        aCep4 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


        // Calc Delta 5
        thisSigma = NextStress  + (7*dSigma1  + 10*dSigma2  + dSigma4)/27;
        thisAlpha = NextAlpha   + (7*dAlpha1  + 10*dAlpha2  + dAlpha4)/27;
        thisFabric = NextFabric + (7*dFabric1 + 10*dFabric2 + dFabric4)/27;
        GetElasticModuli(thisSigma , NextVoidRatio, K, G);
        GetStateDependent(thisSigma, thisAlpha, thisFabric, NextVoidRatio,  alpha_in,  n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma5   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha5   = Macauley(NextDGamma) * two3 * h * b;
        dFabric5   =  -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain5 = NextDGamma * ToCovariant(R);
   
        aCep5 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);


        // Calc Delta 6
        thisSigma = NextStress  + (28*dSigma1  - 125*dSigma2  + 546*dSigma3  + 54*dSigma4  - 378*dSigma5)/625;
        thisAlpha = NextAlpha   + (28*dAlpha1  - 125*dAlpha2  + 546*dAlpha3  + 54*dAlpha4  - 378*dAlpha5)/625;
        thisFabric = NextFabric + (28*dFabric1 - 125*dFabric2 + 546*dFabric3 + 54*dFabric4 - 378*dFabric5)/625;
        GetElasticModuli(thisSigma , NextVoidRatio, K, G);
        GetStateDependent(thisSigma, thisAlpha, thisFabric, NextVoidRatio,  alpha_in,  n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, D, B, C, R);
        r = GetDevPart(NextStress) / p;
        Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
        temp4 = (Kp + 2.0*G*(B-C*GetTrace(SingleDot(n,SingleDot(n,n)))) 
            - K*D*DoubleDot2_2_Contr(n,r));
        if (fabs(temp4) < small) temp4 = small;
        NextDGamma      = (2.0*G*DoubleDot2_2_Mixed(n,dDevStrain) - K*dVolStrain*DoubleDot2_2_Contr(n,r))/temp4;
        dSigma6   = 2.0*G* ToContraviant(dDevStrain) + K*dVolStrain*mI1 - Macauley(NextDGamma)*
             (2.0*G*(B*n-C*(SingleDot(n,n)-1.0/3.0*mI1)) + K*D*mI1);
        dAlpha6   = Macauley(NextDGamma) * two3 * h * b;
        dFabric6   =  -1.0 * Macauley(NextDGamma) * m_cz * Macauley(-1.0*D) * (m_z_max * n + thisFabric);
        dPStrain6 = NextDGamma * ToCovariant(R);

        aCep6 = GetElastoPlasticTangent(thisSigma, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);




        // Update
        dSigma =   ( 14*dSigma1   + 35*dSigma4   + 162*dSigma5   + 125*dSigma6   ) / 336;
        dAlpha =   ( 14*dAlpha1   + 35*dAlpha4   + 162*dAlpha5   + 125*dAlpha6   ) / 336;
        dFabric =  ( 14*dFabric1  + 35*dFabric4  + 162*dFabric5  + 125*dFabric6  ) / 336;
        dPStrain = ( 14*dPStrain1 + 35*dPStrain4 + 162*dPStrain5 + 125*dPStrain6 ) / 336;

        nStress = NextStress + dSigma;
        nAlpha  = NextAlpha  + dAlpha;
        nFabric  = NextFabric  + dFabric;

        // Compute the error 
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

        double curStepError1 = GetNorm_Contr(-42*dSigma1 - 224*dSigma3 - 21*dSigma4 + 162*dSigma5 + 125*dSigma6 )/336;
        if (stressNorm >= 0.5) {curStepError1 /= (2 * stressNorm);}
    
        double curStepError2 = GetNorm_Contr(-42*dAlpha1 - 224*dAlpha3 - 21*dAlpha4 + 162*dAlpha5 + 125*dAlpha6 )/336;
        if (alphaNorm >= 0.5) {curStepError2 /=  (2 * alphaNorm);}
    
        double curStepError = fmax(curStepError1, curStepError2);
    
        if (curStepError > TolE)
        {
            if (debugFlag) 
            {
                opserr << "--- Unsuccessful increment (tag = " << this->getTag() << "): Error =  " << curStepError << endln;
                opserr << "                           T = " << T << ", dT = " << dT << endln;
            }
            
            q = fmax(0.8 * sqrt(TolE / curStepError), 0.1);
            // q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.1);

            if (dT == dT_min) {
                // opserr << "--- ManzariDafalias::RungeKutta45() (tag = " << this->getTag() << "): Warning! Convergence not achieved."<< endln;
                mUseElasticTan = true;

                NextElasticStrain -= dPStrain;// 0.5* (dPStrain1 + dPStrain2);
                NextStress = nStress;
                double eta = sqrt(13.5) * GetNorm_Contr(GetDevPart(NextStress)) / GetTrace(NextStress);
                if (eta > m_Mc)
                    NextStress = one3 * GetTrace(NextStress) * mI1 + m_Mc / eta * GetDevPart(NextStress);
                NextAlpha  = CurAlpha + 3.0 * (GetDevPart(NextStress)/GetTrace(NextStress) - GetDevPart(CurStress)/GetTrace(CurStress));
                
                // ++N_nonconverged;
                // max_error = fmax(max_error, curStepError);

                T += dT;
            }
            dT = fmax(q * dT, dT_min);
        } else {
            
            if (debugFlag) 
                opserr << "+++ Successful increment: T = " << T << ", dT = " << dT <<  "  Error =  " << curStepError << endln;

            NextElasticStrain -= dPStrain;
            NextStress = nStress;
            NextAlpha  = nAlpha;
            NextFabric = nFabric;

            // Stress_Correction(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, NextElasticStrain, NextStress,
            //     NextAlpha, NextFabric, NextDGamma, NextVoidRatio, G, K, aC, aCep, aCep_Consistent);

            T += dT;

            aCep_thisStep =   ( 14*aCep1 + 35*aCep4 + 162*aCep5 + 125*aCep6 ) /336;
            aCep_Consistent = aCep_thisStep * (aD * aCep_Consistent + T * mIImix);
        
            if(curStepError == 0)
                dT = 1 - T;
            else
            {
                // q = fmax(0.8 * sqrt(TolE / curStepError), 0.5);
                q = fmax(0.8 * pow(TolE / curStepError, 0.2), 0.5);
                dT = fmax(q * dT, dT_min);
                dT = fmin(dT, 1 - T);
            }
        }        
    }

    aCep = aCep_thisStep;

    return;
}






int ManzariDafalias::BackwardEuler_CPPM(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma,    double& NextVoidRatio,  double& G, double& K, Matrix& Ce, Matrix& Cep, Matrix& Cep_Consistent, int implicitLevel) 
{
    int errFlag = 1, SchemeControl = 2, mMaxSubStep = 10;
    // errFalg 1 : newton converged and results are fine
    //         0 : newton did not converge in MaxIter number of iterations
    //        -1 : the jacobian is singular or system cannot be solved
    //        -2 : converged stress has p < 0
    //        -3 : max number of sub-stepping reached
    //        -4 : (n:n_tr) < 0

    // check if max number of substepping is reached
    if (implicitLevel > mMaxSubStep){
       if (debugFlag) 
              opserr << "ManzariDafalias (Tag: " << this->getTag() << "): SubStepping did not converge in " << mMaxSubStep << " substeps." << endln;
        return -3;
    }

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
    GetElasticModuli(CurStress, CurVoidRatio, mK, mG); 
    GetElasticModuli(CurStress, CurVoidRatio, K, G); // This is needed. Check why?
    aC          = GetStiffness(K, G);
    TrialStress = CurStress + DoubleDot4_2(aC,(NextElasticStrain - CurElasticStrain));

    // In case of pure elastic response
    NextStress     = TrialStress;
    aCepConsistent = aCep = aC;

    // Trial yield function
    double NextF = GetF(NextStress, NextAlpha);
    double p     = one3 * GetTrace(NextStress);
    if (p < m_Pmin) 
    {
        double NextDLambda = 0;
        if (debugFlag)
             opserr << "The trial stress is in the negative pressure region!" << endln;
    
		// if ((NextF <= mTolF) && (p >= 0))
		// {
		// 			NextDLambda = (m_Pmin - p) / K;
		// 			NextElasticStrain += one3 * NextDLambda * mI1;
		// 			NextStress -= K * NextDLambda * mI1;
		// 			NextDGamma = 0.0;
		// 			return 1;
		// }
		// 
		// 
        if (p < 0) NextStress = m_Pmin * mI1;



        Vector Delta0(20), Delta1(19), InVariants(44), Delta(20);
        Delta1 = SetManzariComponent(NextStress, NextAlpha, NextFabric, NextDGamma);
        for(int ii = 0; ii < 19; ii++)
            Delta0(ii) = Delta1(ii);
        Delta0(19) = NextDLambda;
        InVariants = SetManzariStateInVar(NextStrain, CurStrain, CurStress, CurElasticStrain, CurAlpha, CurFabric,
                        CurVoidRatio, NextVoidRatio, alpha_in);

        // do newton iterations
		// Comment: I have already implemented the tension-cutoff surface. It seems that it's not working properly.
		// Using explicit integrator for the time being.
        //errFlag = NewtonIter2_negP(Delta0, InVariants, Delta, aCepConsistent);
		errFlag = 0;
        // check if newton converged
        if (errFlag == 1)
        {
            NextStress.Extract(Delta, 0, 1.0);
            NextAlpha.Extract(Delta, 6, 1.0);
            NextFabric.Extract(Delta, 12, 1.0);
            NextDGamma = Delta(18);
            NextDLambda = Delta(19);
            // check if the results are acceptible
            errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
        } else {
            if (debugFlag) 
                opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Explicit integration" << endln;
            
            explicit_integrator(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, 
                        NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio, 
                        G, K, aC, aCep, aCepConsistent);
            errFlag = 1;
        }

    } else if (NextF > mTolF) {// elastoplastic response
 
        // Vector n_tr(6);
        // n_tr = GetNormalToYield(NextStress, NextAlpha);

        // if (fabs(DoubleDot2_2_Contr(NextAlpha - alpha_in, n_tr)) < 1.0e-7)
        // {
        //     NextAlpha = GetDevPart(NextStress) / p;
        // }


        Vector Delta0(19), InVariants(44), Delta(19);

        Delta0     = SetManzariComponent(NextStress, NextAlpha, NextFabric, NextDGamma);
        InVariants = SetManzariStateInVar(NextStrain, CurStrain, CurStress, CurElasticStrain, CurAlpha, CurFabric,
                        CurVoidRatio, NextVoidRatio, alpha_in);

        // do newton iterations
        errFlag = NewtonIter2(Delta0, InVariants, Delta, aCepConsistent);
        
        // check if newton converged
        if (errFlag == 1)
        {
            NextStress.Extract(Delta, 0, 1.0);
            NextAlpha.Extract(Delta, 6, 1.0);
            NextFabric.Extract(Delta, 12, 1.0);
            NextDGamma = Delta(18);
            // check if the results are acceptible
            errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
        }

        // try the considerations to get an acceptible solution
        if (mScheme == INT_BackwardEuler)    // I can add the option for pure Backward Euler here.
        {
            // try different approaches and continue until a solution is found
            while(errFlag != 1)
            {
                if (errFlag == -1) SchemeControl = 3; // do an explicit integration
                if (errFlag == -2) SchemeControl = 2; // do sub-stepping

                Vector StrainInc(6), cStress(6), cStrain(6), cAlpha(6), cFabric(6), cAlpha_in(6), cEStrain(6);
                Vector nStrain(6) ,nEStrain(6), nStress(6), nAlpha(6), nFabric(6);
                Matrix nCe(6,6), nCep(6,6), nCepC(6,6);
                double nDGamma, nVoidRatio, nG, nK;
                int numSteps;

                // original strain increment
                StrainInc = NextStrain - CurStrain;

                // create temporary variables
                cStress = CurStress; cStrain = CurStrain; cAlpha = CurAlpha; cFabric = CurFabric;
                cAlpha_in = alpha_in; cEStrain = CurElasticStrain;

                if (SchemeControl == 1)
                {
                    // use numSteps steps of explicit integration as an initial guess to the newton iteration
                    numSteps = 50;
                    if (debugFlag) 
						opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Explicit step as initial guess" << endln;
                    for(int ii=1; ii <= numSteps; ii++)
                    {
                        nStrain = cStrain + StrainInc / numSteps;
                        ForwardEuler(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain,
                                nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, 
                                nG, nK, nCe, nCep, nCepC);
                        //ModifiedEuler(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain,
                        //        nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, 
                        //        nG, nK, nCe, nCep, nCepC);
                        cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; 
                    }
                    // do newton iterations
                    Delta0  = SetManzariComponent(nStress, nAlpha, nFabric, nDGamma);
                    errFlag = NewtonIter2(Delta0, InVariants, Delta, aCepConsistent);
                    // check if newton converged
                    if (errFlag == 1) 
                    {
                        NextStress.Extract(Delta, 0, 1.0);
                        NextAlpha.Extract(Delta, 6, 1.0);
                        NextFabric.Extract(Delta, 12, 1.0);
                        NextDGamma = Delta(18);

                        // check validity of the solution
                        errFlag = Check(TrialStress, NextStress, CurAlpha, NextAlpha);
                        // if not valid do implicit substepping
                        if (errFlag != 1)
                        {
                            SchemeControl += 1;
                            continue;
                        }
                    } else {
                        SchemeControl += 1;
                        continue;
                    }
                // explicit guess did not work, now do implicit substepping
                } else if (SchemeControl == 2) {

                    if (debugFlag) 
                        opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Implicit sub-stepping" << endln;

                    implicitLevel++;
                    nStrain = cStrain + StrainInc / 2;
                    // do a recursive BackwardEuler_CPPM on the first half of the strain increment
                    errFlag = BackwardEuler_CPPM(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain,
                                nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, 
                                nG, nK,nCe, nCep, nCepC, implicitLevel);
                    // check if maximum number of substepping levels is reached
                    if (errFlag == -3){
                        // do explicit integration over this strain increment
                        SchemeControl += 1;
                        continue;
                    }

                    cStress = nStress; cStrain = nStrain; cAlpha = nAlpha; cFabric = nFabric; 
                        
                    nStrain = cStrain + StrainInc / 2;
                    // do a recursive BackwardEuler_CPPM on the second half of the strain increment
                    errFlag = BackwardEuler_CPPM(cStress, cStrain, cEStrain, cAlpha, cFabric, cAlpha_in, nStrain, 
                        nEStrain, nStress, nAlpha, nFabric, nDGamma, nVoidRatio, nG, nK,nCe, nCep, 
                        nCepC, implicitLevel);
                    // check if maximum number of substepping levels is reached
                    if (errFlag == -3){
                        // do explicit integration over this strain increment
                        SchemeControl += 1;
                        continue;
                    }

                    // Update results from substepping
                    if (errFlag == 1)
                    {
                        NextStress = nStress;
                        NextAlpha  = nAlpha;
                        NextFabric = nFabric;
                        NextDGamma = nDGamma;
                        aC = nCe;
                        aCep = nCep;
                        aCepConsistent = nCepC;
                    } else {
                        SchemeControl += 1;
                        continue;
                    }

                } 
                // do explicit integration on this strain increment
                else {
                    if (debugFlag) 
                        opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Explicit integration" << endln;
                    
                    explicit_integrator(CurStress, CurStrain, CurElasticStrain, CurAlpha, CurFabric, alpha_in, NextStrain, 
                                NextElasticStrain, NextStress, NextAlpha, NextFabric, NextDGamma, NextVoidRatio, 
                                G, K, aC, aCep, aCepConsistent);
                    errFlag = 1;
                }
            }
        }


        Vector n(6), d(6), b(6), R(6), dPStrain(6); 
        double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
        GetStateDependent(NextStress, NextAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, 
                alphaDtheta, b0, A, D, B, C, R);
    
        dPStrain          = NextDGamma * ToCovariant(R);
        NextElasticStrain = CurElasticStrain + (NextStrain - CurStrain) - dPStrain;
        // GetElasticModuli(NextStress, CurVoidRatio, NextVoidRatio, NextElasticStrain, CurElasticStrain, K, G);
        // aC                = GetStiffness(K, G);
        aCep              = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
    }

    Ce = aC;
    Cep = aCep;
    Cep_Consistent = aCepConsistent;

    return errFlag;
}


double
ManzariDafalias::IntersectionFactor(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha, 
    double a0, double a1)
{
    double a = a0;
    double G, K, vR, f, f0, f1;
    Vector dSigma(6), dSigma0(6), dSigma1(6), strainInc(6);

    strainInc = NextStrain - CurStrain;

    vR      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + a0 * strainInc);
    GetElasticModuli(CurStress, vR, K, G);
    dSigma0 = a0 * DoubleDot4_2(GetStiffness(K, G), strainInc);
    f0 = GetF(CurStress + dSigma0, CurAlpha);

    vR      = m_e_init - (1 + m_e_init) * GetTrace(CurStrain + a1 * strainInc);
    GetElasticModuli(CurStress, vR, K, G);
    dSigma1 = a1 * DoubleDot4_2(GetStiffness(K, G), strainInc);
    f1 = GetF(CurStress + dSigma1, CurAlpha);

    for (int i = 1; i <= 10; i++)
    {
        a    = a1 - f1 * (a1-a0)/(f1-f0);
        dSigma = a * DoubleDot4_2(GetStiffness(K, G), strainInc);
        f    = GetF(CurStress + dSigma, CurAlpha);
        if (fabs(f) < mTolF) 
        {
            if (debugFlag) opserr << "Found alpha in " << i << " steps" << ", alpha = " << a << endln;
            break;
        }
        if (f * f0 < 0)
        {
            a1 = a;
            f1 = f;
        } else {
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
ManzariDafalias::IntersectionFactor_Unloading(const Vector& CurStress, const Vector& CurStrain, const Vector& NextStrain, const Vector& CurAlpha)
{
    double a = 0.0, a0 = 0.0 , a1 = 1.0, da;
    double G, K, vR, f;
    int nSub = 20;
    Vector dSigma(6), dSigma0(6), dSigma1(6), strainInc(6);

    strainInc = NextStrain - CurStrain;
    
    
    vR    = m_e_init - (1 + m_e_init) * GetTrace(CurStrain ); 
    GetElasticModuli(CurStress, vR, K, G);
    dSigma = DoubleDot4_2(GetStiffness(K, G), strainInc);

    for (int i = 1; i < nSub; i++)
    {
        da = (a1 - a0)/2.0;
        a = a1 - da;
        f    = GetF(CurStress + a * dSigma, CurAlpha);
        if (f > mTolF)
        {
            a1 = a;
        } else if (f < -mTolF) {
            a0 = a;
            break;
        } else {
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


void    
ManzariDafalias::Stress_Correction(const Vector& CurStress, const Vector& CurStrain, const Vector& CurElasticStrain,
        const Vector& CurAlpha, const Vector& CurFabric, const Vector& alpha_in, const Vector& NextStrain,
        Vector& NextElasticStrain, Vector& NextStress, Vector& NextAlpha, Vector& NextFabric,
        double& NextDGamma, double& NextVoidRatio,  double& G, double& K, Matrix& aC, Matrix& aCep, Matrix& aCep_Consistent)
{
    if (!mStressCorrectionInUse) return;

    Vector n(6), d(6), b(6), dPStrain(6), R(6), devStress(6), dSigma(6), dAlpha(6), dSigmaP(6), aBar(6), zBar(6);
    Vector r(6), dfrOverdSigma(6), dfrOverdAlpha(6);
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0;
    double A, B, C, D, p, fr, lambda, NextDLambda;
    int maxIter = 50;

    // see if p < 0
    p = one3 * GetTrace(NextStress) + m_Presidual;
    if (p < m_Pmin + m_Presidual)
    {
        p = m_Pmin + m_Presidual;
        if (false)
        {
        fr = GetF(NextStress, NextAlpha);
        if (fr < mTolF)
        {
            NextDLambda = (m_Pmin - p) / K;
            NextElasticStrain += one3 * NextDLambda * mI1;
            NextStress += K * NextDLambda * mI1;
            NextDGamma = 0.0;
            aCep_Consistent = aCep = aC = GetStiffness(K, G);

        } else {

            // Do Newton iterations to find NextDGamma
            GetStateDependent(NextStress, NextAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);
            R = GetDevPart(R);
            NextDGamma  = 0.0;
            NextDLambda = 0.0;

            Vector N(6); N = GetDevPart(NextStress) - p*NextAlpha;
            double fr1  = GetNorm_Contr(N)-root23*m_m*p;
            double fr2  = m_Pmin - p;
            double J11, J12, J21, J22;

            for (int i = 1; i <= maxIter; i++)
            {
                J11 = DoubleDot2_2_Contr(N/GetNorm_Contr(N),-2.0*G*R+K*D*NextAlpha)+root23*m_m*K*D;
                J12 = DoubleDot2_2_Contr(N/GetNorm_Contr(N),-K*NextAlpha)-root23*m_m*K;
                J21 = K*D;
                J22 = -K;
                
                double det = 1.0 / (J11*J22-J12*J21);

                NextDGamma  -= det * (J22*fr1-J12*fr2);
                NextDLambda -= det * (J11*fr2-J21*fr1);

                N = GetDevPart(NextStress) - p*NextAlpha - 2.0*G*NextDGamma*R + K*(D*NextDGamma-NextDLambda)*NextAlpha;

                fr1  = GetNorm_Contr(N)-root23*m_m*(p-K*(D*NextDGamma-NextDLambda));
                fr2  = m_Pmin - p + K*(D*NextDGamma-NextDLambda);


                if (fabs(fr1) + fabs(fr2)  < mTolF)
                    break;

                if(i == maxIter)
                {
                    if (debugFlag) 
                        opserr << "Still outside with f =  " << fr << endln;
                    NextStress = m_Pmin * mI1;
                    NextAlpha.Zero();
                    return;
                }
                
            }

            p = one3 * GetTrace(NextStress) + m_Presidual;

            Vector dPStrain(6);
            dPStrain = ToCovariant(NextDGamma * R + one3*(NextDGamma*D - NextDLambda) * mI1);
            NextElasticStrain -= dPStrain;
            NextStress -= aC * dPStrain;
        }

    }
        NextStress = p * mI1;
        NextAlpha.Zero();
        return;
    } else {
    
        // See if NextStress is outside yield surface
        fr = GetF(NextStress, NextAlpha);

        if (fabs(fr) < mTolF)
        {
            if (debugFlag) 
                opserr << "ManzariDafalias::StressCorrection() Stress state inside yield surface." << endln;
            return;
        } else {
            Vector nStress = NextStress;
            Vector nAlpha  = NextAlpha;
            for (int i = 1; i <= maxIter; i++)
            {
                if (debugFlag) 
                    opserr << "ManzariDafalias::StressCorrection() Stress state outside yield surface. Correction step =  " << i << ", f = " << fr << endln;
                
                devStress = GetDevPart(nStress);
            
                // do I need to update G and K? check this!
                // GetElasticModuli(CurStress, CurVoidRatio, K, G);

                aC = GetStiffness(K, G);

                GetStateDependent(nStress, nAlpha, NextFabric, NextVoidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);

                dSigmaP = DoubleDot4_2(aC, ToCovariant(R));
                aBar = two3 * h * b;
                r = devStress / p ;
                dfrOverdSigma = n - one3 * DoubleDot2_2_Contr(n, r) * mI1;
                dfrOverdAlpha = - p * n;
                lambda = fr / (DoubleDot2_2_Contr(dfrOverdSigma, dSigmaP)-DoubleDot2_2_Contr(dfrOverdAlpha, aBar));

                if (fabs(GetF(nStress - lambda * dSigmaP, nAlpha + lambda * aBar)) < fabs(fr))
                {
                    nStress -= lambda * dSigmaP;
                    nAlpha  += lambda * aBar;
                } else {
                    lambda = fr / DoubleDot2_2_Contr(dfrOverdSigma, dfrOverdSigma);
                    if (fabs(GetF(nStress - lambda * dfrOverdSigma, nAlpha)) < fabs(fr))
                        nStress -= lambda * dfrOverdSigma;
                    else
                    {
                        if (debugFlag)
                            opserr << "ManzariDafalias::StressCorrection() Couldn't decrease the yield function." << endln;
                        return;
                    }
                }
                
                fr = GetF(nStress, nAlpha);
                if (fabs(fr) < mTolF)
                {
                    NextStress = nStress;
                    NextAlpha  = nAlpha;
                    break;
                }

                if(i == maxIter)
                {
                    if (debugFlag) 
                        opserr << "Still outside with f =  " << fr << endln;
                    if (GetF(CurStress, NextAlpha) < mTolF)
                    {
                        Vector dSigma = NextStress - CurStress;
                        double alpha_up = 1.0;
                        double alpha_mid = 0.5;
                        double alpha_down = 0.0;
                        double fr_old = GetF(CurStress + alpha_mid * dSigma, NextAlpha);
                        for (int jj = 0; jj < maxIter; jj++)
                        {
                            if (fr_old < 0.0)
                            {
                               alpha_down = alpha_mid;
                               alpha_mid = 0.5 * (alpha_up + alpha_mid);
                            } else {
                               alpha_up = alpha_mid;
                               alpha_mid = 0.5 * (alpha_down + alpha_mid);
                            } 

                            fr_old = GetF(CurStress + alpha_mid * dSigma, NextAlpha);

                            if (fabs(fr_old) < mTolF)
                            {
                                NextStress = CurStress + alpha_mid * dSigma;
                                break;
                            }       
                            if(jj == maxIter)
                                //if (debugFlag) 
                                    opserr << "Still outside with f =  " << fr_old << endln;
                        }
                    } else {
                        NextStress = CurStress;
                        NextAlpha  = CurAlpha;
                        NextFabric = CurFabric;
                    }
                }
                
                p = one3 * GetTrace(NextStress) + m_Presidual;
            }
            NextElasticStrain = CurElasticStrain + DoubleDot4_2(GetCompliance(K, G), NextStress - CurStress);
            aCep = GetElastoPlasticTangent(NextStress, NextDGamma, CurStrain, NextStrain, G, K, B, C, D, h, n, d, b);
            aCep_Consistent = aCep;
        }
    }
    return;
}


int
ManzariDafalias::NewtonIter(const Vector& xo, const Vector& inVar, Vector& x, Matrix& aCepPart)
{
    // Newton Iterations, returns 1 : converged 
    //                            0 : did not converge in MaxIter number of iterations
    //                           -1 : the jacobian is singular or system cannot be solved
    // This function uses the full Jacobian matrix
    
    int MaxIter = 50;
    int MaxLS   = 10;
    int ResSize = xo.Size();
    int errFlag = 0;
    bool jacoFlag = true;
    Matrix (ManzariDafalias::*jacoFunc)(const Vector&, const Vector&);
    // Declare variables to be used
    static Vector sol(ResSize);
    static Vector R(ResSize), R2(ResSize);
    static Vector dX(ResSize);
    static Vector norms(ResSize+1);
	static Vector aux;
    static Matrix jaco(ResSize,ResSize);
    static Matrix jInv(ResSize,ResSize);
    double normR1, alpha;
    double aNormR1, aNormR2;

    switch (mJacoType) {
        case 0:
            jacoFunc = &ManzariDafalias::GetFDMJacobian;
            break;
        case 1:
            jacoFunc = &ManzariDafalias::GetJacobian;
            break;
        default :
            jacoFunc = &ManzariDafalias::GetJacobian;
    }

	if (!jacoFlag)
		aux = SetManzariComponent(mSigma_n, mAlpha_n, mFabric_n, mDGamma_n);

    sol = xo;
    alpha = 1.0;

	R = GetResidual(sol, inVar);
	aNormR1 = R.Norm();
	normR1 = 0.0;
	double tolMaterial = 1.0e-4 * aNormR1;

    for(mIter = 1; mIter <= MaxIter; mIter++)
    { 
		if (debugFlag) 
			opserr << "Iteration = " << (int)mIter << " , NewtonDecr = " << normR1 << " (tol = " << mTolR << ")" << ", Actual norm(R) = " << aNormR1 << endln;


		if (aNormR1 < mTolR + tolMaterial)
		{
			errFlag = 1;
			Matrix jInv(19, 19);
			jaco.Invert(jInv);
			aCepPart.Zero();
			aCepPart.Extract(jInv, 0, 0, 1.0);
			break;
		}

        if (jacoFlag)
            jaco = (this->*jacoFunc)(sol, inVar);
        else
            jaco = (this->*jacoFunc)(aux, inVar);
        
        errFlag = jaco.Solve(R, dX);
        if (errFlag != 0) 
        {
            if (debugFlag)
                opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - Jacobian!" << endln;
            errFlag = -1;
            break;
        }
        

        //for (int i = 1; i <= MaxLS; i++)
        //{
        //    R2      = GetResidual(sol - alpha * dX, inVar);
        //    aNormR2 = R2.Norm();
        //    //if (debugFlag) 
        //        opserr << "            LS Iter = " << (int) i << " , alpha = " << alpha << " , norm(R) = " << aNormR2 << endln;
		//
        //    if (aNormR2 < aNormR1)
        //    {
        //        sol -= alpha * dX;
        //        normR1 = fabs(alpha*R2^dX);
        //        aNormR1 = aNormR2;
        //        alpha = 1.0;
        //        break;
        //    } else {
        //        if (i == MaxLS) {
        //            sol -=  dX;
        //            alpha = 1.0;
		//			R = GetResidual(sol , inVar);
		//			aNormR1 = R.Norm();
        //            break;
        //        }
        //        // double alpha_o = alpha;
        //        // alpha = alpha * alpha * aNormR1 / (2.0 * (aNormR2 + alpha * aNormR1 - aNormR1));
        //        // if (alpha < 0.8 * alpha_o) alpha = 0.8 * alpha_o;
		//		alpha *= 0.8;
        //    }
        //}

		sol -= dX;
		R = GetResidual(sol, inVar);
		aNormR1 = R.Norm();
		normR1 = -1.0 * R ^ dX;
		        
    }

    return errFlag;
}

int
ManzariDafalias::NewtonIter2(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart)
{
    // Newton Iterations, returns 1 : converged 
    //                            0 : did not converge in MaxIter number of iterations
    //                           -1 : the jacobian is singular or system cannot be solved

    int MaxIter = 30;
    int MaxLS   = 15;
    int errFlag = 0;
    
    // residuals and incremenets
    Vector delSig(6), delAlph(6), delZ(6);
    Vector del(19), res(19), res2(19);
    double normR1 = 1.0, alpha = 1.0;
    double aNormR1 = 1.0, aNormR2 = 1.0;
    
    sol = xo;
	res.Zero();
	res = NewtonRes(sol, inVar);
	aNormR1 = res.Norm();
	double tolR_loc = aNormR1 * mTolR + mTolR;

    if(debugFlag) 
        opserr << "ManzariDafalias (Tag: " << this->getTag() << ") Newton Iterations:" << endln;

    for (mIter = 1; mIter <= MaxIter; mIter++) {
		if (debugFlag)
			opserr << "Iteration = " << (int)mIter << " , NewtonDecr = " << normR1 << " (tol = " << mTolR << ")" << ", Actual norm(R) = " << aNormR1 << endln;

		if (aNormR1 < tolR_loc)
		{
			errFlag = 1;
			break;
		}

        errFlag        = NewtonSol(sol, inVar, del, aCepPart);
        if (errFlag < 0)
            return errFlag;

        normR1        = res^del;
        
        //if ((normR1 > 0) && fabs(normR1) > 1.0e-4)
        //    del = -1.0 * res;
        
        //for (int i = 1; i <= MaxLS; i++)
        //{
        //    if (alpha * del.Norm() < 1.0e-10)
        //    {
        //        sol += (alpha * del);
        //        alpha = 1.0;
        //        break;
        //    }
        //    
        //    res2    = NewtonRes(sol + (alpha * del), inVar);
        //    aNormR2 = res2.Norm();
		//
        //    if(debugFlag) 
        //        opserr << "            LS Iter = " << (int) i << " , alpha = " << alpha << " , norm(R) = " << aNormR2 << " (normR1 = " << aNormR1 << ")" << endln;
		//
        //    if ((aNormR2 <= aNormR1) || (aNormR2 < mTolR))
        //    {
        //        sol += (alpha * del);
        //        normR1 = alpha*(res2^del);
        //        aNormR1 = aNormR2;
        //        alpha = 1.0;
        //        break;
        //    } else {
        //        // double alpha_o = alpha;
        //        // alpha = alpha * alpha * aNormR1 / (2.0 * (aNormR2 + alpha * aNormR1 - aNormR1));
        //        // if (alpha < 0.8 * alpha_o) alpha = 0.8 * alpha_o;
        //        // 
        //        // alpha = -1.0 * normR2 / (2.0 * (f2 - f1 - normR2));
        //        // alpha = (alpha < 0.3 * alpha_o) ? 0.3 * alpha_o : alpha;
        //        // alpha = (alpha > 1.1 * alpha_o) ? 1.1 * alpha_o : alpha;
		//
        //        alpha *= 0.8;
        //    }
        //    
        //    if (i == MaxLS) {
        //        sol +=  del;
        //        alpha = 1.0;
        //        break;
        //    }
        //}

		sol += del;
		res.Zero();
		res = NewtonRes(sol, inVar);
		aNormR1 = res.Norm();
    }
    
    return errFlag;
}


Vector
ManzariDafalias::NewtonRes(const Vector& x, const Vector& inVar)
{
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6), dEstrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, voidRatio;
    // state dependent variables
    Matrix aD(6,6);
    Vector n(6), d(6), b(6), R(6), devStress(6), r(6), aBar(6), zBar(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
        
    // residuals
    Vector R1(6); Vector R2(6); Vector R3(6); double R4;
    
    // read the trial values
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
    // curVoidRatio = inVar(36);
    voidRatio = inVar(37);
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    aD = GetCompliance(mK, mG);

    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);

    // devStress = GetDevPart(stress);
    // p = one3 * GetTrace(stress);
    // p = p < small ? small : p;
    aBar = two3 * h * b;
    zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
        
    dEstrain = aD * (stress - curStress);
    eStrain = curEStrain + dEstrain;
        
    R1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R);
    R2 = alpha   - curAlpha           - dGamma * aBar;
    R3 = fabric  - curFabric          - dGamma * zBar;
    R4 = GetF(stress, alpha);
    
    // opserr << "res 1 = " << R1.Norm() << ", res 2 = " << R2.Norm() << ", res 3 = " << R3.Norm() << ", f = " << R4 << endln;

    Vector res(19);
    // fill out residual vector
    res.Assemble(R1, 0, 1.0);
    res.Assemble(R2, 6, 1.0);
    res.Assemble(R3, 12, 1.0);
    res(18) = R4;

    return res;
}


int 
ManzariDafalias::NewtonSol(const Vector &xo, const Vector &inVar, Vector& del, Matrix& Cep)
{
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6), dEstrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, voidRatio;
    // state dependent variables
    Matrix aD(6,6), aC(6,6);
    Vector n(6), n2(6), d(6), b(6), R(6), devStress(6), r(6), aBar(6), zBar(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D, p, normR, gc;
        
    // analytical Jacobian
    double AlphaAlphaInDotN;
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
    Vector dAdOverdFabric(6), dDOverdFabric(6), dfOverdSigma(6), dfOverdAlpha(6);

    // Variables needed to solve the system of equations
    Matrix    DAlpha(6,6), DFabric(6,6), DSigma(6,6);
    Matrix    CAlpha(6,6), CFabric(6,6), CSigma(6,6), ASigma(6,6), ZSigma(6,6);
    Vector    ALambda(6), AConstant(6), ZLambda(6), ZConstant(6), LSigma(6), 
            SLambda(6), SConstant(6);
    double    LConstant;

    // Flags to consider the threshold values
    double dpFlag = 1.0, dnFlag = 1.0, dhFlag = 1.0;
    
    // residuals
    Vector R1(6); Vector R2(6); Vector R3(6); double R4;
    
    // read the trial values
    stress.Extract(xo, 0, 1.0);
    alpha.Extract(xo, 6, 1.0);
    fabric.Extract(xo, 12, 1.0);
    dGamma = xo(18);
        
    // current iteration invariants
    strain.Extract(inVar, 0, 1.0);
    curStrain.Extract(inVar, 6, 1.0);
    curStress.Extract(inVar, 12, 1.0);
    curEStrain.Extract(inVar, 18, 1.0);
    curAlpha.Extract(inVar, 24, 1.0);
    curFabric.Extract(inVar, 30, 1.0);
    // curVoidRatio = inVar(36);
    voidRatio = inVar(37);
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    aC = GetStiffness(mK, mG);
    aD = GetCompliance(mK, mG);

    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);
    if (fabs(DoubleDot2_2_Contr(alpha - alpha_in,n)) <= 1.0e-10)
    {
        AlphaAlphaInDotN = 1.0e-10;
        dGamma = 0.0;
        dhFlag = 0.0;
    } else {
        AlphaAlphaInDotN = fabs(DoubleDot2_2_Contr(alpha - alpha_in,n));
        dhFlag = 1.0;
    }

    n2 = SingleDot(n,n);
    devStress = GetDevPart(stress);
    p = one3 * GetTrace(stress);
    p = p < small ? small : p;
    dpFlag = p < small ? 0.0 : 1.0;
    r = devStress - p * alpha;
    normR = GetNorm_Contr(r);
    dnFlag = (normR == 0 ? 0.0 : 1.0);
    gc = g(Cos3Theta, m_c);
    aBar = two3 * h * b;
    zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
        
    dEstrain = aD * (stress - curStress);
    eStrain = curEStrain + dEstrain;
        
    R1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R);
    R2 = alpha   - curAlpha           - dGamma * aBar;
    R3 = fabric  - curFabric          - dGamma * zBar;
    R4 = GetF(stress, alpha);


    // d...OverdSigma : Arranged by order of dependence
    dnOverdSigma          = dnFlag * ( 1.0 / normR * (mIIdevCon - dpFlag*one3*Dyadic2_2(alpha,mI1) - 
        Dyadic2_2(n,n) + dpFlag*one3*DoubleDot2_2_Contr(alpha,n)*Dyadic2_2(n,mI1)));
    dPsiOverdSigma        = dpFlag * one3 * m_ksi * m_lambda_c / m_P_atm * pow(p/m_P_atm, m_ksi-1) * mI1;
    db0OverdSigma         = dpFlag*(-b0 / (6.0*p) * mI1);

    dCos3ThetaOverdSigma  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdSigma));
    dAdOverdSigma         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdSigma));
    dhOverdSigma          = dhFlag * (1.0 / AlphaAlphaInDotN * (db0OverdSigma - 
        h*DoubleDot2_4(alpha-alpha_in, ToCovariant(dnOverdSigma))));

    dgOverdSigma          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdSigma;

    dAlphaDOverdSigma     = m_Mc * exp(m_nd * psi) * (dgOverdSigma + m_nd * gc * dPsiOverdSigma);
    dCOverdSigma          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdSigma;
    dBOverdSigma          = 1.5 * (1.0 - m_c)/m_c * (dgOverdSigma * Cos3Theta + gc * dCos3ThetaOverdSigma);
    dAlphaBOverdSigma     = m_Mc * exp(-1.0*m_nb*psi) * (dgOverdSigma - m_nb * gc * dPsiOverdSigma);

	if (p < 0.05 * m_P_atm)
	{
		double be = 7.2713;
		double temp1 = exp(7.6349 - 7.2713 * p);
		double D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
		dDOverdSigma = D_factor * (dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)))) -
			one3 * Macauley(D) * be * temp1 / pow(1 + temp1, 2) * mI1;
	}
	else {
		dDOverdSigma = dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)));
	}
    dAbarOverdSigma       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdSigma) + 
        root23 * h * (Dyadic2_2(n, dAlphaBOverdSigma)+alphaBtheta * dnOverdSigma));

    dROverdSigma          = B * dnOverdSigma + Dyadic2_2(n, dBOverdSigma) - C * 
        (Trans_SingleDot4T_2(dnOverdSigma,n) + SingleDot2_4(n, dnOverdSigma)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdSigma) + one3 * Dyadic2_2(mI1, dDOverdSigma);
	dZbarOverdSigma       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdSigma)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdSigma;

    // d...OverdAlpha : Arranged by order of dependence
    dnOverdAlpha          = dnFlag * (p / normR * (Dyadic2_2(n,n) - mIIcon));

    dCos3ThetaOverdAlpha  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdAlpha));
    dAdOverdAlpha         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdAlpha));
    dhOverdAlpha          = dhFlag * (-1.0*h / AlphaAlphaInDotN * (n + 
        DoubleDot2_4(alpha-alpha_in,ToCovariant(dnOverdAlpha))));

    dgOverdAlpha          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdAlpha;

    dAlphaDOverdAlpha     = m_Mc * exp(m_nd * psi) * dgOverdAlpha;
    dCOverdAlpha          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdAlpha;
    dBOverdAlpha          = 1.5 * (1.0 - m_c)/m_c * (dgOverdAlpha * Cos3Theta + gc * dCos3ThetaOverdAlpha);
    dAlphaBOverdAlpha     = m_Mc * exp(-1.0*m_nb*psi) * dgOverdAlpha;

    dDOverdAlpha          = dAdOverdAlpha * (root23 * alphaDtheta - 
        DoubleDot2_2_Contr(alpha, n)) + A * (root23 * dAlphaDOverdAlpha -
        n - DoubleDot2_4(alpha, ToCovariant(dnOverdAlpha)));
    dAbarOverdAlpha       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdAlpha) +
        root23 * h * (Dyadic2_2(n, dAlphaBOverdAlpha)+alphaBtheta * dnOverdAlpha) - h * mIIcon);

    dROverdAlpha          = B * dnOverdAlpha + Dyadic2_2(n, dBOverdAlpha) - C * 
        (Trans_SingleDot4T_2(dnOverdAlpha,n) + SingleDot2_4(n, dnOverdAlpha)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdAlpha) + one3 * Dyadic2_2(mI1, dDOverdAlpha);
	dZbarOverdAlpha       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdAlpha)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdAlpha;

    // d...OverdFabric : Arranged by order of dependence
    dAdOverdFabric        = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * n;

    dDOverdFabric         = dAdOverdFabric * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n));
    
    dROverdFabric         = one3 * Dyadic2_2(mI1, dDOverdFabric);

	dZbarOverdFabric      = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdFabric)
		- m_cz * Macauley(-1.0 * D) *  mIIcon;

        
    dfOverdSigma        = n - one3 * (DoubleDot2_2_Contr(n, alpha) + root23 * m_m) * mI1;
    dfOverdAlpha        = -1.0 * p * n;

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
        
    // Jacobian
    Matrix J11(6,6), J12(6,6), J13(6,6); Vector J14(6);
    Matrix J21(6,6), J22(6,6);           Vector J24(6);
    Matrix J31(6,6), J32(6,6), J33(6,6); Vector J34(6);
    Vector J41(6), J42(6);
    
    // inv(J22), inv(J33)
    Matrix J22_1(6,6), J33_1(6,6);
    
    J11        = aD + dGamma * mIIco * ToCovariant(dROverdSigma);
	J12        =      dGamma * mIIco * ToCovariant(dROverdAlpha);
	J13        =      dGamma * mIIco * ToCovariant(dROverdFabric);
    J14        = ToCovariant(R);

    J21        =  -1.0*dGamma * dAbarOverdSigma * mIIco;
    J22        =  mIImix - dGamma * dAbarOverdAlpha * mIIco;
    J24        =  -1.0 *  aBar;

    J31        =  -1.0*dGamma * dZbarOverdSigma * mIIco;
    J32        =  -1.0*dGamma * dZbarOverdAlpha * mIIco;
    J33        =  mIImix - dGamma * dZbarOverdFabric * mIIco;
    J34        =  -1.0 * zBar; 

    J41        = ToCovariant(dfOverdSigma);
    J42        = ToCovariant(dfOverdAlpha);

    
    if (J22.Invert(J22_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CAlpha!" << endln;
        //J22_1 = mIImix;
        return -1;
    }

    if (J33.Invert(J33_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CFabric!" << endln;
        //J33_1 = mIImix;
        return -1;
    }
    

    ASigma         = -1.0 * J22_1 * J21;
    ALambda        = -1.0 * J22_1 * J24;
    AConstant      = -1.0 * J22_1 * R2;

    ZSigma         = -1.0 * J33_1 * (J31 + J32 * ASigma);
    ZLambda        = -1.0 * J33_1 * (J34 + J32 * ALambda);
    ZConstant      = -1.0 * J33_1 * (R3  + J32 * AConstant);

    LSigma         = -1.0 / (J42 ^ ALambda) * (J41 + (ASigma ^ J42))   ;
    LConstant      = -1.0 / (J42 ^ ALambda) * (R4  + (J42 ^ AConstant));

    SConstant      = R1 + J12 * (ALambda * LConstant + AConstant) + 
                    J13 * (ZLambda * LConstant + ZConstant) + J14 * LConstant;
    DSigma         = J11 + J12 * (ASigma + Dyadic2_2(ALambda, LSigma)) + 
                    J13 * (ZSigma + Dyadic2_2(ZLambda, LSigma)) + Dyadic2_2(J14, LSigma);
    

    DSigma        = aC * DSigma;
    if (DSigma.Invert(CSigma) != 0) 
    {
        if (debugFlag) 
			opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - Cep!" << endln;
        //CSigma = aC;
        return -1;
    } else
        CSigma = CSigma * aC;

    Vector delSig(6), delAlph(6), delZ(6);
    double delGamma;
    delSig        = -1.0 *  CSigma * SConstant;
    delGamma    = (LSigma ^ delSig) + LConstant;
    // Check if delGamma is NaN
    if (delGamma != delGamma)
    {
        if (debugFlag)
            opserr << "ManzariDafalias(Tag: " << this->getTag() << "): delGamma is NaN!" << endln;
		return -1;
    } else {
        delZ           = ZSigma * delSig + delGamma * ZLambda + ZConstant;
        delAlph        = ASigma * delSig + delGamma * ALambda + AConstant;
        Cep            = -1.0 * CSigma;
    }
    del                = SetManzariComponent(delSig, delAlph, delZ, delGamma);
//
//       // Check
//       Matrix J(19,19);
//       J.Zero();
//       J.Assemble(J11,0,0 ,1.0);
//       J.Assemble(J12,0,6 ,1.0);
//       J.Assemble(J13,0,12,1.0);
//       J.Assemble(J14,0,18,1.0);
//   
//       J.Assemble(J21,6,0 ,1.0);
//       J.Assemble(J22,6,6 ,1.0);
//       J.Assemble(J24,6,18,1.0);
//   
//       J.Assemble(J31,12,0 ,1.0);
//       J.Assemble(J32,12,6 ,1.0);
//       J.Assemble(J33,12,12,1.0);
//       J.Assemble(J34,12,18,1.0);
//   
//       J.AssembleTranspose(J41,18,0 ,1.0);
//       J.AssembleTranspose(J42,18,6 ,1.0);
//       
//           Vector res(19);
//           // fill out residual vector
//           res.Assemble(R1, 0, 1.0);
//           res.Assemble(R2, 6, 1.0);
//           res.Assemble(R3, 12, 1.0);
//           res(18) = R4;
//   
//       Vector delta_test(19);
//       //J.Solve(res,delta_test);
//       //opserr << "Here is the Full Jacobian Solution: \n" << delta_test ;
//       //opserr << "Here is my Solution: \n" << -1.0*del ;
//       //opserr << "Here is the residual: \n" << res;
//       //opserr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endln << endln;
//   
//   //    del = delta_test;
//
    return 0;    
}











int
ManzariDafalias::NewtonIter3(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart)
{
    // Newton Iterations, returns 1 : converged 
    //                            0 : did not converge in MaxIter number of iterations
    //                           -1 : the jacobian is singular or system cannot be solved

    int MaxIter = 50;
    int MaxLS   = 15;
    int errFlag = 0;
    
    // residuals and incremenets
    Vector delSig(6), delAlph(6), delZ(6);
    Vector del(19), res(19), res2(19), JRes(19), sol2(19);
    double normR1 = 1.0, alpha = 1.0;
    double aNormR1 = 1.0, aNormR2 = 1.0;
    double normDel = 0.0;
    
    sol = xo;

    if(debugFlag) 
        opserr << "ManzariDafalias (Tag: " << this->getTag() << ") Newton Iterations:" << endln;

    for (mIter = 1; mIter <= MaxIter; mIter++) 
    {
        res.Zero(); // don't delete this. required because I'm using Assemble() which adds and not replaces

        errFlag        = NewtonSol2(sol, inVar, res, JRes, del, aCepPart);

        if (errFlag < 0)
            return errFlag;        

        normR1         = JRes^del;
        aNormR1        = res.Norm();
        normDel        = del.Norm();

        //if(debugFlag) 
            opserr << "Iteration = " << (int)mIter << " , NewtonDecr = " << normR1 <<   " (tol = " << mTolR << ")" << ", Actual norm(R) = " << aNormR1 << endln;
        
        if (aNormR1 < mTolR) 
        {
            errFlag = 1;
            break;
        }
        
        for (int i = 1; i <= MaxLS; i++)
        {
            if (alpha * normDel < 1.0e-10)
            {
                //sol = sol2;
                sol += (alpha * del);
                alpha = 1.0;
                break;
            }
            sol2 = sol + (alpha * del);
            res2    = NewtonRes(sol2, inVar);

            aNormR2 = res2.Norm();
            //if(debugFlag) 
                opserr << "            LS Iter = " << (int) i << " , alpha = " << alpha << " , norm(R) = " << aNormR2 << endln;

            //if ((f_new < f_old + alpha * fdd) || (aNormR2 < mTolR))
            if ((aNormR2 < aNormR1) || (aNormR2 < mTolR))
            {
                sol = sol2;
                alpha = 1.0;
                break;
            } else {
                alpha *= 0.8;
            }

            if (i == MaxLS) {
                sol +=  del;
                alpha = 1.0;
                break;
            }
        }
        
    }
    if (debugFlag)
        opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Residual (Norm = " << aNormR1 << "): " << endln << res << endln;
    return errFlag;
}


int 
ManzariDafalias::NewtonSol2(const Vector &xo, const Vector &inVar, Vector& res, Vector& JRes, Vector& del, Matrix& Cep)
{
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6), dEstrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, voidRatio;
    // state dependent variables
    Matrix aD(6,6), aC(6,6);
    Vector n(6), n2(6), d(6), b(6), R(6), devStress(6), r(6), aBar(6), zBar(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D, p, normR, gc;
        
    // analytical Jacobian
    double AlphaAlphaInDotN;
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
    Vector dAdOverdFabric(6), dDOverdFabric(6), dfOverdSigma(6), dfOverdAlpha(6);

    // Variables needed to solve the system of equations
    Matrix    DAlpha(6,6), DFabric(6,6), DSigma(6,6);
    Matrix    CAlpha(6,6), CFabric(6,6), CSigma(6,6), ASigma(6,6), ZSigma(6,6);
    Vector    ALambda(6), AConstant(6), ZLambda(6), ZConstant(6), LSigma(6), 
            SLambda(6), SConstant(6);
    double    LConstant;

    // Flags to consider the threshold values
    double dpFlag = 1.0, dnFlag = 1.0, dhFlag = 1.0;
    
    // residuals
    Vector R1(6); Vector R2(6); Vector R3(6); double R4;
    
    // read the trial values
    stress.Extract(xo, 0, 1.0);
    alpha.Extract(xo, 6, 1.0);
    fabric.Extract(xo, 12, 1.0);
    dGamma = xo(18);

    // current iteration invariants
    strain.Extract(inVar, 0, 1.0);
    curStrain.Extract(inVar, 6, 1.0);
    curStress.Extract(inVar, 12, 1.0);
    curEStrain.Extract(inVar, 18, 1.0);
    curAlpha.Extract(inVar, 24, 1.0);
    curFabric.Extract(inVar, 30, 1.0);
    // curVoidRatio = inVar(36);
    voidRatio = inVar(37);
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    aC = GetStiffness(mK, mG);
    aD = GetCompliance(mK, mG);

    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);
    if (DoubleDot2_2_Contr(alpha - alpha_in,n) <= 1.0e-3)
    {
        AlphaAlphaInDotN = 1.0e-4;
        dhFlag = 0.0;
    } else {
        AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in,n);
        dhFlag = 1.0;
    }

    n2 = SingleDot(n,n);
    devStress = GetDevPart(stress);
    p = one3 * GetTrace(stress);

        
    //p = p < m_Pmin ? m_Pmin : p;
    //dpFlag = p < m_Pmin ? 0.0 : 1.0;
    r = devStress - p * alpha;
    normR = GetNorm_Contr(r);
    dnFlag = normR == 0 ? 0.0 : 1.0;
    gc = g(Cos3Theta, m_c);
    aBar = two3 * h * b;
    zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
        
    dEstrain = aD * (stress - curStress);
    eStrain = curEStrain + dEstrain;

    R1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R);
    R2 = alpha   - curAlpha           - dGamma * aBar;
    R3 = fabric  - curFabric          - dGamma * zBar;
    R4 = GetF(stress, alpha);
    
    // fill out residual vector
    res.Assemble(R1, 0, 1.0);
    res.Assemble(R2, 6, 1.0);
    res.Assemble(R3, 12, 1.0);
    res(18) = R4;

    // d...OverdSigma : Arranged by order of dependence
    dnOverdSigma          = dnFlag * ( 1.0 / normR * (mIIdevCon - dpFlag*one3*Dyadic2_2(alpha,mI1) - 
        Dyadic2_2(n,n) + dpFlag*one3*DoubleDot2_2_Contr(alpha,n)*Dyadic2_2(n,mI1)));
    dPsiOverdSigma        = dpFlag * one3 * m_ksi * m_lambda_c / m_P_atm * pow(p/m_P_atm, m_ksi-1) * mI1;
    db0OverdSigma         = dpFlag*(-b0 / (6.0*p) * mI1);

    dCos3ThetaOverdSigma  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdSigma));
    dAdOverdSigma         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdSigma));
    dhOverdSigma          = dhFlag * (1.0 / AlphaAlphaInDotN * (db0OverdSigma - 
        h*DoubleDot2_4(alpha-alpha_in, ToCovariant(dnOverdSigma))));

    dgOverdSigma          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdSigma;

    dAlphaDOverdSigma     = m_Mc * exp(m_nd * psi) * (dgOverdSigma + m_nd * gc * dPsiOverdSigma);
    dCOverdSigma          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdSigma;
    dBOverdSigma          = 1.5 * (1.0 - m_c)/m_c * (dgOverdSigma * Cos3Theta + gc * dCos3ThetaOverdSigma);
    dAlphaBOverdSigma     = m_Mc * exp(-1.0*m_nb*psi) * (dgOverdSigma - m_nb * gc * dPsiOverdSigma);

    if (p < 0.001 * m_P_atm)
    {
        double be       = 207232.6584 * 2.0 * m_Pmin ;
        double temp1    = exp(20.72326584 - be*p);
        double D_factor = MacauleyIndex(D) / (1+temp1);
        dDOverdSigma    = (D_factor * dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
             A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)))) - 
             one3 * Macauley(D) * be * temp1 / pow(1+temp1,2) * mI1 ;
    } else {
        dDOverdSigma          = dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
             A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)));
    }
    dAbarOverdSigma       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdSigma) + 
        root23 * h * (Dyadic2_2(n, dAlphaBOverdSigma)+alphaBtheta * dnOverdSigma));

    dROverdSigma          = B * dnOverdSigma + Dyadic2_2(n, dBOverdSigma) - C * 
        (Trans_SingleDot4T_2(dnOverdSigma,n) + SingleDot2_4(n, dnOverdSigma)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdSigma) + one3 * Dyadic2_2(mI1, dDOverdSigma);
    dZbarOverdSigma       = -1.0 * m_cz * MacauleyIndex(-1.0*D) * 
        (-1.0*Dyadic2_2(m_z_max*n + fabric, dDOverdSigma) - m_z_max * D * dnOverdSigma);

    // d...OverdAlpha : Arranged by order of dependence
    dnOverdAlpha          = dnFlag * (p / normR * (Dyadic2_2(n,n) - mIIcon));

    dCos3ThetaOverdAlpha  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdAlpha));
    dAdOverdAlpha         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdAlpha));
    dhOverdAlpha          = dhFlag * (-1.0*h / AlphaAlphaInDotN * (n + 
        DoubleDot2_4(alpha-alpha_in,ToCovariant(dnOverdAlpha))));

    dgOverdAlpha          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdAlpha;

    dAlphaDOverdAlpha     = m_Mc * exp(m_nd * psi) * dgOverdAlpha;
    dCOverdAlpha          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdAlpha;
    dBOverdAlpha          = 1.5 * (1.0 - m_c)/m_c * (dgOverdAlpha * Cos3Theta + gc * dCos3ThetaOverdAlpha);
    dAlphaBOverdAlpha     = m_Mc * exp(-1.0*m_nb*psi) * dgOverdAlpha;

    dDOverdAlpha          = dAdOverdAlpha * (root23 * alphaDtheta - 
        DoubleDot2_2_Contr(alpha, n)) + A * (root23 * dAlphaDOverdAlpha -
        n - DoubleDot2_4(alpha, ToCovariant(dnOverdAlpha)));
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
    
    dROverdFabric         = one3 * Dyadic2_2(mI1, dDOverdFabric);

    dZbarOverdFabric      = -1.0*m_cz* MacauleyIndex(-1.0*D) * 
        (-1.0*Dyadic2_2(m_z_max * n + fabric, dDOverdFabric) - D * mIIcon);

        
    dfOverdSigma        = n - one3 * (DoubleDot2_2_Contr(n, alpha) + root23 * m_m) * mI1;
    dfOverdAlpha        = -1.0 * p * n;

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
        
    // Jacobian
    Matrix J11(6,6), J12(6,6), J13(6,6); Vector J14(6);
    Matrix J21(6,6), J22(6,6);           Vector J24(6);
    Matrix J31(6,6), J32(6,6), J33(6,6); Vector J34(6);
    Vector J41(6), J42(6);
    
    // inv(J22), inv(J33)
    Matrix J22_1(6,6), J33_1(6,6);
    
    J11        = aD + dGamma * ToCovariant(dROverdSigma) * mIIco;
    J12        = dGamma * ToCovariant(dROverdAlpha)  * mIIco;
    J13        = dGamma * ToCovariant(dROverdFabric) * mIIco;
    J14        = ToCovariant(R);

    J21        =  -1.0*dGamma * dAbarOverdSigma * mIIco;
    J22        =  mIImix - dGamma * dAbarOverdAlpha * mIIco;
    J24        =  -1.0 *  aBar;

    J31        =  -1.0*dGamma * dZbarOverdSigma * mIIco;
    J32        =  -1.0*dGamma * dZbarOverdAlpha * mIIco;
    J33        =  mIImix - dGamma * dZbarOverdFabric * mIIco;
    J34        =  -1.0 * zBar; 

    J41        = ToCovariant(dfOverdSigma);
    J42        = ToCovariant(dfOverdAlpha);

    // JRes
    Vector temp(6); double temp2;
    temp = (J11^R1) + (J21^R2) + (J31^R3) + R4 * J41;
    JRes.Assemble(temp, 0, 1.0);
    temp = (J12^R1) + (J22^R2) + (J32^R3) + R4 * J42;
    JRes.Assemble(temp, 6, 1.0);
    temp = (J13^R1) + (J33^R3);
    JRes.Assemble(temp, 12, 1.0);
    temp2 = (J14^R1) + (J24^R2) + (J34^R3);
    JRes(18) = temp2;
    
    if (J22.Invert(J22_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CAlpha!" << endln;
        J22_1 = mIImix;
        //return -1;
    }

    if (J33.Invert(J33_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CFabric!" << endln;
        J33_1 = mIImix;
        //return -1;
    }

    ASigma       = -1.0 * J22_1 * ToCovariant(J21);
    ALambda      = -1.0 * J22_1 * ToCovariant(J24);
    AConstant    = -1.0 * J22_1 * ToCovariant(R2);

    ZSigma       = -1.0 * J33_1 * ToCovariant(J31 + J32 * ToContraviant(ASigma));
    ZLambda      = -1.0 * J33_1 * ToCovariant(J34 + J32 * ToContraviant(ALambda));
    ZConstant    = -1.0 * J33_1 * ToCovariant(R3  + J32 * ToContraviant(AConstant));

    LSigma        = -1.0 / (J42 ^ (ToContraviant(ALambda))) * (J41 + (ASigma ^ (ToContraviant(J42))));
    LConstant    = -1.0 / (J42 ^ (ToContraviant(ALambda))) * (R4  + (J42 ^ (ToContraviant(AConstant))));

    SConstant   = R1 + J12 * ToContraviant(ALambda * LConstant + AConstant) + 
                    J13 * ToContraviant(ZLambda * LConstant + ZConstant) + J14 * LConstant;
    DSigma        = J11 + J12 * ToContraviant(ASigma + Dyadic2_2(ALambda, LSigma)) + 
                    J13 * ToContraviant(ZSigma + Dyadic2_2(ZLambda, LSigma)) + Dyadic2_2(J14, LSigma);
    

    DSigma        = aC * DSigma;
    if (DSigma.Invert(CSigma) != 0) 
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - Cep!" << endln;
        CSigma = aC;
        //return -1;
    } else
        CSigma = CSigma * aC;

    Vector delSig(6), delAlph(6), delZ(6);
    double delGamma;
    delSig        = -1.0 *  CSigma * SConstant;
    delGamma    = (LSigma ^ delSig) + LConstant;
    // Check if delGamma is NaN
    if (delGamma != delGamma)
    {
        if (debugFlag)
            opserr << "ManzariDafalias(Tag: " << this->getTag() << "): delGamma is NaN!" << endln;
        delSig.Zero();
        delGamma = 0.0;
        delZ.Zero();
        delAlph.Zero();
        Cep = aC;
    } else {
        delZ        = ToContraviant(ZSigma * delSig + delGamma * ZLambda + ZConstant);
        delAlph     = ToContraviant(ASigma * delSig + delGamma * ALambda + AConstant);
        Cep         = CSigma;
    }
    del            = SetManzariComponent(delSig, delAlph, delZ, delGamma);
    return 0;    
}





int
ManzariDafalias::NewtonIter2_negP(const Vector& xo, const Vector& inVar, Vector& sol, Matrix& aCepPart)
{
    // Newton Iterations, returns 1 : converged 
    //                            0 : did not converge in MaxIter number of iterations
    //                           -1 : the jacobian is singular or system cannot be solved

    int MaxIter = 30;
    int MaxLS   = 15;
    int errFlag = 0;
    
    // residuals and incremenets
    Vector delSig(6), delAlph(6), delZ(6);
    Vector del(20), res(20), res2(20);
    double normR1 = 1.0, alpha = 1.0;
    double aNormR1 = 1.0, aNormR2 = 1.0;

    sol = xo;
	res.Zero();
	res = NewtonRes_negP(sol, inVar);
	aNormR1 = res.Norm();
	double tolR_negP = aNormR1 * mTolR + mTolR;

    if(debugFlag) 
        opserr << "ManzariDafalias (Tag: " << this->getTag() << ") Newton Iterations:" << endln;

    for (mIter = 1; mIter <= MaxIter; mIter++)
    {
		if(debugFlag) 
			opserr << "Iteration = " << (int)mIter << ", Actual norm(R) = " << aNormR1 << " (tol = " << mTolR << ")" << endln;

		if (aNormR1 < tolR_negP)
		{
			errFlag = 1;
			break;
		}

        errFlag        = NewtonSol_negP(sol, inVar, del, aCepPart);

		if (errFlag < 0)
			return errFlag;

		double delNorm = del.Norm();
		if (delNorm < 1.0e-6)
		{
			if(debugFlag)
				opserr << "Not moving further..." << endln;
		}

        normR1        = res^del;
        
        if ((normR1 > 0) && fabs(normR1) > 1.0e-4)
            del = -1.0 * res;

        for (int i = 1; i <= MaxLS; i++)
        {
            if (alpha * del.Norm() < 1.0e-10)
            {
                sol += (alpha * del);
                alpha = 1.0;
                break;
            }
            
            res2    = NewtonRes_negP(sol + (alpha * del), inVar);
            aNormR2 = res2.Norm();
		
            if(debugFlag) 
                opserr << "            LS Iter = " << (int) i << " , alpha = " << alpha << " , norm(R) = " << aNormR2 << " (normR1 = " << aNormR1 << ")" << endln;
		
            if ((aNormR2 <= aNormR1) || (aNormR2 < tolR_negP))
            {
                sol += (alpha * del);
                normR1 = alpha*(res2^del);
                aNormR1 = aNormR2;
				res = res2;
                alpha = 1.0;
                break;
            } else {
                // double alpha_o = alpha;
                // alpha = alpha * alpha * aNormR1 / (2.0 * (aNormR2 + alpha * aNormR1 - aNormR1));
                // if (alpha < 0.8 * alpha_o) alpha = 0.8 * alpha_o;
                // 
                // alpha = -1.0 * aNormR2 / (2.0 * (f2 - f1 - aNormR2));
                // alpha = (alpha < 0.3 * alpha_o) ? 0.3 * alpha_o : alpha;
                // alpha = (alpha > 1.1 * alpha_o) ? 1.1 * alpha_o : alpha;
		
                alpha *= 0.8;
            }
            
            if (i == MaxLS) {
                sol +=  del;
                alpha = 1.0;
                break;
            }
        }

		//sol += del;

		//res.Zero();
		//res = NewtonRes_negP(sol, inVar);
		//aNormR1 = res.Norm();
    }
    
         
    return errFlag;
}



int 
ManzariDafalias::NewtonSol_negP(const Vector &xo, const Vector &inVar, Vector& del, Matrix& Cep)
{
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6), dEstrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, dLambda, voidRatio;
    // state dependent variables
    Matrix aD(6,6), aC(6,6);
    Vector n(6), n2(6), d(6), b(6), R(6), devStress(6), r(6), aBar(6), zBar(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D, p, normR, gc;
        
    // analytical Jacobian
    double AlphaAlphaInDotN;
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
    Vector dAdOverdFabric(6), dDOverdFabric(6), dfOverdSigma(6), dfOverdAlpha(6);

    // Variables needed to solve the system of equations
    Matrix    DAlpha(6,6), DFabric(6,6), DSigma(6,6);
    Matrix    CAlpha(6,6), CFabric(6,6), CSigma(6,6), ASigma(6,6), ZSigma(6,6);
    Vector    ALambda(6), AConstant(6), ZLambda(6), ZConstant(6), LSigma(6), 
            SLambda(6), SConstant(6);
    double    LConstant;

    // Flags to consider the threshold values
    double dpFlag = 1.0, dnFlag = 1.0, dhFlag = 1.0;
    
    // residuals
    Vector R1(6); Vector R2(6); Vector R3(6); double R4, R5;
    
    // read the trial values
    stress.Extract(xo, 0, 1.0);
    alpha.Extract(xo, 6, 1.0);
    fabric.Extract(xo, 12, 1.0);
    dGamma = xo(18);
    dLambda = xo(19);

    // current iteration invariants
    strain.Extract(inVar, 0, 1.0);
    curStrain.Extract(inVar, 6, 1.0);
    curStress.Extract(inVar, 12, 1.0);
    curEStrain.Extract(inVar, 18, 1.0);
    curAlpha.Extract(inVar, 24, 1.0);
    curFabric.Extract(inVar, 30, 1.0);
    // curVoidRatio = inVar(36);
    voidRatio = inVar(37);
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    aC = GetStiffness(mK, mG);
    aD = GetCompliance(mK, mG);

    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);
    if (fabs(DoubleDot2_2_Contr(alpha - alpha_in,n)) <= 1.0e-10)
    {
        AlphaAlphaInDotN = 1.0e-10;
        dGamma = 0.0;
		dLambda = 0.0;
        dhFlag = 0.0;
    } else {
        AlphaAlphaInDotN = fabs(DoubleDot2_2_Contr(alpha - alpha_in,n));
        dhFlag = 1.0;
    }

    n2 = SingleDot(n,n);
    devStress = GetDevPart(stress);
    p = one3 * GetTrace(stress);
    p = p < small ? small : p;
    dpFlag = p < small ? 0.0 : 1.0;
    r = devStress - p * alpha;
    normR = GetNorm_Contr(r);
    dnFlag = (normR == 0 ? 0.0 : 1.0);
    gc = g(Cos3Theta, m_c);
    aBar = two3 * h * b;
    zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
        
    dEstrain = aD * (stress - curStress);
    eStrain = curEStrain + dEstrain;
        
    R1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R) - one3 * dLambda * mI1;
    R2 = alpha   - curAlpha           - dGamma * aBar;
    R3 = fabric  - curFabric          - dGamma * zBar;
    R4 = GetF(stress, alpha);
    R5 = m_Pmin - one3 * GetTrace(stress);


    // d...OverdSigma : Arranged by order of dependence
    dnOverdSigma          = dnFlag * ( 1.0 / normR * (mIIdevCon - dpFlag*one3*Dyadic2_2(alpha,mI1) - 
        Dyadic2_2(n,n) + dpFlag*one3*DoubleDot2_2_Contr(alpha,n)*Dyadic2_2(n,mI1)));
    dPsiOverdSigma        = dpFlag * one3 * m_ksi * m_lambda_c / m_P_atm * pow(p/m_P_atm, m_ksi-1) * mI1;
    db0OverdSigma         = dpFlag*(-b0 / (6.0*p) * mI1);

    dCos3ThetaOverdSigma  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdSigma));
    dAdOverdSigma         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdSigma));
    dhOverdSigma          = dhFlag * (1.0 / AlphaAlphaInDotN * (db0OverdSigma - 
        h*DoubleDot2_4(alpha-alpha_in, ToCovariant(dnOverdSigma))));

    dgOverdSigma          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdSigma;

    dAlphaDOverdSigma     = m_Mc * exp(m_nd * psi) * (dgOverdSigma + m_nd * gc * dPsiOverdSigma);
    dCOverdSigma          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdSigma;
    dBOverdSigma          = 1.5 * (1.0 - m_c)/m_c * (dgOverdSigma * Cos3Theta + gc * dCos3ThetaOverdSigma);
    dAlphaBOverdSigma     = m_Mc * exp(-1.0*m_nb*psi) * (dgOverdSigma - m_nb * gc * dPsiOverdSigma);

	if (p < 0.05 * m_P_atm)
	{
		double be = 7.2713;
		double temp1 = exp(7.6349 - 7.2713 * p);
		double D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
		dDOverdSigma = D_factor * (dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)))) -
			one3 * Macauley(D) * be * temp1 / pow(1 + temp1, 2) * mI1;
	}
	else {
		dDOverdSigma = dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)));
	}
    
    dAbarOverdSigma       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdSigma) + 
        root23 * h * (Dyadic2_2(n, dAlphaBOverdSigma)+alphaBtheta * dnOverdSigma));

    dROverdSigma          = B * dnOverdSigma + Dyadic2_2(n, dBOverdSigma) - C * 
        (Trans_SingleDot4T_2(dnOverdSigma,n) + SingleDot2_4(n, dnOverdSigma)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdSigma) + one3 * Dyadic2_2(mI1, dDOverdSigma);
	dZbarOverdSigma       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdSigma)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdSigma;

    // d...OverdAlpha : Arranged by order of dependence
    dnOverdAlpha          = dnFlag * (p / normR * (Dyadic2_2(n,n) - mIIcon));

    dCos3ThetaOverdAlpha  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdAlpha));
    dAdOverdAlpha         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdAlpha));
    dhOverdAlpha          = dhFlag * (-1.0*h / AlphaAlphaInDotN * (n + 
        DoubleDot2_4(alpha-alpha_in,ToCovariant(dnOverdAlpha))));

    dgOverdAlpha          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdAlpha;

    dAlphaDOverdAlpha     = m_Mc * exp(m_nd * psi) * dgOverdAlpha;
    dCOverdAlpha          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdAlpha;
    dBOverdAlpha          = 1.5 * (1.0 - m_c)/m_c * (dgOverdAlpha * Cos3Theta + gc * dCos3ThetaOverdAlpha);
    dAlphaBOverdAlpha     = m_Mc * exp(-1.0*m_nb*psi) * dgOverdAlpha;

    dDOverdAlpha          = dAdOverdAlpha * (root23 * alphaDtheta - 
        DoubleDot2_2_Contr(alpha, n)) + A * (root23 * dAlphaDOverdAlpha -
        n - DoubleDot2_4(alpha, ToCovariant(dnOverdAlpha)));
    dAbarOverdAlpha       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdAlpha) +
        root23 * h * (Dyadic2_2(n, dAlphaBOverdAlpha)+alphaBtheta * dnOverdAlpha) - h * mIIcon);

    dROverdAlpha          = B * dnOverdAlpha + Dyadic2_2(n, dBOverdAlpha) - C * 
        (Trans_SingleDot4T_2(dnOverdAlpha,n) + SingleDot2_4(n, dnOverdAlpha)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdAlpha) + one3 * Dyadic2_2(mI1, dDOverdAlpha);
	dZbarOverdAlpha       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdAlpha)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdAlpha;

    // d...OverdFabric : Arranged by order of dependence
    dAdOverdFabric        = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * n;

    dDOverdFabric         = dAdOverdFabric * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n));
    
    dROverdFabric         = one3 * Dyadic2_2(mI1, dDOverdFabric);

	dZbarOverdFabric      = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdFabric)
		- m_cz * Macauley(-1.0 * D) *  mIIcon;

        
    dfOverdSigma        = n - one3 * (DoubleDot2_2_Contr(n, alpha) + root23 * m_m) * mI1;
    dfOverdAlpha        = -1.0 * p * n;

    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
    // -------------------------------------------------------------------------
        
    // Jacobian
    Matrix J11(6,6), J12(6,6), J13(6,6); Vector J14(6); Vector J15(6);
    Matrix J21(6,6), J22(6,6);           Vector J24(6);
    Matrix J31(6,6), J32(6,6), J33(6,6); Vector J34(6);
    Vector J41(6), J42(6);
    Vector J51(5);
    
    // inv(J22), inv(J33)
    Matrix J22_1(6,6), J33_1(6,6);
    
    J11        = aD + dGamma * mIIco * ToCovariant(dROverdSigma);
	J12        =      dGamma * mIIco * ToCovariant(dROverdAlpha);
	J13        =      dGamma * mIIco * ToCovariant(dROverdFabric);
    J14        = ToCovariant(R);
    J15        = -one3 * mI1;

    J21        =  -1.0*dGamma * dAbarOverdSigma * mIIco;
    J22        =  mIImix - dGamma * dAbarOverdAlpha * mIIco;
    J24        =  -1.0 *  aBar;

    J31        =  -1.0*dGamma * dZbarOverdSigma * mIIco;
    J32        =  -1.0*dGamma * dZbarOverdAlpha * mIIco;
    J33        =  mIImix - dGamma * dZbarOverdFabric * mIIco;
    J34        =  -1.0 * zBar; 

    J41        = ToCovariant(dfOverdSigma);
    J42        = ToCovariant(dfOverdAlpha);

    J51        = -one3 * mI1;

    
    if (J22.Invert(J22_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CAlpha!" << endln;
        //J22_1 = mIImix;
        return -1;
    }

    if (J33.Invert(J33_1) != 0)
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - CFabric!" << endln;
        //J33_1 = mIImix;
        return -1;
    }

    ASigma         = -1.0 * J22_1 * J21;
    ALambda        = -1.0 * J22_1 * J24;
    AConstant      = -1.0 * J22_1 * R2;

    ZSigma         = -1.0 * J33_1 * (J31 + J32 * ASigma);
    ZLambda        = -1.0 * J33_1 * (J34 + J32 * ALambda);
    ZConstant      = -1.0 * J33_1 * (R3  + J32 * AConstant);

    LSigma         = -1.0 / (J42 ^ ALambda) * (J41 + (ASigma ^ J42))   ;
    LConstant      = -1.0 / (J42 ^ ALambda) * (R4  + (J42 ^ AConstant));

    SConstant      = R1 + J12 * (ALambda * LConstant + AConstant) + 
                    J13 * (ZLambda * LConstant + ZConstant) + J14 * LConstant;
    DSigma         = J11 + J12 * (ASigma + Dyadic2_2(ALambda, LSigma)) + 
                    J13 * (ZSigma + Dyadic2_2(ZLambda, LSigma)) + Dyadic2_2(J14, LSigma);
    

    DSigma        = aC * DSigma;
    if (DSigma.Invert(CSigma) != 0) 
    {
        if (debugFlag) opserr << "ManzariDafalias (Tag: " << this->getTag() << "): Singular Matrix in Newton iterations - Cep!" << endln;
        //CSigma = aC;
        return -1;
    } else
        CSigma = CSigma * aC;

    Vector delSig(6), delAlph(6), delZ(6);
    double delGamma, delLambda;
    delLambda   = (3.0*(mI1^(CSigma * SConstant)) + 9.0*R5) / (mI1^( CSigma * mI1));
    delSig        = CSigma * (one3 * delLambda * mI1 - SConstant);
    delGamma    = (LSigma ^ delSig) + LConstant;
    // Check if delGamma is NaN
    if (delGamma != delGamma)
    {
        if (debugFlag)
            opserr << "ManzariDafalias(Tag: " << this->getTag() << "): delGamma is NaN!" << endln;
		return -1;
    } else {
        delZ           = ZSigma * delSig + delGamma * ZLambda + ZConstant;
        delAlph        = ASigma * delSig + delGamma * ALambda + AConstant;
        Cep            = -1.0 * CSigma;
    }
    Vector del_temp(19);
	del.Zero();
    del_temp       = SetManzariComponent(delSig, delAlph, delZ, delGamma);
    for (int ii = 0; ii < 19; ii++)
        del(ii) = del_temp(ii);
    del(19)        = delLambda;

//       // Check
//       Matrix J(20,20);
//       J.Zero();
//       J.Assemble(J11,0,0 ,1.0);
//       J.Assemble(J12,0,6 ,1.0);
//       J.Assemble(J13,0,12,1.0);
// 	   J.Assemble(J14, 0, 18, 1.0);
// 	   J.Assemble(J15, 0, 19, 1.0);
//   
//       J.Assemble(J21,6,0 ,1.0);
//       J.Assemble(J22,6,6 ,1.0);
//       J.Assemble(J24,6,18,1.0);
//   
//       J.Assemble(J31,12,0 ,1.0);
//       J.Assemble(J32,12,6 ,1.0);
//       J.Assemble(J33,12,12,1.0);
//       J.Assemble(J34,12,18,1.0);
//   
//       J.AssembleTranspose(J41,18,0 ,1.0);
// 	   J.AssembleTranspose(J42, 18, 6, 1.0);
// 	   J.AssembleTranspose(J51, 19, 0, 1.0);
//       
//           Vector res(20);
//           // fill out residual vector
//           res.Assemble(R1, 0, 1.0);
//           res.Assemble(R2, 6, 1.0);
//           res.Assemble(R3, 12, 1.0);
// 		   res(18) = R4;
// 		   res(19) = R5;
//   
// 		   opserr << "Resdiual = " << res + J * del;
//       //Vector delta_test(19);
//       //J.Solve(res,delta_test);
//   
//       //opserr << "Here is the Full Jacobian Solution: \n" << delta_test ;
//       //opserr << "Here is my Solution: \n" << -1.0*del ;
//       //opserr << "Here is the residual: \n" << res;
//       //opserr << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endln << endln;
//   
//   //    del = delta_test;
// 
    return 0;    
}


Vector
ManzariDafalias::NewtonRes_negP(const Vector& x, const Vector& inVar)
{
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6), dEstrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, dLambda, voidRatio;
    // state dependent variables
    Matrix aD(6,6);
    Vector n(6), d(6), b(6), R(6), devStress(6), r(6), aBar(6), zBar(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
        
    // residuals
    Vector R1(6); Vector R2(6); Vector R3(6); double R4, R5;
    
    // read the trial values
    stress.Extract(x, 0, 1.0);
    alpha.Extract(x, 6, 1.0);
    fabric.Extract(x, 12, 1.0);
    dGamma = x(18);
    dLambda = x(19);

    // current iteration invariants
    strain.Extract(inVar, 0, 1.0);
    curStrain.Extract(inVar, 6, 1.0);
    curStress.Extract(inVar, 12, 1.0);
    curEStrain.Extract(inVar, 18, 1.0);
    curAlpha.Extract(inVar, 24, 1.0);
    curFabric.Extract(inVar, 30, 1.0);
    // curVoidRatio = inVar(36);
    voidRatio = inVar(37);
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    aD = GetCompliance(mK, mG);

    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
                    b0, A, D, B, C, R);

    // devStress = GetDevPart(stress);
    // p = one3 * GetTrace(stress);
    // p = p < small ? small : p;
    aBar = two3 * h * b;
    zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);
        
    dEstrain = aD * (stress - curStress);
    eStrain = curEStrain + dEstrain;
        
    R1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R) - dLambda * mI1;
    R2 = alpha   - curAlpha           - dGamma * aBar;
    R3 = fabric  - curFabric          - dGamma * zBar;
    R4 = GetF(stress, alpha);
    R5 = m_Pmin - one3 * GetTrace(stress);
    
                                                                                                                              

    Vector res(20);
    // fill out residual vector
    res.Assemble(R1, 0, 1.0);
    res.Assemble(R2, 6, 1.0);
    res.Assemble(R3, 12, 1.0);
    res(18) = R4;
    res(19) = R5;

    return res;
}




















Vector
ManzariDafalias::GetResidual(const Vector& x, const Vector& inVar)
{
    // This function returns a 19x1 vector containing the residuals (needs to be checked)
    // stress: NextStress, also for other variables

    Vector Res(19);   // Residual Vector
    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6); // Stress and Hardening
    Vector fabric(6), curFabric(6); // Fabric
    double dGamma, voidRatio;

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
    // curVoidRatio = inVar[36];
    voidRatio = inVar[37];
    alpha_in.Extract(inVar,38,1.0);
    

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    
    // state dependent variables
    Vector n(6), d(6), b(6), R(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
        b0, A, D, B, C, R);
    Vector devStress = GetDevPart(stress);
    double p = one3 * GetTrace(stress);
    p = p < small ? small : p;
    Vector aBar(6); aBar = two3 * h * b;
    Vector zBar(6); zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);

    Matrix De = GetCompliance(mK, mG);
    Vector dEstrain(6);
    dEstrain = De * (stress - curStress);
    eStrain = curEStrain + dEstrain;

    // residuals
    Vector g1(6); Vector g2(6); Vector g3(6); double g4;

    g1 = eStrain - TrialElasticStrain + dGamma * ToCovariant(R);
    g2 = alpha   - curAlpha           - dGamma * aBar;
    g3 = fabric  - curFabric          - dGamma * zBar;
    g4 = GetF(stress, alpha);

    // put residuals in a one vector
    Res.Assemble(g1,  0, 1.0);
    Res.Assemble(g2,  6, 1.0);
    Res.Assemble(g3, 12, 1.0);
    Res(18) = g4;
    return Res;
}


Matrix 
ManzariDafalias::GetJacobian(const Vector &x, const Vector &inVar)
{
    // This function returns the full 19x19 Jacobian matrix
    // note: stress: NextStress, also for other variables

    Vector eStrain(6), strain(6), curStrain(6), curEStrain(6), TrialElasticStrain(6); // Strain
    Vector stress(6), alpha(6), curStress(6), curAlpha(6), alpha_in(6);
    Vector fabric(6), curFabric(6);
    double dGamma, voidRatio;
    double AlphaAlphaInDotN;
    
    // Flags to consider the threshold values
    double dpFlag = 1.0, dnFlag = 1.0, dhFlag = 1.0;
    
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
    // curVoidRatio = inVar[36];
    voidRatio = inVar[37];
    alpha_in.Extract(inVar,38,1.0);

    // elastic trial strain
    TrialElasticStrain = curEStrain + (strain - curStrain);
    
    // state dependent variables
    Vector n(6), n2(6), d(6), b(6), R(6); 
    double Cos3Theta, h, psi, alphaBtheta, alphaDtheta, b0, A, B, C, D;
    GetStateDependent(stress, alpha, fabric, voidRatio, alpha_in, n, d, b, Cos3Theta, h, psi, alphaBtheta, alphaDtheta, 
        b0, A, D, B, C, R);
    if (fabs(DoubleDot2_2_Contr(alpha - alpha_in,n)) <= small)
    {
        h = 1.0e10;
        AlphaAlphaInDotN = small;
        dhFlag = 0.0;
    } else {
        AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in,n);
        dhFlag = 1.0;
    }
    
    n2 = SingleDot(n,n);
    Vector devStress = GetDevPart(stress);
    double p = one3 * GetTrace(stress);
    p = p < small ? m_Pmin : p;
    dpFlag = p < small ? 0.0 : 1.0;
    Vector r(6); r = devStress - p * alpha;
    double normR = GetNorm_Contr(r);
    dnFlag = normR == 0 ? 0.0 : 1.0;
    double gc = g(Cos3Theta, m_c);
    Vector aBar(6); aBar = two3 * h * b;
    Vector zBar(6); zBar = -1.0 * m_cz * Macauley(-1.0 * D) * (m_z_max * n + fabric);

    //double G, K;
    //GetElasticModuli(curStress, curVoidRatio, voidRatio, TrialElasticStrain, curEStrain, K, G);
    Matrix aD(6,6);    aD = GetCompliance(mK, mG);
    
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
    dnOverdSigma          = dnFlag * ( 1.0 / normR * (mIIdevCon - dpFlag*one3*Dyadic2_2(alpha,mI1) - 
        Dyadic2_2(n,n) + dpFlag*one3*DoubleDot2_2_Contr(alpha,n)*Dyadic2_2(n,mI1)));
    dPsiOverdSigma        = dpFlag * one3 * m_ksi * m_lambda_c / m_P_atm * pow(p/m_P_atm, m_ksi-1) * mI1;
    db0OverdSigma         = dpFlag*(-b0 / (6.0*p) * mI1);

    dCos3ThetaOverdSigma  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdSigma));
    dAdOverdSigma         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdSigma));
    dhOverdSigma          = dhFlag * (1.0 / AlphaAlphaInDotN * (db0OverdSigma - 
        h*DoubleDot2_4(alpha-alpha_in, ToCovariant(dnOverdSigma))));

    dgOverdSigma          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdSigma;

    dAlphaDOverdSigma     = m_Mc * exp(m_nd * psi) * (dgOverdSigma + m_nd * gc * dPsiOverdSigma);
    dCOverdSigma          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdSigma;
    dBOverdSigma          = 1.5 * (1.0 - m_c)/m_c * (dgOverdSigma * Cos3Theta + gc * dCos3ThetaOverdSigma);
    dAlphaBOverdSigma     = m_Mc * exp(-1.0*m_nb*psi) * (dgOverdSigma - m_nb * gc * dPsiOverdSigma);

	if (p < 0.05 * m_P_atm)
	{
		double be = 7.2713;
		double temp1 = exp(7.6349 - 7.2713 * p);
		double D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
		dDOverdSigma = (D_factor * dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)))) -
			one3 * Macauley(D) * be * temp1 / pow(1 + temp1, 2) * mI1;
	}
	else {
		dDOverdSigma = dAdOverdSigma * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n)) +
			A * (root23 * dAlphaDOverdSigma - DoubleDot2_4(alpha, ToCovariant(dnOverdSigma)));
	}
    dAbarOverdSigma       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdSigma) + 
        root23 * h * (Dyadic2_2(n, dAlphaBOverdSigma)+alphaBtheta * dnOverdSigma));

    dROverdSigma          = B * dnOverdSigma + Dyadic2_2(n, dBOverdSigma) - C * 
        (Trans_SingleDot4T_2(dnOverdSigma,n) + SingleDot2_4(n, dnOverdSigma)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdSigma) + one3 * Dyadic2_2(mI1, dDOverdSigma);
	dZbarOverdSigma       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdSigma)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdSigma;

    // d...OverdAlpha : Arranged by order of dependence
    dnOverdAlpha          = dnFlag * (p / normR * (Dyadic2_2(n,n) - mIIcon));

    dCos3ThetaOverdAlpha  = 3.0 * sqrt(6.0) * DoubleDot2_4(n2, ToCovariant(dnOverdAlpha));
    dAdOverdAlpha         = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * 
        DoubleDot2_4(fabric, ToCovariant(dnOverdAlpha));
    dhOverdAlpha          = dhFlag * (-1.0*h / AlphaAlphaInDotN * (n + 
        DoubleDot2_4(alpha-alpha_in,ToCovariant(dnOverdAlpha))));

    dgOverdAlpha          = pow(gc,2.0) * (1.0-m_c)/(2.0*m_c) * dCos3ThetaOverdAlpha;

    dAlphaDOverdAlpha     = m_Mc * exp(m_nd * psi) * dgOverdAlpha;
    dCOverdAlpha          = 3.0 * sqrt(1.5) * (1 - m_c)/m_c * dgOverdAlpha;
    dBOverdAlpha          = 1.5 * (1.0 - m_c)/m_c * (dgOverdAlpha * Cos3Theta + gc * dCos3ThetaOverdAlpha);
    dAlphaBOverdAlpha     = m_Mc * exp(-1.0*m_nb*psi) * dgOverdAlpha;

    dDOverdAlpha          = dAdOverdAlpha * (root23 * alphaDtheta - 
        DoubleDot2_2_Contr(alpha, n)) + A * (root23 * dAlphaDOverdAlpha -
        n - DoubleDot2_4(alpha, ToCovariant(dnOverdAlpha)));
    dAbarOverdAlpha       = two3 * (Dyadic2_2(root23*alphaBtheta*n-alpha, dhOverdAlpha) +
        root23 * h * (Dyadic2_2(n, dAlphaBOverdAlpha)+alphaBtheta * dnOverdAlpha) - h * mIIcon);

    dROverdAlpha          = B * dnOverdAlpha + Dyadic2_2(n, dBOverdAlpha) - C * 
        (Trans_SingleDot4T_2(dnOverdAlpha,n) + SingleDot2_4(n, dnOverdAlpha)) -
        Dyadic2_2((n2 - one3 * mI1),dCOverdAlpha) + one3 * Dyadic2_2(mI1, dDOverdAlpha);
	dZbarOverdAlpha       = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdAlpha)
		- m_cz * Macauley(-1.0 * D) * m_z_max * dnOverdAlpha;

    // d...OverdFabric : Arranged by order of dependence
    dAdOverdFabric        = m_A0 * MacauleyIndex(DoubleDot2_2_Contr(fabric, n)) * n;

    dDOverdFabric         = dAdOverdFabric * (root23 * alphaDtheta - DoubleDot2_2_Contr(alpha, n));
    
    dROverdFabric         = one3 * Dyadic2_2(mI1, dDOverdFabric);

    dZbarOverdFabric      = m_cz * MacauleyIndex(-1.0 * D) * Dyadic2_2(m_z_max*n + fabric, dDOverdFabric)
		- m_cz * Macauley(-1.0 * D) *  mIIcon;

    // Derivatives of residuals
    Matrix dR1OverdSigma(6,6), dR2OverdSigma(6,6), dR3OverdSigma(6,6); Vector dR1OverdDGamma(6);
    Matrix dR1OverdAlpha(6,6), dR2OverdAlpha(6,6), dR3OverdAlpha(6,6); Vector dR2OverdDGamma(6);
    Matrix dR1OverdFabric(6,6), dR2OverdFabric(6,6), dR3OverdFabric(6,6); Vector dR3OverdDGamma(6);
    Vector dR4OverdSigma(6), dR4OverdAlpha(6), dR4OverdFabric(6);
    double dR4OverdDGamma;
    
    dR1OverdSigma         = aD + dGamma * mIIco * ToCovariant(dROverdSigma);
    dR1OverdAlpha         =      dGamma * mIIco * ToCovariant(dROverdAlpha);
	dR1OverdFabric        =      dGamma * mIIco * ToCovariant(dROverdFabric);
    dR1OverdDGamma        = ToCovariant(R);

    dR2OverdSigma         =  -1.0*dGamma * dAbarOverdSigma * mIIco;
    dR2OverdAlpha         =  mIImix - dGamma * dAbarOverdAlpha * mIIco;
    dR2OverdFabric.Zero();
    dR2OverdDGamma        =  -1.0 *  aBar;

    dR3OverdSigma         =  -1.0*dGamma * dZbarOverdSigma * mIIco;
    dR3OverdAlpha         =  -1.0*dGamma * dZbarOverdAlpha * mIIco;
    dR3OverdFabric        =  mIImix - dGamma * dZbarOverdFabric * mIIco;
    dR3OverdDGamma        =  -1.0 * zBar;

    dR4OverdSigma         = ToCovariant(n - one3 * (DoubleDot2_2_Contr(n,alpha) + root23 * m_m) * mI1);
    dR4OverdAlpha         = -1.0 * ToCovariant(p * n);
    dR4OverdFabric.Zero();  
    dR4OverdDGamma         = 0;


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



Matrix
ManzariDafalias::GetFDMJacobian(const Vector& delta, const Vector& inVar)
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
        h = sqrt(2.0 * mEPS);
        if (h == 0.0) h = mEPS;
        x(i) = temp + h;
        h = x(i) - temp;
        f = GetResidual(x, inVar);
        x(i) = temp;
        j.Assemble(((f - fn) / h),0,i,1.0);
    }
    return j;
}



Vector 
ManzariDafalias::SetManzariComponent(const Vector& stress, const Vector& alpha,
                             const Vector& fabric, const double& dGamma)
{
    // flush the all data field
    // mSize = 19;
    // Caution: Vector::Assemble() adds the number to current values
    Vector result(19);
    result.Assemble(stress, 0);        // Stress
    result.Assemble(alpha, 6);        // Alpha
    result.Assemble(fabric, 12);    // Fabric
    result(18) = dGamma;            // DGamma
    return result;
}


Vector 
ManzariDafalias::SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, const Vector& cEStrain, 
                const Vector& cAlpha, const Vector& cFabric, const double& cVoidRatio, const double& nVoidRatio, 
                const Vector& Alpha_in)
{
    // flush the all data field
    // mSize = 44;
    // Caution: Vector::Assemble() adds the number to current values
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



double
ManzariDafalias::machineEPS()
{
    double eps = 1.0;
    while ( ((double) 1.0 + eps) > ((double) 1.0) )
        eps /= 2.0;
    return eps;
}


double ManzariDafalias::Macauley(double x)
{
    // Macauley bracket
    return (x > 0 ? x : 0.0);
}


double ManzariDafalias::MacauleyIndex(double x)
{
    // Macauley index
    return (x > 0 ? 1.0 : 0.0);
}


double 
ManzariDafalias::g(const double cos3theta, const double c)
{
    return 2 * c / ((1 + c) - (1 - c) * cos3theta);
}


double 
ManzariDafalias::GetF(const Vector& nStress, const Vector& nAlpha)
{
    // Manzari's yield function
    Vector s(6); s = GetDevPart(nStress);
    double p = one3 * GetTrace(nStress) + m_Presidual;
    s -= p * nAlpha;
    return GetNorm_Contr(s) - root23 * m_m * p;
}


double 
ManzariDafalias::GetPSI(const double& e, const double& p)
{
    return e - (m_e0 - m_lambda_c * pow((p / m_P_atm),m_ksi));
}


double 
ManzariDafalias::GetLodeAngle(const Vector& n)
// Returns cos(3*theta)
{
    double Cos3Theta = sqrt(6.0) * GetTrace(SingleDot(n,SingleDot(n,n)));
    Cos3Theta = Cos3Theta > 1 ? 1 : Cos3Theta;
    Cos3Theta = Cos3Theta < -1 ? -1 : Cos3Theta;
    return Cos3Theta;
}


void
ManzariDafalias::GetElasticModuli(const Vector& sigma, const double& en, const double& en1, const Vector& nEStrain, 
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
    if (mElastFlag == 0) 
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init);
    else
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init) * sqrt(pn / m_P_atm);
    K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


void
ManzariDafalias::GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G, const double& D)
// Calculates G, K
{
    double pn = one3 * GetTrace(sigma);
    pn = (pn <= m_Pmin) ? m_Pmin : pn;

    if (mElastFlag == 0) 
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init);
    else
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init) * sqrt(pn / m_P_atm);
    K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


void
ManzariDafalias::GetElasticModuli(const Vector& sigma, const double& en, double &K, double &G)
// Calculates G, K
{
    double pn = one3 * GetTrace(sigma);
    pn = (pn <= m_Pmin) ? m_Pmin : pn;

    if (mElastFlag == 0) 
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init);
    else
        G = m_G0 * m_P_atm * pow((2.97 - m_e_init),2) / (1 + m_e_init) * sqrt(pn / m_P_atm);
    K = two3 * (1 + m_nu) / (1 - 2 * m_nu) * G;
}


Matrix
ManzariDafalias::GetStiffness(const double& K, const double& G)
// returns the stiffness matrix in its contravarinat-contravariant form
{
    Matrix C(6,6);
    double a = K + 4.0*one3 * G;
    double b = K - 2.0*one3 * G;
    C(0,0) = C(1,1) = C(2,2) = a;
    C(3,3) = C(4,4) = C(5,5) = G;
    C(0,1) = C(0,2) = C(1,2) = b;
    C(1,0) = C(2,0) = C(2,1) = b;
    return C;
}


Matrix
ManzariDafalias::GetCompliance(const double& K, const double& G)
// returns the compliance matrix in its covariant-covariant form
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


Matrix
ManzariDafalias::GetElastoPlasticTangent(const Vector& NextStress, const double& NextDGamma, 
                    const Vector& CurStrain, const Vector& NextStrain,
                    const double& G, const double& K, const double& B, 
                    const double& C,const double& D, const double& h, 
                    const Vector& n, const Vector& d, const Vector& b) 
{    
    double p = one3 * GetTrace(NextStress) + m_Presidual;
    p = (p < small + m_Presidual) ? small + m_Presidual : p;
    // Vector r = GetDevPart(NextStress) / p;
	Vector r = GetDevPart(NextStress); r /= p;
    double Kp = two3 * p * h * DoubleDot2_2_Contr(b, n);
    
    Matrix aC(6,6), aCep(6,6);
    Vector temp0(6), temp1(6), temp2(6), R(6);
    double temp3;

    aC  = GetStiffness(K, G);
    // R = ToCovariant((B * n ) - (C * (SingleDot(n,n)-one3*mI1)) + (one3 * D * mI1));
	temp0 = n; temp0 *= B;
	temp1 = mI1; temp1 *= (-1.0 * one3); temp1 += SingleDot(n, n); temp1 *= C;
	temp2 = mI1; temp2 *= (one3 * D);
	temp0 -= temp1; temp0 += temp2;
	R = ToCovariant(temp0);

    temp1 = DoubleDot4_2(aC, ToCovariant(R));
    // temp2 = DoubleDot2_4(ToCovariant(n - one3 * DoubleDot2_2_Contr(n,r) * mI1), aC);
	temp0 = mI1; temp0 *= (-1.0 * one3 * DoubleDot2_2_Contr(n, r)); temp0 += n;
	temp0 = ToCovariant(temp0);
	temp2 = DoubleDot2_4(temp0, aC);
    temp3 = DoubleDot2_2_Contr(temp2, R) + Kp;
    if (fabs(temp3) < small) return aC;
    
    // aCep = (aC - (MacauleyIndex(NextDGamma) / temp3 * (Dyadic2_2(temp1, temp2))));
	aCep = Dyadic2_2(temp1, temp2);
	aCep *= (-1.0 * MacauleyIndex(NextDGamma) / temp3);
	aCep += aC;
    return aCep;
}


Vector
ManzariDafalias::GetNormalToYield(const Vector &stress, const Vector &alpha)
{
    // Vector devStress(6); devStress = GetDevPart(stress);

    double p = one3 * GetTrace(stress) + m_Presidual;

    Vector n(6); 
    if (fabs(p) < small)
    {
        n.Zero();
    } else {
        // n = devStress - p * alpha;
        // double normN = GetNorm_Contr(n);
        // normN = (normN < small) ? 1.0 : normN;
        // n = n / normN;
		n = alpha; n *= (-p);
		n += GetDevPart(stress);
		double normN = GetNorm_Contr(n);
		normN = (normN < small) ? 1.0 : normN;
		n /= normN;
    }

    return n;
}


int
ManzariDafalias::Check(const Vector& TrialStress, const Vector& stress, const Vector& CurAlpha, const Vector& NextAlpha)
// Check if the solution of implicit integration makes sense
{
    int result = 1;
    
    if (GetTrace(stress) < 0) 
    {
        if(debugFlag)
            opserr << "p < 0 !!!!!! Check this : ManzariDafalias::Check()" << endln;
        //result = -2;
    }
    
    Vector n(6);    n    = GetNormalToYield(stress, CurAlpha);
    Vector n_tr(6); n_tr = GetNormalToYield(TrialStress, CurAlpha);
    
    // check the direction of stress and trial stress
    if (DoubleDot2_2_Contr(n, n_tr) < 0) 
    {
        if(debugFlag)        
            opserr << "Direction of n and n_tr are more than 90 degrees apart!" << endln;
        result = -4;
    }
    
    // add any other checks here
    
    return result;
}


void 
ManzariDafalias::GetStateDependent(const Vector &stress, const Vector &alpha, const Vector &fabric
                , const double &e, const Vector &alpha_in, Vector &n, Vector &d, Vector &b
                , double &cos3Theta, double &h, double &psi, double &alphaBtheta
                , double &alphaDtheta, double &b0, double& A, double& D, double& B
                , double& C, Vector& R)
{
	Vector tmp0(6), tmp1(6);
    double D_factor = 1.0;
    double p = one3 * GetTrace(stress) + m_Presidual;
    p = (p < small) ? small : p;

    n = GetNormalToYield(stress, alpha);

    double AlphaAlphaInDotN;
    // AlphaAlphaInDotN = DoubleDot2_2_Contr(alpha - alpha_in,n);
	tmp0 = alpha; tmp0 -= alpha_in;
	AlphaAlphaInDotN = DoubleDot2_2_Contr(tmp0, n);

    psi = GetPSI(e, p);

    cos3Theta = GetLodeAngle(n);

    alphaBtheta = g(cos3Theta, m_c) * m_Mc * exp(-1.0 * m_nb * psi) - m_m;
    
    alphaDtheta = g(cos3Theta, m_c) * m_Mc * exp(m_nd * psi) - m_m;

    b0 = m_G0 * m_h0 * (1.0 - m_ch * e) / sqrt(p / m_P_atm);
    
    // d    = root23 * alphaDtheta * n - alpha;
	d = n; d *= (root23 * alphaDtheta); d -= alpha;
    // b    = root23 * alphaBtheta * n - alpha;
	b = n; b *= (root23 * alphaBtheta); b -= alpha;
        
	if (fabs(AlphaAlphaInDotN) < small)
		h = 1.0e10;
	else
		h = b0 / AlphaAlphaInDotN;

    A = m_A0 * (1 + Macauley(DoubleDot2_2_Contr(fabric, n)));

    D = A * DoubleDot2_2_Contr(d, n);

    // Apply a factor to D so it doesn't go very big when p is small
    if (p < 0.05 * m_P_atm)
    {
        D_factor = 1.0 / (1.0 + (exp(7.6349 - 7.2713 * p)));
    } else {
        D_factor = 1.0;
    }

    D *= D_factor;

    B = 1.0 + 1.5 * (1 - m_c)/ m_c * g(cos3Theta, m_c) * cos3Theta;

    C = 3.0 * sqrt(1.5) * (1 - m_c)/ m_c * g(cos3Theta, m_c);

    // R = B * n - C * (SingleDot(n,n) - one3 * mI1) + one3 * D * mI1;
	R = n; R *= B;
	tmp0 = mI1; tmp0 *= (-1.0 * one3); tmp0 += SingleDot(n, n); tmp0 *= C;
	tmp1 = mI1; tmp1 *= (one3 * D);
	R -= tmp0; R += tmp1;
}

int
ManzariDafalias::Elastic2Plastic()
{
    if ((GetTrace(mSigma) < 3.0 * m_Pmin) || (GetTrace(mSigma_n) < 3.0 * m_Pmin))
    {
        mSigma = mSigma_n = m_Pmin * mI1;
        mAlpha.Zero();
        mAlpha_n.Zero();
        return 0;
    }
    double curM = sqrt(1.5) * GetNorm_Contr(GetDevPart(mSigma));
    double p = one3 * GetTrace(mSigma) + m_Presidual;
    curM = curM / p;
    if (curM > m_Mc)
    {  
         m_Mc = 1.1 * curM;
    //     opserr << "Outside Bounding!" << endln;
    //     opserr << "Before = " << mSigma;
    //     mAlpha = mAlpha_n = (m_Mc - m_m) / curM * GetDevPart(mSigma) / p;
    //     mSigma = mSigma_n = one3 * GetTrace(mSigma) * mI1 + (m_Mc / curM) * GetDevPart(mSigma);
    //     mAlpha_in = mAlpha_in_n = mAlpha;
    //     opserr << "After = " << mSigma << endln;
    }
    // Vector n(6), d(6), b(6), R(6);
    // double cos3Theta, h, psi, aB, aD, b0, A, D, B, C;


    // GetStateDependent(mSigma, mAlpha, mFabric, mVoidRatio, mAlpha_in, 
    //           n, d, b, cos3Theta, h, psi, aB, aD, b0, A, D, B, C, R);
    //  double q = sqrt(1.5 * DoubleDot2_2_Contr(GetDevPart(mSigma), GetDevPart(mSigma)));
    //  opserr << "Committed stress (tag = " << this->getTag() << ") = " << mSigma << "Yield = " << GetF(mSigma, mAlpha) << endln << "p = " << p << ", q = " << q << ", eta = " << q/p << endln;

    //  opserr << "psi = " << psi << endln;
    //  opserr << "alpha_b = " << aB << endln;
    //  opserr << "alpha_d = " << aD << endln;
    //  opserr << "b0 = " << b0 << endln;
    //  opserr << "d = " << d;
    //  opserr << "b = " << b;
    //  opserr << "h = " << h << endln;
    //  opserr << "A = " << A << endln;
    //  opserr << "D = " << D << endln;
    //  opserr << "B = " << B << endln;
    //  opserr << "C = " << C << endln;
    //  opserr << "R = " << R;
    //  opserr << "n = " << n;
    //  opserr << endln;
     return 0;
}













/*************************************************************/
/*************************************************************/
//            SYMMETRIC TENSOR OPERATIONS                    //
/*************************************************************/
/*************************************************************/
// In all the functions below, contravariant means a stress-like tensor
// and covariant means a strain-like tensor

double
ManzariDafalias::GetTrace(const Vector& v) 
// computes the trace of the input argument
{
    if (v.Size() != 6)
        opserr << "\n ERROR! ManzariDafalias::GetTrace requires vector of size(6)!" << endln;

    return (v(0) + v(1) + v(2));
}

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

Matrix
ManzariDafalias::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
// caution: second coordinate of the first matrix should be in opposite 
// variant form of the first coordinate of second matrix
{
    if ((m1.noCols() != 6) || (m1.noRows() != 6) || (m2.noCols() != 6) || (m2.noRows() != 6)) 
        opserr << "\n ERROR! ManzariDafalias::DoubleDot4_4 requires 6-by-6 matrices " << endln;

    return m1*m2;
}

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

double ManzariDafalias::Det(const Vector& aV)
{
    if (aV.Size() != 6)
        opserr << "\n ERROR! ManzariDafalias::Det requires vector of size(6)!" << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    return (     aV[0] * aV[1] * aV[2] 
             + 2*aV[3] * aV[4] * aV[5] 
             -   aV[0] * aV[5] * aV[5] 
             -   aV[2] * aV[3] * aV[3] 
             -   aV[1] * aV[4] * aV[4]);
}

Vector ManzariDafalias::Inv(const Vector& aV)
{
    if (aV.Size() != 6)
        opserr << "\n ERROR! ManzariDafalias::Inv requires vector of size(6)!" << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    double det = Det(aV);
    if (det == 0) 
    {
        opserr << "\n Error! ManzariDafalias::Inv - Singular tensor - return 0 tensor" << endln;
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

Vector ManzariDafalias::ToContraviant(const Vector& v1)
{
    if (v1.Size() != 6)
        opserr << "\n ERROR! ManzariDafalias::ToContraviant requires vector of size(6)!" << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    Vector res = v1;
    res(3) *= 0.5;
    res(4) *= 0.5;
    res(5) *= 0.5;
    
    return res;
}

Vector ManzariDafalias::ToCovariant(const Vector& v1)
{
    if (v1.Size() != 6)
        opserr << "\n ERROR! ManzariDafalias::ToCovariant requires vector of size(6)!" << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    Vector res = v1;
    res(3) *= 2.0;
    res(4) *= 2.0;
    res(5) *= 2.0;
    
    return res;
}

Matrix ManzariDafalias::ToContraviant(const Matrix& m1)
{
    if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
        opserr << "\n ERROR! ManzariDafalias::ToContraviant requires 6-by-6 matrix " << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    Matrix res = m1;
    for (int ii = 0; ii < 6; ii++)
    {
        res(3,ii) *= 0.5;
        res(4,ii) *= 0.5;
        res(5,ii) *= 0.5;
    }
    
    return res;
}

Matrix ManzariDafalias::ToCovariant(const Matrix& m1)
{
    if ((m1.noCols() != 6) || (m1.noRows() != 6)) 
        opserr << "\n ERROR! ManzariDafalias::ToCovariant requires 6-by-6 matrix " << endln;
    // aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
    Matrix res = m1;
    for (int ii = 0; ii < 6; ii++)
    {
        res(3,ii) *= 2.0;
        res(4,ii) *= 2.0;
        res(5,ii) *= 2.0;
    }
    
    return res;
}

// send back the strain
const Vector& 
ManzariDafalias::getEStrain() 
{
    opserr << "ManzariDafalias::getEStrain - base class function called. This is an error\n ";
    return mEpsilonE; 
} 

const Vector& 
ManzariDafalias::getPStrain() 
{
    opserr << "ManzariDafalias::getPStrain - base class function called. This is an error\n ";
    static Vector result(6);
    result = mEpsilon - mEpsilonE;
    return result; 
} 


