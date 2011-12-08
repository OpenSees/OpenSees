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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the implementation of the ManzariDafalias class.

#include <ManzariDafalias.h>
#include <ManzariDafalias3D.h>
#include <ManzariDafaliasPlaneStrain.h>

#include <Information.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

const double ManzariDafalias::one3   = 1.0/3.0 ;
const double ManzariDafalias::two3   = 2.0/3.0;
const double ManzariDafalias::root23 = sqrt(2.0/3.0);
const double ManzariDafalias::PI = 4.0 * atan(1.0);

#include <elementAPI.h>

static int numManzariDafaliasMaterials = 0;

void *
OPS_NewManzariDafaliasMaterial(void)
{
  if (numManzariDafaliasMaterials == 0) {
    numManzariDafaliasMaterials++;
    opserr << "ManzariDafalias nDmaterial - Written: P.Arduino, C.McGann, U.Washington\n";
  }

  NDMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 23) {
    opserr << "Want: nDMaterial ManzariDafalias tag? " << endln;
    return 0;	
  }
  
  int tag;
  double dData[22];

  int numData = 1;
  if (OPS_GetInt(&numData, &tag) != 0) {
    opserr << "WARNING invalid nDMaterial ManzariDafalias material  tag" << endln;
    return 0;
  }
  numData = 22;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING invalid material data for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
    return 0;
  }

  theMaterial = new ManzariDafalias(tag, 0, dData[0] , dData[1] , dData[2] , dData[3] , dData[4] , dData[5] , 
                                            dData[6] , dData[7] , dData[8] , dData[9] , dData[10], dData[11],
                                            dData[12], dData[13], dData[14], dData[15], dData[16], dData[17],
											dData[18], dData[19], dData[20], dData[21]);
  
  if (theMaterial == 0) {
    opserr << "WARNING ran out of memory for nDMaterial ManzariDafalias material  with tag: " << tag << endln;
  }

  return theMaterial;
}

// full constructor
ManzariDafalias::ManzariDafalias(int tag, int classTag, double Ko, double Go, double v, double b, double Patm,
	                                                    double Ao, double ho, double Cm, double Me, double Mc,
														double kBE, double kBC, double kDE, double kDC, double ecRef,
														double lambda, double Pref, double m, double Fmax, double Cf,
														double eo, double mDen)
  : NDMaterial(tag,ND_TAG_ManzariDafalias),
  mParam(21),
  mEpsilonE(6),
  mAlpha(6),
  mFabric(6),
  mCe(6,6),
  mCep(6,6),
  mState(7),
  mEpsilon(6),
  mEpsilon_n(6),
  mSigma(6),
  mSigma_n(6),
  mI1(6),
  mIIco(6,6),
  mIIcon(6,6),
  mIImix(6,6),
  mIIvol(6,6),
  mIIdevCon(6,6),
  mIIdevMix(6,6)
{
	massDen = mDen;
	//Add all other member parameters
	mParam(0)  = Ko;      // Ko
	mParam(1)  = Go;      // Go
	mParam(2)  = v;       // v
	mParam(3)  = b;       // b
	mParam(4)  = Patm;    // Patm
	mParam(5)  = Ao;      // Ao
	mParam(6)  = ho;      // ho
	mParam(7)  = Cm;      // Cm
	mParam(8)  = Me;      // Me
	mParam(9)  = Mc;      // Mc
	mParam(10) = kBE;     // kBE
	mParam(11) = kBC;     // KBC
	mParam(12) = kDE;     // KDE
	mParam(13) = kDC;     // KDC
	mParam(14) = ecRef;   // ecRef
	mParam(15) = lambda;  // lamda
	mParam(16) = Pref;    // Pref
	mParam(17) = m;       // m
	mParam(18) = Fmax;    // Fmax
	mParam(19) = Cf;      // Cf
	mParam(20) = eo;      // eo

	/*
	// Possible fixed values
	mParam(0)  = 31400;   // Ko
	mParam(1)  = 23550;   // Go
	mParam(2)  = 0.2;     // v
	mParam(3)  = 0.5;     // b
	mParam(4)  = 100;     // Patm
	mParam(5)  = 2.64;    // Ao
	mParam(6)  = 1200;    // ho
	mParam(7)  = 0.0;     // Cm
	mParam(8)  = 1.14;    // Me
	mParam(9)  = 1.3;     // Mc
	mParam(10) = 2.0;     // kBE
	mParam(11) = 3.975;   // KBC
	mParam(12) = 0.07;    // KDE
	mParam(13) = 4.2;     // KDC
	mParam(14) = 0.8;     // ecRef
	mParam(15) = 0.025;   // lamda
	mParam(16) = 160.0;   // Pref
	mParam(17) = 0.05;    // m
	mParam(18) = 100;     // Fmax
	mParam(19) = 100;     // Cf
	mParam(20) = 0.81;    // eo
	*/
	this->initialize();
}

   
// null constructor
ManzariDafalias ::ManzariDafalias() 
  : NDMaterial()
{
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
		clone = new ManzariDafaliasPlaneStrain(this->getTag(), mParam(0),  mParam(1),  mParam(2),  mParam(3),  mParam(4),
															   mParam(5),  mParam(6),  mParam(7),  mParam(8),  mParam(9),
															  mParam(10), mParam(11), mParam(12), mParam(13), mParam(14),
															  mParam(15), mParam(16), mParam(17), mParam(18), mParam(19),
															  mParam(20), massDen);
		return clone;
	} else if (strcmp(type,"ThreeDimensional")==0 || strcmp(type,"3D") ==0) {
		ManzariDafalias3D *clone;
     	clone = new ManzariDafalias3D(this->getTag(), mParam(0),  mParam(1),  mParam(2),  mParam(3),  mParam(4),
													  mParam(5),  mParam(6),  mParam(7),  mParam(8),  mParam(9),
													 mParam(10), mParam(11), mParam(12), mParam(13), mParam(14),
													 mParam(15), mParam(16), mParam(17), mParam(18), mParam(19),
													 mParam(20), massDen);
	 	return clone;
  	} else {
	  	opserr << "ManzariDafalias::getCopy failed to get copy: " << type << endln;
	  	return 0;
  	}
}

int 
ManzariDafalias::commitState(void)
{
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

/*************************************************************/
// 2nd Order Tensor Operations
/*************************************************************/
double
ManzariDafalias::GetContraNorm(const Vector& v)
// computes contravariant (stress-type) norm of input 6x1 tensor
{
	double result=0.0;	
	for (int i = 0; i < 3; i++) {
		result += v(i)*v(i);
	}
	for (int i = 3; i < 6; i++) {
		result += 2.0*v(i)*v(i);
	}

	return sqrt(result);
}
		
double
ManzariDafalias::GetCovariantNorm(const Vector& v) 
// computes the norm of the input argument (for strain-type storage)
{
	if (v.Size() != 6) {
		opserr << "\n ERROR! BoundingCamClay::NormEngStrain requires vector of size(6)!" << endln;
	}
    double result=0.0;
	for (int i = 0; i < 3; i++) {
		result += v(i)*v(i);
	}
	for (int i = 3; i < 6; i++) {
		result += 0.5*v(i)*v(i);
	}

    return sqrt(result);
}

double
ManzariDafalias::GetTrace(const Vector& v) 
// computes the trace of the input argument
{
	if (v.Size() != 6)
		opserr << "\n ERROR! BoundingCamClay::GetTrace requires vector of size(6)!" << endln;

	return (v(0) + v(1) + v(2));
}

//================================================
//  GetDevPart()
//================================================
Vector ManzariDafalias::GetDevPart(const Vector& aV)
{
	if (aV.Size() != 6)
	opserr << "\n ERROR! ManzariDafalias::GetDevPart requires vector of size(6)!" << endln;
	Vector result(6);
	double p = GetTrace(aV);
	result = aV;
	result(0) -= p/3.0;
	result(1) -= p/3.0;
	result(2) -= p/3.0;

	return result;
}

//================================================
//  GetJ2() const
//================================================
double ManzariDafalias::GetJ2(const Vector& aV)
{
	Vector devTensor(6);
	devTensor = GetDevPart(aV);
	double n = 	GetContraNorm(devTensor);
	return 0.5 * n * n;
}

//================================================
//  GetJ3() const
//================================================
double ManzariDafalias::GetJ3(const Vector& aV)
{
	Vector devTensor(6);
	devTensor = GetDevPart(aV);
	return (Det(devTensor));
}

//================================================
//  GetLodeAngle() const
//================================================
double ManzariDafalias::GetLodeAngle(const Vector& aV)
{
	double term = 1.5 * sqrt(3.0) * GetJ3(aV) / pow(GetJ2(aV), 1.5);
	if (term >1.0)	term = 1.0;
	if (term <-1.0) term = -1.0;
	double threeTheta = acos(term);
	return threeTheta / 3.0;
}

//================================================
//  Det
//================================================
double ManzariDafalias::Det(const Vector& aV)
{
	// aV(i) -> T(i,j) 1 = 11, 2=22, 3=33, 4=12, 5=23, 6=13
	return( -   aV[6] * aV[6] * aV[2] 
		    + 2*aV[4] * aV[6] * aV[5] 
			-   aV[1] * aV[5] * aV[5] 
			-   aV[4] * aV[4] * aV[3] 
			+ 	aV[1] * aV[2] * aV[3]);
}

//================================================
//  g - elliptical function
//================================================
//--- Counter-clockwise 60 dgree rotation of Willam
double ManzariDafalias::g(const double Theta, const double c)
{
	double aNumer, aDeno;
	double term1;
	term1  = 4.0 * (1.0 - c*c) * cos(Theta - PI/3.0) * cos(Theta - PI/3.0) + 5.0 * c * c - 4.0 * c;
	aNumer = 2.0 * (1.0 - c*c) * cos(Theta - PI/3.0) + (2.0*c - 1.0) * sqrt(term1);
	aDeno  = 4.0 * (1.0 - c*c) * cos(Theta - PI/3.0) * cos(Theta - PI/3.0) + (1.0 - 2.0*c) * (1.0 - 2.0*c);
	return aNumer / aDeno;
}

/*************************************************************/
//            g  (Manzari)                                   //
/*************************************************************/
/*double 
ManzariDafalias::g(const double cos3theta, const double c)
{
	return 2 * c / ((1 + c) - (1 - c) * cos3theta);
}*/


/*************************************************************/
//            Macauley                                       //
/*************************************************************/
double ManzariDafalias::Macauley(double x)
{
	// Macauley bracket
	double result;
	if (x > 0.0)
		result = x;
	else
		result = 0;
	return result;
}

/*************************************************************/
//            MacauleyIndex                                  //
/*************************************************************/
int ManzariDafalias::MacauleyIndex(double x)
{
	int result;
	result = 0;
	if (x > 0.0)
	    result = 1;
	else
		result = 0;
	return result;
}

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

/*************************************************************/
// Tensor Products
/*************************************************************/
Vector ManzariDafalias::Dot(const Vector& v1)
{
    Vector result(6);
	for (int i = 0; i < v1.Size(); i++) {
		result(i)= v1(i) * v1(i);
	}
	return result;
}

double
ManzariDafalias::DoubleDot2_2(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments
{
	double result = 0.0;
	
	if (v1.Size() != v2.Size()) {
		opserr << "\n ERROR! BoundingCamClay::DoubleDot2_2 function requires vectors of equal size!" << endln;
	}

	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i);
	}

	return result;
}

Vector
ManzariDafalias::DoubleDot2_4(const Vector& v1, const Matrix& m1)
// computes doubledot product for vector-matrix arguments
{
	Vector result(6);
	result.Zero();

	if (v1.Size() != m1.noRows()) {
		opserr << "\n ERROR! BoundingCamClay::DoubleDot2_4 function requires Size(v1) = noRows(m1) " << endln;
	}

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(j) += v1(i) * m1(i,j);
	}

	return result;
}

Vector
ManzariDafalias::DoubleDot4_2(const Matrix& m1, const Vector& v1)
// computes doubledot product for matrix-vector arguments
{
	Vector result(6);
	result.Zero();

	if (m1.noCols() != v1.Size()) {
		opserr << "\n ERROR! BoundingCamClay::DoubleDot4_2 function requires noCols(m1) = Size(v1) " << endln;
	}

	for (int i = 0; i < m1.noRows(); i++) {
		for (int j = 0; j < m1.noCols(); j++) 
			result(i) += m1(i,j) * v1(j);
	}

	return result;
}

Matrix
ManzariDafalias::DoubleDot4_4(const Matrix& m1, const Matrix& m2)
// computes doubledot product for matrix-matrix arguments
{
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < m1.noRows(); i ++) {
		for (int j = 0; j < m2.noCols(); j++) {
			for (int k = 0; k<m1.noRows(); k++) 
				result(i,j) += m1(i,k) * m2(k,j);
		}
	}
	return result;
}

Matrix 
ManzariDafalias::Dyadic2_2(const Vector& v1, const Vector& v2)
// computes dyadic product for two vector-storage arguments
{
	Matrix result(6,6);
	result.Zero();

	for (int i = 0; i < v1.Size(); i++) {
		for (int j = 0; j < v2.Size(); j++) 
			result(i,j) = v1(i) * v2(j);
	}
	
	return result;
}

// -------------------------------------------------------------------------------------------------------
/*************************************************************/
// Plastic Integrator
/*************************************************************/
void ManzariDafalias::plastic_integrator() 
{
	Vector CurStress(6);
	Vector NextStress(6);
	Vector CurStrain(6);
	Vector NextStrain(6);
	Vector CurElasticStrain(6);
	Vector NextElasticStrain(6);
	Matrix aC(6,6);
	Matrix aCep(6,6), aCepPart(6,6);
	Vector NextAlpha(6);
	Vector NextFabric(6);
	double NextM;
	double NextDGamma;

	NextStrain	= mEpsilon;
	CurStrain	= mEpsilon_n;
	CurElasticStrain  = mEpsilonE;
	NextElasticStrain  = CurElasticStrain + (NextStrain - CurStrain);
	CurStress	= mSigma_n;
	NextAlpha	= mAlpha;
	NextFabric	= mFabric;
	NextM		= mM;
	NextDGamma	= 0.0;
	//Trial Elastic Stress
	NextStress = HypoElastic(CurStress, NextElasticStrain, CurElasticStrain, aC);
	//In case of pure elastic response
	aCep = aC;
	mCep = aCep;
	//Trial yield function
	double NextF = F(NextStress, NextAlpha, NextM);
	
	if (NextF > mTolF)
	{
		Vector Delta0(20), InVariants(37), Delta(20);
		Delta0 = SetManzariComponent(NextElasticStrain, NextAlpha, NextM, NextFabric, NextDGamma);
		InVariants = SetManzariStateInVar(NextStrain, CurStrain, CurStress, CurElasticStrain, NextAlpha, NextM, NextFabric);
	    Delta = NewtonSolve(Delta0, InVariants, aCepPart);

		NextElasticStrain.Extract(Delta, 0, 1.0);
		NextAlpha.Extract(Delta, 6, 1.0);
		NextM = Delta(12);
		NextFabric.Extract(Delta, 13, 1.0);
		NextDGamma = Delta(19);
		// Converged ElasticModulus and Stress
		NextStress = HypoElastic(CurStress, NextElasticStrain, CurElasticStrain, aC);
		// Old Sign
		//*(aNextState.Cep) = -DoubleDot(aCepPart, aC);
		// New sign
		//*(aNextState.Cep) = DoubleDot(aCepPart, aC);
	}

	//update State variables
	mEpsilonE = NextElasticStrain;
	mEpsilon_n = NextStrain;
	mSigma_n = NextStress;
	mSigma = NextStress;
	mAlpha = NextAlpha;
	mFabric = NextFabric;
	mM = NextM;
	mDGamma = NextDGamma;
	mCe  = aC;
	//mCep = DoubleDot4_4(aCepPart, aC);
	mCep = aC;
	return;
}


/*************************************************************/
//            HypoElastic                                    //
/*************************************************************/
Vector 
ManzariDafalias::HypoElastic(const Vector& cStress, const Vector& nEStrain, 
					 const Vector& cEStrain, Matrix& C)
{
	double K, G;
	// Material Parameters
	double Ko = mParam(0), Go = mParam(1), v = mParam(2), b = mParam(3), Patm = mParam(4);
	double curP = GetTrace(cStress)/3.0;
	if (curP < Patm/100) curP=Patm/100; 
	Vector curDevStress(6); curDevStress = cStress  - curP*mI1;
	double volEStrain = GetTrace(nEStrain) - GetTrace(cEStrain);
	if (fabs(volEStrain) < 1e-13)
		K = Ko * (pow(curP, b)) / (pow(Patm, b));
	else
	{
		double term1 = pow(curP, (1 - b));
		double term2 = (1 - b) * Ko * volEStrain / (pow(Patm , b));
		double p  = pow((term1 + term2), 1 / (1 - b));
		K = (p -curP) / volEStrain;
	}
	G = 1.5 * K * (1.0 - 2.0 * v) / (1.0 + v);
	//Set Elastic Components
	mK = K;
	mG = G;
	C = 3.0*K*mIIvol + 2.0*G*mIIdevCon;

	Vector result;
	result = cStress + DoubleDot4_2(C, nEStrain - cEStrain);
	//result = cStress +(3.0*K*volEStrain*mI1+2.0*G*mIIdevCon*(nEStrain - cEStrain);
	return result;
}

 /*************************************************************/
//            GetF											  //
/*************************************************************/
double 
ManzariDafalias::F(Vector& nextStress, Vector& nAlpha, double m)
{
	// Manzari's yield function
	Vector s(6); s = GetDevPart(nextStress);
	double p = GetTrace(nextStress)/3.0;
	Vector r(6); r = s - p * nAlpha;
	double f = GetContraNorm(r) - root23 * m * p;
	return f;
}

// zero internal variables
void ManzariDafalias::initialize()
{
	// strain and stress terms
	mEpsilon.Zero();
	mEpsilon_n.Zero();
	mSigma.Zero();
	mSigma_n.Zero();

	mEpsilonE.Zero();
	mAlpha.Zero();
	mDGamma = 0.0;
	mFabric.Zero();
	mState.Zero();
	mCe.Zero();
	mCep.Zero();
	mM = mParam(17);
	mEPS = machineEPS();

	mTolF = 1.0E-5;
	mTolR = 1.0E-5;

	initializeState = true;

	// 2nd-order Identity Tensor
	mI1.Zero();
	mI1(0) = 1.0;
	mI1(1) = 1.0;
	mI1(2) = 1.0;

	// 4th-order mixed variant identity
	mIImix.Zero();
	for (int i = 0; i<6; i++) {
		mIImix(i,i) = 1.0;
	}

	// 4th-order covariant identity
	mIIco = mIImix;
	mIIco(3,3) = 2.0;
	mIIco(4,4) = 2.0;
	mIIco(5,5) = 2.0;

	// 4th-order contravariant identity
	mIIcon = mIImix;
	mIIcon(3,3) = 0.5;
	mIIcon(4,4) = 0.5;
	mIIcon(5,5) = 0.5;

	// 4th-order Volumetric Tensor, IIvol = I1 tensor I1
	mIIvol.Zero();
	for (int i = 0; i<3; i++) {
		mIIvol(i,0) = 1.0;
		mIIvol(i,1) = 1.0;
		mIIvol(i,2) = 1.0;
	}

	// 4th-order contravariant deviatoric tensor
	mIIdevCon = mIIcon - one3*mIIvol;

	// 4th-order mixed variant deviatoric tensor
	mIIdevMix = mIImix - one3*mIIvol;

}


/*************************************************************/
// SetComponent for Manzari Model
/*************************************************************/
Vector 
ManzariDafalias::SetManzariComponent(const Vector& eStrain, const Vector& alpha, const double& m,
							 const Vector& fabric, const double& dGamma)
{
	// flush the all data field
	//mSize = 20;
	Vector result(20);
	result.Zero();
	// Assign new data
	for (int i=0; i < 6; i++)      // Elastic Strain
		result[i] = eStrain[i];
	for (int j=6; j < 12; j++)     // Alpha
		result[j] = alpha[j-6];
	result[12] = m;                //m
	for (int k=13; k < 19; k++)    // Fabric
		result[k] = fabric[k-13];
	result[19] = dGamma;           // DGamma
	return result;
}

/*************************************************************/
// Set Manzari State Invariants Component for Manzari Model
/*************************************************************/
Vector 
ManzariDafalias::SetManzariStateInVar(const Vector& nStrain, const Vector& cStrain, const Vector& cStress, 
									  const Vector& cEStrain, const Vector& cAlpha, const double& cM,
									  const Vector& cFabric)
{
	// flush the all data field
	//mSize = 37;
	Vector result(37);
	result.Zero();
	// Assign new data
	for (int i=0; i < 6; i++)      // Next Strain
		result[i] = nStrain[i];
	for (int j=6; j < 12; j++)     // Cur Strain
		result[j] = cStrain[j-6];
	for (int k=12; k < 18; k++)    // Cur Stress
		result[k] = cStress[k-12];
	for (int l=18; l < 24; l++)    // Cur Elastic Strain
		result[l] = cEStrain[l-18];
	for (int m=24; m < 30; m++)    // Cur Alpha
		result[m] = cAlpha[m-24];
    result[30] = cM;               //Cur m
	for (int n=31; n < 37; n++)    // Cur Fabric
		result[n] = cFabric[n-31];
	return result;
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
	Vector sol(20), R(20), dX(20);
	Matrix  jaco(20,20), jInv(20,20);
	sol = xo;
	R = GetResidual(sol, inVar);
	while ((R.Norm() > mTolR) && (mIter <= MaxIter))
	{
		mIter++;
		jaco = GetJacobian(sol, inVar);
		jaco.Invert(jInv);
		dX   = -1.0*jInv * R;
		sol += dX;
		R    = GetResidual(sol, inVar);
	}
	aCepPart.Extract(jInv, 0, 0, 1.0);
	return sol;
}


Vector
ManzariDafalias::GetResidual(const Vector& x, const Vector& inVar)
{
	// Inside this function, if head of variable doesn't have 'Cur', it indicate 'Next' State
	// Initialization of Variables to be used
	// Material Parameters
	double Ao = mParam(5), ho = mParam(6), Cm = mParam(7), Fmax = mParam(18), Cf = mParam(19), eo = mParam(20);

	Vector g(20);   // Residual Vector
	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6); // Strain
	Vector stress(6), alpha(6), curStress(6), curAlpha(6);
	Vector fabric(6), curFabric(6);
	double m, dGamma, curM;
	// Next Iteration variable setting
	eStrain.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	m = x(12);
	fabric.Extract(x, 13, 1.0);
	dGamma = x(19);
	// Current Iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curM = inVar(30);
	curFabric.Extract(inVar, 31, 1.0);
	// Hypo elastic Stress calculation
	Matrix aC(6,6);  // Actually not use in here but for function argument compatibility
	stress = HypoElastic(curStress, eStrain, curEStrain, aC);
	// State Dependent variables
	Vector n(6), d(6), b(6); 
	double bref, psi, Theta, alphaBtheta;
	StateDepend(stress, alpha, m, strain, n, d, b, bref, psi, Theta, alphaBtheta);
	double A = Ao * (1 + Macauley(DoubleDot2_2(fabric, n)));
	double D = A * DoubleDot2_2(d, n);
	Vector R(6); R = n + D * mI1 / 3.0;
	double h = ho * fabs(DoubleDot2_2(b,n)) / (bref - fabs(DoubleDot2_2(b, n)));
	Vector aBar(6); aBar = h * b;
	double mBar = Cm * (1 + eo) * D;
	// Get Residual
	Vector TrialEStrain(6); TrialEStrain = curEStrain + (strain - curStrain);
	Vector g1(6); Vector g2(6); Vector g4(6); double g3, g5;

	g1 = eStrain - TrialEStrain + dGamma * R;
	g2 = alpha   - curAlpha     - dGamma * aBar;
	g3 = m       - curM         - dGamma * mBar;
	g4 = fabric  - curFabric    + dGamma * Cf * Macauley(-D) * (Fmax * n + fabric);
	g5 = F(stress, alpha, m);

	// Arrange Residual as Row vector
	g.Assemble(g1,  0, 1.0);
	g.Assemble(g2,  6, 1.0);
	g(12) = g3;
	g.Assemble(g4, 13, 1.0);
	g(19) = g5;
	return g;
}

/*************************************************************/
//            SetStateDepend                                 //
/*************************************************************/
void 
ManzariDafalias::StateDepend( const Vector &stress, const Vector &alpha, const double &m
							, const Vector &strain, Vector &n, Vector &d, Vector &b
							, double &bref, double &psi, double &Theta, double &alphaBtheta)
{
	double Me = mParam(8), Mc = mParam(9), kBE = mParam(10), kBC = mParam(11), kDE = mParam(12), kDC = mParam(13);
	Vector devStress(6); devStress = GetDevPart(stress);
	double p = GetTrace(stress)/3.0;
	Vector r(6); r = devStress - p * alpha;
	Vector rBar(6); rBar = r / p;
	psi = PSI(strain, p);
	double cos3Theta = cos(3 * GetLodeAngle(rBar));
	Theta = acos(cos3Theta) / 3.0;
	//cosTheta = cos(aTheta);
	double c  = Me / Mc;
	double cb = kBE / kBC;
	double cd = kDE / kDC;
	alphaBtheta = g(Theta, c) * Mc + g(Theta, cb) * kBC * Macauley(-psi) - m;
	double alphaDtheta = g(Theta, c) * Mc + g(Theta, cd) * kDC * psi - m;
	// obtain Return variables
	n    = r / GetContraNorm(r);
	d    = root23 * alphaDtheta * n - alpha;
	b    = root23 * alphaBtheta * n - alpha;
	bref = 2.0 * root23 * alphaBtheta;
}

/*************************************************************/
//            GetPSI                                         //
/*************************************************************/
double 
ManzariDafalias::PSI(const Vector& strain, const double& p)
{
	double ecRef = mParam(14), lamda = mParam(15), Pref = mParam(16), eo = mParam(20);

	double ecs = ecRef - lamda * log(p / Pref);
	double volStrain = GetTrace(strain);
	double e = eo - (1 + eo) * volStrain;
	double debugPSI = e - ecs;
	return (e - ecs);
}

/////////////////////////////////////////////////////////
///   From Here Jacobain Related Functions         //////
/////////////////////////////////////////////////////////
/*************************************************************/
//            GetJacobian                                     //
/*************************************************************/
Matrix 
ManzariDafalias::GetJacobian(const Vector &x, const Vector &inVar)
{
	// Next Iteration Variable setting
	double Ao = mParam(5), ho = mParam(6), Cm = mParam(7), Fmax = mParam(18), Cf = mParam(19), eo = mParam(20);
	double Me = mParam(8), Mc = mParam(9), kBE = mParam(10), kBC = mParam(11), kDE = mParam(12), kDC = mParam(13);
	Vector eStrain(6), strain(6), curStrain(6), curEStrain(6); // Strain
	Vector stress(6) , alpha(6) , curStress(6), curAlpha(6);     // Stress Quantities
	Vector  fabric(6), curFabric(6);
	double m, dGamma, curM;
	// Next Iteration variable setting
	eStrain.Extract(x, 0, 1.0);
	alpha.Extract(x, 6, 1.0);
	m = x(12);
	fabric.Extract(x, 13, 1.0);
	dGamma = x(19);
	// Current Iteration invariants
	strain.Extract(inVar, 0, 1.0);
	curStrain.Extract(inVar, 6, 1.0);
	curStress.Extract(inVar, 12, 1.0);
	curEStrain.Extract(inVar, 18, 1.0);
	curAlpha.Extract(inVar, 24, 1.0);
	curM = inVar(30);
	curFabric.Extract(inVar, 31, 1.0);
	// Hypo elastic Stress calculation
	Matrix aC(6,6);  // Actually not use in here but for function argument compatibility
	stress = HypoElastic(curStress, eStrain, curEStrain, aC);
	// State Dependent variables
	Vector n(6), d(6), b(6); double bref, psi, Theta, alphaBtheta;
	StateDepend(stress, alpha, m, strain, n, d, b, bref, psi, Theta, alphaBtheta);
	Vector devStress(6); devStress = GetDevPart(stress);
	double p = GetTrace(stress)/3.0;
	Vector r(6); r = devStress - p * alpha;
	double c  = Me / Mc;
	double cb = kBE / kBC;
	double cd = kDE / kDC;
	double A = Ao * (1 + Macauley(DoubleDot2_2(fabric, n)));
	double D = A * DoubleDot2_2(d, n);
	Vector R(6); R = n + D * mI1 / 3.0;
	double h = ho * fabs(DoubleDot2_2(b,n)) / (bref - fabs(DoubleDot2_2(b, n)));
	// Up to here same as GetStateDepend
	double K = mK, G = mG;
	// 4th level
	Vector dCosThetaOverdEE(6), dCosThetaOverdAlpha(6);
	dCosThetaOverdEE = getdCosThetaOverdEE(r, p, devStress, K, G, Theta);
	dCosThetaOverdAlpha = getdCosThetaOverdAlpha( r/p, Theta);
	// 3rd Level no dependence
	Vector dAlphaDOverdEE(6), dAlphaDOverdAlpha(6), dAlphaBOverdEE(6), dAlphaBOverdAlpha(6), dbrefOverdEE(6), dbrefOverdAlpha(6);
	dAlphaDOverdEE = getdAlphaDOverdEE(c, Theta, dCosThetaOverdEE, psi, cd, K, p);
	dAlphaDOverdAlpha = getdAlphaDOverdAlpha(c, Theta, dCosThetaOverdAlpha, psi, cd);
	double dAlphaDOverdM = -1.0;                                                                         
	dAlphaBOverdEE = getdAlphaBOverdEE(c, Theta, dCosThetaOverdEE, psi, cb, K, p);
	dAlphaBOverdAlpha = getdAlphaBOverdAlpha(c, Theta, dCosThetaOverdAlpha, psi, cb, K, p);
	double dAlphaBOverdM = -1.0;
	dbrefOverdEE = 2.0 * root23 * dAlphaBOverdEE;
	dbrefOverdAlpha = 2.0 * root23 * dAlphaBOverdAlpha;
	double dbrefOverdM = -2.0 * root23;

	// 2nd level
	Matrix dnOverdEE(6,6), dbOverdEE(6,6), dbOverdAlpha(6,6), dnOverdAlpha(6,6);
	Vector dnOverdM(6), dDOverdEE(6), dDOverdAlpha(6), dbOverdM(6);
	dnOverdEE = -(1/GetContraNorm(r))
				* (K * (Dyadic2_2(alpha, mI1) - DoubleDot2_2(n, alpha) * Dyadic2_2(n, mI1))
				- 2 * G * ( mIIcon - Dyadic2_2(mI1, mI1) * one3 - Dyadic2_2(n,n)));
	dnOverdAlpha = (-p / GetContraNorm(r)) * (mIIcon - Dyadic2_2(n,n));
	dnOverdM.Zero();
	dDOverdEE = Ao * MacauleyIndex(DoubleDot2_2(fabric, n)) * DoubleDot2_4(fabric, dnOverdEE) * DoubleDot2_2(d, n)
		      + A * (sqrt(two3) * dAlphaDOverdEE - DoubleDot2_4(alpha, dnOverdEE));
/*	dDOverdEE = Ao * MacauleyIndex(DoubleDot(fabric, n)) * DoubleDot(fabric, dnOverdEE) * DoubleDot(d, n)
		      + A * (sqrt(TwoOverThree) * dAlphaDOverdEE + DoubleDot(alpha, dnOverdEE));*/ //--Changed in June 11, 2003
	dDOverdAlpha = Ao * MacauleyIndex(DoubleDot2_2(fabric, n)) *  DoubleDot2_4(fabric, dnOverdAlpha)* DoubleDot2_2(d, n)
              + A * (sqrt(two3) * dAlphaDOverdAlpha - n - DoubleDot4_2(dnOverdAlpha, alpha)); 
	double dDOverdM = A * root23 * dAlphaDOverdM;
	dbOverdEE = root23 * (Dyadic2_2(n, dAlphaBOverdEE) + alphaBtheta * dnOverdEE);
	dbOverdAlpha = root23 * Dyadic2_2(n, dAlphaBOverdAlpha) + root23 * alphaBtheta * dnOverdAlpha - mIIcon;
	dbOverdM = root23 * dAlphaBOverdM * n ;
	/////////////////////////////////////
	Vector dDoubleDotBNOverdEE(6), dhOverdEE(6), dhOverdAlphaTmp1(6), dhOverdAlphaTmp2(6), dhOverdAlpha(6);
	double debugSS1 = DoubleDot2_2(b, n);
	dDoubleDotBNOverdEE = sign(DoubleDot2_2(b, n)) * (DoubleDot4_2(dbOverdEE, n) + DoubleDot2_4(b, dnOverdEE));
	dhOverdEE = ho * (bref * dDoubleDotBNOverdEE - fabs(DoubleDot2_2(b, n)) * dbrefOverdEE) / pow((bref - fabs(DoubleDot2_2(b,n))), 2);
	dhOverdAlphaTmp1 = sign(DoubleDot2_2(b,n)) * (DoubleDot4_2(dbOverdAlpha, n) + DoubleDot2_4(b, dnOverdAlpha)) * bref;
	dhOverdAlphaTmp2 = fabs(DoubleDot2_2(b,n)) * dbrefOverdAlpha;
	dhOverdAlpha = ho *  (dhOverdAlphaTmp1 - dhOverdAlphaTmp2) / pow((bref - fabs(DoubleDot2_2(b,n))), 2);
	double dhOverdMTmp1 = bref * sign(DoubleDot2_2(b,n)) * DoubleDot2_2(dbOverdM, n);
	double dhOverdMTmp2 = fabs(DoubleDot2_2(b,n)) * dbrefOverdM;
	double dhOverdM = ho *  (dhOverdMTmp1 - dhOverdMTmp2) / pow((bref - fabs(DoubleDot2_2(b,n))), 2);
	// For Fabric Tensor
	Vector dDOverdFabric(6), dMauDOverdEE(6), dMauDOverdAlpha(6), dMauDOverdFabric(6);
//**dDOverdFabric    =  Ao * MacauleyIndex(DoubleDot(fabric, n)) * DoubleDot(d, n) * DoubleDot(I4_2, n);
	dDOverdFabric    =  Ao * MacauleyIndex(DoubleDot2_2(fabric, n)) * DoubleDot2_2(d, n) * n;
	dMauDOverdEE     = -MacauleyIndex(-D) * dDOverdEE;
	dMauDOverdAlpha  = -MacauleyIndex(-D) * dDOverdAlpha;
	double dMauDOverdM    = -MacauleyIndex(-D) * dDOverdM;
	dMauDOverdFabric = -MacauleyIndex(-D) * dDOverdEE;

	// First Level
	Matrix dR1OverdEE(6,6), dR1OverdAlpha(6,6), dR1OverdFabric(6,6);
	Vector dR1OverdM(6), dR1OverdDGamma(6);
	dR1OverdEE     = mIIcon + dGamma * (dnOverdEE + Dyadic2_2(mI1, dDOverdEE) * one3);
	dR1OverdAlpha  = dGamma * (dnOverdAlpha + Dyadic2_2(mI1, dDOverdAlpha) * one3);
	dR1OverdM      = dGamma * (dnOverdM + dDOverdM * mI1 / 3.0);
	dR1OverdFabric = dGamma * Dyadic2_2(mI1, dDOverdFabric) * one3; 
	dR1OverdDGamma = R;

	Matrix dR2OverdEE(6,6), dR2OverdAlpha(6,6), dR2OverdFabric(6,6);
	Vector dR2OverdM(6), dR2OverdDGamma(6);
	dR2OverdEE     =  -dGamma * (Dyadic2_2(b, dhOverdEE) + h * dbOverdEE);
	dR2OverdAlpha  =  mIIcon - dGamma * (Dyadic2_2(b, dhOverdAlpha) + h * dbOverdAlpha);
	dR2OverdM      =  -dGamma * (dhOverdM * b + h * dbOverdM);
	dR2OverdFabric.Zero();
	dR2OverdDGamma =  -h * b;

	Vector dR3OverdEE(6), dR3OverdAlpha(6), dR3OverdFabric(6);
	double dR3OverdM, dR3OverdDGamma;
	dR3OverdEE     =  -dGamma * Cm * (1 + eo) * dDOverdEE;
	dR3OverdAlpha  =  -dGamma * Cm * (1 + eo) * dDOverdAlpha;
	dR3OverdM      =  1 - dGamma * Cm * (1 + eo) * dDOverdM;
	dR3OverdFabric =  -dGamma * Cm * (1 + eo) * dDOverdFabric;
	dR3OverdDGamma =  -Cm * (1 + eo ) * D;

	Matrix dR4OverdEE(6,6), dR4OverdAlpha(6,6), dR4OverdFabric(6,6);
	Vector dR4OverdM(6), dR4OverdDGamma(6);
	dR4OverdEE     = dGamma * Cf * (Dyadic2_2((Fmax * n + fabric), dMauDOverdEE) + Macauley(-D) * Fmax * dnOverdEE);
	dR4OverdAlpha  = dGamma * Cf * (Dyadic2_2((Fmax * n + fabric), dMauDOverdAlpha) + Macauley(-D) * Fmax * dnOverdAlpha);
	dR4OverdM      = dGamma * Cf * ((Fmax * n + fabric) * dMauDOverdM + Macauley(-D) * Fmax * dnOverdM);
	dR4OverdFabric = mIIcon + dGamma * Cf * (Dyadic2_2((Fmax * n + fabric), dMauDOverdFabric) + Macauley(-D) * mIIcon);
	dR4OverdDGamma = Cf * Macauley(-D) * (Fmax * n + fabric);

	Vector dR5OverdEE(6), dR5OverdAlpha(6), dR5OverdFabric(6);
	double dR5OverdM, dR5OverdDGamma;
//**dR5OverdEE     =  2 * G * DoubleDot2_4(n, mIIcon) - sqrt(two3) * m * K * mI1 - K * DoubleDot2_2(alpha, n) * mI1;
	dR5OverdEE     =  2 * G * n - sqrt(two3) * m * K * mI1 - K * DoubleDot2_2(alpha, n) * mI1;
//**dR5OverdAlpha  = -p * DoubleDot2_4(n, mIIcon);
	dR5OverdAlpha  = -p * n;
	dR5OverdM      = -root23 * p;
	dR5OverdFabric.Zero();  
	dR5OverdDGamma = 0.0;

	// Arrange Jacobian
	Matrix j(20, 20);
	j.Zero();
    j.Assemble(dR1OverdEE    , 0,  0, 1.0);
	j.Assemble(dR1OverdAlpha , 0,  6, 1.0);
	j.Assemble(dR1OverdM     , 0, 12, 1.0);
	j.Assemble(dR1OverdFabric, 0, 13, 1.0);
	j.Assemble(dR1OverdDGamma, 0, 19, 1.0);

	j.Assemble(dR2OverdEE    , 6,  0, 1.0);
	j.Assemble(dR2OverdAlpha , 6,  6, 1.0);
	j.Assemble(dR2OverdM     , 6, 12, 1.0);
	j.Assemble(dR2OverdFabric, 6, 13, 1.0);
	j.Assemble(dR2OverdDGamma, 6, 19, 1.0);

	j.AssembleTranspose(dR3OverdEE    , 12,  0, 1.0);
	j.AssembleTranspose(dR3OverdAlpha , 12,  6, 1.0);
	j(12,12) = dR3OverdM;
	j.AssembleTranspose(dR3OverdFabric, 12, 13, 1.0);
	j(12,19) = dR3OverdDGamma;

	j.Assemble(dR4OverdEE    , 13,  0, 1.0);
	j.Assemble(dR4OverdAlpha , 13,  6, 1.0);
	j.Assemble(dR4OverdM     , 13, 12, 1.0);
	j.Assemble(dR4OverdFabric, 13, 13, 1.0);
	j.Assemble(dR4OverdDGamma, 13, 19, 1.0);

	j.AssembleTranspose(dR5OverdEE    , 19,  0, 1.0);
	j.AssembleTranspose(dR5OverdAlpha , 19,  6, 1.0);
	j(19,12) = dR5OverdM;
	j.AssembleTranspose(dR5OverdFabric, 19, 13, 1.0);
	j(19,19) = dR5OverdDGamma;

	return j;
}

Vector 
ManzariDafalias::getdCosThetaOverdEE(const Vector& r, const double &p, const Vector& s, double K, double G, double theta)
{
	//Scalar theta = acos(aCosTheta);
	Vector result(6);
	Vector  rBar(6); rBar = r / p;
	double j2 = GetJ2(rBar);
	double j3 = GetJ3(rBar);
	Vector term1(6);
	Matrix term2(6,6);
	term1 = (pow(j2,1.5) * Dot(rBar) - 1.5 * sqrt(j2) * j3 * rBar) / (j2 * j2* j2);
	term2 = (2.0 * G * p * (mIIcon - Dyadic2_2(mI1, mI1)/3) - K * Dyadic2_2(s, mI1)) / (p * p);

	double aTheta; //-- Non-continuous Willam's equation correction at theta = 0 and pi/3
	if (fabs(theta) <= sqrt(mEPS))
		aTheta = sqrt(mEPS);
	else if ((fabs(theta)-PI/3.0) <= sqrt(mEPS))
		aTheta = theta - sqrt(mEPS);
	else
		aTheta = theta;

	result = (sin(aTheta) / (3.0*sin(3.0*aTheta))) * 1.5 * sqrt(3.0) * DoubleDot2_4(term1, term2);
	return result;
}

Vector
ManzariDafalias::getdCosThetaOverdAlpha(const Vector& rBar, double theta)
{
	//Scalar theta = acos(aCosTheta);
	Vector result(6);
	double j2 = GetJ2(rBar);
	double j3 = GetJ3(rBar);
	Vector term1(6); term1 = Dot(rBar);
	Vector term(6);
	term = -pow(j2, 1.5) * term1 + 1.5 * sqrt(j2) * j3 * rBar;

	double aTheta; //-- Non-continuous Willam's equation correction at theta = 0 and pi/3
	if (fabs(theta) <= sqrt(mEPS))
		aTheta = sqrt(mEPS);
	else if ((fabs(theta)-PI/3.0) <= sqrt(mEPS))
		aTheta = theta - sqrt(mEPS);
	else
		aTheta = theta;

	result = (sin(aTheta) / (3.0*sin(3.0*aTheta))) * 1.5 * sqrt(3.0) * term / (j2 * j2 * j2);
	return result;
}

Vector 
ManzariDafalias::getdAlphaDOverdEE(double c, double theta, const Vector& dCosThetaOverdEE
								  ,double psi, double cd, double K, double p)
{
	Vector result(6);
	double Mc = mParam(9), kDE = mParam(12), kDC = mParam(13), lamda = mParam(15);
	double dgCOverdCosTheta  = getdgOverdCosTheta(theta, c); 
	double dgCdOverdCosTheta = getdgOverdCosTheta(theta, cd);
	result = Mc * dgCOverdCosTheta * dCosThetaOverdEE +
			 kDC * psi * dgCdOverdCosTheta * dCosThetaOverdEE +
			 kDC * g(theta, cd) * lamda * K * mI1 / p;
	return result;
}

Vector 
ManzariDafalias::getdAlphaDOverdAlpha(double c, double theta, const Vector& dCosThetaOverdAlpha
									 ,double psi, double cd)
{
	Vector result(6);
	double Mc = mParam(9), kDC = mParam(13);
	double dgCOverdCosTheta  = getdgOverdCosTheta(theta, c); 
	double dgCdOverdCosTheta = getdgOverdCosTheta(theta, cd);
	result = Mc * dgCOverdCosTheta * dCosThetaOverdAlpha +
			 kDC * psi * dgCdOverdCosTheta * dCosThetaOverdAlpha;
	return result;
}

Vector
ManzariDafalias::getdAlphaBOverdEE(double c, double theta, const Vector& dCosThetaOverdEE
								  ,double psi, double cb, double K, double p)
{
	Vector result(6);
	double Mc = mParam(9), kBE = mParam(10), kBC = mParam(11), lamda = mParam(15);
	double dgCOverdCosTheta  = getdgOverdCosTheta(theta, c); 
	double dgCbOverdCosTheta = getdgOverdCosTheta(theta, cb);
	result = Mc * dgCOverdCosTheta * dCosThetaOverdEE
			 + kBC * Macauley(-psi) * dgCbOverdCosTheta * dCosThetaOverdEE
			 - MacauleyIndex(-psi) * kBC * g(theta, cb) * lamda * K * mI1 / p;
	return result;
}

Vector
ManzariDafalias::getdAlphaBOverdAlpha(double c, double theta, const Vector dCosThetaOverdAlpha
									 ,double psi, double cb, double K, double p)
{
	Vector result(6);
	double Mc = mParam(9), kBC = mParam(11); 
	double dgCOverdCosTheta  = getdgOverdCosTheta(theta, c); 
	double dgCbOverdCosTheta = getdgOverdCosTheta(theta, cb);
	result = Mc * dgCOverdCosTheta * dCosThetaOverdAlpha
			 + kBC * Macauley(-psi) * dgCbOverdCosTheta * dCosThetaOverdAlpha;
	return result;
}

/*   No 60 degree rotation
double
ManzariDafalias::getdgOverdCosTheta(const double theta, const double c)
{
	double g1, dg1, g2, dg2, result;
	double term1 = 4.0 * (1.0 - c*c) * cos(theta) * cos(theta) + 5.0 * c * c - 4.0 * c;
	g1  = 2.0 * (1.0 - c*c) * cos(theta) + (2.0*c - 1.0) * sqrt(term1);
	g2  = 4.0 * (1.0 - c*c) * cos(theta) * cos(theta) + (1.0 - 2.0*c) * (1.0 - 2.0*c);
	double term2 = 4.0 * (2.0*c - 1.0) * (1.0 - c*c) * cos(theta) * sin(theta);
	dg1 = -2.0 * (1.0 - c*c) * sin(theta) - term2 / sqrt(term1);
	dg2 = -8.0 * (1.0 - c*c) * cos(theta) * sin(theta);
	result = (dg1*g2 - g1*dg2) / (g2 * g2);
	return result;
}
*/

//--- Counter-clockwise 60 dgree rotation of Willam
double
ManzariDafalias::getdgOverdCosTheta(const double theta, const double c)
{
	double g1, dg1, g2, dg2, result;
	double term1 = 4.0 * (1.0 - c*c) * cos(theta - PI/3.0) * cos(theta - PI/3.0) + 5.0 * c * c - 4.0 * c;
	g1  = 2.0 * (1.0 - c*c) * cos(theta - PI/3.0) + (2.0*c - 1.0) * sqrt(term1);
	g2  = 4.0 * (1.0 - c*c) * cos(theta - PI/3.0) * cos(theta - PI/3.0) + (1.0 - 2.0*c) * (1.0 - 2.0*c);
	double term2 = 4.0 * (2.0*c - 1.0) * (1.0 - c*c) * cos(theta - PI/3.0) * sin(theta - PI/3.0);
	dg1 = -2.0 * (1.0 - c*c) * sin(theta - PI/3.0) - term2 / sqrt(term1);
	dg2 = -8.0 * (1.0 - c*c) * cos(theta - PI/3.0) * sin(theta - PI/3.0);
	result = (dg1*g2 - g1*dg2) / (g2 * g2);
	return result;
}

int ManzariDafalias::sign(double x)
{
	int result = 1;
	if (x < 0)
		result = -1;
	return result;
}

// ---------------------------------------------------------------------------------------------------------

Vector
ManzariDafalias::GetState()
// returns vector of state parameters for recorders
{
	return mState;
}

Response*
ManzariDafalias::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	if (strcmp(argv[0],"stress") == 0 || strcmp(argv[0],"stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());
	else if (strcmp(argv[0],"strain") == 0 || strcmp(argv[0],"strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());
	else if (strcmp(argv[0], "state") == 0)
		return new MaterialResponse(this, 3, this->GetState());
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
				*(matInfo.theVector) = GetState();
			return 0;
		default:
			return -1;
	}
}

int
ManzariDafalias::sendSelf(int commitTag, Channel &theChannel)
{
  // we place all the data needed to define material and it's state into a vector object
  static Vector data(8);
  int cnt = 0;

  // send the vector object to the channel
  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ManzariDafalias::sendSelf - failed to send vector to channel\n";
    return -1;
  }
  

  return 0;
 
}

int 
ManzariDafalias::recvSelf(int commitTag, Channel &theChannel, 
                                         FEM_ObjectBroker &theBroker)    
{
  // recv the vector object from the channel which defines material param and state
  static Vector data(7);
  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
    opserr << "ManzariDafalias::recvSelf - failed to recv vector from channel\n";
    return -1;
  }

  // set the material parameters and state variables
  int cnt = 0;
  this->setTag(data(cnt++));

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
	
	return 0;
}