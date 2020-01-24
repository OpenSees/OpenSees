//
// Created by Alex Hartloper on 10.07.18.
//

#include "UVCplanestress.h"

#include <cmath>
#include <iostream>

#include <elementAPI.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

static int numUVCplanestress = 0;

// NOTE: Do not use the OPS_GetNumRemainingInputArgs() function or the
// OPS_GetString() function: causes crash with .dll
void* OPS_UVCplanestress(void) {
  if (numUVCplanestress == 0) {
    opserr << "Using the UVCplanestress material, see "
      "https://www.epfl.ch/labs/resslab/resslab-tools/" << endln;
    numUVCplanestress++;
  }
  NDMaterial* theMaterial = 0;

  // Parameters for parsing
  const unsigned int N_TAGS = 1;
  const unsigned int N_BASIC_PROPERTIES = 5;
  const unsigned int N_UPDATED_PROPERTIES = 2;
  const unsigned int N_PARAM_PER_BACK = 2;
  const unsigned int MAX_BACKSTRESSES = 8;
  const unsigned int BACKSTRESS_SPACE = MAX_BACKSTRESSES * N_PARAM_PER_BACK;

  std::string inputInstructions = "Invalid args, want:\n"
    "nDMaterial UVCplanestress "
    "tag? E? nu? fy? QInf? b? DInf? a? "
    "N? C1? gamma1? <C2? gamma2? C3? gamma3? ... C8? gamma8?>\n"
    "Note: to neglect the updated model, set DInf = 0.0";

  // Containers for the inputs
  int nInputsToRead;
  int nBackstresses[1];  // for N
  int materialTag[N_TAGS];  // for the tag
  double basicProps[N_BASIC_PROPERTIES];  // holds E, nu, fy, QInf, b
  double updProps[N_UPDATED_PROPERTIES];  // holds DInf, a
  double backstressProps[BACKSTRESS_SPACE];  // holds C's and gamma's
  std::vector<double> cK;
  std::vector<double> gammaK;

  // Get the material tag
  nInputsToRead = N_TAGS;
  if (OPS_GetIntInput(&nInputsToRead, materialTag) != 0) {
    opserr << "WARNING invalid nDMaterial UVCplanestress tag" << endln;
    return 0;
  }

  // Get E, nu, fy, qInf, b
  nInputsToRead = N_BASIC_PROPERTIES;
  if (OPS_GetDoubleInput(&nInputsToRead, basicProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }

  // Read in the updated model paramters
  nInputsToRead = N_UPDATED_PROPERTIES;
  if (OPS_GetDoubleInput(&nInputsToRead, updProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }

  // Get the number of backstresses
  nInputsToRead = 1;
  if (OPS_GetIntInput(&nInputsToRead, nBackstresses) != 0) {
    opserr << "WARNING N must be an integer" <<
      inputInstructions.c_str() << endln;
    return 0;
  }
  if (nBackstresses[0] > MAX_BACKSTRESSES) {
    opserr << "WARNING: Too many backstresses defined, maximum is: " <<
      MAX_BACKSTRESSES << endln <<
      inputInstructions.c_str() << endln;
    return 0;
  }

  // Get the backstress parameters
  nInputsToRead = 2 * nBackstresses[0];
  if (OPS_GetDoubleInput(&nInputsToRead, backstressProps) != 0) {
    opserr << inputInstructions.c_str() << endln;
    return 0;
  }
  // cK's alternate with gammaK's
  for (unsigned int i = 0; i < nBackstresses[0]; ++i) {
    cK.push_back(backstressProps[2 * i]);
    gammaK.push_back(backstressProps[1 + 2 * i]);
  }

  // Allocate the material
  theMaterial = new UVCplanestress(materialTag[0],
    basicProps[0], basicProps[1], basicProps[2],
    basicProps[3], basicProps[4],
    updProps[0], updProps[1],
    cK, gammaK);

  return theMaterial;
}


/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param tag Material tag
* @param E elastic modulus
* @param poissonRatio Poisson's Ratio
* @param sy0 initial yield stress
* @param qInf increase in yield surface size at isotropic saturation (>= 0)
* @param b controls the saturation rate of the isotropic hardening, must be
nonzero
* @param dInf difference in initial and steady state yield stress (>=0), set
*              this to zero for the original Chaboche model
* @param a controls the saturation of dInf, must be nonzero
* @param cK backstress kinematic hardening moduli
* @param gammaK controls the saturation rate of the kinematic hardening
*/
UVCplanestress::UVCplanestress(int tag, double E, double poissonRatio,
  double sy0, double qInf, double b,
  double dInf, double a,
  std::vector<double> cK, std::vector<double> gammaK)
  : NDMaterial(tag, ND_TAG_UVCplanestress),
  elasticModulus(E),
  poissonRatio(poissonRatio),
  initialYield(sy0),
  qInf(qInf),
  bIso(b),
  dInf(dInf),
  aIso(a),
  cK(cK),
  gammaK(gammaK),
  shearModulus(E / (2. * (1. + poissonRatio))),
  bulkModulus(E / (3. * (1. - 2. * poissonRatio))),
  strainConverged(N_DIMS),
  strainTrial(N_DIMS),
  strainPlasticConverged(N_DIMS),
  strainPlasticTrial(N_DIMS),
  strainPEqConverged(0.),
  strainPEqTrial(0.),
  stressConverged(N_DIMS),
  stressTrial(N_DIMS),
  plasticLoading(false),
  elasticMatrix(Matrix(N_DIMS, N_DIMS)),
  stiffnessInitial(Matrix(N_DIMS, N_DIMS)),
  stiffnessConverged(Matrix(N_DIMS, N_DIMS)),
  stiffnessTrial(Matrix(N_DIMS, N_DIMS)),
  pMat(Matrix(N_DIMS, N_DIMS)),
  qMat(Matrix(N_DIMS, N_DIMS)),
  qMatT(Matrix(N_DIMS, N_DIMS)),
  lambdaC(N_DIMS),
  lambdaP(N_DIMS)
{
  // Set the number of backstresses
  nBackstresses = cK.size();
  for (unsigned int i = 0; i < nBackstresses; ++i) {
    alphaKTrial.push_back(Vector(N_DIMS));
    alphaKConverged.push_back(Vector(N_DIMS));
  }

  // Zero all the vectors and matrices
  revertToStart();

  // Set the Eigendecomposition matrices
  initializeEigendecompositions();

  // Set elastic parameters and elastic stiffness matrix
  calculateElasticStiffness();
  stiffnessInitial = elasticMatrix;
  stiffnessTrial = elasticMatrix;
  stiffnessConverged = elasticMatrix;
};

/* ----------------------------------------------------------------------------------------------------------------- */

UVCplanestress::UVCplanestress()
  : NDMaterial(0, ND_TAG_UVCplanestress),
  elasticModulus(0.),
  poissonRatio(0.),
  initialYield(0.),
  qInf(0.),
  bIso(0.),
  dInf(0.),
  aIso(0.),
  cK(0.),
  gammaK(0.),
  shearModulus(0. / (2. * (1. + poissonRatio))),
  bulkModulus(0. / (3. * (1. - 2. * poissonRatio))),
  strainConverged(N_DIMS),
  strainTrial(N_DIMS),
  strainPlasticConverged(N_DIMS),
  strainPlasticTrial(N_DIMS),
  strainPEqConverged(0.),
  strainPEqTrial(0.),
  stressConverged(N_DIMS),
  stressTrial(N_DIMS),
  plasticLoading(false),
  elasticMatrix(Matrix(N_DIMS, N_DIMS)),
  stiffnessInitial(Matrix(N_DIMS, N_DIMS)),
  stiffnessConverged(Matrix(N_DIMS, N_DIMS)),
  stiffnessTrial(Matrix(N_DIMS, N_DIMS)),
  pMat(Matrix(N_DIMS, N_DIMS)),
  qMat(Matrix(N_DIMS, N_DIMS)),
  qMatT(Matrix(N_DIMS, N_DIMS)),
  lambdaC(N_DIMS),
  lambdaP(N_DIMS)
{
  // Set the number of backstresses
  // todo: this probably wont work with parallel processing?
  nBackstresses = cK.size();
  for (unsigned int i = 0; i < nBackstresses; ++i) {
    alphaKTrial.push_back(Vector(N_DIMS));
    alphaKConverged.push_back(Vector(N_DIMS));
  }

  // Zero all the vectors and matrices
  revertToStart();

  // Set the Eigendecomposition matrices
  initializeEigendecompositions();

  // Set elastic stiffness matrix
  calculateElasticStiffness();
  stiffnessInitial = elasticMatrix;
  stiffnessTrial = elasticMatrix;
  stiffnessConverged = elasticMatrix;
}

/* ----------------------------------------------------------------------------------------------------------------- */

UVCplanestress::~UVCplanestress() {

}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @return 0 if successful
*/
int UVCplanestress::returnMapping() {
  // Initialize all the variables
  int retVal = 0;
  bool converged = true;
  unsigned int iterationNumber = 0;
  Vector alpha = Vector(N_DIMS);
  Vector stressRelative = Vector(N_DIMS);
  Vector etaTrial = Vector(N_DIMS);
  Vector etaTilde = Vector(N_DIMS);
  Vector eta = Vector(N_DIMS);
  Vector alphaTilde = Vector(N_DIMS);
  Vector alphaTildePrime = Vector(N_DIMS);
  Vector gammaDiag = Vector(N_DIMS);
  Vector gammaDiagPrime = Vector(N_DIMS);
  double yieldFunction = 0., yieldStress = 0., consistParam = 0., isotropicModulus = 0., eK = 0.,
    consistDenom = 0., fBarSquared = 0., fBar = 0., beta = 0., gammaDenom = 0., betaPrime = 0.;

  // Elastic trial step
  alpha.Zero();
  for (unsigned int i = 0; i < nBackstresses; ++i)
    alpha = alpha + alphaKConverged[i];
  stressTrial = elasticMatrix * (strainTrial - strainPlasticConverged);
  etaTrial = qMatT * (stressTrial - alpha);
  eta = etaTrial;

  // Yield condition
  yieldStress = calculateYieldStress();
  fBarSquared = 1. / 3. * pow(eta(0), 2) + pow(eta(1), 2) + 2. * pow(eta(2), 2);
  yieldFunction = 1. / 2. * fBarSquared - 1. / 3. * pow(yieldStress, 2);
  if (yieldFunction > RETURN_MAP_TOL) {
    converged = false;
  }

  // Do the return mapping if plastic loading
  while (!converged && iterationNumber < MAXIMUM_ITERATIONS) {
    iterationNumber++;

    // Isotropic hardening parameters
    yieldStress = calculateYieldStress();
    isotropicModulus = calculateIsotropicModulus();
    // Kinematic hardening parameters
    beta = 0.;
    alphaTilde.Zero();
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      beta += cK[i] / gammaK[i] * (1. - eK);
      alphaTilde += alphaKConverged[i] * eK;
    }
    alphaTilde = alpha - alphaTilde;
    beta = 1. + beta / yieldStress;

    // Update the relative stress and eta
    gammaDenom = beta + 2. * shearModulus * consistParam;
    gammaDiag(0) = 1. / (beta + consistParam * elasticModulus / (3. * (1. - poissonRatio)));
    gammaDiag(1) = 1. / gammaDenom;
    gammaDiag(2) = 1. / gammaDenom;


    etaTilde = etaTrial + qMatT * alphaTilde;
    eta = vecMult3(etaTilde, gammaDiag);
    fBarSquared = 1. / 3. * pow(eta(0), 2) + pow(eta(1), 2) + 2. * pow(eta(2), 2);
    fBar = sqrt(fBarSquared);

    // Calculate Newton denominator
    betaPrime = 0.;
    alphaTildePrime.Zero();
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      betaPrime = betaPrime - cK[i] * isotropicModulus / (gammaK[i] * pow(yieldStress, 2)) * (1. - eK)
        + cK[i] * eK / yieldStress;
      alphaTildePrime = alphaTildePrime + gammaK[i] * eK * alphaKConverged[i];
    }
    betaPrime = betaPrime * sqrt(2. / 3.) * fBar;
    alphaTildePrime = alphaTildePrime * sqrt(2. / 3.) * fBar;
    for (unsigned int i = 0; i < N_DIMS; ++i)
      gammaDiagPrime(i) = -pow(gammaDiag(i), 2) * (betaPrime + lambdaP(i) * lambdaC(i));

    consistDenom = dotprod3(vecMult3(lambdaP, eta),
      vecMult3(gammaDiagPrime, etaTilde) + vecMult3(gammaDiag, qMatT * alphaTildePrime))
      - sqrt(2. / 3.) * 2. / 3. * yieldStress * isotropicModulus * fBar;

    // Newton step
    yieldFunction = 1. / 2. * fBarSquared - 1. / 3. * pow(yieldStress, 2);
    consistParam = consistParam - yieldFunction / (consistDenom + RETURN_MAP_TOL);
    strainPEqTrial = strainPEqConverged + sqrt(2. / 3.) * consistParam * fBar;

    // Check convergence
    if (fabs(yieldFunction) < RETURN_MAP_TOL) {
      converged = true;
    }
  }

  // Condition for plastic loading is whether or not iterations were performed
  if (iterationNumber == 0) {
    plasticLoading = false;
  }
  else {
    // Update the variables
    plasticLoading = true;
    etaTilde = etaTrial + qMatT * alphaTilde;
    eta = vecMult3(gammaDiag, etaTilde);
    stressRelative = qMat * eta;
    yieldStress = calculateYieldStress();
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      alphaKTrial[i] = alphaKConverged[i] * eK + stressRelative / yieldStress * cK[i] / gammaK[i] * (1. - eK);
    }
    strainPlasticTrial = strainPlasticConverged + consistParam * pMat * stressRelative;
    stressTrial = elasticMatrix * (strainTrial - strainPlasticTrial);
  }

  // Update the stiffness
  calculateStiffness(consistParam, fBar, stressRelative);

  // Warn the user if the algorithm did not converge and return -1
  if (iterationNumber >= MAXIMUM_ITERATIONS && fabs(yieldFunction) > RETURN_MAP_TOL) {
    opserr << "UVCplanestress::returnMapping return mapping in UVCplanestress did not converge!" << endln;
    opserr << "\tDelta epsilon 11 = " << strainTrial[0] - strainConverged[0] << endln;
    opserr << "\tDelta epsilon 22 = " << strainTrial[1] - strainConverged[1] << endln;
    opserr << "\tDelta epsilon 12 = " << strainTrial[3] - strainConverged[3] << endln;
    opserr << "\tExiting with yield function = " << yieldFunction << " > " << RETURN_MAP_TOL << endln;
    retVal = -1;
  }

  return retVal;
}

/* ----------------------------------------------------------------------------------------------------------------- */

void UVCplanestress::calculateStiffness(double consistParam, double fBar, const Vector& stressRelative) {
  if (!plasticLoading) {
    stiffnessTrial = elasticMatrix;
  }
  else  // plastic loading
  {
    double yieldStress = 0., isotropicModulus = 0., eK = 0., beta = 0., theta_2 = 0., theta_1 = 0.;
    Vector hPrime = Vector(N_DIMS), nHat = Vector(N_DIMS), nTilde = Vector(N_DIMS), hTilde = Vector(N_DIMS);
    Matrix complianceMatrix = Matrix(N_DIMS, N_DIMS), hOutN = Matrix(N_DIMS, N_DIMS), iD3 = Matrix(N_DIMS, N_DIMS),
      aMat = Matrix(N_DIMS, N_DIMS), xiTilde = Matrix(N_DIMS, N_DIMS), xiTildeA = Matrix(N_DIMS, N_DIMS),
      nOutN = Matrix(N_DIMS, N_DIMS);

    iD3.Zero();
    iD3(0, 0) = iD3(1, 1) = iD3(2, 2) = 1.;
    complianceMatrix = calculateComplianceMatrix();

    // Isotropic hardening parameters
    yieldStress = calculateYieldStress();
    isotropicModulus = calculateIsotropicModulus();

    // Kinematic hardening related parameters
    nHat = stressRelative / fBar;
    for (unsigned int i = 0; i < nBackstresses; ++i)
      beta += cK[i] / gammaK[i] * (1. - eK);
    beta = 1. + beta / yieldStress;
    hPrime = -(beta - 1.) * isotropicModulus * stressRelative / yieldStress;
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      hPrime += cK[i] * eK / yieldStress * stressRelative - gammaK[i] * eK * alphaKConverged[i];
    }
    hPrime *= sqrt(2. / 3.);
    hOutN = hPrime % nHat;
    aMat = matinv3(beta * iD3 + consistParam * hOutN * pMat);

    nTilde = nHat - consistParam * aMat * hPrime;
    xiTilde = matinv3(complianceMatrix + consistParam * pMat * aMat);
    xiTildeA = aMat * xiTilde;

    theta_2 = 1. - 2. / 3. * isotropicModulus * consistParam;
    hTilde = hPrime + xiTilde * (pMat * nTilde);
    theta_1 = 2. / 3. * isotropicModulus + theta_2 * dotprod3(nHat, pMat * (aMat * hTilde));
    nOutN = nTilde % nHat;
    stiffnessTrial.Zero();
    stiffnessTrial = xiTilde - theta_2 / theta_1 * xiTilde * pMat * nOutN * pMat * xiTildeA;

    // Take the symmetric approximation
    stiffnessTrial.addMatrixTranspose(0.5, stiffnessTrial, 0.5);
  }
  return;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v new total strain vector.
* @return 0 if successful, -1 if return mapping did not converge.
*/
int UVCplanestress::setTrialStrain(const Vector& v) {

  int rm_convergence;
  // Reset the trial state
  revertToLastCommit();

  // Set the trial strain
  strainTrial = v;

  // Do the return mapping and calculate the tangent modulus
  rm_convergence = returnMapping();

  return rm_convergence;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v new total strain vector
* @param v new strain rate vector - unused
* @return 0 if successful
*
* Note that this material model is rate independent.
*/
int UVCplanestress::setTrialStrain(const Vector& v, const Vector& r) {

  // Reset the trial state
  revertToLastCommit();

  // Set the trial strain
  strainTrial = v;

  // Do the return mapping and calculate the tangent modulus
  returnMapping();

  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v strain increment vector
* @return 0 if successful
*/
int UVCplanestress::setTrialStrainIncr(const Vector& v) {

  // Reset the trial state
  revertToLastCommit();

  // Set the trial strain
  strainTrial += v;

  // Do the return mapping and calculate the tangent modulus
  returnMapping();

  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v strain increment vector
* @param v strain rate vector - unused
* @return 0 if successful
*
* Note that this material model is rate independent.
*/
int UVCplanestress::setTrialStrainIncr(const Vector& v, const Vector& r) {

  // Reset the trial state
  revertToLastCommit();

  // Set the trial strain
  strainTrial += v;

  // Do the return mapping and calculate the tangent modulus
  returnMapping();

  return 0;
}


/* ----------------------------------------------------------------------------------------------------------------- */

const Vector& UVCplanestress::getStrain() {
  return strainTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Vector& UVCplanestress::getStress() {
  return stressTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Matrix& UVCplanestress::getTangent() {
  return stiffnessTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Matrix& UVCplanestress::getInitialTangent() {
  // todo: can make more efficient by changing this to elasticMatrix and removing stiffnessInitial as a variable
  return stiffnessInitial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @return 0 if successful
*/
int UVCplanestress::commitState() {
  strainConverged = strainTrial;
  strainPlasticConverged = strainPlasticTrial;
  strainPEqConverged = strainPEqTrial;
  stressConverged = stressTrial;
  alphaKConverged = alphaKTrial;
  stiffnessConverged = stiffnessTrial;
  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @return 0 if successful
*/
int UVCplanestress::revertToLastCommit() {
  strainTrial = strainConverged;
  strainPlasticTrial = strainPlasticConverged;
  strainPEqTrial = strainPEqConverged;
  stressTrial = stressConverged;
  alphaKTrial = alphaKConverged;
  stiffnessTrial = stiffnessConverged;
  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @return 0 if successful
*/
int UVCplanestress::revertToStart() {
  strainConverged.Zero();
  strainPlasticConverged.Zero();
  strainPEqConverged = 0.;
  stressConverged.Zero();
  plasticLoading = false;
  stiffnessConverged.Zero();
  for (unsigned int i = 0; i < nBackstresses; ++i) {
    alphaKConverged[i].Zero();
  }
  revertToLastCommit();
  return 0;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
* Returns a new NDMaterial with all the internal values copied.
* @return a to pointer to the copy
*
* This is called by GenericSectionXD
*/
NDMaterial* UVCplanestress::getCopy() {

  UVCplanestress* theCopy;
  theCopy = new UVCplanestress(this->getTag(), elasticModulus, poissonRatio,
    initialYield, qInf, bIso,
    dInf, aIso,
    cK, gammaK);

  // Copy all the internals
  theCopy->strainConverged = strainConverged;
  theCopy->strainTrial = strainTrial;
  theCopy->strainPlasticConverged = strainPlasticConverged;
  theCopy->strainPlasticTrial = strainPlasticTrial;
  theCopy->strainPEqConverged = strainPEqConverged;
  theCopy->strainPEqTrial = strainPEqTrial;
  theCopy->stressConverged = stressConverged;
  theCopy->stressTrial = stressTrial;
  theCopy->alphaKConverged = alphaKConverged;
  theCopy->alphaKTrial = alphaKTrial;
  theCopy->stiffnessConverged = stiffnessConverged;
  theCopy->stiffnessTrial = stiffnessTrial;
  theCopy->plasticLoading = plasticLoading;

  return theCopy;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
* Returns a new NDMaterial if the code matches the type specification.
* @param code the type specification of the material copy requested
* @return a to pointer to the copy
*
* This is called by the continuum elements.
*/
NDMaterial* UVCplanestress::getCopy(const char* code) {
  if (strcmp(code, getType()) == 0) {
    UVCplanestress* theCopy;
    theCopy = new UVCplanestress(this->getTag(), elasticModulus, poissonRatio,
      initialYield, qInf, bIso,
      dInf, aIso,
      cK, gammaK);
    return theCopy;
  }
  else {
    // todo: output to opserr
    opserr << "UVCplanestress::getCopy invalid NDMaterial type, expecting " << code << endln;
    return 0;
  }
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
* todo: fill out
* @param commitTag
* @param theChannel
* @return 0 if successful
*/
int UVCplanestress::sendSelf(int commitTag, Channel& theChannel) {

  /*
  static Vector data(26);  // enough space for 4 backstresses
  // Material properties
  data(0) = elasticModulus;
  data(1) = initialYield;
  data(2) = qInf;
  data(3) = bIso;
  data(4) = dInf;
  data(5) = aIso;
  data(6) = stiffnessInitial;

  // Internal variables
  data(7) = strainConverged;
  data(8) = strainPEqConverged;
  data(9) = stressConverged;
  data(10) = stiffnessConverged;
  data(11) = flowDirection;
  data(12) = plasticLoading;

  // Kinematic hardening related, 12 total spaces required
  unsigned int cKStart = 13;  // starts at the 13th space
  unsigned int gammaKStart = cKStart + nBackstresses;
  unsigned int alpha_k_start = gammaKStart + nBackstresses;
  for (unsigned int i = 0; i < nBackstresses; ++i) {
  data(cKStart + i) = cK[i];
  data(gammaKStart + i) = gammaK[i];
  data(alpha_k_start + i) = alphaKConverged[i];
  }

  data(25) = this->getTag();

  if (theChannel.sendVector(this->getDbTag(), commitTag, data) < 0) {
  opserr << "UVCplanestress::sendSelf() - failed to sendSelf\n";
  return -1;
  }
  */

  opserr << "Fatal: Paralleliziation for UVCplanestress is not implemented yet!" << endln;
  return -1;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
* todo: fill out
* @param commitTag
* @param theChannel
* @param theBroker
* @return 0 if successful
*/
int UVCplanestress::recvSelf(int commitTag, Channel& theChannel,
  FEM_ObjectBroker& theBroker) {
  /*
  static Vector data(26);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
  opserr << "UVCplanestress::recvSelf() - failed to recvSelf\n";
  return -1;
  }

  // Material properties
  elasticModulus = data(0);
  yieldStress = data(1);
  qInf = data(2);
  bIso = data(3);
  dInf = data(4);
  aIso = data(5);
  stiffnessInitial = data(6);

  // Internal variables
  strainConverged = data(7);
  strainPEqConverged = data(8);
  stressConverged = data(9);
  stiffnessConverged = data(10);
  flowDirection = data(11);
  plasticLoading = bool(data(12));

  // Kinematic hardening related, 12 total spaces required
  unsigned int cKStart = 13;  // starts at the 13th space
  unsigned int gammaKStart = cKStart + nBackstresses;
  unsigned int alpha_k_start = gammaKStart + nBackstresses;
  for (unsigned int i = 0; i < nBackstresses; ++i) {
  cK[i] = (cKStart + i);
  gammaK[i] = (gammaKStart + i);
  alphaKConverged[i] = (alpha_k_start + i);
  }

  this->setTag(int(data(25)));

  // Set the trial to the converged values
  revertToLastCommit();
  */

  opserr << "Fatal: Paralleliziation for UVCplanestress is not implemented yet!" << endln;
  return -1;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param s the opensees output stream
* @param flag is 2 for standard output, 25000 for JSON output
(see OPS_Globals.h)
*/
void UVCplanestress::Print(OPS_Stream& s, int flag) {

  // todo: change these back when not only .dll
  // if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
  if (flag == 2) {
    s << "UVCplanestress tag: " << this->getTag() << endln;
    s << "   E: " << elasticModulus << " ";
    s << "  fy: " << initialYield << " ";
    s << "   Q: " << qInf << " ";
    s << "   b: " << bIso << " ";
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      s << "  C" << (i + 1) << ": " << cK[i] << " ";
      s << "gam" << (i + 1) << ": " << gammaK[i] << " ";
    }
  }

  // if (flag == OPS_PRINT_PRINTMODEL_JSON) {
  if (flag == 25000) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"UVCplanestress\", ";
    s << "\"E\": " << elasticModulus << ", ";
    s << "\"fy\": " << initialYield << ", ";
    s << "\"Q\": " << qInf << ", ";
    s << "\"b\": " << bIso << ", ";
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      s << "\"C\": " << cK[i] << ", ";
      s << "\"gam\": " << gammaK[i] << ", ";
    }
  }

}

/* ----------------------------------------------------------------------------------------------------------------- */

void UVCplanestress::calculateElasticStiffness() {
  double eDenom = elasticModulus / (1. - pow(poissonRatio, 2));
  elasticMatrix.Zero();
  elasticMatrix(0, 0) = elasticMatrix(1, 1) = 1.0 * eDenom;
  elasticMatrix(1, 0) = elasticMatrix(0, 1) = poissonRatio * eDenom;
  elasticMatrix(2, 2) = (1. - poissonRatio) / 2. * eDenom;
}

/* ----------------------------------------------------------------------------------------------------------------- */

Matrix UVCplanestress::calculateComplianceMatrix() {
  double eDenom = elasticModulus;
  Matrix complianceMatrix = Matrix(N_DIMS, N_DIMS);
  complianceMatrix.Zero();
  complianceMatrix(0, 0) = complianceMatrix(1, 1) = 1.0 / eDenom;
  complianceMatrix(1, 0) = complianceMatrix(0, 1) = -poissonRatio / eDenom;
  complianceMatrix(2, 2) = 2. * (1. + poissonRatio) / eDenom;
  return complianceMatrix;
}

/* ----------------------------------------------------------------------------------------------------------------- */

void UVCplanestress::initializeEigendecompositions() {

  // Orthogonal matrix (eigenvectors)
  double qDenom = sqrt(2.);
  qMat.Zero();
  qMat(0, 0) = 1. / qDenom;  qMat(0, 1) = -1. / qDenom; qMat(0, 2) = 0;
  qMat(1, 0) = 1. / qDenom;  qMat(1, 1) = 1. / qDenom; qMat(1, 2) = 0;
  qMat(2, 0) = 0;  qMat(2, 1) = 0; qMat(2, 2) = 1.;
  // Transpose
  qMatT.Zero();
  qMatT.addMatrixTranspose(0., qMat, 1.0);

  // Projection matrix
  pMat.Zero();
  pMat(0, 0) = pMat(1, 1) = 2. / 3.;
  pMat(1, 0) = pMat(0, 1) = -1. / 3.;
  pMat(2, 2) = 2.;
  lambdaP.Zero();
  lambdaP(0) = 1. / 3.;
  lambdaP(1) = 1.;
  lambdaP(2) = 2.;

  // Elastic matrix, diagonal
  lambdaC.Zero();
  lambdaC(0) = elasticModulus / (1. - poissonRatio);
  lambdaC(1) = 2. * shearModulus;
  lambdaC(2) = shearModulus;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v1 length 3 vector
* @param v2 length 3 vector
* @return the dot product of the two vectors
*
*/
double UVCplanestress::dotprod3(const Vector& v1, const Vector& v2) {
  double res = 0.;
  for (unsigned int i = 0; i < N_DIMS; ++i)
    res += v1(i) * v2(i);
  return res;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v1 length 3 vector
* @param v2 length 3 vector
* @return vector containing the component-wise multiplication of the two vectors
*
*/
Vector UVCplanestress::vecMult3(const Vector& v1, const Vector& v2) {
  Vector res = Vector(N_DIMS);
  for (unsigned int i = 0; i < N_DIMS; ++i)
    res(i) = v1(i) * v2(i);
  return res;
}

/* ----------------------------------------------------------------------------------------------------------------- */

double UVCplanestress::calculateYieldStress() {
  double sigmaY1, sigmaY2;
  sigmaY1 = qInf * (1. - exp(-bIso * strainPEqTrial));
  sigmaY2 = dInf * (1. - exp(-aIso * strainPEqTrial));
  return initialYield + sigmaY1 - sigmaY2;
}


/* ----------------------------------------------------------------------------------------------------------------- */

double UVCplanestress::calculateIsotropicModulus() {
  double sigmaY1, sigmaY2;
  sigmaY1 = qInf * (1. - exp(-bIso * strainPEqTrial));
  sigmaY2 = dInf * (1. - exp(-aIso * strainPEqTrial));
  return bIso * (qInf - sigmaY1) - aIso * (dInf - sigmaY2);
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param i the i'th backstress
* @return computed eK factor
*
*/
double UVCplanestress::calculateEk(unsigned int i) {
  return exp(-gammaK[i] * (strainPEqTrial - strainPEqConverged));
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param A a 3x3 invertible matrix
* @return the inverse of A
*
*/
Matrix UVCplanestress::matinv3(const Matrix& A) {
  double detInv;
  Matrix B = Matrix(3, 3);

  // Calculate the determinant
  detInv = 1. /
    (A(0, 0) * A(1, 1) * A(2, 2) - A(0, 0) * A(1, 2) * A(2, 1)
      - A(0, 1) * A(1, 0) * A(2, 2) + A(0, 1) * A(1, 2) * A(2, 0)
      + A(0, 2) * A(1, 0) * A(2, 1) - A(0, 2) * A(1, 1) * A(2, 0));

  // Calculate the inverse
  B(0, 0) = +detInv * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1));
  B(1, 0) = -detInv * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0));
  B(2, 0) = +detInv * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
  B(0, 1) = -detInv * (A(0, 1) * A(2, 2) - A(0, 2) * A(2, 1));
  B(1, 1) = +detInv * (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0));
  B(2, 1) = -detInv * (A(0, 0) * A(2, 1) - A(0, 1) * A(2, 0));
  B(0, 2) = +detInv * (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1));
  B(1, 2) = -detInv * (A(0, 0) * A(1, 2) - A(0, 2) * A(1, 0));
  B(2, 2) = +detInv * (A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0));

  return B;
}

/* ----------------------------------------------------------------------------------------------------------------- */
