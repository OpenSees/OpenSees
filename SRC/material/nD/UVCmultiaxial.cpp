//
// Created by Alex Hartloper on 09.07.18.
//

#include "UVCmultiaxial.h"

#include <cmath>
#include <iostream>

#include <elementAPI.h>
#include <Channel.h>
#include <Information.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <OPS_Globals.h>

#include "classTags.h"
static int numUVCmultiaxial = 0;

// NOTE: Do not use the OPS_GetNumRemainingInputArgs() function or the
// OPS_GetString() function: causes crash with .dll
void* OPS_UVCmultiaxial(void) {
  if (numUVCmultiaxial == 0) {
    //std::cout << "Using the UVCmultiaxial material, see "
    //  "https://www.epfl.ch/labs/resslab/resslab-tools/" << std::endl;
    numUVCmultiaxial++;
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
    "nDMaterial UVCmultiaxial "
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
    opserr << "WARNING invalid nDMaterial UVCmultiaxial tag" << endln;
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
  theMaterial = new UVCmultiaxial(materialTag[0],
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
UVCmultiaxial::UVCmultiaxial(int tag, double E, double poissonRatio,
  double sy0, double qInf, double b,
  double dInf, double a,
  std::vector<double> cK, std::vector<double> gammaK)
  : NDMaterial(tag, ND_TAG_UVCmultiaxial),
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
  flowNormal(N_DIMS),
  plasticLoading(false),
  elasticMatrix(Matrix(N_DIMS, N_DIMS)),
  stiffnessInitial(Matrix(N_DIMS, N_DIMS)),
  stiffnessConverged(Matrix(N_DIMS, N_DIMS)),
  stiffnessTrial(Matrix(N_DIMS, N_DIMS))
{
  // Set the number of backstresses
  nBackstresses = cK.size();
  for (unsigned int i = 0; i < nBackstresses; ++i) {
    alphaKTrial.push_back(Vector(N_DIMS));
    alphaKConverged.push_back(Vector(N_DIMS));
  }

  // Zero all the vectors and matrices
  revertToStart();

  // Set elastic parameters and elastic stiffness matrix
  calculateElasticStiffness();
  stiffnessInitial = elasticMatrix;
  stiffnessTrial = elasticMatrix;
  stiffnessConverged = elasticMatrix;
};

/* ----------------------------------------------------------------------------------------------------------------- */

UVCmultiaxial::UVCmultiaxial()
  : NDMaterial(0, ND_TAG_UVCmultiaxial),
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
  strainConverged(Vector(N_DIMS)),
  strainTrial(Vector(N_DIMS)),
  strainPlasticConverged(Vector(N_DIMS)),
  strainPlasticTrial(Vector(N_DIMS)),
  strainPEqConverged(0.),
  strainPEqTrial(0.),
  stressConverged(Vector(N_DIMS)),
  stressTrial(Vector(N_DIMS)),
  flowNormal(Vector(N_DIMS)),
  plasticLoading(false),
  elasticMatrix(Matrix(N_DIMS, N_DIMS)),
  stiffnessInitial(Matrix(N_DIMS, N_DIMS)),
  stiffnessConverged(Matrix(N_DIMS, N_DIMS)),
  stiffnessTrial(Matrix(N_DIMS, N_DIMS))
{
  // Set the number of backstresses
  // todo: this probably wont work with parallel processing?
  nBackstresses = cK.size();
  for (unsigned int i = 0; i < nBackstresses; ++i) {
    alphaKTrial.push_back(Vector(N_DIMS));
    alphaKConverged.push_back(Vector(N_DIMS));
  }

  // Set elastic parameters and elastic stiffness matrix
  calculateElasticStiffness();
  stiffnessInitial = elasticMatrix;
  stiffnessTrial = elasticMatrix;
  stiffnessConverged = elasticMatrix;
}

/* ----------------------------------------------------------------------------------------------------------------- */

UVCmultiaxial::~UVCmultiaxial() {

}

/* ----------------------------------------------------------------------------------------------------------------- */

int UVCmultiaxial::returnMapping(){
  // Initialize all the variables
  int ret_val = 0;
  bool converged = true;
  unsigned int iterationNumber = 0;
  double sigmaY1 = 0.;
  double sigmaY2 = 0.;

  Vector alpha = Vector(N_DIMS);
  Vector alphaUpd = Vector(N_DIMS);
  Vector stressRelative = Vector(N_DIMS);
  Vector stressDeviatoric = Vector(N_DIMS);
  Vector alphaDiff = Vector(N_DIMS);
  double yieldFunction;
  double yieldStress;
  double stressRelativeNorm;
  double stressHydro;
  double consistParam = 0.;
  double isotropicModulus;
  double eK;
  double kinematicModulus;
  double aDotN;
  double pMultNumer, pMultDenom;

  // Elastic trial step
  alpha.Zero();
  for (unsigned int i = 0; i < nBackstresses; ++i)
    alpha = alpha + alphaKConverged[i];
  stressTrial = elasticMatrix * (strainTrial - strainPlasticConverged);
  stressHydro = (stressTrial(0) + stressTrial(1) + stressTrial(2)) / 3.;
  stressDeviatoric = stressTrial;
  for (unsigned int i = 0; i < N_DIRECT; ++i)
    stressDeviatoric[i] = stressTrial[i] - stressHydro;
  stressRelative = stressDeviatoric - alpha;
  stressRelativeNorm = sqrt(dotprod6(stressRelative, stressRelative));
  flowNormal = stressRelative / (RETURN_MAP_TOL + stressRelativeNorm);

  // Yield condition
  yieldStress = calculateYieldStress();
  isotropicModulus = calculateIsotropicModulus();
  yieldFunction = stressRelativeNorm - sqrt(2. / 3.) * yieldStress;
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
    kinematicModulus = 0.;
    alphaUpd.Zero();
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      kinematicModulus += cK[i] * eK - sqrt(2. / 3.) * gammaK[i] * eK * dotprod6(flowNormal, alphaKConverged[i]);
      alphaUpd += eK * alphaKConverged[i] + sqrt(2. / 3.) * cK[i] / gammaK[i] * (1. - eK) * flowNormal;
    }
    aDotN = dotprod6(alphaUpd - alpha, flowNormal);

    // Local Newton step
    pMultNumer = stressRelativeNorm - (2. * shearModulus * consistParam + sqrt(2. / 3.) * yieldStress + aDotN);
    pMultDenom = -2.0 * shearModulus * (1. + (kinematicModulus + isotropicModulus) / (3. * shearModulus));
    consistParam = consistParam - pMultNumer / pMultDenom;
    strainPEqTrial = strainPEqConverged + sqrt(2. / 3.) * consistParam;

    // Check convergence
    if (abs(pMultNumer) < RETURN_MAP_TOL) {
      converged = true;
    }
  }

  // Condition for plastic loading is whether or not iterations were performed
  if (iterationNumber == 0) {
    plasticLoading = false;
  }
  else {
    plasticLoading = true;
    // Update the internal variables
    for (unsigned int i = 0; i < N_DIRECT; ++i)
      strainPlasticTrial(i) = strainPlasticConverged(i) + consistParam * flowNormal(i);
    for (unsigned int i = N_DIRECT; i < N_DIMS; ++i)  // 2x since engineering strain definition
      strainPlasticTrial(i) = strainPlasticConverged(i) + 2. * consistParam * flowNormal(i);
    stressTrial = elasticMatrix * (strainTrial - strainPlasticTrial);
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      alphaKTrial[i] = eK * alphaKConverged[i] + sqrt(2. / 3.) * cK[i] / gammaK[i] * (1. - eK) * flowNormal;
    }
    alphaDiff = alphaUpd - alpha;
  }

  // Update the stiffness
  calculateStiffness(consistParam, stressRelativeNorm, alphaDiff);

  // Warn the user if the algorithm did not converge
  if (iterationNumber >= MAXIMUM_ITERATIONS - 1) {
    std::cerr << "UVCmultiaxial::returnMapping return mapping in UVCmultiaxial did not converge!" << endln;
    std::cerr << "\tDelta epsilon 11 = " << strainTrial[0] - strainConverged[0] << std::endl;
    std::cerr << "\tDelta epsilon 22 = " << strainTrial[1] - strainConverged[1] << std::endl;
    std::cerr << "\tDelta epsilon 12 = " << strainTrial[3] - strainConverged[3] << std::endl;
    std::cerr << "\tExiting with yield function = " << pMultNumer << " > " << RETURN_MAP_TOL << std::endl;
    ret_val = -1;
  }

  return ret_val;
}

/* ----------------------------------------------------------------------------------------------------------------- */

void UVCmultiaxial::calculateStiffness(double consistParam, double stressRelativeNorm, Vector alphaDiff) {
  if (!plasticLoading) {
    stiffnessTrial = elasticMatrix;
  }
  else  // plastic loading
  {
    double yieldStress, isotropicModulus, kinematicModulus, eK, beta, theta_1, theta_2, theta_3,
      id2OutId2, nOutN, alphaOutN;
    // 2nd order identity tensor
    std::vector<double> id2 = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
    // Symmetric 4th order identity tensor
    Matrix id4 = Matrix(N_DIMS, N_DIMS);
    for (unsigned int i = 0; i < N_DIRECT; ++i)
      id4(i, i) = 1.0;
    for (unsigned int i = N_DIRECT; i < N_DIMS; ++i)
      id4(i, i) = 1.0 / 2.0;

    // Isotropic hardening parameters
    yieldStress = calculateYieldStress();
    isotropicModulus = calculateIsotropicModulus();
    // Kinematic hardening parameters
    kinematicModulus = 0.;
    for (unsigned int i = 0; i < nBackstresses; ++i) {
      eK = calculateEk(i);
      kinematicModulus += cK[i] * eK - sqrt(2. / 3.) * gammaK[i] * eK * dotprod6(flowNormal, alphaKConverged[i]);
    }

    beta = 1.0 + (kinematicModulus + isotropicModulus) / (3.0 * shearModulus);
    theta_1 = 1.0 - 2.0 * shearModulus * consistParam / stressRelativeNorm;
    theta_3 = 1.0 / (beta * stressRelativeNorm);
    theta_2 = 1.0 / beta + (dotprod6(flowNormal, alphaDiff)) * theta_3 - (1.0 - theta_1);
    stiffnessTrial.Zero();
    for (unsigned int i = 0; i < N_DIMS; ++i) {
      for (unsigned int j = 0; j < N_DIMS; ++j) {
        id2OutId2 = id2[i] * id2[j];
        nOutN = flowNormal[i] * flowNormal[j];
        alphaOutN = alphaDiff[i] * flowNormal[j];
        stiffnessTrial(i, j) = bulkModulus * id2OutId2
          + 2. * shearModulus * theta_1 * (id4(i, j) - 1. / 3. * id2OutId2)
          - 2. * shearModulus * theta_2 * nOutN
          + 2. * shearModulus * theta_3 * alphaOutN;
      }
    }
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
int UVCmultiaxial::setTrialStrain(const Vector &v) {
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
int UVCmultiaxial::setTrialStrain(const Vector &v, const Vector &r) {

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
int UVCmultiaxial::setTrialStrainIncr(const Vector &v) {

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
int UVCmultiaxial::setTrialStrainIncr(const Vector &v, const Vector &r) {

  // Reset the trial state
  revertToLastCommit();

  // Set the trial strain
  strainTrial += v;

  // Do the return mapping and calculate the tangent modulus
  returnMapping();

  return 0;
}


/* ----------------------------------------------------------------------------------------------------------------- */

const Vector &UVCmultiaxial::getStrain() {
  return strainTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Vector &UVCmultiaxial::getStress() {
  return stressTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Matrix &UVCmultiaxial::getTangent() {
  return stiffnessTrial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

const Matrix &UVCmultiaxial::getInitialTangent() {
  // todo: can make more efficient by changing this to elasticMatrix and removing stiffnessInitial as a variable
  return stiffnessInitial;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @return 0 if successful
*/
int UVCmultiaxial::commitState() {
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
int UVCmultiaxial::revertToLastCommit() {
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
int UVCmultiaxial::revertToStart() {
  strainConverged.Zero();
  strainPlasticConverged.Zero();
  strainPEqConverged = 0.;
  stressConverged.Zero();
  flowNormal.Zero();
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
NDMaterial* UVCmultiaxial::getCopy() {

  UVCmultiaxial* theCopy;
  theCopy = new UVCmultiaxial(this->getTag(), elasticModulus, poissonRatio,
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
  theCopy->flowNormal = flowNormal;
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
NDMaterial* UVCmultiaxial::getCopy(const char *code) {
  if (strcmp(code, getType()) == 0) {
    UVCmultiaxial* theCopy;
    theCopy = new UVCmultiaxial(this->getTag(), elasticModulus, poissonRatio,
      initialYield, qInf, bIso,
      dInf, aIso,
      cK, gammaK);
    return theCopy;
  }
  else {
    // todo: change to opserr
    std::cerr << "UVCmultiaxial::getCopy invalid NDMaterial type, expecting " << code << std::endl;
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
int UVCmultiaxial::sendSelf(int commitTag, Channel &theChannel) {

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
  opserr << "UVCmultiaxial::sendSelf() - failed to sendSelf\n";
  return -1;
  }
  */

  opserr << "Fatal: Paralleliziation for UVCmultiaxial is not implemented yet!" << endln;
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
int UVCmultiaxial::recvSelf(int commitTag, Channel &theChannel,
  FEM_ObjectBroker &theBroker) {
  /*
  static Vector data(26);

  if (theChannel.recvVector(this->getDbTag(), commitTag, data) < 0) {
  opserr << "UVCmultiaxial::recvSelf() - failed to recvSelf\n";
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

  opserr << "Fatal: Paralleliziation for UVCmultiaxial is not implemented yet!" << endln;
  return -1;
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param s the opensees output stream
* @param flag is 2 for standard output, 25000 for JSON output
(see OPS_Globals.h)
*/
void UVCmultiaxial::Print(OPS_Stream &s, int flag) {

  // todo: change these back when not only .dll
  // if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
  if (flag == 2) {
    s << "UVCmultiaxial tag: " << this->getTag() << endln;
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
    s << "\"type\": \"UVCmultiaxial\", ";
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

void UVCmultiaxial::calculateElasticStiffness() {
  double id2OutId2;
  // 2nd order identity tensor
  //std::vector<double> id2 = { 1.0, 1.0, 1.0, 0.0, 0.0, 0.0 };
  std::vector<double> id2(6);
  id2[0] = id2[1] = id2[2] = 1.0;
  id2[3] = id2[4] = id2[5] = 0.0;
  // Symmetric 4th order identity tensor
  Matrix id4 = Matrix(N_DIMS, N_DIMS);
  for (unsigned int i = 0; i < N_DIRECT; ++i)
    id4(i, i) = 1.0;
  for (unsigned int i = N_DIRECT; i < N_DIMS; ++i)
    id4(i, i) = 1.0 / 2.0;
  for (unsigned int i = 0; i < N_DIMS; ++i) {
    for (unsigned int j = 0; j < N_DIMS; ++j) {
      id2OutId2 = id2[i] * id2[j];
      elasticMatrix(i, j) = id2OutId2 * bulkModulus + 2. * shearModulus * (id4(i, j) - 1. / 3. * id2OutId2);
    }
  }
}

/* ----------------------------------------------------------------------------------------------------------------- */

/**
*
* @param v1 vector of length 6 that represents a symmetric 2nd order tensor
* @param v2 vector of length 6 that represents a symmetric 2nd order tensor
* @return the dot product of v1 and v2
*
* This function assumes that the last 3 components represent the symmetric terms.
*/
double UVCmultiaxial::dotprod6(Vector v1, Vector v2) {
  double val = 0.;
  for (unsigned int i = 0; i < N_DIRECT; ++i)
    val += v1[i] * v2[i];
  for (unsigned int i = N_DIRECT; i < N_DIMS; ++i)
    val += 2.0 * (v1[i] * v2[i]);
  return val;
}

/* ----------------------------------------------------------------------------------------------------------------- */

// todo: these methods should inherit from a mother class for multiaxial and plane stress -> wont work for dll?
double UVCmultiaxial::calculateYieldStress() {
  double sigmaY1, sigmaY2;
  sigmaY1 = qInf * (1. - exp(-bIso * strainPEqTrial));
  sigmaY2 = dInf * (1. - exp(-aIso * strainPEqTrial));
  return initialYield + sigmaY1 - sigmaY2;
}

/* ----------------------------------------------------------------------------------------------------------------- */

double UVCmultiaxial::calculateIsotropicModulus(){
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
double UVCmultiaxial::calculateEk(unsigned int i) {
  return exp(-gammaK[i] * (strainPEqTrial - strainPEqConverged));
}

/* ----------------------------------------------------------------------------------------------------------------- */
