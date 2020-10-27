//
// Created by Alex Hartloper on 09.07.18.
//

#ifndef CPP_UVC_MA_H
#define CPP_UVC_MA_H

#include <vector>
#include "NDMaterial.h"
#include <OPS_Globals.h>
#include <elementAPI.h>
#include "Matrix.h"
#include "Vector.h"


/* ------------------------------------------------------------------------ */

class UVCmultiaxial : public NDMaterial
{

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  //! Constructor, called by clients
  UVCmultiaxial(int tag, double E, double poissonRatio, double sy0, 
    double qInf, double b, double dInf, double a,
    std::vector<double> cK, std::vector<double> gammaK);

  //! Constructor, parallel processing
  UVCmultiaxial(void);

  //! Destructor
  ~UVCmultiaxial();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  //! Returns the class type
  const char *getClassType(void) const { return "UVCmultiaxial"; };

  //! Returns the type of ND material
  const char* getType() const { return "ThreeDimensional"; };

  //! Returns the number of vector components 
  int getOrder() const { return 6; };

  //! Calculates the trial strain and stress, provided the total strain
  int setTrialStrain(const Vector &v);
  int setTrialStrain(const Vector &v, const Vector &r);

  //! Calculates the trial strain and stress, provided the strain increment
  int setTrialStrainIncr(const Vector &v);
  int setTrialStrainIncr(const Vector &v, const Vector &r);


  //! Returns the trial strain
  const Vector &getStrain(void);

  //! Returns the trial stress
  const Vector &getStress(void);

  //! Returns the trial elastoplastic tangent modulus
  const Matrix &getTangent(void);

  //! Returns the tangent modulus in the undeformed configuration
  const Matrix &getInitialTangent(void);

  //! Returns the mass density of the material - zero mass assumed
  double getRho(void) { return 0.; };

  //! Sets the converged state to be the current trial state
  int commitState(void);

  //! Sets the trial state to be the converged state
  int revertToLastCommit(void);

  //! Sets the converged state to the undeformed configuration
  int revertToStart(void);

  //! Returns a copy of the material in the current state
  NDMaterial *getCopy(void); 

  //! Returns a copy of the material without copying the state variables
  NDMaterial *getCopy(const char *code);

  //! todo: fill out
  int sendSelf(int commitTag, Channel &theChannel);

  //! todo: fill out
  int recvSelf(int commitTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker);

  //! Adds the print information to the stream
  void Print(OPS_Stream &s, int flag = 0);

private:
  //! Determines the trial stress for the given strain increment
  int returnMapping();

  //! Sets the elastoplastic tangent modulus based on the trial state
  void calculateStiffness(double plasticMultiplier, double stressRelativeNorm, 
    Vector alphaDiff);

  //! Returns the equivalent dot product of a 2nd order symmetric tensor
  double dotprod6(Vector v1, Vector v2);

  //! Calculates the elastic stiffness matrix
  void calculateElasticStiffness(void);

  //! Returns the current yield stress
  double calculateYieldStress(void);

  //! Returns the isotropic hardening modulus
  double calculateIsotropicModulus(void);

  //! Returns the current eK value
  // todo: inline this?
  double calculateEk(unsigned int i);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */

private:
  // Parameters
  const unsigned int N_BASIC_PARAMS = 5;
  const unsigned int N_PARAM_PER_BACK = 2;
  const double RETURN_MAP_TOL = 1.0e-10;
  const unsigned int MAXIMUM_ITERATIONS = 1000;
  const unsigned int N_DIRECT = 3;
  const unsigned int N_DIMS = 6;

  // Material properties, set by the constructor
  double elasticModulus;
  double shearModulus;
  double bulkModulus;
  double poissonRatio;
  double initialYield;
  double qInf;
  double bIso;
  double dInf;
  double aIso;
  Matrix stiffnessInitial;
  Matrix elasticMatrix;
  std::vector<double> cK;
  std::vector<double> gammaK;
  unsigned int nBackstresses;

  // Internal variables
  Vector strainConverged;
  Vector strainTrial;
  Vector strainPlasticConverged;
  Vector strainPlasticTrial;
  double strainPEqConverged;  // Equivalent plastic strain
  double strainPEqTrial;
  Vector stressConverged;
  Vector stressTrial;
  std::vector<Vector> alphaKConverged;
  std::vector<Vector> alphaKTrial;
  Matrix stiffnessConverged;
  Matrix stiffnessTrial;
  Vector flowNormal;
  bool plasticLoading;
};

/* ------------------------------------------------------------------------ */

#endif //CPP_UVC_MA_H
