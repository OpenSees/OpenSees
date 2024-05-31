//
// Created by Alex Hartloper on 17.04.18.
//

#ifndef CPP_UVC_UA_H
#define CPP_UVC_UA_H

#include <vector>
#include "UniaxialMaterial.h"
#include <OPS_Globals.h>
#include <elementAPI.h>

/* ------------------------------------------------------------------------ */

class UVCuniaxial : public UniaxialMaterial
{

  /* ------------------------------------------------------------------------ */
  /* Constructors/Destructors                                                 */
  /* ------------------------------------------------------------------------ */

public:
  //! Constructor
  UVCuniaxial(int tag, double E, double sy0, double qInf, double b,
    double dInf, double a,
    std::vector<double> cK, std::vector<double> gammaK);
  UVCuniaxial();
  //! Destructor
  ~UVCuniaxial();

  /* ------------------------------------------------------------------------ */
  /* Methods                                                                  */
  /* ------------------------------------------------------------------------ */

public:
  //! Returns the class type
  const char *getClassType(void) const { return "UVCuniaxial"; };

  //! Calculates the trail strain and stress
  int setTrialStrain(double strain, double strainRate = 0);

  //! Returns the trial strain
  double getStrain(void);

  //! Returns the trial stress
  double getStress(void);

  //! Returns the trial elastoplastic tangent modulus
  double getTangent(void);

  //! Returns the tangent modulus in the undeformed configuration
  double getInitialTangent(void);

  //! Sets the converged state to be the current trial state
  int commitState(void);

  //! Sets the trial state to be the converged state
  int revertToLastCommit(void);

  //! Sets the converged state to the undeformed configuration
  int revertToStart(void);

  //! Returns a copy of the material in the current state
  UniaxialMaterial *getCopy(void);

  //! todo: fill out
  int sendSelf(int commitTag, Channel &theChannel);

  //! todo: fill out
  int recvSelf(int commitTag, Channel &theChannel,
    FEM_ObjectBroker &theBroker);

  //! Adds the print information to the stream
  void Print(OPS_Stream &s, int flag = 0);

private:
  //! Determines the trial stress for the given strain increment
  void returnMapping(double strain_inc);

  //! Sets the elastoplastic tangent modulus based on the trial state
  void calculateStiffness();

  //! Returns the sign of val
  template <typename T> int sgn(T val);

  /* ------------------------------------------------------------------------ */
  /* Members                                                                  */
  /* ------------------------------------------------------------------------ */

private:
  // Parameters
  //const int N_BASIC_PARAMS = 4;
  //const int N_PARAM_PER_BACK = 2;
  //const double RETURN_MAP_TOL = 10.0e-10;
  //const int MAXIMUM_ITERATIONS = 1000;
  enum {MAX_BACKSTRESSES = 8};
  
  // Material properties, set by the constructor
  double elasticModulus;
  double yieldStress;
  double qInf;
  double bIso;
  double dInf;
  double aIso;
  double stiffnessInitial;
  std::vector<double> cK;
  std::vector<double> gammaK;
  int nBackstresses;

  // Internal variables
  double strainConverged;
  double strainTrial;
  double strainPEqConverged;  // Equivalent plastic strain
  double strainPEqTrial;
  double stressConverged;
  double stressTrial;
  std::vector<double> alphaKConverged;
  std::vector<double> alphaKTrial;
  double stiffnessConverged;
  double stiffnessTrial;
  double flowDirection;
  bool plasticLoading;
};

/* ------------------------------------------------------------------------ */

#endif //CPP_UVC_UA_H
