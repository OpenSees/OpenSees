// Written: Quan Gu and Zhijian Qiu
// Created: 2013.7
//
// Reference:JP Conte, MK. Jagannath, Seismic relibility analysis of concrete
// gravity dams, A Report on Research, Rice University, 1995. 
//           EA de Souza Neto, D Peri´c, DRJ Owen, Computational methods for 
// plasticity, Theory and applications (see pages 357 to 366), 2008.
// 
// 3D J2 plasticity model with linear isotropic and kinematic hardening
//  
// -------------------

#ifndef PlaneStressSimplifiedJ2_h
#define PlaneStressSimplifiedJ2_h

#include <NDMaterial.h>

#include <T2Vector.h> 
#include <Matrix.h>
#include <Vector.h>
 
class PlaneStressSimplifiedJ2 : public NDMaterial
{
 public:

  PlaneStressSimplifiedJ2 (int tag, 
			   int nd,
			   NDMaterial &the3DMaterial);
  
  
  virtual ~PlaneStressSimplifiedJ2 ();
  
  const char *getClassType(void) const {return "PlaneStressSimplifiedJ2";};
  const char *getType(void) const{ return "PlaneStress";};
  int setTrialStrain (const Vector &strain);
  int setTrialStrain(const Vector &v, const Vector &r);
  int setTrialStrainIncr(const Vector &v);
  int setTrialStrainIncr(const Vector &v, const Vector &r);
  
  // Calculates current tangent stiffness.
  const Matrix &getTangent (void);
  const Matrix &getInitialTangent (void);
  
  // Calculates the corresponding stress increment (rate), for a given strain increment. 
  const Vector &getStress (void);
  const Vector &getStrain (void);
  const Vector &getCommittedStress (void);
  const Vector &getCommittedStrain (void);
  
  /*
    int setTrialStrain (const Tensor &v) {return 0;}
    int setTrialStrain (const Tensor &v, const Tensor &r) {return 0;}
    int setTrialStrainIncr (const Tensor &v) {return 0;}
    int setTrialStrainIncr (const Tensor &v, const Tensor &r) {return 0;}
  */
  
  int commitState (void);
  int revertToLastCommit (void);
  int revertToStart(void);
  
  NDMaterial *getCopy (void);
  NDMaterial *getCopy (const char *type);
  
  int plastIntegrator();
  
  
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel, 
	       FEM_ObjectBroker &theBroker);    
  
  Response *setResponse (const char **argv, int argc, OPS_Stream &s);
  int getResponse (int responseID, Information &matInformation);
  void Print(OPS_Stream &s, int flag =0);
  
  int setParameter(const char **argv, int argc, Parameter &param);
  int updateParameter(int responseID, Information &eleInformation);	
  
 protected:
  
 private:
  
  int ndm;
  
  NDMaterial * the3DMaterial;
  
  // --- trial variables
  Vector stress;
  Vector strain;

  // ---- committed variables
  Vector Cstress;
  Vector Cstrain;
  
  // --
  Matrix theTangent; 
  double savedStrain33;
  double CsavedStrain33;
  // ---  define classwide variables

  static Vector tmpVector;
  static Matrix tmpMatrix;
};
#endif

