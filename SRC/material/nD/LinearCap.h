// Written: Quan Gu and Zhijian Qiu  
// Created: 2015/01/25 

//------------------------------------------


#ifndef LinearCap_h
#define LinearCap_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>


class LinearCap : public NDMaterial {

//-------------------Declarations-------------------------------

  public :

 
  LinearCap( int    tag,
			     double G,
			     double K,
			     double rho,
			     double theta,
				 double alpha,
				 double T, 
				 int ndm,
				 double pTol_k
			     ) ;
  LinearCap( const LinearCap & a);
  
      ~LinearCap( ) ;

      double getRho(void);
      int setTrialStrain(const Vector &v);
      int setTrialStrain(const Vector &v, const Vector &r);
      int setTrialStrainIncr(const Vector &v);
      int setTrialStrainIncr(const Vector &v, const Vector &r);
      const Matrix &getTangent(void);
      const Matrix &getInitialTangent(void) ;

      const Vector &getStress(void);
      const Vector &getStrain(void);

      int commitState(void) ;
      int revertToLastCommit(void) ;
      int revertToStart(void) ;

      NDMaterial *getCopy(void);
      NDMaterial *getCopy(const char *code) ;

      const char *getType(void) const;
      int getOrder(void) const ;
	  
	  int sendSelf(int commitTag, Channel &theChannel) ;
	  int recvSelf(int commitTag, Channel &theChannel,
				   FEM_ObjectBroker &theBroker ) ;

      Response *setResponse (const char **argv, int argc, OPS_Stream &output);
      int getResponse (int responseID, Information &matInformation);
	  
	  void Print(OPS_Stream &s, int flag = 0) ;	
  
  private:
 
	  double CapBoundL(double k);
	  double CapBoundX(double k); 
	  double failureEnvelop(double I1);                           //Fe
	  double CapSurface(double normS,double I1, double k);         //Fc
	  double failureEnvelopDeriv(double I); 

	  int    findMode(double normS, double I1);

// ----------------------------  for consistent tangent modulus ---------------

	  int computeConsistentTangent(double gammar1, double gammar2, double gammar3, int mode);


  private:


  //int    tag;
  double shearModulus;
  double bulkModulus;
  double rho;
  double X;
  double D;
  double W;
  double R;
  double lambda;
  double theta;
  double beta;
  double alpha;
  double T;
  double deltPlastStrainI1;


  int theMode; // mode =1-6:  
  /* (1) the tension curoff mode: f3
     (2) the tension corner
	 (3) the cap mode:  f2
	 (4) the compressive corner mode
	 (5) the failure mode:   f1
	 (6) the elastic region
  */

  int ndm; 
  double tol_k;
  int flag; // show whether Newton is converged in mode =3

  double stressI1;     // trace of stress at n+1 step
  Vector stressDev;    // stress at n+1 step

  Matrix theTangent; 

// ---- history variables ---

  Vector CStrain;
  Vector CPlastStrain;
  Vector CStress;
  double CHardening_k;  // hardening parameter

// --- trial variables ---

  Vector strain;
  Vector plastStrain;
  Vector stress; 
  double hardening_k;

  int debug;

// ---  classwide variables --

  static Matrix tempMatrix;
  static Vector tempVector;

// ---------------------sensitivity -----------------------
public:
    int            setParameter             (const char **argv, int argc, Parameter &param);
    int            updateParameter          (int parameterID, Information &info);
	int            activateParameter        (int parameterID);
	const Vector & getStressSensitivity     (int gradNumber, bool conditional);
	int            commitSensitivity        (const Vector & strainGradient, int gradNumber, int numGrads);

private:	
 
	Matrix *SHVs;
	int parameterID;

};
#endif
