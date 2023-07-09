/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright for this material model is by Quan Gu  Xiamen Unv.   **
**                                                                    **
** Commercial use of this program without express permission of the   **
** Authors is strictly prohibited.                                    **
**                                                                    **
** Developed by:                                                      **
**   Quan Gu (quangu@xmu.edu.cn)                                      **
**                                                                    **
** reference: Hofstetter, G., J.C. Simo and R.L. Taylor, \A modified  **
** cap-model: Closest point solution algorithms\, Computers &         **
** Structures, 46 (1993),203-214.                                     **
**                                                                    **
** ****************************************************************** */



#ifndef CapPlasticity_h
#define CapPlasticity_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <NDMaterial.h>


class CapPlasticity : public NDMaterial {
  
  //-------------------Declarations-------------------------------
  
 public :
  
  
  CapPlasticity( int    tag,
		 double G,
		 double K,
		 double rho,
		 double X,
		 double D,
		 double W,
		 double R,
		 double lambda,
		 double theta,
		 double beta,
		 double alpha,
		 double T, 
		 int ndm,
		 double pTol_k
		 ) ;
  CapPlasticity( const CapPlasticity & a);
  
  ~CapPlasticity( ) ;
  
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
  
  Response*setResponse (const char **argv, int argc, OPS_Stream &output);

  int getResponse (int responseID, Information &matInformation);
  
  void Print(OPS_Stream &s, int flag = 0) ;	
  
 private:
  
  double CapBoundL(double k);
  double CapBoundX(double k); 
  double failureEnvelop(double I1);                           //Fe
  double CapSurface(double normS,double I1, double k);         //Fc
  double failureEnvelopDeriv(double I);
  double hardeningParameter_H(double k1, double k2); 
  double Newton_k (double tol, int mode /*, double normS, double I1, double k*/); // mode = 1,or 2, or 5-2)
  double Newton_I1(double tol, int mode, double normTS, double I1_trial ); // mode =5-1, or 3
  double Bisection(double tol, double normTS, double I1_trial ); // mode = 3
  
  int    findMode(double normS, double I1, double k);
  
  // ----------------------------  for consistent tangent modulus ---------------
  Matrix & dF2dSigma ( int mode);
  Vector & dFdSigma (int mode);
  double dFdk (int OrderOfDerivative);  // order one or two, depending on OrderOfDerivative
  Vector dF2dSigmadk ( void);
  double dHdk(double k2); 
  double dFdIdk(void);
  
  int computeConsistentTangent(double gammar1, double gammar2, double gammar3, int mode);
  double tripleTensorProduct (Vector &A, Matrix &B, Vector &C);  // result = A:B:C
  
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
  
  // ---- history variables ---
  
  Vector CStrain;
  Vector CPlastStrain;
  Vector CStress;
  double CHardening_k;  // hardening parameter
  
  // --- trial variables ---
  
  Vector strain;
  Vector plastStrain;
  Vector stress;

  double stressI1;     // trace of stress at n+1 step
  Vector stressDev;    // stress at n+1 step
  Matrix theTangent;   
  
  double hardening_k;

  int debug;

// ---  classwide variables --

  static Matrix tempMatrix;
  static Vector tempVector;


	////////////////////add sensitivity ////////////////////////
public:
  int            setParameter             (const char **argv, int argc, Information &info);
  int            updateParameter          (int parameterID, Information &info);
  int            activateParameter        (int parameterID);
  const Vector & getStressSensitivity     (int gradNumber, bool conditional);
  //	int            commitSensitivity        (Vector & strainGradient, int gradNumber, int numGrads);
  //    const Matrix & getInitialTangentSensitivity(int gradNum);
  
  
private:

	Matrix *SHVs;
	int parameterID;
	bool isKAdjusted; // 6-14-2013    0: not adjusted; 1: adjusted
	

} ; 

#endif
