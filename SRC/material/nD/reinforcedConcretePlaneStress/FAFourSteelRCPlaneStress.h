#ifndef FAFourSteelRCPlaneStress_h
#define FAFourSteelRCPlaneStress_h

// File: FAFourSteelRCPlaneStress.h
//
// Written: JZhong
// Created: 2004.11
//
// Description: This file contains the class definition for 
// FAFourSteelRCPlaneStress material.
// Hsu's Model 2002
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class FAFourSteelRCPlaneStress : public NDMaterial
{
  public:
    FAFourSteelRCPlaneStress ( int      tag, 
		                           double   RHO,
								   UniaxialMaterial *s1,
								   UniaxialMaterial *s2,
								   UniaxialMaterial *s3,
								   UniaxialMaterial *s4,
								   UniaxialMaterial *c1,
								   UniaxialMaterial *c2,
								   double   ANGLE1,
								   double   ANGLE2,
								   double   ANGLE3,
								   double   ANGLE4,
								   double   ROU1,
								   double   ROU2,
								   double   ROU3,
								   double   ROU4,
								   double   FPC,
								   double   FY,
								   double   E,
								   double   EPSC0);					  
	FAFourSteelRCPlaneStress();
	~FAFourSteelRCPlaneStress();				  
								  
    double getRho(void);

    int setTrialStrain(const Vector &v); // really used one
    int setTrialStrain(const Vector &v, const Vector &r);
    int setTrialStrainIncr(const Vector &v);
    int setTrialStrainIncr(const Vector &v, const Vector &r);
    const Matrix &getTangent(void);
    const Matrix &getInitialTangent(void) {return this->getTangent();};

    const Vector &getStress(void);
    const Vector &getStrain(void);
    
    const Vector &getCommittedStress(void);
    const Vector &getCommittedStrain(void);    
    
    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);
    
    NDMaterial *getCopy(void);
    NDMaterial *getCopy(const char *type);

    void Print(OPS_Stream &s, int flag = 0);
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    const char *getType(void) const { return "PlaneStress"; };
    int getOrder(void) const { return 3;};

  protected:
    
  private:

    double   rho; 
    UniaxialMaterial **theMaterial; // pointer of the materials 
    Response **theResponses; // pointer to material responses needed for Concrete

    double   angle1;    // angel of the first steel layer to x coordinate 
    double   angle2;    // angel of the second steel layer to x coordinate
    double   angle3;    // angel of the third steel layer to x coordinate 
    double   angle4;    // angel of the forth steel layer to x coordinate

    double   rou1;      // steel ratio of the first steel layer
    double   rou2;      // steel ratio of the second steel layer
    double   rou3;      // steel ratio of the third steel layer
    double   rou4;      // steel ratio of the forth steel layer

    double   fpc;       // compressive strength of the concrete
    double   fy;        // yield stress of the bare steel bar
    double   E0;        // young's modulus of the steel
    double   epsc0;     // compressive strain of the concrete corresponding to fpc

    double   Tstress[3];  // Trial stresses
    double   lastStress[3];  // Last committed stresses, added for x, k
    
    int      steelStatus;  // check if steel yield, 0 not yield, 1 yield
    int      dirStatus;    // check if principle direction has exceed 90 degree, 1 yes, 0 no
    
    double   citaStrain;      // principle strain direction
    double   citaStress;     // principle stress direction
    double   miu12;        // Hsu/Zhu ratio
    double   miu21;        // Hsu/Zhu ratio
    double   G12;          // Shear Modulus
    
    
    // for damgage factor D=1-0.4*epslonC'/epslon0; epslon0=0.002
	
    // Trial values
    int TOneReverseStatus;         // Trial reverse status for concrete One, 1 reverse, 0 no
    double TOneNowMaxComStrain;
    double TOneLastMaxComStrain;
    
    int TTwoReverseStatus;         // Trial reverse status for concrete Two, 1 reverse, 0 no
    double TTwoNowMaxComStrain;
    double TTwoLastMaxComStrain;
    
    // Converged values
    int COneReverseStatus;         // Converged reverse status for concrete One, 1 reverse, 0 no
    double COneNowMaxComStrain;
    double COneLastMaxComStrain;
    
    int CTwoReverseStatus;         // Converged reverse status for concrete Two, 1 reverse, 0 no
    double CTwoNowMaxComStrain;
    double CTwoLastMaxComStrain;
    
    double DDOne; // damage factor for concrete One
    double DDTwo; // damage factor for concrete Two
    
    
    Vector strain_vec;
    Vector stress_vec;	
    Matrix tangent_matrix;


    int determineTrialStress(void);
    double getPrincipalStressAngle(double inputAngle);
    double getAngleError(double inputCita);
 
};

#endif
