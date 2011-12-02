#ifndef PrestressedConcretePlaneStress_h
#define PrestressedConcretePlaneStress_h

// File: PrestressedConcretePlaneStress.h
//
// Written: ALaskar
// Created: 2007.09
//
// Description: This file contains the class definition for 
// PrestressedConcretePlaneStress material.
// Hsu's Model 2002
// For Detailed explanation of the model, please refer to the book
// entitled "Unified Theory of Concrete Structures,"
// by Thomas T.C. Hsu and Y.L. Mo, John Wiley & Sons, April 2010.

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class PrestressedConcretePlaneStress : public NDMaterial
{
  public:
    PrestressedConcretePlaneStress ( int      tag, 
				     double   RHO,
				     UniaxialMaterial *t1,
				     UniaxialMaterial *s1,
				     UniaxialMaterial *c1,
				     UniaxialMaterial *c2,
				     double   ANGLE1,
				     double   ANGLE2,
				     double   ROU1,
				     double   ROU2,
				     double	PSTRAIN,
				     double   FPC,
				     double   FY1,
				     double	FY2,
				     double   E,
				     double   EPSC0);					  
    PrestressedConcretePlaneStress();
    ~PrestressedConcretePlaneStress();				  
    
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
    Response **theResponses;

    double   angle1;    // angel of the tendon layer to x coordinate 
    double   angle2;    // angel of the steel layer to x coordinate
    double   rou1;      // steel ratio of the first steel layer
    double   rou2;      // steel ratio of the second steel layer
    double	 pstrain;	//initial strain in tendons
    double   fpc;       // compressive strength of the concrete
    double   fy1;       // yield stress of the bare tendons
    double	 fy2;		// yield stress of bare steel bar
    double   E0;        // young's modulus of the steel
    double   epsc0;     // compressive strain of the concrete
    double   Tstress[3];  // Trial stresses
    double   lastStress[3];  // Last committed stresses, added for x, k
    
    int      steelStatus;  // check if steel yield, 0 not yield, 1 yield
    int      dirStatus;    // check if principle direction has exceed 90 degree, 1 yes, 0 no
    
    double   citaStrain;      // principle strain direction
    double   citaStress;     // principle stress direction
    double   miu12;        // Hsu/Zhu ratio
    double   miu21;        // Hsu/Zhu ratio
    double   G12;
    

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
    
    double tt1;
    double tt2;
    double xxx;
    double kkk;
    
    double beta1;
    double beta2;
    
    Vector strain_vec;
    Vector stress_vec;	
    Matrix tangent_matrix;
    
    int determineTrialStress(void);
 
};

#endif
