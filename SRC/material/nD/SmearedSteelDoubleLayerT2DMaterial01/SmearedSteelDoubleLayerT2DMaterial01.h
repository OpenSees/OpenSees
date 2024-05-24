// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 03/2023
// 
// Description: This file contains the class definition for SmearedSteelDoubleLayerT2DMaterial01. 
// A SmearedSteelDoubleLayerT2DMaterial01 is a subclass of the class NDMaterial and corresponds to the abstract representation
// of a double perpendicular Smeared Steel layers (plane stress) 2D Material with a tangent formulation. Each layer works only 
// in the direction of the bars, so a uniaxial constitutive model is used to represent the behavior of reinforcing steel bars in each direction.
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\nD\SmearedSteelDoubleLayerT2DMaterial01
//
// Rev: 1.0

#ifndef SmearedSteelDoubleLayerT2DMaterial01_h
#define SmearedSteelDoubleLayerT2DMaterial01_h

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class SmearedSteelDoubleLayerT2DMaterial01 : public NDMaterial
{
  public:

	  SmearedSteelDoubleLayerT2DMaterial01(int tag,					// nDMaterial tag
		  UniaxialMaterial* steelULayer1Obj,						// steel layer 1
		  UniaxialMaterial* steelULayer2Obj,						// steel layer 2
		  double ratioSteelLayer1,									// ratio of smeared steel layer 1
		  double ratioSteelLayer2,									// ratio of smeared steel layer 2
		  double OrientationEmbeddedSteel);							// orientation of the smeared steel layer

	  SmearedSteelDoubleLayerT2DMaterial01();

	  ~SmearedSteelDoubleLayerT2DMaterial01();

	  const char* getClassType(void) const { return "SmearedSteelDoubleLayerT2DMaterial01";};
	  int setTrialStrain(const Vector& v);
	  int setTrialStrain(const Vector& v, const Vector& r);
	  int setTrialStrainIncr(const Vector& v);
	  int setTrialStrainIncr(const Vector& v, const Vector& r);
	  const Matrix& getTangent(void);
	  const Matrix& getInitialTangent(void);

	  Response* setResponse(const char** argv, int argc, OPS_Stream& theOutputStream);
	  int getResponse(int responseID, Information& matInformation);

	  const Vector& getStress(void);
	  const Vector& getStrain(void);

	  // public methods to set the state of the material
	  int commitState(void);
	  int revertToLastCommit(void);
	  int revertToStart(void);

	  // Create a copy of material parameters and state variables
	  NDMaterial* getCopy(void);
	  // Create a copy of just the material parameters
	  NDMaterial* getCopy(const char* type);
	  
	  // Return a string indicating the type of material model
	  const char* getType(void) const { return "PlaneStress"; };
	  int getOrder(void) const { return 3; };

	  int sendSelf(int commitTag, Channel& theChannel);
	  int recvSelf(int commitTag, Channel& theChannel,
		  FEM_ObjectBroker& theBroker);

	  void Print(OPS_Stream& s, int flag = 0);

  protected:

  private:
	  
	  int setTrialStrainPrincipalDirection(const Vector& v);									// calculate trial stress and tangent 
	  void calculateStrainPrincipalDirections01(void);											// calculate the principal direction for the strains (11, 22, 12) using the calculateAngle01 method
	  void calculateAngle01(double cosTheta, double sinTheta, double& theta);					// calculate the theta angle [-pi,pi] from the cos(theta) and sin(theta)
	  void calculateStrainTransformationMatrix(double* pTmatrixStrain, double theta);			// calculate the Strain Transformation Matrix that goes from the Local Coord System to the orientation of the Principal Direction of strain

	  // Functions for recorders
	  Vector getStrainStressSteel1(void);
	  Vector getStrainStressSteel2(void);

	  UniaxialMaterial** theMaterial;		// pointer to a material

	  double ratioLayer1;					// store the reinforcing ratio of the smeared steel layer 1
	  double ratioLayer2;					// store the reinforcing ratio of the smeared steel layer 1
	  double thetaSmearedSteel;				// store the orientation of the smeared steel layers
	  // Material history variables..............................................................................
	  double thetaPrincipalDirection;	    // store the orientation of the principal direction in the material
	  Vector strainPrincipalDirection;	    // store the principal strains
	  // Committed state variables...............................................................................
	  Vector CstrainLayer;					// store the Committed Strain for the Layer
	  Vector Cstrain;						// store the Committed Strain
	  Vector Cstress;						// store the Committed Stress
	  Vector Ctangent;						// store the Committed Tangent 
	  // Trial state variables...................................................................................
	  Vector stressNDM;				        // store the Trial Stress
	  Matrix tangentNDM;			        // store the Trial Tangent
	  Matrix initialTangentNDM;		

	  Vector TstrainLayer;			        // store the Trial Strain for the Layer
	  Vector Tstrain;				        // store the Trial Strain
	  Vector Tstress;				        // store the Trial Stress
	  Vector Ttangent;				        // store the Trial Tangent

	  const double pi;
};

#endif
