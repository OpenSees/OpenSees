// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 12/2022
// 
// Description: This file contains the OrthotropicRotatingAngleConcreteT2DMaterial01 class definition.
// An OrthotropicRotatingAngleConcreteT2DMaterial01 is a subclass of the class NDMaterial and corresponds to the abstract representation
// of an Orthotropic Concrete Layer (plane stress) 2D Material with the Rotating Angle and Tangent Formulation for Cycling or Reversed
// Loading with damage that is used in Finite Element Method or Structural Analysis.
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massones, L. M. (2016). A nonlinear quadrilateral layered membrane with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\nD\OrthotropicRotatingAngleConcreteT2DMaterial01
//
// Rev: 1.0

#ifndef OrthotropicRotatingAngleConcreteT2DMaterial01_h
#define OrthotropicRotatingAngleConcreteT2DMaterial01_h

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <ID.h>

class OrthotropicRotatingAngleConcreteT2DMaterial01 : public NDMaterial
{
public:

	OrthotropicRotatingAngleConcreteT2DMaterial01(int tag,		// nDMaterial tag
		UniaxialMaterial* concreteUniaxialMatObject1,			// concrete uniaxial material objects in a principal direction
		UniaxialMaterial* concreteUniaxialMatObject2,			// concrete uniaxial material objects in the other principal direction 
		double strainAtFcr,										// strain at tension cracking of the concrete
		double strainAtFc,										// strain at the compresion strength of the concrete
		double rho,											    // density
		double damageCte1 = 0.14,								// damage constant 1
		double damageCte2 = 0.6);								// damage constant 2
										
	OrthotropicRotatingAngleConcreteT2DMaterial01();

	~OrthotropicRotatingAngleConcreteT2DMaterial01();

	double getRho(void);

	int setTrialStrain(const Vector& v);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrianIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);

	Response* setResponse(const char** argv, int argc, OPS_Stream& theOutputStream);
	int getResponse(int responseID, Information& matInformation);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);
	NDMaterial* getCopy(const char* type);

	void Print(OPS_Stream& s, int flag = 0);
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	const char* getType(void) const { return "PlaneStress"; };
	int getOrder(void) const { return 3; };

	const Vector& getStress(void);
	const Vector& getStrain(void);

protected:

private:

	int setTrialStrainPrincipalDirection(const Vector& v);									// calculate trial stress and tangent
	void calculateStrainPrincipalDirections01(void);										// calculate the principal direction for the strains (11, 22, 12) using the calculateAngle01 method
	void calculateAngle01(double cosTheta, double sinTheta, double& theta);					// calculate the theta angle [-pi,pi] from the cos(theta) and sin(theta)
	void calculatePoissonRatios(double e1, double e2);										// calculate the Vecchio Poisson Ratios
	void calculateStrainTransformationMatrix(double* pTmatrixStrain, double theta);			// calculate the Strain Transformation Matrix that goes from the Local Coord System to the orientation of the Principal Direction of strain
	void calculateStressTransformationMatrix(double* pTmatrixStress, double theta);			// calculate the Stress Transformation Matrix that goes from the orientation of the Principal Direction T(-thetaPD) to the Local Coord System

	// Function added for MEFI3D
	Vector getEc(void);		                // return input parameters

	// Pointers to material arrays
	UniaxialMaterial** theMaterial;         // pointer of the materials

	Vector poissonRatios;					// store the Poisson Ratio (nu12 and nu21) of the Orthotropic Concrete Layer
	double thetaPrincipalDirection;			// store the orientation of the Principal Direction in the material
	Vector strainPrincipalDirection;		// store the principal strains
	double ecr;								// store the strain at the tension cracking of the concrete
	double ec;								// store the strain at the compresion strength of the concrete
	bool isConcreteCracked;					// store a flag indicating if the concrete has cracking
	double rhoNDM;							// density

	// Committed Sate Variables........................................................
	Vector Cstrain;							// store the committed strain vector in the orthotropic concrete layer
	Vector CMaxMinStrainRec;				// store the maximun and minimun strain recorded in the concrete layer

	// Trial State Variables...........................................................
	Vector strain_vec;					    // store the strain vector at the current trial in the Orthotropic Concrete Layer
	Vector stress_vec;						// store the stress vector in the Orthotropic Concrete Layer
	Matrix tangent_matrix;					// store the tangent matrix in the Orthotropic Concrete Layer
	Matrix initialTangentNDM;				// store the initial tangent matrix in the Orthotropic Concrete Layer
	Vector TMaxMinStrainRec;				// store the maximun and minimun strain recorded in the concrete layer
	double Eo;								// store the concrete Young's modulus
	double Gmin;							// store the minimun shear tangent in the concrete layer 
	double damageConstant1;					// store the constant to determinate the damage factor
	double damageConstant2;					// store the constant to determinate the damage factor

	const double pi;
};

#endif 
