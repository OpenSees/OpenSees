// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 04/2023
// 
// Description: This file contains the ReinforcedConcreteLayeredMembraneSection class definition
// A ReinforcedConcreteLayeredMembraneSection is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a reinforced concrete layered membrane element in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massone, L. M. (2016). A nonlinear quadrilateral layered membrane element with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\LayeredMembraneSection
//
// Rev: 1.0

#ifndef ReinforcedConcreteLayeredMembraneSection_h
#define ReinforcedConcreteLayeredMembraneSection_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>

class ReinforcedConcreteLayeredMembraneSection : public SectionForceDeformation {
public:

	ReinforcedConcreteLayeredMembraneSection(int tag,			// section tag
		int nSteelLayer,										// number of reinforced steel layers
		int nConcLayer,											// number of concrete layers
		NDMaterial** reinforcedSteelMaterialObjects,			// array of nDMaterial reinforced steel tags for each layer
		NDMaterial** concrete2DMaterialObjects,					// array of nDMaterial concrete tags for each layer
		double* concThickness);									// array of concrete layers thickness

	ReinforcedConcreteLayeredMembraneSection();

	~ReinforcedConcreteLayeredMembraneSection();

	int setTrialSectionDeformation(const Vector& newTrialSectionStrain);

	// Public methods to obtain strain, stress, tangent and residual information
	const Vector& getSectionDeformation(void);
	const Vector& getStressResultant(void);
	const Matrix& getSectionTangent(void);
	const Matrix& getInitialTangent(void);

	// Public methods to obtain a copy of section and other information
	SectionForceDeformation* getCopy(void);
	const char* getClassType(void) const { return "ReinforcedConcreteLayeredMembraneSection"; };
	const ID& getType(void);
	int getOrder(void) const;
	
	// Public methods to set the state of the section
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// Public methods for output
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag);
	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& info);

	// Functions used for recorders
	const Vector& getCommittedStrain(void);
	const Vector& getCommittedStress(void);

	// Function used by MEFI3D
	double getRho(void);															// Return the concrete density

private:
	
	void calculateStrainPrincipalDirections01(void);								// Calculate the principal direction for the strains (11, 22, 12) using the calculateAngle01 method
	void calculateAngle01(double cosTheta, double sinTheta, double& theta);			// Calculate the theta angle [-pi,pi] from the cos(theta) and sin(theta)
	void calculatePoissonRatios(double e1, double e2);								// Calculate the Vecchio Poisson Ratios
	void setCrackPattern(void);														// Obtain the crack pattern for the PlaneElementSection
	
	// Functions used by MEFI3D
	double getEcAvg(void);															// Return the average young's modulus of concrete
	Vector getBendingParameters(void);												// Return bending parameters
	
	// Functions used for recorders
	const Vector& getCrackPattern(void);											// Return the crack pattern
	double getThetaPDAngle(void);													// Return the principal strain direction 
	Vector getSectionStressAvg(void);											    // Return the average section stress

	// Private attributes
	NDMaterial** TheConcrete2DMaterial;												// Array of ND concrete materials
	NDMaterial** TheReinforcedSteel2DMaterial;										// Array of ND reinforced steel materials

	double thetaPrincipalDirection;													// Store the orientation of the Principal Direction in the material
	Vector strainPrincipalDirection;												// Store the principal strains
	Vector poissonRatios;															// Store the Poisson Ratios (nu12 and nu21) 
	double ecr;																		// Store the strain at the tension cracking of the concrete
	double ec;																		// Store the strain at the compresion strength of the concrete
	bool isConcreteCracked;															// Store a flag indicating if the concrete has cracking

	int numberConcreteLayers;														// Store the number of Concrete layers
	int numberReinforcedSteelLayers;												// Store the number of Reinforced Steel layers
	
	// Material History Variables
	Vector crackPattern;															// Store the propierties of the Crack Pattern
	
	// Committed State Variables
	Vector CSectionStrain;															// Store the commit strains of the membrane section
	Vector CSectionStress;															// Store the commit stress of the membrane section
	Matrix CSectionTangent;															// Store the commit tangent of the membrane section
	// Trial State Variables
	Vector TSectionStrain;															// Store the trial strains of the membrane section
	Vector TSectionStress;															// Store the trial stress of the membrane section
	Matrix TSectionTangent;															// Store the trial tangent of the membrane section
	Matrix InitialTangent;															// Store the initial tangent of the membrane section

	// Section array
	double* t;																		// Concrete layer thickness

	// Calculated section parameters
	double h;																		// Store the total thickness of the section 

	const double pi;

	static ID array;
};

#endif // !ReinforcedConcreteLayeredMembraneSection_h
