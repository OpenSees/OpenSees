// Code written/implemented by: Fabian Rojas B.
//								Maria Jose Nunez
//
// Created: 07/2023
// 
// Description: This file contains the ReinforcedConcreteLayerMembraneSection02 class definition
// A ReinforcedConcreteLayerMembraneSection02 is a subclass of the sectionForceDeformation class and corresponds to the abstract representation
// for the stress-strain behavior for a Reinforced Concrete Layer Membrane Element in the Finite Element Method or Structural Analysis. 
//
// Reference:
// 1. Rojas, F., Anderson, J. C., Massones, L. M. (2016). A nonlinear quadrilateral layered membrane with drilling degrees of freedom for 
// the modeling of reinforced concrete walls. Engineering Structures, 124, 521-538.
//
// Source: \OpenSees\SRC\material\section\ReinforcedConcreteLayerMembraneSection
//
// Rev: 1.0

#ifndef ReinforcedConcreteLayerMembraneSection02_h
#define ReinforcedConcreteLayerMembraneSection02_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>

class ReinforcedConcreteLayerMembraneSection02 : public SectionForceDeformation {
public:

	ReinforcedConcreteLayerMembraneSection02(int tag,			// section tag
		NDMaterial* RCMaterialObject,							// reinforced concrete nDMaterial tag
		double Thickness);										// section thickness

	ReinforcedConcreteLayerMembraneSection02();

	~ReinforcedConcreteLayerMembraneSection02();

	int setTrialSectionDeformation(const Vector& newTrialSectionStrain);

	// Public methods to obtain strain, stress, tangent and residual information
	const Vector& getSectionDeformation(void);
	const Vector& getStressResultant(void);
	const Matrix& getSectionTangent(void);
	const Matrix& getInitialTangent(void);

	// Public methods to obtain a copy of section and other information
	SectionForceDeformation* getCopy(void);
	const char* getClassType(void) const { return "ReinforcedConcreteLayerMembrane"; };
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

	// Function used by MEFI3D
	double getEcAvg(void);															// Return the average young's modulus of RC material
	Vector getBendingParameters(void);												// Return input parameters

	// Private attributes
	NDMaterial* TheRCMaterial;														// Pointer to a nd material
	
	// Committed State Variables
	Vector CSectionStrain;															// Store the commit strains of the membrane section
	Vector CSectionStress;															// Store the commit stress of the membrane section
	Matrix CSectionTangent;															// Store the commit tangent of the membrane section
	// Trial State Variables
	Vector TSectionStrain;															// Store the trial strains of the membrane section
	Vector TSectionStress;															// Store the trial stress of the membrane section
	Matrix TSectionTangent;															// Store the trial tangent of the membrane section
	Matrix InitialTangent;															// Store the initial tangent of the membrane section

	double t;																		// Store the thickness of the section 

	static ID array;
};

#endif // !ReinforcedConcreteLayerMembraneSection02_h
