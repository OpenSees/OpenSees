// Code written/implemented by:	Kristijan Kolozvari (kkolozvari@fullerton.edu)
//								California State University, Fullerton 
//								Kutay Orakcal
//								Bogazici University, Istanbul, Turkey
//								John Wallace
//								University of California, Los Angeles
//
// Created: 07/2015
//
// Description: This file contains the class definition for Fixed-Strut Angle
// Model (FSAM, Ulugtekin, 2010; Orakcal et al., 2012) which is a plane-stress
// constitutive model for simulating the behavior of RC panel elements under
// generalized, in-plane, reversed-cyclic loading conditions. The model assumes
// perfect bond assumption between concrete and reinforcing steel bars. The
// reinforcing steel bars develop uniaxial stresses under strains in their
// longitudinal direction, the behavior of concrete is defined using stress-strain
// relationships in biaxial directions, the orientation of which is governed by
// the state of cracking in concrete, and also incorporates biaxial softening
// effects including compression softening and biaxial damage. For transfer of
// shear stresses across the cracks, a friction-based elasto-plastic shear aggregate
// interlock model is adopted, together with a linear elastic model for representing
// dowel action on the reinforcing steel bars (Kolozvari, 2013).
//
// References:
// 1) Orakcal, K., Massone L.M., Ulugtekin, D.,"Constitutive Modeling of Reinforced Concrete
// Panel Behavior under Cyclic Loading", Proceedings of the 15th World Conference on
// Earthquake Engineering, Lisbon, Portugal, 2012.
// 2) Ulugtekin, D., "Analytical Modeling of Reinforced Concrete Panel Elements under
// Reversed Cyclic Loadings", M.S. Thesis, Bogazici University, Istanbul, Turkey, 2010.
// 3) Kolozvari K. (2013). "Analytical Modeling of Cyclic Shear-Flexure Interaction in
// Reinforced Concrete Structural Walls", PhD Dissertation, University of California, Los Angeles.
//
// Source: /usr/local/cvs/OpenSees/SRC/material/nD/reinforcedConcretePlaneStress/FSAM.h
//
// Rev: 1


#ifndef FSAM_h
#define FSAM_h

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class FSAM : public NDMaterial
{
public:

	// Constructor
	FSAM (int tag,              // nDMaterial tag
		double RHO,             // density
		UniaxialMaterial *s1,   // steel X
		UniaxialMaterial *s2,   // steel Y
		UniaxialMaterial *c1,   // concrete 1.1 - uncracked
		UniaxialMaterial *c2,   // concrete 1.2 - uncracked
		UniaxialMaterial *cA1,  // concrete A1 - 1st crack, 2nd crack
		UniaxialMaterial *cA2,  // concrete A2 - 1st crack, 2nd crack
		UniaxialMaterial *cB1,  // concrete B1 - 2nd crack
		UniaxialMaterial *cB2,  // concrete B2 - 2nd crack
		double ROUX,            // Reinforcing ratio of Steel 1
		double ROUY,			// Reinforcing ratio of Steel 2
		double NU,				// Friction coefficient of shear aggregate interlock
		double ALFADOW);		// Stiffness parameter of dowel action

	// Blank constructor
	FSAM();
	
	// Destructor
	~FSAM();				  

	double getRho(void);

	int setTrialStrain(const Vector &v);
	int setTrialStrain(const Vector &v, const Vector &r);
	int setTrialStrainIncr(const Vector &v);
	int setTrialStrainIncr(const Vector &v, const Vector &r);
	const Matrix &getTangent(void);
	const Matrix &getInitialTangent(void);
	
	Response *setResponse (const char **argv, int argc, OPS_Stream &theOutputStream);
	int getResponse (int responseID, Information &matInformation);

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

	// Functions used for recorders
	const Vector &getStress(void);
	const Vector &getStrain(void);

	const Vector &getCommittedStress(void); 
	const Vector &getCommittedStrain(void);  

protected:

private:

	// Determine trial stress and tangent
	int determineTrialStressAndTangent(void);

	// Stages of panel constitutive behavior
	void Stage1(double &ex, double &ey, double &gamma); // Uncracked behavior
	void Stage2(double &ex, double &ey, double &gamma); // 1st Crack Formed
	void Stage3(double &ex, double &ey, double &gamma); // 2nd Crack Formed

	// Additional constitive functions
	void betaf4(double &eo, double &epc, double &fc, double &epscmax); // biaxial damage
	void InterLocker_improved(double &e_cr_normal, double &f_cr_normal, double &e_cr_parallel, double &e_cr_parallel_old, double &epc, double &Ec, double &Tau_Interlock_old); // shear aggregate interlock
	void dowel_action(double &gama, double &Es); // reinforcement dowel action
	Vector getInputParameters(void); // return input parameters

	// Pointers to material arrays
	UniaxialMaterial **theMaterial; // pointer of the materials 
	Response **theResponses;		// pointer to material responses needed for Concrete

	// Input variables
	double   rho;					// density
	double	 TeTaSt;				// angel of the first steel layer to x coordinate
	double   roux;					// steel ratio of the first steel layer
	double   rouy;					// steel ratio of the second steel layer
	double   FyX;					// yield stress of the bare steel bar in X direction
	double   FyY;					// yield stress of the bare steel bar in Y direction
	double   E0x;					// steel Young's modulus - X Direction
	double   E0y;					// steel Young's modulus - Y Direction
	double   epcc;					// concrete strain at compressive stress
	double   et;					// monotonic Cracking Strain (et = 0.00008)
	double   fpc;					// concrete compressive strength
	double   Ec;					// concrete Young's modulus
	double   nu;					// friction Coefficient of Shear Aggregate Interlock 
	double   alfadow;				// stiffness Coefficient of Dowel Action 
	Vector ConcreteInput;			// vector for storing concrete input variables

	// Biaxial Damage Parameters ......................
	// Variables for getting value from betaf4 function
	double beta;
	double delbeta;
	double epsiloncmax; 

	double Tepscmax1;
	double Tepscmax2;
	double Cepscmax1;
	double Cepscmax2;

	double TepscmaxA1;
	double TepscmaxA2;
	double CepscmaxA1;
	double CepscmaxA2;

	double TepscmaxB1;
	double TepscmaxB2;
	double CepscmaxB1;
	double CepscmaxB2;

	// Shear Aggregate Interlock ........................
	// "Transfer" variables
	double Tau_Interlock;
	double dTau_de12;
	double dTau_dfcnormal;

	// Interlock Stress
	double Ttau_Interlock_A;
	double Ttau_Interlock_B;
	double Ctau_Interlock_A; 
	double Ctau_Interlock_B; 

	// Crack Slip
	double TeA12; 
	double TeB12; 
	double CeA12; 
	double CeB12; 

	// for 2nd crack criteria
	double TepsA2; 
	double CepsA2;

	// Dowel Action .....................................
	// "Transfer" variables
	double Tau_Dowel;
	double dTau_dgamma;

	// State Variables ..................................
	int crackA;
	int crackB;
	
	// Current ..........................................
	Vector strain_vec; // Strain
	Vector stress_vec; // Stress
	Matrix tangent_matrix; // Tangent

	// Strut Orientation ................................
	double   alpha_strain;			// principle strain direction
	double   alfa_crackA;			// direction of 1st strut
	double   alfa_crackB;			// direction of 2nd strut

	// Principal Strains ................................
	double	 Tprstrain1;			// principal strain 1 (uncracked)
	double   Tprstrain2;			// principal strain 2 (uncracked)
	double	 Cprstrain1;			// principal strain 1 (uncracked)
	double   Cprstrain2;			// principal strain 2 (uncracked)

	// Committed .........................................
	// Panel Strains and Stresses
	Vector   CStress;  // Committed stesses
	Vector   CStrain;  // Committed strain

	const double pi; 

	// For recorders ....................................
	// Variables
	Vector TPanelConcStress;
	Vector TPanelSteelStress;
	Vector TStrainStressSteel1;
	Vector TStrainStressSteel2;
	Vector TStrainStressConc1;
	Vector TStrainStressConc2;
	Vector TStrainStressInterlock1; 
	Vector TStrainStressInterlock2;
	
	Vector CPanelConcStress;
	Vector CPanelSteelStress;
	Vector CStrainStressSteel1;
	Vector CStrainStressSteel2;
	Vector CStrainStressConc1;
	Vector CStrainStressConc2;
	Vector CStrainStressInterlock1; 
	Vector CStrainStressInterlock2;
	Vector CCrackingAngles;
	
	// Functions
	Vector getPanelStressConcrete(void);
	Vector getPanelStressSteel(void);
	Vector getStrainStressSteel1(void);
	Vector getStrainStressSteel2(void);
	Vector getStrainStressConcrete1(void);
	Vector getStrainStressConcrete2(void);
	Vector getStrainStressInterlock1(void);
	Vector getStrainStressInterlock2(void);
	Vector getCrackingAngles(void);

};

#endif
