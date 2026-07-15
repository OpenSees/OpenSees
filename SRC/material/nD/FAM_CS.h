// Code written/implemented by:	Shaohui Zhang
//								Tsinghua University, Beijing, China
//								Xiaodong Ji (jixd@mail.tsinghua.edu.cn)
//								Tsinghua University, Beijing, China
//								Yue Yu
//								Tsinghua University, Beijing, China 
//
// Created: 12/2023
//
// Description: This file contains the definition for
// FAM_CS
// For Detailed explanation of the model, please refer to the references.
//
// References:
// 1) Zhang S., Ji X., Sun L., et al. New OpenSees material model for simulating reinforced concrete shear walls subjected to 
//    coupled axial tension and cyclic lateral loads. Eng Struct 2024;318:118774. https://doi.org/10.1016/j.engstruct.2024.118774.
// 2) Zhang S. Seismic shear force and shear design of reinforced concrete core walls [Ph.D. Thesis].
//    Beijing: Tsinghua University; 2026. (in Chinese)

// modified from FSAM.cpp    author: Kristijan Kolozvari, Kutay Orakcal, John Wallace

#ifndef FAM_CS_h
#define FAM_CS_h

#include <NDMaterial.h>
#include <UniaxialMaterial.h>

#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

class FAM_CS : public NDMaterial
{
public:

	// Constructor
	FAM_CS (int tag,              // nDMaterial tag
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
		double DY,                 // diameter of Steel 2
		double GAMAX,            // Maximum size of aggregate
	    double LM0,             //Claculated average crack spacing
	    double SH);		    //Spacing of horizontal rebar

	// Blank constructor
	FAM_CS();
	
	// Destructor
	~FAM_CS();				  

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
	void InterLocker_improved(double& fpc, double& crack_spacing, double& e_cr_normal_true, double&e_cr_parallel, double&Gc, double&Beta_pmax, double&Tau0_pmax, double&Beta_nmax, double&Tau0_nmax); // shear aggregate interlock
	void Bisection(double& fpc, double& crack_spacing, double& e_cr_normal_true, double& gamma, double& Gc, double& Beta_pmax, double& Tau0_pmax, double& Beta_nmax, double& Tau0_nmax); // Bisection method
	void dowel_action_0(double &gamma, double &Es); // reinforcement dowel action
	void dowel_action(double& fpc, double& dY, double& Es, double& lm, double& crack_angle, double& gamma, double& e_cr_normal_true, double& DI_p, double& DI_n);
	Vector getInputParameters(void); // return input parameters



	// Pointers to material arrays
	UniaxialMaterial **theMaterial; // pointer of the materials 
	Response **theResponses;		// pointer to material responses needed for Concrete

	// Input variables
	double   rho;					// density
	double	 TeTaSt;				// angel of the first steel layer to x coordinate
	double   roux;					// steel ratio of the first steel layer
	double   rouy;					// steel ratio of the second steel layer
	double   E0x;					// steel Young's modulus - X Direction
	double   E0y;					// steel Young's modulus - Y Direction
	double   epcc;					// concrete strain at compressive stress
	double   et;					// monotonic Cracking Strain (et = 0.00008)
	double   fpc;					// concrete compressive strength
	double   Ec;					// concrete Young's modulus 
	double   dY;				     // diameter of Steel 2
	double   Gamax;                  // Maximum size of aggregate
	double   lm0;				    ///Claculated average crack spacing
	double   sh;				    //Spacing of horizontal rebar

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
	double gamma_cr;
	double Tau_Interlock;
	double Tau_Interlock0;
	double dTau_de12_cr;
	double dTau_denormal;
	
	//Maekawa model add by Z
	double Beta;

	// Interlock Stress
	double Tgamma_cr_A;
	double Tgamma_cr_B;
	double Ttau_Interlock_A;
	double Ttau_Interlock_B;

	//Maekawa model add by Z
	double Ttau0_pmaxA;
	double Ttau0_nmaxA;
	double Ttau0_pmaxB;
	double Ttau0_nmaxB;
	double Ctau0_pmaxA;	
	double Ctau0_nmaxA;	
	double Ctau0_pmaxB;	
	double Ctau0_nmaxB;

	// Crack Slip
	double TeA12; 
	double TeB12; 
	
	//Maekawa model add by Z
	double Tbeta_A;
	double Tbeta_pmaxA;	
	double Tbeta_nmaxA;
	double Tbeta_B;
	double Tbeta_pmaxB;
	double Tbeta_nmaxB;
	double Cbeta_pmaxA;
	double Cbeta_nmaxA;
	double Cbeta_pmaxB;
	double Cbeta_nmaxB;

	// for 2nd crack criteria
	double TepsA2; 
	double CepsA2;

	// Dowel Action .....................................
	// "Transfer" variables
	double Tau_Dowel_0;
	double dTau_dgamma_0;
	double Tau_Dowel;
	double dTau_dgamma;
	double dTau_deps;

	//Damage index
	double TDI;
	double TDI_A_p;
	double CDI_A_p;
	double TDI_A_n;
	double CDI_A_n;
	double TDI_B_p;
	double CDI_B_p;
	double TDI_B_n;
	double CDI_B_n;

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
