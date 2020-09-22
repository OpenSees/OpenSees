/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: csasj 
// $Revision: 1.1 $
// $Date: 15/09/2020 $
//
// Description: This file contains the class implementation for MtSandPISA material. 
//				Provide moment distributed spring for sands according to PISA project. 

#ifndef MtSandPISA_h
#define MtSandPISA_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class MtSandPISA : public UniaxialMaterial
{
public:
	// Full constructor 
	MtSandPISA(int tag, double z, double D, double L, double sigma, double Go, double P, double gs, double A0, double Au);
	// Null constructor
	MtSandPISA();
	// Destructor
	~MtSandPISA();

	// OpenSees methods 
	int setTrialStrain(double newy, double yRate);
	double getStrain(void);
	double getStress(void);
	double getTangent(void);

	double getInitialTangent(void);

	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	UniaxialMaterial* getCopy(void);

	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	void Print(OPS_Stream& s, int flag = 0);

protected:
	// Input parameters
	double depth;				// depth below the original seafloor (m)
	double diameter;			// pile outside diameter (m)
	double embedded_length;		// embedded pile length (m)
	double sigma_vo_eff;		// vertical effective soil stress (kPa)
	double gmax;				// small-strain shear modulus (kPa)
	double Pvalue;				// Current value of the local distributed load
	double grid;				// grid spacing (m) 
	double A0;					// Scaling factor for initial stiffness
	double Au;					// Scaling factor for ultimate stress 

private:
	// Functions to get normalised parameters for each soil reaction curve
	void Mt_conic_function(double depth, double diameter, double embedded_length, double A0, double Au);

	// Function to evaluate conic function from PISA project
	double conic_function(double epsilon_ratio, double k0, double n, double Epsilon_u, double Sigma_u);

	// Normalised parameters for conic function
	double epsilon_u;			// Ultimate strain 
	double initial_stiffness;	// Initial stiffness
	double curvature;			// Curvature 
	double sigma_u;				// Ultimate reaction 

	// Vectors to stock p-y curves 
	Matrix data;
	int numSlope;
	int tSlope;

	// Trial values
	double tStrain;				// current t strain
	double tStress;				// current t stress
	double tTangent;			// current t tangent
	// Commited values  
	double cStrain;				// last ced strain
	double cStress;				// last ced stress
	double cTangent;			// last cted  tangent
};
#endif