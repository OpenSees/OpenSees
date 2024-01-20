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
// $Revision: 1.0 $
// $Date: 26/07/2023 $
//
// Description: This file contains the class implementation for TzSandCPT material. 
//              Calculates the load-transfer curves that can be used to represent 
//				the axial load-displacement relationship between the pile and the 
//				soil following the unified CPT method for sands (Lehane et al. 2020).
//
// References
//				Lehane, B. M., Li, L., & Bittar, E. J. (2020). Cone penetration test-based load-transfer formulations for driven piles in sand. 
//				Geotechnique Letters, 10(4), 568-574.

#ifndef TzSandCPT_h
#define TzSandCPT_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class TzSandCPT : public UniaxialMaterial
{
public:
	// Full constructor 
	TzSandCPT(int tag, double qc, double Sv, double D, double t, double h, double dz);
	// Null constructor
	TzSandCPT();
	// Destructor
	~TzSandCPT();

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
	double q_c;				    // cone resistance (kPa)
	double sigma_vo_eff;		// vertical effective soil stress (kPa)
	double diameter;			// pile outside diameter (m)
	double wall_thickness;	    // pile wall thickness (m)
	double h_dist;              // distance to the pile tip (m)
	double delta_h;		        // local pile height (m)

private:
	// Function to calculate the ultimate shaft capacity for sands
	void ultimate_capacity(double qc, double Sv, double D, double t, double h);

	// t-z curve parameters
	double tau_f;               // ultimate shaft friction (kPa)
	double w_f;			        // ultimate settlement (m)

	// Vectors to stock t-z curves 
	Matrix data_c;              // compression
	Matrix data_t;              // tension
	int numSlope;               // number of points
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