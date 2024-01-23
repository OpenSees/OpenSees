/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** This file contains the class implementation for TzSandCPT          **
** material. This uniaxial material material represents the           **
** shaft-load transfer curve (t-z spring) according to the new Unifi- **
** ed CPT-based method for driven piles in sands. The formulation     ** 
** incorporates the maximum skin friction and end-bearing CPT values  **
** based on the new ISO-19901-4.                                      **
**                                                                    **
**                                                                    **
** Written:                                                           **
**   csasj  (carlos.sastre.jurado@vub.be)                             **
**                                                                    **
**                                                                    **
** References:                                                        **
**   - LEHANE, Barry M., et al. A new'unified'CPT-based axial pile    **
**   capacity design method for driven piles in sand. In: 4th Inter-  **
**   national Symposium on Frontiers in Offshore Geotechnics          **
**   (postponed). 2020. p. 462-477.                                   **
**   - LEHANE, B. M.; LI, Lin; BITTAR, E. J. Cone penetration         **
**   test-based load-transfer formulations for driven piles in sand.  **
**   Geotechnique Letters, 2020, 10.4: 568-574.                       **
** ****************************************************************** */

// $Revision: 1.0    $
// $Date: 23/01/2024 $

#ifndef TzSandCPT_h
#define TzSandCPT_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

#define IFA_DEFAULT   29.        // interface friction angle=29deg

class TzSandCPT : public UniaxialMaterial
{
public:
	// Full constructor 
	TzSandCPT(int tag, double qc, double Sv, double D, double t, double h, double dz, 
	double dcpt, double pa, double d_f = IFA_DEFAULT);
	// Null constructor
	TzSandCPT();
	// Destructor
	~TzSandCPT();

	const char *getClassType(void) const {return "TzSandCPT";};

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

	int sendSelf(int cTag, Channel &theChannel);
	int recvSelf(int cTag, Channel &theChannel,
		FEM_ObjectBroker &theBroker);

	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	/***     Input parameters      ***/
	double q_c;				    // cone resistance 
	double sigma_vo_eff;		// vertical effective soil stress 
	double diameter;			// pile outside diameter 
	double wall_thickness;	    // pile wall thickness 
	double h_dist;              // distance to the pile tip 
	double delta_h;		        // local pile height 
	double d_cpt;               // diameter of the standard CPT probe (35.7mm)
	double p_a;                 // atmospheric pressure
	double delta_f;             // interface friction angle (degrees)

	/***     Axial calculations    ***/
	double tau_f;               // ultimate shaft friction 
	double w_f;			        // peak settlement 
	void ultimate_capacity(
		double qc, double Sv, double D, double t, double h, 
		double dcpt, double pa, double d_f);

	/*** CONVERGED State Variables ***/  
	double cStrain;				// last ced strain
	double cStress;				// last ced stress
	double cTangent;			// last cted tangent

	/*** TRIAL State Variables     ***/
	double tStrain;				// current t strain
	double tStress;				// current t stress
	double tTangent;			// current t tangent

	// Vectors to stock t-z curves 
	Matrix data_c;              // compression
	Matrix data_t;              // tension
	int numSlope;               // number of points
	int tSlope;
};
#endif