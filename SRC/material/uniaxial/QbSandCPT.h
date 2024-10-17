/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** This file contains the class implementation for QbSandCPT          **
** material. This uniaxial material material represents the           **
** shaft-load transfer curve (q-z spring) according to the new Unifi- **
** ed CPT-based method for driven piles in sands. The formulation     **
** incorporates the maximum skin friction and end-bearing CPT values  **
** based on the new ISO-19901-4.                                      **
**                                                                    **
**                                                                    **
** Written:                                                           **
**   csasj (carlos.sastre.jurado@vub.be)                              **
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

#ifndef QbSandCPT_h
#define QbSandCPT_h

#include <UniaxialMaterial.h>
#include <Matrix.h>

class QbSandCPT : public UniaxialMaterial
{
public:
	// Full constructor 
	QbSandCPT(int tag, double qp, double D, double t, double dcpt);
	// Null constructor
	QbSandCPT();
	// Destructor
	~QbSandCPT();

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

private:
	/***     Input parameters      ***/
	double q_p;				    // end bearing capacity mobilised at large displacements at the pile tip (kPa)
	double diameter;			// pile outside diameter (m)
	double wall_thickness;	    // pile wall thickness (m)
	double d_cpt;               // diameter of the standard CPT probe (35.7mm)

	// Function to calculate the ultimate shaft capacity for sands
	void ultimate_capacity(double qp, double D, double t, double dcpt);

	// qb-wb curve parameters
	double q_b01;               // base resistance at 0.1D displacement 

	// Vectors to stock qb-wb curves 
	Matrix data;                // compression
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