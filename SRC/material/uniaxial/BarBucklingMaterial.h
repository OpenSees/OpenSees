/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for a mechanical bar-buckling model
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#ifndef BarBucklingMaterial_h
#define BarBucklingMaterial_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class BarBucklingMaterial : public UniaxialMaterial
{
public:
	// Constructors
	BarBucklingMaterial(int tag, UniaxialMaterial* mat, double d, double sl, double imp = 0.001, double eps = 0.0, int iter = 50, double tol = 1.0e-4);
	BarBucklingMaterial();

	// Destructor
	~BarBucklingMaterial();

	// Method to Get Class Type
	const char* getClassType(void) const { return "BarBucklingMaterial"; };

	// Methods to Set the State of Material
	int setTrialStrain(double strain, double strainRate = 0.0);
	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	// Methods to Get Response Variables
	double getStrain(void);
	double getStress(void);
	double getTangent(void);
	double getInitialTangent(void);

	// Method to Copy Material
	UniaxialMaterial* getCopy(void);

	// Methods to Do Parallel Processing
	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

	// Methods to Obtain Information
	Response* setResponse(const char** argv, int argc,
		OPS_Stream& theOutputStream);
	int getResponse(int responseID, Information& matInformation);
	void Print(OPS_Stream& s, int flag = 0);

protected:

private:
	// Constant Variables
	double LD;			// bar slenderness ratio (unbraced length to diameter)
	double D;			// bar diameter
	double L ;			// bar length
	double IMP;			// initial imperfection
	double barA;		// bar area
	double minBcklEps;	// minimimum strain before buckling can occur

	int maxIterNo;		// maximum iteration number
	double maxTol;		// maximum acceptable tolerance

	double secX[10];		// array of normalized locations of IPs along bar length (x/L)
	double secL[10];		// array of weights of IPs along bar length
	double secY[10];		// array of locations of IPs over bar cross section
	double secA[10];		// array of areas of IPs over bar cross section

	Matrix J_o;			// initial Jacobian matrix

	const int N = 10;		// number of IPs along bar length and over cross section area
	const double PI = 3.14159265359;

	// State Variables
	double eps_T;		// trial average strain
	double eps_C;		// committed average strain
	double sig_T;		// trial average stress
	double sig_C;		// committed average stress

	double eps_min;		// minimum reached strain

	double eps_eT;		// trial axial strain at bar end
	double eps_eC;		// committed axial strain at bar end
	double eps_qT;		// trial axial strain at bar quarter-length
	double eps_qC;		// committed axial strain at bar quarter-length
	double phi_eT;		// trial curvature at bar end
	double phi_eC;		// committed curvature at bar end

	double N_e;			// end section axial force
	double N_q;			// quarter-length section axial force
	double M_e;			// end section moment

	double w_m;
	double E_tanT;
	double E_tanC;

	// Material Model Pointers
	UniaxialMaterial** mats;	// an array of pointers to the UniaxialMaterials representing bar's axial stress vs. axial strain response

	// Private Methods
	int updateSectionStrains(void);		// update trial section strains
	void updateSectionStresses(void);	// get trial section stresses
	void getSectionTangents(double& N2eps_e, double& N2phi_e, double& N2eps_q, double& M2eps_e, double& M2phi_e);		// get trial section tangents
};


#endif

