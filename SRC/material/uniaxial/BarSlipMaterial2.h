#ifndef BarSlipMaterial2_h
#define BarSlipMaterial2_h

#include <UniaxialMaterial.h>
#include <Vector.h>
#include <Matrix.h>

class BarSlipMaterial2 : public UniaxialMaterial
{
public:
	// Constructors
	BarSlipMaterial2(int tag, UniaxialMaterial** aMats, UniaxialMaterial** bMats, double ab, double cb, double l, int n, double Lc = 0.0, int iter = 50, double tol = 1.0e-4, bool dbl = false);
	BarSlipMaterial2();

	// Destructor
	~BarSlipMaterial2();

	// Method to Get Class Type
	const char* getClassType(void) const { return "BarSlipMaterial2"; };

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
	double Ab;			// bar area
	double Cb;			// bar circumference
	double L;			// embedded bar length
	double lambda;		// integration constant
	double dx;			// constant spacing of IPs
	double lc;			// characteristic length
	int N;				// number of IPs
	int maxIterNo;		// maximum iteration number
	double maxTol;		// maximum acceptable tolerance
	bool doubleSided;	// true if there is bar slip at both sides of joint
	Matrix H_inv;		// inverse of [H] matrix for slip nonlocality relation

	// State Variables
	double s_L_T;		// trial slip at bar end (at joint surface; material model's strain input!)
	double s_L_C;		// committed slip at bar end (at joint surface; material model's strain input!)
	double tau_L;		// bond stress at end IP (at joint surface)
	double sig_L_T;		// trial axial stress at end IP (model output)
	double sig_L_C;		// committed axial stress at end IP
	double E_tan_T;		// tangent
	double E_tan_C;		// committed tangent

	double* s_u_T;		// trial unknown slips (at all but end IPs)
	double* s_u_C;		// committed unknown slips (at all but end IPs)
	double* s_ub;		// trial unknown nonlocal slips (at all but end IPs)
	double* eps_T;		// trial axial strains (at all IPs)
	double* eps_C;		// committed axial strains (at all IPs)
	double* tau_u;		// trial bond stress corresponding to unknown slip values
	double* sig;		// trial axial stresses (at all IPs)

	Matrix J_T;			// trial Jacobian matrix
	Matrix J_o;			// initial Jacobian matrix
	Matrix J_C;			// committed Jacobian matrix

	// Material Model Pointers
	UniaxialMaterial** axialMats;	// an array of pointers to the UniaxialMaterials representing bar's axial stress vs. axial strain response
	UniaxialMaterial** bondMats;	// an array of pointers to the UniaxialMaterials representing bar's bond stress vs. slip response

	// Private Methods
	Vector getAxialStrains(void);
	Vector getSlips(void);
	Vector getAxialStresses(void);
	Vector getBondStresses(void);
};


#endif

