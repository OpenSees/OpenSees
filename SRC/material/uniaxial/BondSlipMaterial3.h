/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for BondSlip material model used in model BarSlip2
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#ifndef BondSlipMaterial3_h
#define BondSlipMaterial3_h

// Directives
#include <UniaxialMaterial.h>
#include <Channel.h>
#include <float.h>

class BondSlipMaterial3 : public UniaxialMaterial
{
public:
	// Constructors
	BondSlipMaterial3(int tag, double so, double s1, double to, double t1, double ko);
	BondSlipMaterial3();

	// Destructor
	~BondSlipMaterial3();

	// Method to Get Class Type
	const char *getClassType(void) const { return "BondSlipMaterial3"; };

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
	UniaxialMaterial *getCopy(void);

	// Methods to Do Parallel Processing
	int sendSelf(int commitTag, Channel &theChannel);
	int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

	// Methods to Obtain Information
	void Print(OPS_Stream &s, int flag = 0);

protected:

private:
	// Input Variables
	double s_o;			// slip corresponding to maximum bond stress under monotonic loading
	double s_1;			// clear spacing between bar ribs
	double tau_o;		// initial maximum bond stress
	double tau_1;		// initial residual bond stress

	// Constants
	double k_init;		// initial stiffness
	double k_o;			// unloading/reloading stiffness for bearing and friction components
	double PI;			// pi value
	double E_o;			// monotonic energy dissipation
	double E_of;		// monotonic friction energy

	// State Variables
	double s_commit;		// committed slip
	double s_trial;			// trial slip

	double tau_commit;		// committed bond stress
	double tau_trial;		// trial bond stress

	double k_tan;			// trial tangent stiffness

	double tau_or_pos;		// reduced maximum bond stress in positive direction
	double tau_1r_pos;		// reduced positive residual resistance
	double tau_f_pos;		// friction resistance after first reversal (updated if maximum slip changes)
	double tau_fr_pos;		// reduced positive friction resistance

	double tau_or_neg;		// reduced maximum bond stress in negative direction
	double tau_1r_neg;		// reduced negative residual resistance
	double tau_f_neg;		// friction resistance after first reversal (updated if maximum slip changes)
	double tau_fr_neg;		// reduced negative friction resistance
	
	double s_max;			// maximum reached slip
	double s_max_commit;	// committed maximum reached slip
	double s_min;			// minimum reached slip
	double s_min_commit;	// committed minimum reached slip

	double s_max_env;		// maximum reached slip used to determine envelope
	double s_min_env;		// minimum reached slip used to determine envelope

	double E_d;				// total dissipated energy
	double E_f_pos;			// total friction energy for computing friction damage factor in positive direction
	double E_f_neg;			// total friction energy for computing friction damage factor in negative direction
};

#endif