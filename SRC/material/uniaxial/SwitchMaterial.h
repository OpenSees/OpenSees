/* Written by: Mohammad Salehi (mohammad.salehi@rice.edu)
** Created: 2020
** Description: The source code for Switch material model used to simulate repaired RC members
**
**
** Reference:
**
** Mohammad Salehi, Petros Sideris, and Reginald DesRoches (2022)
** “Numerical modeling of repaired reinforced concrete bridge columns”
** Engineering Structures, 253: 113801
*/

#ifndef SwitchMaterial_h
#define SwitchMaterial_h

#include <UniaxialMaterial.h>

class SwitchMaterial : public UniaxialMaterial
{
public:
	SwitchMaterial(int tag, UniaxialMaterial** theMaterialModels, double ts, double dts, int st);
	SwitchMaterial();
	~SwitchMaterial();

	const char* getClassType(void) const { return "SwitchMaterial"; };

	int setTrialStrain(double strain, double strainRate = 0.0);
	double getStrain(void);
	double getStrainRate(void);
	double getStress(void);
	double getTangent(void);
	double getDampTangent(void);
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
	double trialStrain;
	double trialStrainRate;

	double t_s;				// switch start time
	double dt_s;			// switch transition time
	double strain_s;		// new material model's strain shift
	double strain_res;		// strain at which first model's stress was compressive last time (for modification mode only)

	int sType;		// switch type: 1 for replacement and 2 for modification

	bool switched;	// true if mateiral switch has been initiated

	UniaxialMaterial** theModels; // an array of pointers to the UniaxialMaterials
};


#endif

