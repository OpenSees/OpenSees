#ifndef RemoveMaterial_h
#define RemoveMaterial_h

#include <UniaxialMaterial.h>

class RemoveMaterial : public UniaxialMaterial
{
public:
	RemoveMaterial(int tag, UniaxialMaterial* theMaterialModel, double ts, double dts);
	RemoveMaterial();
	~RemoveMaterial();

	const char* getClassType(void) const { return "RemoveMaterial"; };

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

	double t_s;			// switch start time
	double dt_s;		// switch transition time

	UniaxialMaterial* theModel; // pointer to the UniaxialMaterial
};


#endif

