// by Quan Gu UCSD
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_RESPSURFACESIMULATION_H)
#define AFX_RESPSURFACESIMULATION_H
#include "ReliabilityDomain.h"
#include "RandomNumberGenerator.h"
#include "PrincipalAxis.h"
#include "SurfaceDesign.h"


class RespSurfaceSimulation 
{
public:
	void setNormOfGradientU( double normU);
	double getCov();

	double getFailureProbability();

	RespSurfaceSimulation( 	ReliabilityDomain * pReliabilityDomain,					
		                Vector * pDesignPoint,
						char * pFileName,
						SurfaceDesign * pSurfaceDesign,
						RandomNumberGenerator * pRandomNumberGenerator,
                        double pTargetCOV,
						int pNumberOfSimulations,
						bool isTimeVariant = false,
						double littleDtOverDt = 0.05,
						double pNormOfGradient = 0.0);



	virtual ~RespSurfaceSimulation();

	int runSimulationAnalysis();
private:
	double normOfGradient;
	double littleDtOverDt;
	bool isTimeVariant;
	double theCOV;
	ReliabilityDomain * theReliabilityDomain;

 
	Vector * theDesignPoint;
	char fileName[30];
	SurfaceDesign * theSurfaceDesign;
	RandomNumberGenerator * theRandomNumberGenerator;
	double targetCOV;
	double failureProbability;
	int numOfSimulations;
};

#endif // !defined(AFX_RESPSURFACESIMULATION_H)
