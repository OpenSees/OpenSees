/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 2001, The Regents of the University of California    **
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
** Reliability module developed by:                                   **
**   Terje Haukaas (haukaas@ce.berkeley.edu)                          **
**   Armen Der Kiureghian (adk@ce.berkeley.edu)                       **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/SelectLoadInitialStaticAnalysis.h,v $
                                                                        
#ifndef SelectLoadStaticAnalysis_h
#define SelectLoadStaticAnalysis_h

#include <fstream>
#include <iomanip>
#include <iostream>

using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;

#include<TaggedObjectStorage.h>
#include<ArrayOfTaggedObjects.h>
#include<StaticAnalysis.h>
#include<AnalysisModel.h>
#include<EquiSolnAlgo.h>
#include<ConstraintHandler.h>
#include<DOF_Numberer.h>
#include<LinearSOE.h>
#include<SensitivityAlgorithm.h>
#include<SensitivityIntegrator.h>//Abbas
#include<StaticIntegrator.h>
#include<ConvergenceTest.h>
#include<CTestNormUnbalance.h>
#include<CTestNormDispIncr.h>
#include<NewtonRaphson.h>
#include<PlainHandler.h>
#include<RCM.h>
#include<LoadControl.h>
#include<ProfileSPDLinDirectSolver.h>
#include<ProfileSPDLinSolver.h>
#include<ProfileSPDLinSOE.h>
#include<BandGenLinLapackSolver.h>
#include<BandGenLinSOE.h>
#include<ReliabilityStaticAnalysis.h>
#include<NewStaticSensitivityIntegrator.h>//Abbas
#include<NewSensitivityAlgorithm.h>//Abbas
#include <InitialStaticAnalysis.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Node.h>
#include <NodeIter.h>
#include <RandomVariablePositioner.h>
#include <RandomVariable.h>
#include<Integrator.h>



class SelectLoadInitialStaticAnalysis : public InitialStaticAnalysis
{
public:
	SelectLoadInitialStaticAnalysis(ReliabilityDomain* theReliabilityDomain,
							 Domain* theStructuralDomain,
							 int nstep,
							 int numLoadPatterns,
							 int* StaticLoadPatterns,
							 bool print);
	~SelectLoadInitialStaticAnalysis();
	void activateSensitivity(void);
	void inactivateSensitivity();
	void analyze(Vector x);
	void analyzeMean(void);
	void recoverLoads();
	void constLoads(double);
	void resetconstLoads(double);
	void constandrecoverLoads(double);

protected:

private:
	int Nstep;
	int NumLoadPatterns;
	int* StaticLoadPatterns;
	LoadPatternIter* theOrgPatternIter;
	TaggedObjectStorage* theOrgPatterns;
	void saveLoads();
	void createStaticAnalysis(void);
	void reset();
	void modifyLoads();
	bool modified;

	ReliabilityStaticAnalysis *theReliabilityStaticAnalysis;
	AnalysisModel *theAnalysisModel;
	EquiSolnAlgo *theAlgorithm;
	ConstraintHandler *theHandler;
	DOF_Numberer *theNumberer;
	LinearSOE *theSOE;
	StaticIntegrator *theStaticIntegrator;
	ConvergenceTest *theTest;
	
	bool activeSensitivity;
//	SensitivityAlgorithm *theSensitivityAlgorithm;
           Integrator        *theSensitivityAlgorithm;
	SensitivityIntegrator *theSensitivityIntegrator;//Abbas
};

#endif
