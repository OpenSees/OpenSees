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
// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/RandomVibrationSimulation.h,v $

#ifndef RandomVibrationSimulation_h
#define RandomVibrationSimulation_h

#include <ReliabilityDomain.h>
#include <Domain.h>
#include <FunctionEvaluator.h>
#include <RandomProcess.h>
#include <RandomNumberGenerator.h>
#include <GeneralRandGenerator.h>
#include <ProbabilityTransformation.h>
#include <tcl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ofstream;
using std::ios;
using std::setw;
using std::setprecision;
using std::setiosflags;
using std::scientific;
using std::showpoint;

class RandomVibrationSimulation
{
 public:
  RandomVibrationSimulation(ReliabilityDomain* passedReliabilityDomain,
			    Domain* passedDomain,
			    FunctionEvaluator* passedGFunEvaluator,
			    ProbabilityTransformation* passedTransformation,
			    double passedFragMin,
			    double passedFragInt,
			    int passednFrag,
			    bool passedtwoside,
			    bool passedsystem,
			    int passedmaxSim,
			    int passedcheckinterval,
			    double passedeps,
			    int passedinstantaneous,
			    int passedfirstpassage,
			    TCL_Char *passedFileName,
			    char* passedFileBinary,
			    Tcl_Interp *passedTclInterp);
  
  virtual ~RandomVibrationSimulation();
  virtual void analyze();
  virtual void generateRV();
  virtual void crudeInstantaneousSimulation()=0;
  virtual void crudeFisrtpassageSimulation()=0;
  virtual void samplingInstantaneousSimulation()=0;
  virtual void samplingFisrtpassageSimulation()=0;
  
  protected: 

	ReliabilityDomain* theReliabilityDomain;
	Domain* theDomain;
	FunctionEvaluator* theGFunEvaluator;
	RandomProcess* theRandomProcess;
	RandomNumberGenerator* theRandomNumberGenerator;
	ProbabilityTransformation* theTransformation;
	int NumTotalPulse;
	double delta_Pulse;
	double delta;
	bool twoside;
	bool system;

	int MaxSim;
	double eps;
	int checkinterval;
	int instantaneous;
	int firstpassage;
	int	numFragility;
	Vector* Fragility;
	int numRV;
	int numRVPos;
	int numLsf;
	char* fileName;
	char* fileBinary; 
    Tcl_Interp* theTclInterp;
	Vector* uRV;
	Vector* xRV;

  private:
};

#endif
