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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/StatRandomVibrationSimulation.h,v $


#ifndef StatRandomVibrationSimulation_h
#define StatRandomVibrationSimulation_h

#include <RandomVibrationSimulation.h>
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

class StatRandomVibrationSimulation : public RandomVibrationSimulation

{
  public:
	StatRandomVibrationSimulation(ReliabilityDomain* passedReliabilityDomain,
								  Domain* passedDomain,
						          FunctionEvaluator* passedGFunEvaluator,
							      ProbabilityTransformation* passedTransformation,
						          double passedStartTime,
						          double passedEndTime,
						          double passedTimeInterval,
						          double passedFragMin,
						          double passedFragInt,
						          int passednFrag,
								  double passedstationarytime,
								  bool passedtwoside,
								  bool passedsystem,
						          int passedmaxSim,
						          int passedcheckinterval,
						          double passedeps,
						          int passedinstantaneous,
						          int passedfirstpassage,
	  				              TCL_Char *passedFileName,
						          char* passedFileBinary,
	                              Tcl_Interp *passedTclInterp,
								  bool passedprint,
								  double passedsampleAmp,
								  double passedsampleTime);
//								  Vector *pStartPoint);

	~StatRandomVibrationSimulation();
	void crudeInstantaneousSimulation();
	void crudeFisrtpassageSimulation();
	void samplingInstantaneousSimulation();
	void samplingFisrtpassageSimulation();

  protected: 
    
  private:
	int numTimePoints;
	Vector* timepoints;
	int* anaSteps;
	double stationaryTime;
	int stationaryStep;
	bool print;
	bool sample;
	int ifsample;
	double sampleAmp;
	int sampleStep;
	ofstream output;
	ofstream outputFile;
	int checkTimePoints();
	Vector* samplExc;
	Vector* samplResp;
	ofstream samplRec;
	ifstream samplRead;
	int nsample;
	Matrix* respMat;
	Matrix* excMat;
//	Vector *startPoint;
	int sampletofile(int nstep, ofstream& record);
	long samplefromfile(int nstep,ifstream& record);

};

#endif



