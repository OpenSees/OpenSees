// $Revision: 1.1 $
// $Date: 2008-02-29 19:43:53 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NonStatRandomVibrationSimulation.h,v $

#ifndef NonStatRandomVibrationSimulation_h
#define NonStatRandomVibrationSimulation_h

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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NonStatRandomVibrationSimulation.h,v $
                                                                        
// Written: fmk 
// Created: 11/96
// Revision: A
//
// Purpose: This file contains the class interface for Node.
// A Node provides the abstraction for a node in the FEM.
// Nodes have original position, trial displacement, velocity and 
// acceleration and committed displacement, velocity and acceleration 
// (the last committed() trial quantities).
//
// What: "@(#) Node.h, revA"
class NonStatRandomVibrationSimulation : public RandomVibrationSimulation

{
 public:
  NonStatRandomVibrationSimulation(ReliabilityDomain* passedReliabilityDomain,
				   Domain* passedDomain,
				   FunctionEvaluator* passedGFunEvaluator,
				   ProbabilityTransformation* passedTransformation,
				   double passedStartTime,
				   double passedEndTime,
				   double passedTimeInterval,
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
				   Tcl_Interp *passedTclInterp,
				   bool passedprint);
  
  virtual ~NonStatRandomVibrationSimulation();
  void crudeInstantaneousSimulation();
  void crudeFisrtpassageSimulation();
  void samplingInstantaneousSimulation();
  void samplingFisrtpassageSimulation();

  protected: 
    
  private:
  int numTimePoints;
  Vector* timepoints;
  int* anaSteps;
  bool print;
  ofstream output;
  ofstream outputFile;
  int checkTimePoints();
};

#endif



