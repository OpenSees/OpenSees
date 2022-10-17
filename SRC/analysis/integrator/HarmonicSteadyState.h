
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

#ifndef HARMONICSTEADYSTATE_H
#define HARMONICSTEADYSTATE_H

// Written: Seweryn Kokot, Opole University of Technology, Poland
// Created: 2021
//
// based on LoadControl.cpp
// written: fmk
//
// Description: This file contains the class definition for HarmonicSteadyState.
// HarmonicSteadyState is an algorithmic class for performing a quasi-static harmonic
// steady-state analysis.

#include "StaticIntegrator.h"
#include <classTags.h>

class LinearSOE;
class AnalysisModel;
class FE_Element;
class Vector;
class EquiSolnAlgo;
class ReliabilityDomain;

class HarmonicSteadyState : public StaticIntegrator
{
  public:
    HarmonicSteadyState(double deltaLambda, double loadPeriod, int numIncr,
			   double minLambda, double maxlambda,
			   int classtag=INTEGRATOR_TAGS_HarmonicSteadyState);

    ~HarmonicSteadyState();

    int newStep(void);
    int update(const Vector &deltaU);
    int setDeltaLambda(double newDeltaLambda);

    // Public methods for Output
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
			 FEM_ObjectBroker &theBroker);

    void Print(OPS_Stream &s, int flag =0);

    int formEleTangent(FE_Element *theEle);
    int formEleResidual(FE_Element *theEle);

    // Adding sensitivity
    int formSensitivityRHS(int gradNum);
    int formIndependentSensitivityRHS();
    int saveSensitivity(const Vector &v, int gradNum, int numGrads);
    int commitSensitivity(int gradNum, int numGrads);
    int computeSensitivities(void);//Abbas
    bool computeSensitivityAtEachIteration();


    ///////////////////////

protected:

  private:
    double deltaLambda;  // dlambda at step (i-1)
    double loadPeriod; // load period in seconds (p = 2pi/T)

    double specNumIncrStep, numIncrLastStep; // Jd & J(i-1)
    double dLambdaMin, dLambdaMax; // min & max values for dlambda at step (i)

    // Adding sensitivity
    int gradNumber;
    int sensitivityFlag;
   // EquiSolnAlgo *theAlgorithm;
    ReliabilityDomain *theDomain;
    ////////////////////
};

#endif
