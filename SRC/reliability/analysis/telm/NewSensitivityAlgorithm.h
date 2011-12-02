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
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/analysis/telm/NewSensitivityAlgorithm.h,v $

#ifndef NewSensitivityAlgorithm_h
#define NewSensitivityAlgorithm_h

#include <SensitivityAlgorithm.h>
#include <ReliabilityDomain.h>
#include <EquiSolnAlgo.h>
#include <SensitivityIntegrator.h>
#include <fstream>
using std::ofstream;


class NewSensitivityAlgorithm: public SensitivityAlgorithm
{
  public:
  NewSensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
			  EquiSolnAlgo *passedAlgorithm,
			  SensitivityIntegrator *passedSensitivityIntegrator,
			  int analysisTypeTag);
  ~NewSensitivityAlgorithm();
  //	int computeSensitivities(void);
  bool shouldComputeAtEachStep(void);
  //////////////////////////////////////////////////////////////////
  //// added by K Fujimura /////////////////////////////////////////
  //////////////////////////////////////////////////////////////////
  int computeSensitivities(bool fromFEM=false);
  int sensitivityDomainChanged(void);
  bool newAlgorithm(void);
  
 protected:
  
 private:
  ReliabilityDomain *theReliabilityDomain;
  EquiSolnAlgo *theAlgorithm;
  SensitivityIntegrator *theSensitivityIntegrator;
  int analysisTypeTag; 
  Vector* zeroVector;
  LinearSOE *theSOE;
  IncrementalIntegrator *theIncInt;
  int numGrads;
  int numPos;
  //	int** gradPositioner;
  RandomVariablePositioner** gradRVPositioner;
  ParameterPositioner** gradParaPositioner;
  int** idGradPositioner;
  int* numGradPositioner;
  ofstream output;
};

#endif

