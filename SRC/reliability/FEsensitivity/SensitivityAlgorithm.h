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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-02-29 19:47:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityAlgorithm.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SensitivityAlgorithm_h
#define SensitivityAlgorithm_h

class ReliabilityDomain;
class EquiSolnAlgo;
class SensitivityIntegrator;

class SensitivityAlgorithm
{
 public:
  SensitivityAlgorithm(ReliabilityDomain *passedReliabilityDomain,
		       EquiSolnAlgo *passedAlgorithm,
		       SensitivityIntegrator *passedSensitivityIntegrator,
		       int analysisTypeTag);

  ~SensitivityAlgorithm();
  int computeSensitivities(void);
  bool shouldComputeAtEachStep(void);
  int sensitivityDomainChanged(void) {return 0;}
  bool newAlgorithm(void) {return true;};
  
 protected:
  
 private:
    ReliabilityDomain *theReliabilityDomain;
    EquiSolnAlgo *theAlgorithm;
    SensitivityIntegrator *theSensitivityIntegrator;
    int analysisTypeTag; 
};

#endif

