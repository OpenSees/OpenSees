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
                                                                        
// $Revision: 1.7 $
// $Date: 2008-12-03 23:47:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/reliability/FEsensitivity/SensitivityAlgorithm.h,v $


//
// Written by Terje Haukaas (haukaas@ce.berkeley.edu)
//

#ifndef SensitivityAlgorithm_h
#define SensitivityAlgorithm_h

class Domain;
class ReliabilityDomain;
class EquiSolnAlgo;
class Integrator;
//class SensitivityIntegrator;
class SensitivityAlgorithm
{
 public:

    SensitivityAlgorithm(Domain *passedDomain,
		       EquiSolnAlgo *passedAlgorithm,
		       Integrator *passedSensitivityIntegrator,
                       		      
		       int analysisTypeTag);

  ~SensitivityAlgorithm();
//  int computeSensitivities(void);
//  bool shouldComputeAtEachStep(void);
//  int sensitivityDomainChanged(void) {return 0;}
  // This method needs to go -- MHS
//  bool newAlgorithm(void) {return true;};
  
 protected:
 //int gradNumber;//Abbas......... 
 private:
    Domain *theDomain;
    ReliabilityDomain *theReliabilityDomain;
    EquiSolnAlgo *theAlgorithm;
    //  SensitivityIntegrator *theSensitivityIntegrator;
    int analysisTypeTag; 
};

#endif

