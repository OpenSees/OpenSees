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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/EigenAnalysis.h,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/eigenAnalysis/EigenAnalysis.h
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenAnalysis.
// EigenAnalysis is a subclass of Analysis, it is used to perform the 
// eigen value analysis on the FE_Model.
//
// This class is inheritanted from the base class of Analysis
// which was created by fmk (Frank).


#ifndef EigenAnalysis_h
#define EigenAnalysis_h

#include <Analysis.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class EigenAlgorithm;
class EigenIntegrator; 
class EigenSOE;

class EigenAnalysis : public Analysis
{
  public:
     EigenAnalysis(Domain &theDomain,
		   ConstraintHandler &theHandler,
		   DOF_Numberer &theNumberer,
		   AnalysisModel &theModel,
		   EigenAlgorithm &theAlgo,
		   EigenSOE &theSOE,
		   EigenIntegrator &theIntegrator);
     
     virtual ~EigenAnalysis();
     
     virtual int analyze(int numModes);
     void clearAll(void);	         
     virtual int domainChanged();
     
     virtual int setAlgorithm(EigenAlgorithm &theAlgo);
     virtual int setIntegrator(EigenIntegrator &theIntegrator);
     virtual int setEigenSOE(EigenSOE &theSOE);
  
  protected:
     ConstraintHandler	*getConstraintHandlerPtr() const;
     DOF_Numberer	*getDOF_NumbererPtr() const;
     AnalysisModel	*getAnalysisModelPtr() const;
     EigenAlgorithm	*getEigenAlgorithm() const;
     EigenSOE 		*getEigenSOE() const;
     EigenIntegrator	*getEigenIntegrator() const;
  
  private:
     ConstraintHandler 	*theConstraintHandler;
     DOF_Numberer	*theDOF_Numberer;
     AnalysisModel 	*theAnalysisModel;
     EigenAlgorithm	*theAlgorithm;
     EigenSOE		*theSOE;
     EigenIntegrator	*theIntegrator;
     int domainStamp;
}; 

#endif

