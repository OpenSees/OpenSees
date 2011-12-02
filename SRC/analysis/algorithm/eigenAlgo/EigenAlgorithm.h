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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/EigenAlgorithm.h,v $
                                                                        
                                                                        
// File: ~/analysis/algorithm/eigenAlgo/EigenAlgorithm.h
//
// Written: Jun Peng
// Created: Wed Jan 27, 1999
// Revision: A
//
// Description: This file contains the class definition of EigenAlgorithm.
// EigenAlgorithm is a class which performs a eigen solution algorithm
// to solve the equations. 
//
// This class is inheritanted from the base class of SolutionAlgorithm
// which was created by fmk (Frank).


#ifndef EigenAlgorithm_h
#define EigenAlgorithm_h

#include <SolutionAlgorithm.h>

class AnalysisModel;
class EigenSOE;
class EigenIntegrator;

class EigenAlgorithm : public SolutionAlgorithm
{
  public:
     EigenAlgorithm(int classTag);
     virtual ~EigenAlgorithm();
     
     // public functions defined for subclasses
     virtual void setLinks(AnalysisModel &theModel,
		   EigenIntegrator &theIntegrator,
		   EigenSOE &theSOE);
     
     // pure virtural functions
     virtual int solveCurrentStep(int numModes) = 0;
     virtual void Print(ostream &s, int flag=0) = 0;
     
     virtual AnalysisModel *getAnalysisModelPtr() const;
     virtual EigenIntegrator *getEigenIntegratorPtr() const;
     virtual EigenSOE *getEigenSOEptr() const;
 
  protected:
  
  private:
     AnalysisModel  *theModel;
     EigenIntegrator *theIntegrator;
     EigenSOE *theSOE;
     
};

#endif
