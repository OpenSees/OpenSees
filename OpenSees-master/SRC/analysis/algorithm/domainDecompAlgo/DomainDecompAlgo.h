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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/domainDecompAlgo/DomainDecompAlgo.h,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/DomainDecompAlgo.h 
// 
// Written: fmk 
// Created: 10/96 
// Revision: A 
//

// Description: This file contains the class definition for 
// DomainDecompAlgo. DomainDecompAlgo is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses deifine
// the sequence of operations to be performed in the analysis by static
// equilibrium of a finite element model.  
// 
// What: "@(#)DomainDecompAlgo.h, revA"

#ifndef DomainDecompAlgo_h
#define DomainDecompAlgo_h

#include <SolutionAlgorithm.h>

class AnalysisModel;
class IncrementalIntegrator;
class LinearSOE;
class DomainSolver;
class DomainDecompAnalysis;
class Subdomain;

class DomainDecompAlgo: public SolutionAlgorithm
{
  public:
    DomainDecompAlgo();    
    ~DomainDecompAlgo();

    // public functions defined for subclasses
    int solveCurrentStep(void);

    void setLinks(AnalysisModel         &theModel, 
		  IncrementalIntegrator &theIntegrator,
		  LinearSOE             &theSOE,
		  DomainSolver          &theSolver,		     
		  Subdomain             &theSubdomain);	

    int sendSelf(int cTag, Channel &theChannel);
    int recvSelf(int cTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker);
    
  protected:
    
  private:
    AnalysisModel	  *theModel;
    IncrementalIntegrator *theIntegrator;    
    LinearSOE             *theLinearSOE;    
    DomainSolver          *theSolver;    
    Subdomain             *theSubdomain;
};

#endif


