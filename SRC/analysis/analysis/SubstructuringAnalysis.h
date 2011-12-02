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
                                                                        
// $Revision: 1.3 $
// $Date: 2007-04-02 23:43:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/SubstructuringAnalysis.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SubstructuringAnalysis. SubstructuringAnalysis is a subclass 
// of AnalysisAnalysis, it is used when performing a domain decomposition
// analysis. It provides methods which can be invoked by a subdomain to
// perform the numerical computations required.
//
// What: "@(#) SubstructuringAnalysis.h, revA"

#ifndef SubstructuringAnalysis_h
#define SubstructuringAnalysis_h

#include <DomainDecompositionAnalysis.h>
#include <Matrix.h>
#include <Vector.h>

class ConstraintHandler;
class DOF_Numberer;
class AnalysisModel;
class IncrementalIntegrator;
class LinearSOE;
class DomainSolver;
class DomainDecompAlgo;
class Subdomain;
class Vector;

class SubstructuringAnalysis: public DomainDecompositionAnalysis
{
  public:
    SubstructuringAnalysis(Subdomain &theDomain,
			   ConstraintHandler &theHandler,
			   DOF_Numberer &theNumberer,
			   AnalysisModel &theModel,
			   DomainDecompAlgo &theSolnAlgo,		   
			   IncrementalIntegrator &theIntegrator,	
			   LinearSOE &theSOE,
			   DomainSolver &theSolver,
			   ConvergenceTest *theTest);

    virtual int analyze(void);
    virtual ~SubstructuringAnalysis();
    
  protected: 
    
  private:
};

#endif











