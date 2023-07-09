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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-04-02 23:43:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/SubstructuringAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/Analysis/SubstructuringAnalysis.C
// 
// Written: fmk 
// Created: Tue Sept 17 16:34:47: 1996
// Revision: A
//
// Description: This file contains the class definition for 
// SubstructuringAnalysis. SubstructuringAnalysis is a subclass 
// of DomainDecompositionAnalysis, it is used to perform the operations
// required of a domain decomposition substructuring method.
//
// What: "@(#) SubstructuringAnalysis.C, revA"

#include <SubstructuringAnalysis.h>
#include <ConstraintHandler.h>
#include <DOF_Numberer.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <DomainDecompAlgo.h>
#include <DomainSolver.h>

#include <IncrementalIntegrator.h>
#include <Subdomain.h>

#include <Matrix.h>
#include <ID.h>


// Constructor
//    sets theModel and theSysOFEqn to 0 and the Algorithm to the one supplied

SubstructuringAnalysis::SubstructuringAnalysis(Subdomain &the_Domain,
					       ConstraintHandler &handler,
					       DOF_Numberer &numberer,
					       AnalysisModel &model,
					       DomainDecompAlgo &theSolnAlgo,
					       IncrementalIntegrator &integrator,
					       LinearSOE &theLinSOE,
					       DomainSolver &theDDSolver,
					       ConvergenceTest *theTest)
:DomainDecompositionAnalysis(the_Domain,
			     handler,
			     numberer,
			     model,
			     theSolnAlgo,
			     integrator,
			     theLinSOE,
			     theDDSolver,
			     theTest)
{

}    


SubstructuringAnalysis::~SubstructuringAnalysis()
{

}    

int 
SubstructuringAnalysis::analyze(void)
{
    opserr << "SubstructuringAnalysis::analyze(void)";
    opserr << "does nothing and should not have been called\n";
    return -1;
}

