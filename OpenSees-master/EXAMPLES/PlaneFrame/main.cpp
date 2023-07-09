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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-02-19 22:12:16 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/PlaneFrame/main.cpp,v $
                                                                        
// File: ~/model/main.C
//
// Written: fmk 08/99
//
// Purpose: this file contains a C++ main procedure to perform the analysis
// of example1 (found in most documents). In the main() procedure:
// 	1) each object of the domain, i.e. Nodes, Elements, Constraints,
//	   and LoadPattern objects are created and then added to the Domain.
//	2) the components of the analysis object are constructed and then
//	   the Analysis object is created.
//	3) the analysis is performed.
//	4) the results are printed - here the contents of Domain and end of
//	   the analysis operation.


// standard C++ includes
#include <stdlib.h>
#include <OPS_Globals.h>
#include <StandardStream.h>

// includes for the domain classes
#include <Domain.h>
#include <PlaneFrame.h>

// includes for the analysis classes
#include <StaticAnalysis.h>
#include <AnalysisModel.h>
#include <Linear.h>
#include <PenaltyConstraintHandler.h>
#include <DOF_Numberer.h>
#include <RCM.h>
#include <LoadControl.h>
#include <BandSPDLinSOE.h>
#include <BandSPDLinLapackSolver.h>


// init the global variabled defined in OPS_Globals.h
StandardStream sserr;
OPS_Stream &opserr = sserr;

double        ops_Dt = 0;
Domain       *ops_TheActiveDomain = 0;
Element      *ops_TheActiveElement = 0;


int main(int argc, char **argv)
{
    
    //
    //	now create a domain and a modelbuilder
    //  and build the model
    //

    Domain *theDomain = new Domain();
    
    // create the Modelbuilder and invoke buildFE_Model() on it to create
    // the domain components.

    ModelBuilder *theModelBuilder = new PlaneFrame(*theDomain);
    theModelBuilder->buildFE_Model();
    
    // create an Analysis object to perform a static analysis of the model
    //  - constructs:
    //    AnalysisModel of type AnalysisModel,
    //	  EquiSolnAlgo of type Linear
    //	  StaticIntegrator of type LoadControl
    //	  ConstraintHandler of type Penalty
    //    DOF_Numberer which uses RCM
    //    LinearSOE of type Band SPD
    // and then the StaticAnalysis object
    
    AnalysisModel     *theModel = new AnalysisModel();
    EquiSolnAlgo      *theSolnAlgo = new Linear();
    StaticIntegrator  *theIntegrator = new LoadControl(1.0, 1, 1.0, 1.0);
    ConstraintHandler *theHandler = new PenaltyConstraintHandler(1.0e8,1.0e8);
    RCM               *theRCM = new RCM();
    DOF_Numberer      *theNumberer = new DOF_Numberer(*theRCM);    
    BandSPDLinSolver  *theSolver = new BandSPDLinLapackSolver();       
    LinearSOE         *theSOE = new BandSPDLinSOE(*theSolver);        

    StaticAnalysis    theAnalysis(*theDomain,
				  *theHandler,
				  *theNumberer,
				  *theModel,
				  *theSolnAlgo,
				  *theSOE,
				  *theIntegrator);

    // perform the analysis & print out the results for the domain
    theAnalysis.analyze(1);
    opserr << *theDomain;

    exit(0);
}	
	
