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
// $Date: 2003-02-14 23:00:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/domainDecompAlgo/DomainDecompAlgo.cpp,v $
                                                                        
                                                                        
// File: ~/OOP/analysis/algorithm/DomainDecompAlgo.C
// 
// Written: fmk 
// Created: Sun Sept 15 15:06:47: 1996 
// Revision: A 
//

// Description: This file contains the class definition for 
// DomainDecompAlgo. DomainDecompAlgo is an abstract base class, 
// i.e. no objects of it's type can be created.  Its subclasses deifine
// the sequence of operations to be performed in the analysis by static
// equilibrium of a finite element model.  
// 
// What: "@(#)DomainDecompAlgo.h, revA"

#include <DomainDecompAlgo.h>
#include <AnalysisModel.h>
#include <IncrementalIntegrator.h>
#include <LinearSOE.h>
#include <DomainSolver.h>
#include <Subdomain.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

DomainDecompAlgo::DomainDecompAlgo()
:SolutionAlgorithm(DomDecompALGORITHM_TAGS_DomainDecompAlgo),
 theModel(0), theIntegrator(0), theLinearSOE(0), theSolver(0),
 theSubdomain(0)
{

}


DomainDecompAlgo::~DomainDecompAlgo()
{

}

int
DomainDecompAlgo::solveCurrentStep(void)
{
    if (theModel == 0 || theIntegrator == 0 || theLinearSOE == 0 ||
	theSolver == 0 || theSubdomain != 0 ) {

	const Vector &extResponse = 
	    theSubdomain->getLastExternalSysResponse();

	theSolver->setComputedXext(extResponse);
	theSolver->solveXint();

	theIntegrator->update(theLinearSOE->getX());
	
	return 0;
    }
    else {
	opserr << "DomainDecompAlgo::solveCurrentStep() ";
	opserr << "no links have been set\n";
	return -1;
    }
}

void 
DomainDecompAlgo::setLinks(AnalysisModel &theAnaModel, 
				  IncrementalIntegrator &theInteg,
				  LinearSOE &theSOE,
				  DomainSolver &theDomainSolver,
				  Subdomain &theSub)
{
    theModel = &theAnaModel;
    theIntegrator = &theInteg;
    theLinearSOE = &theSOE;
    theSolver = &theDomainSolver;
    theSubdomain = &theSub;
}
    


int
DomainDecompAlgo::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int
DomainDecompAlgo::recvSelf(int ctag, Channel &theChannel, 
			   FEM_ObjectBroker &theBroker)
{
    return 0;
}

