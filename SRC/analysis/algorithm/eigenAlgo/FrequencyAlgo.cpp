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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/FrequencyAlgo.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/algorithm/eigenAlgo/FrequencyAlgo.C
//
// Written: Jun Peng
// Created: Mon Feb. 8, 1999
// Revision: A
//
// Description: This file contains the class definition of FrequencyAlgo.
// FrequencyAlgo is a class which performs a eigen solution algorithm
// to solve the Generalized eigen equations. It is not expected that 
// this class will have subclasses.
//
// This class is inheritanted from the base class of SolutionAlgorithm
// which was created by fmk (Frank).


#include <FrequencyAlgo.h>
#include <AnalysisModel.h>
#include <EigenAnalysis.h>
#include <EigenIntegrator.h>
#include <EigenSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

FrequencyAlgo::FrequencyAlgo()
  :EigenAlgorithm(EigenALGORITHM_TAGS_Frequency)
{
    // do nothing here.
}

FrequencyAlgo::~FrequencyAlgo()
{
    // do nothing here.
}

int 
FrequencyAlgo::solveCurrentStep(int numModes)
{
    AnalysisModel *theModel = this->getAnalysisModelPtr();
    EigenSOE *theSOE = this->getEigenSOEptr();
    EigenIntegrator *theIntegrator = this->getEigenIntegratorPtr();

    if ((theModel == 0) || (theIntegrator == 0) || (theSOE == 0)) {
       opserr << "WARNING FrequencyAlgo::solverCurrentStep() - ";
       opserr << "setLinks() has not been called. \n";
       return -1;
    }

    if (theIntegrator->formK() < 0) {
       opserr << "WARNING FrequencyAlgo::solverCurrentStep() - ";
       opserr << "the Integrator failed in formK().\n";
       return -2;
    }

    if (theIntegrator->formM() < 0) {
       opserr << "WARNING FrequencyAlgo::solverCurrentStep() - ";
       opserr << "the Integrator failed in formK().\n";
       return -3;
    }

    if (theSOE->solve(numModes) < 0) {
       opserr << "Warning FrequencyAlgo::solveCurrentStep() - ";
       opserr << "the EigenSOE failed in solve().\n";
       return -4;
    }

    // now set the eigenvalues and eigenvectors in the model
    theModel->setNumEigenvectors(numModes);
    Vector theEigenvalues(numModes);
    for (int i=1; i<=numModes; i++) {
	theEigenvalues[i-1] = theSOE->getEigenvalue(i);
	theModel->setEigenvector(i, theSOE->getEigenvector(i));
    }    
    theModel->setEigenvalues(theEigenvalues);
    
    return 0;
}

int 
FrequencyAlgo::sendSelf(int cTag, Channel &theChannel)
{
    return 0;
}

int 
FrequencyAlgo::recvSelf(int cTag, Channel &theChannel,
			FEM_ObjectBroker &theBroker)
{
    return 0;
}

void 
FrequencyAlgo::Print(OPS_Stream &s, int flag)
{
    s << "\t Eigen Algorithm \n";
}


