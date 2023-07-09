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
// $Date: 2009-07-29 21:57:42 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/StandardEigenAlgo.cpp,v $
                                                                        
// Written: MHS
// Created: Oct 2001
//
// Description: This file contains the class definition of StandardEigenAlgo.
// StandardEigenAlgo is a class which performs a eigen solution algorithm
// to solve standard eigenvalue equations. It is not expected that 
// this class will have subclasses.

#include <StandardEigenAlgo.h>
#include <AnalysisModel.h>
#include <EigenAnalysis.h>
#include <EigenIntegrator.h>
#include <EigenSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Timer.h>

StandardEigenAlgo::StandardEigenAlgo()
  :EigenAlgorithm(EigenALGORITHM_TAGS_Standard)
{
  // do nothing here.
}

StandardEigenAlgo::~StandardEigenAlgo()
{
  // do nothing here.
}

int 
StandardEigenAlgo::solveCurrentStep(int numModes)
{
  AnalysisModel *theModel = this->getAnalysisModelPtr();
  EigenSOE *theSOE = this->getEigenSOEptr();
  EigenIntegrator *theIntegrator = this->getEigenIntegratorPtr();
  
  if ((theModel == 0) || (theIntegrator == 0) || (theSOE == 0)) {

    opserr << "StandardEigenAlgo::solverCurrentStep() -- setLinks() has not been called\n";
    return -1;
  }
  
  if (theIntegrator->formK() < 0) {
    opserr << "StandardEigenAlgo::solverCurrentStep() -- the Integrator failed in formK()\n";
    return -2;
  }
  
  if (theSOE->solve(numModes, false) < 0) {
    opserr << "StandardEigenAlgo::solverCurrentStep() -- the EigenSOE failed in solve()\n";
    return -4;
  }
  
  // now set the eigenvalues and eigenvectors in the model
  theModel->setNumEigenvectors(numModes);
  Vector theEigenvalues(numModes);
  for (int i = 1; i <= numModes; i++) {
    theEigenvalues[i-1] = theSOE->getEigenvalue(i);
    theModel->setEigenvector(i, theSOE->getEigenvector(i));
  }    
  theModel->setEigenvalues(theEigenvalues);
  
  return 0;
}

int 
StandardEigenAlgo::sendSelf(int cTag, Channel &theChannel)
{
  return 0;
}

int 
StandardEigenAlgo::recvSelf(int cTag, Channel &theChannel,
			  FEM_ObjectBroker &theBroker)
{
  return 0;
}

void 
StandardEigenAlgo::Print(OPS_Stream &s, int flag)
{
  s << "\tStandardEigenAlgo\n";
}
