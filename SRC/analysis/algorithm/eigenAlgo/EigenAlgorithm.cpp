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
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/eigenAlgo/EigenAlgorithm.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/algorithm/eigenAlgo/EigenAlgorithm.C
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


#include <EigenAlgorithm.h>
#include <AnalysisModel.h>
#include <EigenIntegrator.h>
#include <EigenSOE.h>

EigenAlgorithm::EigenAlgorithm(int classTag)
  :SolutionAlgorithm(classTag),
   theModel(0), theIntegrator(0), theSOE(0)
{
    // need do nothing here.
}


EigenAlgorithm::~EigenAlgorithm()
{
    // do nothing here.
}

void 
EigenAlgorithm::setLinks(AnalysisModel &theNewModel,
			 EigenIntegrator &theNewIntegrator,
			 EigenSOE &theNewSOE)
{
    theModel = &theNewModel;
    theIntegrator = &theNewIntegrator;
    theSOE = &theNewSOE;
}

AnalysisModel *
EigenAlgorithm::getAnalysisModelPtr() const
{
    return theModel;
}

EigenIntegrator *
EigenAlgorithm::getEigenIntegratorPtr() const
{
    return theIntegrator;
}

EigenSOE *
EigenAlgorithm::getEigenSOEptr() const
{
    return theSOE;
}

