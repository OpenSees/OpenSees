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

// $Revision: 1.6 $
// $Date: 2007-04-02 23:42:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/integrator/StagedLoadControl.cpp,v $


//
// Written: fmk
// Created: 07/98
// Revision: A
//
// Description: This file contains the class definition for StagedLoadControl.
// StagedLoadControl is an algorithmic class for performing a static analysis
// using a load control integration scheme.
//
// What: "@(#) StagedLoadControl.h, revA"



#include <StagedLoadControl.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <Vector.h>
#include <Channel.h>
#include <Element.h>
#include <FE_Element.h>
#include <FE_EleIter.h>
#include <Node.h>
#include <DOF_Group.h>
#include <DOF_GrpIter.h>
#include <LoadPattern.h>
#include <LoadPatternIter.h>
#include <Domain.h>
#include <Parameter.h>
#include <ParameterIter.h>
#include <EquiSolnAlgo.h>
#include <elementAPI.h>
#include <iostream>

#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#endif

// static double doubleone = 1.0;
// static Matrix one(&doubleone,1,1);
// static ID dofid(1);

#define SEQUENTIAL_SECTION_BEGIN  for (int proc = 0; proc < nproc; ++proc){ if (rank == proc) {
#define SEQUENTIAL_SECTION_END } MPI_Barrier(MPI_COMM_WORLD);}


void* OPS_StagedLoadControlIntegrator()
{
    if (OPS_GetNumRemainingInputArgs() < 1) {
        opserr << "insufficient arguments\n";
        return 0;
    }

    double lambda;
    int numData = 1;
    if (OPS_GetDoubleInput(&numData, &lambda) < 0) {
        opserr << "WARNING failed to read double lambda\n";
        return 0;
    }

    int numIter = 1;
    double mLambda[2] = {lambda, lambda};
    if (OPS_GetNumRemainingInputArgs() > 2) {
        if (OPS_GetIntInput(&numData, &numIter) < 0) {
            opserr << "WARNING failed to read int numIter\n";
            return 0;
        }
        numData = 2;
        if (OPS_GetDoubleInput(&numData, &mLambda[0]) < 0) {
            opserr << "WARNING failed to read double min and max\n";
            return 0;
        }
    }

    return new StagedLoadControl(lambda, numIter, mLambda[0], mLambda[1]);
}


StagedLoadControl::StagedLoadControl()
    : LoadControl(0, 0, 0, 0, INTEGRATOR_TAGS_StagedLoadControl)
{
}



StagedLoadControl::StagedLoadControl(double dLambda, int numIncr, double min, double max)
    : LoadControl(dLambda, numIncr, min, max, INTEGRATOR_TAGS_StagedLoadControl)
{
}


int StagedLoadControl::formTangent(int statFlag)
{
    int rank = 0;
    int nproc = 1;

    #ifdef _PARALLEL_PROCESSING
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    #endif

    // Run a typocal LoadControl formTangent call
    int errflag = this->IncrementalIntegrator::formTangent(statFlag);

    if (errflag < 0)
    {
        return errflag;
    }

    // Now detect inactive nodes and add 1 to the tangent diagonal there

    AnalysisModel *theAnalysisModel = this->getAnalysisModel();
    Domain *theDomain = theAnalysisModel->getDomainPtr();
    LinearSOE *theSOE = this->getLinearSOE();
    int numEqn = theSOE->getNumEqn();

    int * nodedofs = new int[numEqn + 1];
    #ifdef _PARALLEL_PROCESSING
    int * allnodedofs = new int[numEqn + 1];
    #endif

    for (int i = 0; i < numEqn; ++i)
    {
        nodedofs[i] = 0;
    #ifdef _PARALLEL_PROCESSING
        allnodedofs[i] = 0;
    #endif
    }

    FE_Element *elePtr = 0;

    FE_EleIter &theEles = theAnalysisModel->getFEs();

    while ((elePtr = theEles()) != 0) {
        const ID& elenodedofs = elePtr->getID();

        for (int i = 0; i < elenodedofs.Size(); ++i)
        {
            int dof = elenodedofs(i);
            if (dof > numEqn)
            {
                std::cout << "i = " << i << std::endl;
                std::cout << "numEqn = " << numEqn << std::endl;
                std::cout << "elenodedofs(i) = " << dof << std::endl;
                exit(-1);
            }
            // std::cout << "i = " << i << " numEqn = " << numEqn << " dof = " << dof << std::endl;
            if (dof >= 0 && elePtr->isActive())
            {

                nodedofs[dof] = (int) 1;
            }
            
        }
    }


    #ifdef _PARALLEL_PROCESSING
    MPI_Allreduce(nodedofs, allnodedofs, numEqn, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    #endif


    for (int i = 0; i < numEqn; ++i)
    {
        #ifdef _PARALLEL_PROCESSING
        bool is_lonely_dof = allnodedofs[i] == 0;
        #else
        bool is_lonely_dof = nodedofs[i] == 0;
        #endif

        if (is_lonely_dof)
        {
            // opserr << "i = " << i << " nodedofs(i) = " << nodedofs[i] << endln;
            double uno = 1.0;
            static ID dofid(1);
            static Matrix one(1, 1);
            one(0, 0) = uno;
            dofid(0) = i;
            theSOE->addA(one, dofid);
        }
    }

    delete [] nodedofs;

    #ifdef _PARALLEL_PROCESSING
    delete [] allnodedofs;
    #endif

    return errflag;
}
