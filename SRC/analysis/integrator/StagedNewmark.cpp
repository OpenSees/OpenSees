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


#include <StagedNewmark.h>
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



static bool converged = false;
static int count = 0;


#ifdef _PARALLEL_PROCESSING
#include <mpi.h>
#endif




void *
OPS_StagedNewmark(void)
{
  // Pointer to a uniaxial material that will be returned
  TransientIntegrator *theIntegrator = 0;

  int argc = OPS_GetNumRemainingInputArgs();
  if (argc != 2 && argc != 4) {
    opserr << "WARNING - incorrect number of args want StagedNewmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }

  bool dispFlag = true;
  double dData[2];
  int numData = 2;
  if (OPS_GetDouble(&numData, dData) != 0) {
    opserr << "WARNING - invalid args want StagedNewmark $gamma $beta <-form $typeUnknown>\n";
    return 0;
  }
  
  if (argc == 2)
    theIntegrator = new StagedNewmark(dData[0], dData[1]);
  else {
    //    char nextString[10];
    const char *nextString = OPS_GetString();
    //    OPS_GetString(nextString, 10);
    if (strcmp(nextString,"-form") == 0) {
      //      OPS_GetString(nextString, 10);
      nextString = OPS_GetString();
      if ((nextString[0] == 'D') || (nextString[0] == 'd')) 
	dispFlag = true;
      else if ((nextString[0] == 'A') || (nextString[0] == 'a')) 
	dispFlag = false;      
    }    
    theIntegrator = new StagedNewmark(dData[0], dData[1], dispFlag);
  }

  if (theIntegrator == 0)
    opserr << "WARNING - out of memory creating StagedNewmark integrator\n";

  return theIntegrator;
}



StagedNewmark::StagedNewmark()
    : Newmark( 0,  0,  0,  0,  INTEGRATOR_TAGS_StagedNewmark)
{
    
}



StagedNewmark::StagedNewmark(double _gamma, double _beta, bool dispFlag, bool aflag)
    : Newmark( _gamma,  _beta,  dispFlag,  aflag,  INTEGRATOR_TAGS_StagedNewmark)
{
    
}




int StagedNewmark::formTangent(int statFlag)
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
