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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-11-30 13:01:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/PFEMAnalysis.h,v $
                                                                        
                                                                        
// Written: Minjie  
// Created: Nov 2012
// Revision: A
//
// Description: This file contains the class definition for 
// PFEMAnalysis. PFEMAnalysis is a 
// subclass of DirectIntegrationAnalysis. It is used to perform a 
// dynamic analysis for PFEM using a direct integration scheme.
//
// What: "@(#) PFEMAnalysis.h, revA"

#include <PFEMAnalysis.h>
#include <EquiSolnAlgo.h>
#include <AnalysisModel.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <ConvergenceTest.h>
#include <TransientIntegrator.h>
#include <Domain.h>
#include <Node.h>
#include <SP_Constraint.h>
#include <Pressure_Constraint.h>
#include <Element.h>
#include <SP_ConstraintIter.h>
#include <Pressure_ConstraintIter.h>
#include <ElementIter.h>
#include <map>

#ifdef _PARALLEL_INTERPRETERS
#include <mpi.h>
#endif

PFEMAnalysis::PFEMAnalysis(Domain& theDomain, 
                           ConstraintHandler& theHandler,
                           DOF_Numberer& theNumberer,
                           AnalysisModel& theModel,
                           EquiSolnAlgo& theSolnAlgo,		   
                           LinearSOE& theSOE,
                           TransientIntegrator& theIntegrator,
                           ConvergenceTest* theTest,
                           double max, double min, double g, double r)
    :DirectIntegrationAnalysis(theDomain,theHandler,theNumberer,theModel,theSolnAlgo,
                               theSOE,theIntegrator,theTest),
     dtmax(max), dtmin(min), gravity(g), ratio(r),
     dt(max), next(max), newstep(true)
{
    
}

PFEMAnalysis::~PFEMAnalysis()
{
}

int 
PFEMAnalysis::analyze()
{
    Domain* theDomain = this->getDomainPtr();
    double current = theDomain->getCurrentTime();
    if(newstep) {
        next = current + dtmax;
    }
    bool instep = false;
	
    while(true) {

	// identify domain
	if (this->identify() < 0) {
	    opserr<<"WARNING: failed to identify domain -- PFEMAnalysis\n";
	    return -1;
	}

        // analyze
#ifdef _PARALLEL_INTERPRETERS
	int myid;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if(myid == 0) {
	    opserr<<"\n\nTime = "<<current<<", dt = "<<dt<<"\n\n";
	}
#else
        opserr<<"\n\nTime = "<<current<<", dt = "<<dt<<"\n\n";
#endif

	// analysis
        int converged  = DirectIntegrationAnalysis::analyze(1, dt);

        // if failed
        if(converged < 0) {
            dt *= ratio;
            if(dt < dtmin) {
                return -1;
            }
            instep = true;
            newstep = false;
            continue;
        }

        // if converged
        if(instep) {
            current = theDomain->getCurrentTime();
            dt = next - current;
        } else {
            newstep = true;
            dt = dtmax;
        }
	
	break;
    }

    return 0;
}

int
PFEMAnalysis::identify()
{
    Domain* theDomain = this->getDomainPtr();

    if (theDomain->getNumPCs() == 0) return 0;

    // disconnect pcs from structural elements
    Pressure_ConstraintIter& thePCs = theDomain->getPCs();
    Pressure_Constraint* thePC = 0;
    while((thePC = thePCs()) != 0) {
        thePC->disconnect();
    }

    // to connect pc to all elements
    Element *elePtr;
    ElementIter &theElemIter = theDomain->getElements();    
    while((elePtr = theElemIter()) != 0) {

        const ID& nodeids = elePtr->getExternalNodes();
        for(int i=0; i<nodeids.Size(); i++) {
            Pressure_Constraint* thePC = theDomain->getPressure_Constraint(nodeids(i));
            if(thePC != 0) {
                thePC->connect(elePtr->getTag(), false);
            }
        }
    }

    // get all SPs and MPs
    SP_ConstraintIter& theSPs = theDomain->getDomainAndLoadPatternSPs();
    SP_Constraint* theSP = 0;
    std::map<int,ID> sps;
    while((theSP = theSPs()) != 0) {
        sps[theSP->getNodeTag()].insert(theSP->getDOF_Number());
    }

    // get number of parameters
    int numParams = theDomain->getNumParameters();

    // move the isolated nodes
    Pressure_ConstraintIter& thePCs1 = theDomain->getPCs();
    while((thePC = thePCs1()) != 0) {

        if(thePC->isIsolated()) {
            int tag = thePC->getTag();
            Node* node = theDomain->getNode(tag);
            if(node == 0) {
                opserr<<"WARNING: node "<<tag;
                opserr<<" does not exist -- PFEMAnalysis::identify\n";
		return -1;
            }
            std::map<int,ID>::iterator it = sps.find(tag);

            // responses
            const Vector& disp = node->getDisp();
            const Vector& vel = node->getVel();
            const Vector& accel = node->getAccel();
            Vector ndisp = disp;
            Vector nvel = vel;
            Vector naccel = accel;
            naccel.Zero();
            naccel(1) = gravity;
            ndisp.addVector(1.0, vel, ops_Dt);
            ndisp.addVector(1.0, accel, 0.5*ops_Dt*ops_Dt);
            nvel.addVector(1.0, accel, ops_Dt);

            if(it != sps.end()) {
                const ID& dofs = it->second;
                for(int i=0; i<dofs.Size(); i++) {
                    ndisp(dofs(i)) = disp(dofs(i));
                    nvel(dofs(i)) = vel(dofs(i));
                    naccel(dofs(i)) = accel(dofs(i));
                }
            }


            node->setTrialDisp(ndisp);
            node->setTrialVel(nvel);
            node->setTrialAccel(naccel);
            node->commitState();

            // sensitivity
            for(int grad=0; grad<numParams; grad++) {
                Vector sensdisp(disp.Size());
                Vector sensvel(vel.Size());
                Vector sensaccel(accel.Size());

                for(int dof=0; dof<disp.Size(); dof++) {
                    bool fixed = false;
                    if(it != sps.end()) {
                        const ID& dofs = it->second;
                        for(int i=0; i<dofs.Size(); i++) {
                            if(dof == dofs(i)) {
                                fixed = true;
                                break;
                            }
                        }
                    }
                    if(!fixed) {
                        sensvel(dof) = node->getVelSensitivity(dof,grad);
                        sensdisp(dof) = node->getDispSensitivity(dof,grad)+ops_Dt*sensvel(dof);
                    }
                }

                node->saveDispSensitivity(sensdisp,grad,numParams);
                node->saveVelSensitivity(sensvel,grad,numParams);
                node->saveAccelSensitivity(sensaccel,grad,numParams);
            }
        }
    }

    return 0;
}
