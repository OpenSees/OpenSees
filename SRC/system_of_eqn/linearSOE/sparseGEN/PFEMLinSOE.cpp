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
// $Date: 2012-08-31 11:36:36 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMLinSOE.h,v $

// Written: Minjie Zhu
// Created: August 2012
//
// Description: This file contains the class definition for PFEMLinSOE
// PFEMLinSOE is a subclass of SparseGenColLinSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. It solves the equations using the Fractional Step Method in PFEM.
//
// What: "@(#) PFEMLinSOE.h, revA"

#include "PFEMLinSOE.h"
#include <SparseGenColLinSolver.h>
#include <AnalysisModel.h>
#include <Pressure_ConstraintIter.h>
#include <Pressure_Constraint.h>
#include <Node.h>
#include <DOF_Group.h>
#include <Domain.h>


PFEMLinSOE::PFEMLinSOE(SparseGenColLinSolver& solver)
    :SparseGenColLinSOE(solver, LinSOE_TAGS_PFEMLinSOE), theModel(0), mID(0), pID(0), piID(0), 
     mIDall(0), Mid(0), Mhatid(0), Gid(0), Gtid(0), Lid(0), Qtid(0), G(0), Gt(0), L(0), Qt(0)
{
    solver.setLinearSOE(*this);
}


PFEMLinSOE::PFEMLinSOE(int classTag)
    :SparseGenColLinSOE(classTag), theModel(0), mID(0), pID(0), piID(0), mIDall(0), 
     Mid(0), Mhatid(0), Gid(0), Gtid(0), Lid(0), Qtid(0), G(0), Gt(0), L(0), Qt(0)
{
}

PFEMLinSOE::PFEMLinSOE(SparseGenColLinSolver &theSolver, int classTag)
    :SparseGenColLinSOE(theSolver, classTag), theModel(0), mID(0), pID(0), piID(0), mIDall(0), 
     Mid(0), Mhatid(0), Gid(0), Gtid(0), Lid(0), Qtid(0), G(0), Gt(0), L(0), Qt(0)
{
}

PFEMLinSOE::PFEMLinSOE()
    :SparseGenColLinSOE(LinSOE_TAGS_PFEMLinSOE), theModel(0), mID(0), pID(0), piID(0), mIDall(0), 
     Mid(0), Mhatid(0), Gid(0), Gtid(0), Lid(0), Qtid(0), G(0), Gt(0), L(0), Qt(0)
{
}

PFEMLinSOE::~PFEMLinSOE()
{
    if(mID != 0) delete mID;
    if(pID != 0) delete pID;
    if(piID != 0) delete piID;
    if(mIDall != 0) delete mIDall;
    if(Mid != 0) delete Mid;
    if(Mhatid != 0) delete Mhatid;
    if(Gid != 0) delete Gid;
    if(Gtid != 0) delete Gtid;
    if(Lid != 0) delete Lid;
    if(Qtid != 0) delete Qtid;
    if(G != 0) cs_spfree(G);
    if(Gt != 0) cs_spfree(Gt);
    if(L != 0) cs_spfree(L);
    if(Qt != 0) cs_spfree(Qt);
}

int PFEMLinSOE::solve()
{
    // fill matrices
    for(int ind=0; ind<Gid->Size(); ind++) {   // local index
        G->x[ind] = -A[(*Gid)(ind)];            // global index
    }
    for(int ind=0; ind<Gtid->Size(); ind++) {   // local index
        Gt->x[ind] = A[(*Gtid)(ind)];           // global index
    }
    for(int ind=0; ind<Lid->Size(); ind++) {   // local index
        L->x[ind] = A[(*Lid)(ind)];            // global index
    }
    for(int ind=0; ind<Qtid->Size(); ind++) {   // local index
        Qt->x[ind] = A[(*Qtid)(ind)];           // global index
    }

    // predictor : delta U* = Md^{-1} * rm
    double* deltaUs = new double[mID->Size()];
    for(int j=0; j<mID->Size(); j++) {   // local column id
        int col = (*mID)(j);             // global column id
        int ind = (*Mid)(j);             // global index
        if(A[ind] == 0.0) {
            opserr<<"Zero mass at location "<<j<<" ";
            opserr<<" - PFEMLinSOE::solve()\n";
            return -1;
        }
        deltaUs[j] = B[col]/A[ind];          // Md^{-1}*rm
    }

    // pressure : (L+Gt*Md^{-1}*G) * delta P = rp - Gt*delta U*
    cs* GtMd = cs_spalloc(Gt->m, Gt->n, Gt->nzmax, 1, 0);
    GtMd->p[0] = 0;
    for(int j=0; j<mID->Size(); j++) {        // local column id
        int ind = (*Mid)(j);                  // global index
        GtMd->p[j+1] = Gt->p[j+1];            // copy column start
        for(int Gtind=Gt->p[j]; Gtind<Gt->p[j+1]; Gtind++) { // local index
            GtMd->x[Gtind] = Gt->x[Gtind] / A[ind];     // Gt*Md^{-1}
            GtMd->i[Gtind] = Gt->i[Gtind];              // copy row
        }
    } 

    cs* S1 = cs_multiply(GtMd, G);            // Gt*Md^{-1}*G
    cs* S = cs_add(L, S1, 1.0, 1.0);          // L + Gt*Md^{-1}*G
    cs_spfree(S1);
    cs_spfree(GtMd);

    // right hand side : rp - Gt*delta U*
    double* deltaP = new double[pID->Size()];
    for(int i=0; i<pID->Size(); i++) deltaP[i]=0.0;
    cs_gaxpy(Gt, deltaUs, deltaP);           // Gt*delta U*
    for(int i=0; i<pID->Size(); i++) {       // local row id
        int row = (*pID)(i);                 // global row id
        deltaP[i] = B[row] - deltaP[i];      // rp - Gt*delta U*
    }

    // solve pressure : delta P
    cs_qrsol(3, S, deltaP);
    cs_spfree(S);
    

    // corrector
    double* deltaU = new double[mID->Size()];
    for(int i=0; i<mID->Size(); i++) deltaU[i] = 0.0;
    cs_gaxpy(G, deltaP, deltaU);        // G*deltaP
    for(int i=0; i<mID->Size(); i++) {  // local row id
        int ind = (*Mid)(i);            // global index
        deltaU[i] /= A[ind];            // Md^{-1}*G*deltaP
        deltaU[i] += deltaUs[i];        // deltaU* + Md^{-1}*G*deltaP
    }

    // pressure gradient
    double* deltaPi = new double[piID->Size()];
    for(int i=0; i<piID->Size(); i++) deltaPi[i] = 0.0;
    cs_gaxpy(Qt, deltaP, deltaPi);       // Qt*deltaP
    for(int i=0; i<piID->Size(); i++) {  // local row id
        int row = (*piID)(i);            // global row id
        int ind = (*Mhatid)(i);          // global index
        deltaPi[i] = (B[row] - deltaPi[i])/A[ind];  // Mhatd^{-1}*(rpi - Qt*deltaP)
    }

    // copy to X
    for(int i=0; i<mID->Size(); i++) {
        X[(*mID)(i)] = deltaU[i];
    }
    for(int i=0; i<pID->Size(); i++) {
        X[(*pID)(i)] = deltaP[i];
    }
    for(int i=0; i<piID->Size(); i++) {
        X[(*piID)(i)] = deltaPi[i];
    }

    delete [] deltaUs;
    delete [] deltaU;
    delete [] deltaP;
    delete [] deltaPi;

    return 0;
}

int PFEMLinSOE::setLinks(AnalysisModel& model)
{
    theModel = &model;
    return 0;
}

int PFEMLinSOE::setSize(Graph& theGraph)
{
    // this object set size
    this->SparseGenColLinSOE::setSize(theGraph);

    // set Dof IDs
    this->setDofIDs();

    // set matrix IDs
    cs* S = this->setMatIDs();
    if(S == 0) {
        opserr<<"can't create S - PFEMLinSOE::setSize()\n";
        return -1;
    }

    // free memory
    cs_spfree(S);

    return 0;
}

AnalysisModel* PFEMLinSOE::getAnalysisModel()
{
    return theModel;
}

int PFEMLinSOE::sendSelf(int commitTag, Channel &theChannel)
{
    return 0;
}

int PFEMLinSOE::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    return 0;
}

int PFEMLinSOE::setDofIDs()
{
    if(theModel == 0) {
        opserr << "Analysis model has not been linked - PFEMLinSOE::setDofIDs()\n";
        return -1;
    }

    if(mID != 0) delete mID;
    if(pID != 0) delete pID;
    if(piID != 0) delete piID;
    if(mIDall != 0) delete mIDall;
    mID = new ID(0,128);
    pID = new ID(0,128);
    piID = new ID(0,128);
    mIDall = new ID(0,128);

    Domain* domain = theModel->getDomainPtr();
    Pressure_ConstraintIter& thePCs = domain->getPCs();
    Pressure_Constraint* thePC = 0;
    while((thePC = thePCs()) != 0) {
        int nodeTag = thePC->getNodeConstrained();
        Node* theNode = domain->getNode(nodeTag);
        DOF_Group* theDof = theNode->getDOF_GroupPtr();
        const ID& id = theDof->getID();
        for(int i=0; i<id.Size(); i++) {
            if(id(i) >= 0) {
                if(i<2) {
                    mID->insert(id(i));
                } else if(i<3) {
                    pID->insert(id(i));
                } else if(i<5) {
                    piID->insert(id(i));
                } 
            } 
        }
    }

    for(int col=0; col<size; col++) {   // loop all columns
        if(pID->getLocationOrdered(col) >= 0) {  // pressure columns
            continue;
        }
        if(piID->getLocationOrdered(col) >= 0) {  // pressure gradient columns
            continue;
        }
        
        // general momentum columns
        mIDall->insert(col);
    }
    
    return 0;
}

cs* PFEMLinSOE::setMatIDs()
{
    // allocate Mid
    if(Mid != 0 && Mid->Size() != mID->Size()) {
        delete Mid;
        Mid = 0;
    }
    if(Mid == 0) Mid = new ID(mID->Size());

    // allocate Gt1
    cs* Gt1 = cs_spalloc(pID->Size(), mID->Size(), 1, 1, 1);

    // allocate G1
    cs* G1 = cs_spalloc(mID->Size(), pID->Size(), 1, 1, 1);

    // allocate L1
    cs* L1 = cs_spalloc(pID->Size(), pID->Size(), 1, 1, 1);

    // allocate Qt1
    cs* Qt1 = cs_spalloc(piID->Size(), pID->Size(), 1, 1, 1);

    // allocate Mhatid
    if(Mhatid != 0 && Mhatid->Size() != piID->Size()) {
        delete Mhatid;
        Mhatid = 0;
    }
    if(Mhatid == 0) Mhatid = new ID(piID->Size());

    // loop momentum columns
    for(int j=0; j<mID->Size(); j++) {                              // local column id
        int col = (*mID)(j);                                        // global column id
        (*Mid)(j) = -1;

        int i = -1;                                                 // local row id
        for(int ind=colStartA[col]; ind<colStartA[col+1]; ind++) {  // global index
            int row = rowA[ind];                                    // global row id
            if(row == col) {                                        // Mid
                (*Mid)(j) = ind;
            } else if((i=pID->getLocationOrdered(row)) >= 0) {      // Gt
                cs_entry(Gt1, i, j, ind);
            }
        }
        if((*Mid)(j) == -1) {
            opserr<<"WARNING: can't find Mass for PFEM node - PFEMLinSOE::setMatIDs()\n";
            return 0;
        }
    }

    // loop pressure columns;
    for(int j=0; j<pID->Size(); j++) {      // local column id
        int col = (*pID)(j);                // global column id

        int i = -1;                         // local row id
        for(int ind=colStartA[col]; ind<colStartA[col+1]; ind++) {  // global index
            int row = rowA[ind];            // global row id
            if((i=mID->getLocationOrdered(row)) >= 0) {           // G
                cs_entry(G1, i, j, ind);
            } else if((i=pID->getLocationOrdered(row)) >= 0) {    // L
                cs_entry(L1, i, j, ind);
            } else if((i=piID->getLocationOrdered(row)) >= 0) {   // Qt
                cs_entry(Qt1, i, j, ind);
            }
        }
    }

    // loop pressure gradient columns;
    for(int j=0; j<piID->Size(); j++) {    // local column id
        int col = (*piID)(j);              // global column id
        (*Mhatid)(j) = -1;

        for(int ind=colStartA[col]; ind<colStartA[col+1]; ind++) {  // global index
            int row = rowA[ind];           // global row id
            if(row == col) {    // Mhat
                (*Mhatid)(j) = ind;
            }
        }
        if((*Mhatid)(j) == -1) {
            opserr<<"WARNING: can't find Mhat for PFEM node - PFEMLinSOE::setMatIDs()\n";
            return 0;
        }
    }

    // convert to compressed format
    if(G != 0) cs_spfree(G);
    G = cs_compress(G1);
    cs_spfree(G1);

    if(Gt != 0) cs_spfree(Gt);
    Gt = cs_compress(Gt1);
    cs_spfree(Gt1);

    if(L != 0) cs_spfree(L);
    L = cs_compress(L1);
    cs_spfree(L1);

    if(Qt != 0) cs_spfree(Qt);
    Qt = cs_compress(Qt1);
    cs_spfree(Qt1);

    // copy the global index
    if(Gid != 0 && Gid->Size() != G->nzmax) {
        delete Gid;
        Gid = 0;
    }
    if(Gid == 0) Gid = new ID(G->nzmax);
    for(int i=0; i<G->nzmax; i++) {           // local index
        int ind = static_cast<int>(G->x[i]);  // global index
        (*Gid)(i) = ind;                      // copy global index
        G->x[i] = 0.0;                        // zero entry
    }

    if(Gtid != 0 && Gtid->Size() != Gt->nzmax) {
        delete Gtid;
        Gtid = 0;
    }
    if(Gtid == 0) Gtid = new ID(Gt->nzmax);
    for(int i=0; i<Gt->nzmax; i++) {           // local index
        int ind = static_cast<int>(Gt->x[i]);  // global index
        (*Gtid)(i) = ind;                      // copy global index
        Gt->x[i] = 0.0;                        // zero entry
    }

    if(Lid != 0 && Lid->Size() != L->nzmax) {
        delete Lid;
        Lid = 0;
    }
    if(Lid == 0) Lid = new ID(L->nzmax);
    for(int i=0; i<L->nzmax; i++) {           // local index
        int ind = static_cast<int>(L->x[i]);  // global index
        (*Lid)(i) = ind;                      // copy global index
        L->x[i] = 0.0;                        // zero entry
    }

    if(Qtid != 0 && Qtid->Size() != Qt->nzmax) {
        delete Qtid;
        Qtid = 0;
    }
    if(Qtid == 0) Qtid = new ID(Qt->nzmax);
    for(int i=0; i<Qt->nzmax; i++) {           // local index
        int ind = static_cast<int>(Qt->x[i]);  // global index
        (*Qtid)(i) = ind;                      // copy global index
        Qt->x[i] = 0.0;                        // zero entry
    }

    // calculate S = L+Gt*M*G
    cs* S1 = cs_multiply(Gt, G);
    cs* S = cs_add(L, S1, 1.0, 1.0);
    cs_spfree(S1);

    return S;
}
