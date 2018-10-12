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

// $Revision $
// $Date$

// Written: Minjie Zhu
//
// Description: This file contains the class definition for PFEMDiaLinSOE
// PFEMDiaLinSOE is a subclass of SparseGenColLinSOE.
// It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the
// matrix A. It solves the equations using the Diagonal
// Fractional Step Method in PFEM.
//

#include <PFEMDiaSolver.h>
#include <PFEMDiaLinSOE.h>
#include <iostream>
#include <cmath>
#include <Timer.h>

void* OPS_PFEMDiaSolver()
{
    PFEMDiaSolver* theSolver = new PFEMDiaSolver();
    return new PFEMDiaLinSOE(*theSolver);
}

PFEMDiaSolver::PFEMDiaSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMDiaSolver), theSOE(0)
{
}

PFEMDiaSolver::~PFEMDiaSolver()
{
}

int
PFEMDiaSolver::solve()
{
    // Timer timer;
    // timer.start();

    // Gt and -G
    cs* Gt = theSOE->Gt;
    cs* G = theSOE->G;

    // Mp and M
    Vector& Mp = theSOE->Mp;
    Vector& M = theSOE->M;

    // size
    const Vector& B = theSOE->getB();
    const ID& dofType = theSOE->getDofType();
    const ID& dofID = theSOE->getDofID();

    int Vsize = M.Size();
    int Psize = Mp.Size();
    int size = B.Size();

    // tolv, tolp, numiter
    double tolv = 1e-6;
    double tolp = 1e-4;
    int numiter = 1000;

    // timer.pause();
    // timer.start();

    // initial dv and dp
    Vector dP(Psize);
    Vector dV(Vsize);
    double* dP_ptr = 0;
    double* dV_ptr = 0;
    if (Psize > 0) {
	dP_ptr = &dP(0);
    }
    if (Vsize > 0) {
	dV_ptr = &dV(0);
    }

    bool converged = false;
    int numi = 0;
    for (int j=0; j<numiter; ++j) {

	// (-G)*dp
	if (Psize > 0 && Vsize > 0) {
	    cs_gaxpy(G, dP_ptr, dV_ptr);
	}

	// dv
	for(int i=0; i<size; i++) {        // row

	    int rowtype = dofType(i);      // row type
	    int rowid = dofID(i);          // row id
	    if(rowtype < 3 && rowtype>=0) {
		if (M(rowid) == 0) {
		    opserr << "M("<<rowid<<") = 0\n";
		    return -1;
		}
		dV(rowid) = (B(i)-dV(rowid))/M(rowid);
	    }
	}

	// Gt*dv
	if (Psize > 0 && Vsize > 0) {
	    cs_gaxpy(Gt, dV_ptr, dP_ptr);
	}

	// dp
	for(int i=0; i<size; i++) {        // row

	    int rowtype = dofType(i);      // row type
	    int rowid = dofID(i);          // row id
	    if(rowtype == 3) {
		if (Mp(rowid) == 0) {
		    opserr << "Mp("<<rowid<<") = 0\n";
		    return -1;
		}
		dP(rowid) = (B(i)-dP(rowid))/Mp(rowid);
	    }
	}

	if (dV.Norm()<tolv && dP.Norm()<tolp) {
	    numi = j+1;
	    converged = true;
	    break;
	}
    }

    if (converged) {
	opserr<<"Converged in "<<numi+1<<"iterations\n";
    } else {
	opserr<<"Failed to converged, norm(dv) = "<<dV.Norm()<<", norm(dp) = "<<dP.Norm()<<"\n";
	return -1;
    }

    // copy to X
    Vector X(size);
    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i);
        if(rowtype>=0 && rowtype<3) {
            X(i) = dV(rowid);
        } else if(rowtype == 3) {
            X(i) = dP(rowid);
        }

    }
    theSOE->setX(X);

    // timer.pause();
    // opserr<<"solving time for PFEMDiasolver = "<<timer.getReal()<<"\n";

    return 0;
}

int PFEMDiaSolver::setSize()
{
    return 0;
}

int
PFEMDiaSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMDiaSolver::recvSelf(int ctag,
		  Channel &theChannel,
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int
PFEMDiaSolver::setLinearSOE(PFEMDiaLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
