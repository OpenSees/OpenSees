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
// $Date: 2012-09-17 10:51:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver.h
//
// Written: Minjie 
// Created: Sep 17 2012
//

#include <PFEMCompressibleSolver.h>
#include <PFEMCompressibleLinSOE.h>
#include <iostream>
#include <cmath>
#include <Timer.h>

void* OPS_PFEMCompressibleSolver()
{
    PFEMCompressibleSolver* theSolver = new PFEMCompressibleSolver();
    return new PFEMCompressibleLinSOE(*theSolver);
}

PFEMCompressibleSolver::PFEMCompressibleSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMCompressibleSolver), theSOE(0)
{
    // set default control parameters
    umfpack_di_defaults(Control);
}

PFEMCompressibleSolver::~PFEMCompressibleSolver()
{
    
}

int
PFEMCompressibleSolver::solve()
{
    // Timer timer;
    // timer.start();
    cs* M = theSOE->M;
    cs* Gt = theSOE->Gt;
    cs* G = theSOE->G;
    const Vector& B = theSOE->getB();
    Vector& Mp = theSOE->Mp;
    const ID& dofType = theSOE->getDofType();
    const ID& dofID = theSOE->getDofID();
    
    int Vsize = M->n;
    int Psize = Mp.Size();
    int size = B.Size();

    if(Vsize<=0 || Psize<=0) {
        opserr<<"WARNING: Fsize or Psize or Pisize <= 0 -- ";
        opserr<<"PFEMCompressibleSolver::solve\n";
        return -1;
    }
    // timer.pause();
    // timer.start();

    // get velocity tangent and rhs
    cs* GMp = G;
    for(int j=0; j<Psize; j++) {
        if(Mp(j) == 0.0) {
            opserr<<"WARNING: Mp is zero at "<<j<<"\n";
            return -1;
        }
        for(int k=GMp->p[j]; k<GMp->p[j+1]; k++) {
            double& x = GMp->x[k];
            x /= Mp(j);
        }
    }

    cs* GMpGt = cs_multiply(GMp, Gt);
    cs* H = cs_add(M, GMpGt, 1.0, -1.0);
    cs_spfree(GMpGt);
    
    Vector dP(Psize);
    double* dP_ptr = &dP(0);
    for(int i=0; i<size; i++) {        // row

        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype == 3) {
            dP(rowid) = B(i);          // rp
        }
    }

    Vector dV(Vsize);
    double* dV_ptr = &dV(0);
    cs_gaxpy(GMp, dP_ptr, dV_ptr);     // GMp*rp
    for(int i=0; i<size; i++) {        // row
        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype < 3 && rowtype>=0) {
            dV(rowid) = B(i)-dV(rowid);         // rm-GMp*rp
        }
    }
    // timer.pause();
    // opserr<<"velocity tangent time = "<<timer.getReal()<<"\n";
    // timer.start();

    // reorder rows of H for umfpack
    int* Ap = H->p;
    int* Ai = H->i;
    double* Ax = H->x;

    for (int j=0; j<Vsize; j++) {
	ID col(0, Ap[j+1]-Ap[j]);
	Vector colval(Ap[j+1]-Ap[j]);
	ID col0(colval.Size());
	int index = 0;
	for (int k=Ap[j]; k<Ap[j+1]; k++) {
	    col.insert(Ai[k]);
	    col0(index) = Ai[k];
	    colval(index++) = Ax[k];
	}

	index = 0;
	for (int k=Ap[j]; k<Ap[j+1]; k++) {
	    Ai[k] = col[index++];
	    Ax[k] = colval(col0.getLocation(Ai[k]));
	}
    }

    
    // symbolic analysis of H
    void* Symbolic = 0;
    int status = umfpack_di_symbolic(Vsize,Vsize,Ap,Ai,Ax,&Symbolic,Control,Info);

    // check errors
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: symbolic analysis returns "<<status<<" -- PFEMCompressibleSolver::solve\n";
	return -1;
    }
    // timer.pause();
    // opserr<<"analysis time = "<<timer.getReal()<<"\n";
    // timer.start();

    // numerical analysis
    void* Numeric = 0;
    status = umfpack_di_numeric(Ap,Ai,Ax,Symbolic,&Numeric,Control,Info);

    // check error
    if (status!=UMFPACK_OK) {
    	opserr<<"WARNING: numeric analysis returns "<<status<<" -- PFEMCompressibleSolver::solve\n";
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
    	return -1;
    }

    // timer.pause();
    // opserr<<"numeric time = "<<timer.getReal()<<"\n";
    // timer.start();
    
    // solve velocity
    Vector x(Vsize);
    double* x_ptr = &x(0);
    status = umfpack_di_solve(UMFPACK_A,Ap,Ai,Ax,x_ptr,dV_ptr,Numeric,Control,Info);
    dV = x;
    
    // check error
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: solving returns "<<status<<" -- PFEMCompressibleSolver::solve\n";
	if (Numeric != 0) {
	    umfpack_di_free_numeric(&Numeric);
	}
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
	return -1;
    }
	
    // timer.pause();
    // opserr<<"solve velocity time = "<<timer.getReal()<<"\n";
    // timer.start();
    
    // solve pressure
    dP *= -1;                          // -rp
    cs_gaxpy(Gt, dV_ptr, dP_ptr);      // Gt*dV - rp
    for(int i=0; i<Psize; i++) {       // (rp - Gt*dV) / Mp
        dP(i) /= -Mp(i);
    }
    // timer.pause();
    // opserr<<"solve pressure time = "<<timer.getReal()<<"\n";

    // clean up
    if (Numeric != 0) {
	umfpack_di_free_numeric(&Numeric);
    }
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
    }
    cs_spfree(H);


    // copy to X
    Vector X(size);
    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i); 
        if(rowtype>=0 && rowtype<3) {
            X(i) = dV(rowid);
        } else if(rowtype == 3) {
            X(i) = dP(rowid);
        } else if(rowtype == 4) {
            //X(i) = dPi(rowid);            
        }

    }
    theSOE->setX(X);

    // timer.pause();
    // opserr<<"solving time for PFEMCompressiblesolver = "<<timer.getReal()<<"\n";

    return 0;
}

int PFEMCompressibleSolver::setSize()
{
    return 0;
}

int
PFEMCompressibleSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMCompressibleSolver::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMCompressibleSolver::setLinearSOE(PFEMCompressibleLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}
