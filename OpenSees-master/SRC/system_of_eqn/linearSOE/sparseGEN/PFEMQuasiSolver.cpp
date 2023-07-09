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
                                                                        
// $Revision$
// $Date$
// $URL$

// Written: Minjie 
// Description: Solve PFEM with quasi-incompressible form, no global matrix operations
//

#include <PFEMQuasiSolver.h>
#include <PFEMQuasiLinSOE.h>
#include <iostream>
#include <cmath>

void* OPS_PFEMQuasiSolver()
{
    PFEMQuasiSolver* theSolver = new PFEMQuasiSolver();
    return new PFEMQuasiLinSOE(*theSolver);
}

PFEMQuasiSolver::PFEMQuasiSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMQuasiSolver), theSOE(0)
{
    // set default control parameters
    umfpack_di_defaults(Control);
}

PFEMQuasiSolver::~PFEMQuasiSolver()
{
    
}

int
PFEMQuasiSolver::solve()
{
    // M and Mp
    cs* M = theSOE->M;
    cs* Mp = theSOE->Mp;
    cs* Gt = theSOE->Gt;

    // dofs
    const Vector& B = theSOE->getB();
    const ID& dofType = theSOE->getDofType();
    const ID& dofID = theSOE->getDofID();
    
    int Vsize = M->n;
    int Psize = Mp->n;
    int size = B.Size();

    if(Vsize<=0 || Psize<=0) {
        opserr<<"WARNING: Fsize or Psize or Pisize <= 0 -- ";
        opserr<<"PFEMQuasiSolver::solve\n";
        return -1;
    }

    // rv
    Vector dV(Vsize);
    double* dV_ptr = &dV(0);
    for(int i=0; i<size; i++) {        // row
	
        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype < 3 && rowtype>=0) {
            dV(rowid) = B(i);         // rv
        }
    }

    // M
    reorder(M);

    // symbolic analysis of M
    void* Symbolic = 0;
    int status = 0;
    status = umfpack_di_symbolic(Vsize,Vsize,M->p,M->i,M->x,
				 &Symbolic,Control,Info);
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: M symbolic analysis returns "<<status;
	opserr <<" -- PFEMQuasiSolver::solve\n";
	return -1;
    }

    // numerical analysis of M
    void* Numeric = 0;
    status = umfpack_di_numeric(M->p,M->i,M->x,Symbolic,&Numeric,Control,Info);
    if (status!=UMFPACK_OK) {
    	opserr<<"WARNING: M numeric analysis returns "<<status;
	opserr<<" -- PFEMQuasiSolver::solve\n";
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
    	return -1;
    }

    // solve velocity
    Vector x(Vsize);
    double* x_ptr = &x(0);
    status = umfpack_di_solve(UMFPACK_A,M->p,M->i,M->x,x_ptr,dV_ptr,
			      Numeric,Control,Info);
    dV = x;
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: M solving returns "<<status;
	opserr<<" -- PFEMQuasiSolver::solve\n";
	if (Numeric != 0) {
	    umfpack_di_free_numeric(&Numeric);
	}
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
	return -1;
    }

    // clean up
    if (Numeric != 0) {
	umfpack_di_free_numeric(&Numeric);
	Numeric = 0;
    }
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
	Symbolic = 0;
    }

    // rp
    Vector dP(Psize);
    double* dP_ptr = &dP(0);
    for(int i=0; i<size; i++) {        // row

        int rowtype = dofType(i);      // row type
        int rowid = dofID(i);          // row id
        if(rowtype == 3) {
            dP(rowid) = -B(i);          // -rp
        }
    }
    cs_gaxpy(Gt, dV_ptr, dP_ptr);   // Gt*dV - rp

    // Mp
    reorder(Mp);
    
    // Symoblic of Mp
    status = umfpack_di_symbolic(Psize,Psize,Mp->p,Mp->i,Mp->x,
				 &Symbolic,Control,Info);
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: Mp symbolic analysis returns "<<status;
	opserr <<" -- PFEMQuasiSolver::solve\n";
	return -1;
    }

    // Numeric of Mp
    status = umfpack_di_numeric(Mp->p,Mp->i,Mp->x,Symbolic,&Numeric,
				Control,Info);
    if (status!=UMFPACK_OK) {
    	opserr<<"WARNING: Mp numeric analysis returns "<<status;
	opserr<<" -- PFEMQuasiSolver::solve\n";
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
    	return -1;
    }

    // solve presssure: -dP
    Vector y(Psize);
    double* y_ptr = &y(0);
    status = umfpack_di_solve(UMFPACK_A,Mp->p,Mp->i,Mp->x,y_ptr,dP_ptr,
			      Numeric,Control,Info);
    dP = y;
    dP *= -1;
    if (status!=UMFPACK_OK) {
	opserr<<"WARNING: Mp solving returns "<<status;
	opserr<<" -- PFEMQuasiSolver::solve\n";
	if (Numeric != 0) {
	    umfpack_di_free_numeric(&Numeric);
	}
	if (Symbolic != 0) {
	    umfpack_di_free_symbolic(&Symbolic);
	}
	return -1;
    }

    // clean up
    if (Numeric != 0) {
	umfpack_di_free_numeric(&Numeric);
	Numeric = 0;
    }
    if (Symbolic != 0) {
	umfpack_di_free_symbolic(&Symbolic);
	Symbolic = 0;
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

    return 0;
}

int PFEMQuasiSolver::setSize()
{
    return 0;
}

int
PFEMQuasiSolver::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMQuasiSolver::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMQuasiSolver::setLinearSOE(PFEMQuasiLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}

void
PFEMQuasiSolver::reorder(cs* A)
{
    int* Ap = A->p;
    int* Ai = A->i;
    double* Ax = A->x;

    for (int j=0; j<A->n; j++) {
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
}
