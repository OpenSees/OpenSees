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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver_Mumps.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleSolver_Mumps.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
#include <PFEMCompressibleSolver_Mumps.h>
#include <PFEMCompressibleLinSOE.h>
#include <iostream>
#include <cmath>
#include <Timer.h>
#include <mpi.h>

PFEMCompressibleSolver_Mumps::PFEMCompressibleSolver_Mumps(int r, int e, int s)
    :PFEMCompressibleSolver(), theSOE(0),
     sid(), myid(0), relax(r), err(e)
{
    // mumps id
    sid.job = JOB_INIT;
    sid.par = 1;
    sid.sym = s;
    sid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&sid);

    // get process id
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
}

PFEMCompressibleSolver_Mumps::~PFEMCompressibleSolver_Mumps()
{
    sid.job = JOB_END;
    dmumps_c(&sid);
}

int
PFEMCompressibleSolver_Mumps::solve()
{
    // Timer timer;
    if (myid == 0) {
	// timer.start();	
    }

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

    MPI_Bcast(&Vsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Psize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(Vsize<=0 || Psize<=0) {
        opserr<<"WARNING: Fsize or Psize or Pisize <= 0 -- ";
        opserr<<"PFEMCompressibleSolver_Mumps::solve\n";
        return -1;
    }
    if (myid == 0) {
	// timer.pause();
	// opserr<<"setup time = "<<timer.getReal()<<"\n";
	// timer.start();
    }

    // velocity full tangent
    cs* H = 0;
    Vector dP, dV;
    if (myid == 0) {

	// G*Mp^{-1}*Gt
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
	H = cs_add(M, GMpGt, 1.0, -1.0);

	// G*Mp^{-1}*dp
	dP.resize(Psize);
	dP.Zero();
	double* dP_ptr = &dP(0);
	for(int i=0; i<size; i++) {        // row

	    int rowtype = dofType(i);      // row type
	    int rowid = dofID(i);          // row id
	    if(rowtype == 3) {
		dP(rowid) = B(i);          // rp
	    }
	}

	dV.resize(Vsize);
	dV.Zero();
	double* dV_ptr = &dV(0);
	cs_gaxpy(GMp, dP_ptr, dV_ptr);     // GMp*rp
	for(int i=0; i<size; i++) {        // row
	    int rowtype = dofType(i);      // row type
	    int rowid = dofID(i);          // row id
	    if(rowtype < 3 && rowtype>=0) {
		dV(rowid) = B(i)-dV(rowid);         // rm-GMp*rp
	    }
	}

	// clean
	cs_spfree(GMpGt);

	// timer.pause();
	// opserr<<"velocity tangent time = "<<timer.getReal()<<"\n";
	// timer.start();
    }
    
    // solve velocity
    if (myid == 0) {
	// timer.start();	
    }
    
    // assembled format
    ICNTL(sid,5,0);

    // input matrix is centralized on the host
    ICNTL(sid,18,0);

    // workspace relaxation: 20%
    if(relax <= 0) {
	relax = 20;
    }
    ICNTL(sid,14,relax);

    // dense right hand side
    ICNTL(sid,20,0);

    // centralized right hand side
    ICNTL(sid,21,0);

    // No error messages
    if(err < 0) err = 0;
    ICNTL(sid,1,err);
    ICNTL(sid,2,err);
    ICNTL(sid,3,err);
    ICNTL(sid,4,err);

    // host
    if(myid == 0) {

	// lhs
	sid.n = H->n;
	sid.nz = H->nzmax;
	sid.a = H->x;
	sid.irn = new int[sid.nz];
	sid.jcn = new int[sid.nz];

	for(int j=0; j<sid.n; j++) {
	    for(int k=H->p[j]; k<H->p[j+1]; k++) {
		sid.irn[k] = H->i[k]+1;
		sid.jcn[k] = j+1;
	    }
	}

	// rhs
	double* dV_ptr = &dV(0);
	sid.rhs = dV_ptr;
	sid.nrhs = 1;
    }

    // call mumps
    sid.job = JOB_SOLVE;
    dmumps_c(&sid);

    // error
    if (myid == 0) {
	if(sid.info[0] != 0) {
	    opserr<<"WARNING: failed to solve velocity with error code "<<sid.info[0];
	    opserr<<" -- PFEMCompressibleSolver_Mumps::setSize\n";
	    return -1;
	}

	cs_spfree(H);
	if(sid.irn != 0) delete [] sid.irn;
	if(sid.jcn != 0) delete [] sid.jcn;
	// timer.pause();
	// opserr<<"dV = "<<dV.Norm()<<"\n";
	// opserr<<"solve velocity time = "<<timer.getReal()<<"\n";
	// timer.start();
    }
 
    // solve pressure
    if (myid == 0) {
	dP *= -1;                          // -rp
	double* dV_ptr = &dV(0);
	double* dP_ptr = &dP(0);
	cs_gaxpy(Gt, dV_ptr, dP_ptr);      // Gt*dV - rp
	for(int i=0; i<Psize; i++) {       // (rp - Gt*dV) / Mp
	    dP(i) /= -Mp(i);
	}
	// timer.pause();
	// opserr<<"solve pressure time = "<<timer.getReal()<<"\n";
    }

    // copy to X
    if (myid == 0) {
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
    }
 

    return 0;
}

int PFEMCompressibleSolver_Mumps::setSize()
{
    return 0;
}

int
PFEMCompressibleSolver_Mumps::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMCompressibleSolver_Mumps::recvSelf(int ctag,
				       Channel &theChannel, 
				       FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMCompressibleSolver_Mumps::setLinearSOE(PFEMCompressibleLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}

void
PFEMCompressibleSolver_Mumps::ICNTL(DMUMPS_STRUC_C& id, int I, int val)
{
    id.icntl[I-1] = val;
}
