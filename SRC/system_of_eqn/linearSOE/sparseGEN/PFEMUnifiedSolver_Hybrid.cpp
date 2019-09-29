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
// $Date: 2014-09-06 13:53:05 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMUnifiedSolver_Hybrid.h,v $

// Written: Minjie Zhu

#include <PFEMUnifiedSolver_Hybrid.h>
#include <PFEMGeneralLinSOE.h>
#include <iostream>
#include <cmath>
#include <vector>
#include <Timer.h>
#include <Matrix.h>
#include <mpi.h>

PFEMUnifiedSolver_Hybrid::PFEMUnifiedSolver_Hybrid(int r, int e, int h, int s)
    :LinearSOESolver(SOLVER_TAGS_PFEMUnifiedSolver_Hybrid), theSOE(0), id(), myid(0), handle(),
     relax(r), err(e), host(h), sym(s)
{
    // mumps id
    id.job = JOB_INIT;
    id.par = 1;
    id.sym = sym;
    id.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&id);

    // get process id
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    // create cusolver handle
    if(cusolverDnCreate(&handle) != CUSOLVER_STATUS_SUCCESS) {
	opserr<<"WARNING: failed to create cusolver handle--PFEMUnifiedSolver_Hybrid\n";
	exit(-1);
    }
}

PFEMUnifiedSolver_Hybrid::~PFEMUnifiedSolver_Hybrid()
{
    id.job = JOB_END;
    dmumps_c(&id);
    cusolverDnDestroy(handle);
}

int
PFEMUnifiedSolver_Hybrid::solve()
{
    Timer timer;
    timer.start();
    
    // factorization
    Value schur;
    Index psc, irn, jcn;
    if(myid == 0) {
	Index& rowInd = theSOE->rowInd;
	Index& colPtr = theSOE->colPtr;
	Value& A = theSOE->A;
	const ID& dofType = theSOE->getDofType();
	const ID& dofID = theSOE->newDofID;
	
	// matrix
	id.n = theSOE->N;
	id.nz = A.size();
	if(id.n<=0 || id.nz<=0) return 0;
	id.a = &A[0];
	irn.resize(id.nz,1);
	jcn.resize(id.nz,1);
	for(int j=0; j<id.n; j++) {
	    for(int k=colPtr[j]; k<colPtr[j+1]; k++) {
		// increment row index for mumps fortran indexing
		irn[k] = rowInd[k]+1;
		jcn[k] = j+1;
	    }
	}
	id.irn = &irn[0];
	id.jcn = &jcn[0];

	// schur
	for(int i=0; i<dofType.Size(); i++) {
	    if(dofType(i) == 3) {
		psc.push_back(dofID(i)+1);
	    }
	}
	id.size_schur = psc.size();
	id.listvar_schur = &psc[0];
	id.nprow = 1;
	id.npcol = 1;
	id.mblock = 100;
	id.nblock = 100;
	id.schur_lld = id.size_schur;
	schur.resize(id.size_schur*id.size_schur,0.0);
	id.schur = &schur[0];
    }

    // call mumps for factorization
    id.job = JOB_FACTORIZATION;
    dmumps_c(&id);

    if(id.info[0] != 0) {
	opserr<<"WARNING: failed to factorize - return "<<id.info[0];
	opserr<<" -- PFEMUnifiedSolver_Hybrid::solve ";
	opserr<<" from process "<<myid<<"\n";
	return -1;
    }

    timer.pause();
    if(myid == 0) {
	opserr<<"factorization time = "<<timer.getReal()<<"\n";
    }
    timer.start();

    // predictor
    Vector dV;
    if(myid == 0) {
	const ID& dofID = theSOE->newDofID;
	const Vector& B = theSOE->getB();
	const ID& dofType = theSOE->getDofType();
	dV.resize(id.n);
	for(int i=0; i<dofID.Size(); i++) {
	    int ind = dofID(i);
	    if(dofType(i) >= 0) {
		dV(ind) = B(i);
	    }
	}
	double* dV_ptr = &dV(0);
	id.rhs = dV_ptr;
	id.nrhs = 1;
    }

    ICNTL(26,0);
    id.job = JOB_SOLUTION;
    dmumps_c(&id);
    if(id.info[0] != 0) {
	opserr<<"WARNING: failed to solve predictor -- PFEMUnifiedSolver_Hybrid::solve\n";
	return -1;
    }

    timer.pause();
    if(myid == 0) {
	opserr<<"dV1 = "<<dV.Norm()<<"\n";
	opserr<<"predictor time = "<<timer.getReal()<<"\n";
    }
    timer.start();

    // pressure rhs
    Vector dP;
    if(myid == 0) {
	dP.resize(id.size_schur);
	dP.Zero();
	double* dP_ptr = &dP(0);
	id.redrhs = dP_ptr;
	id.lredrhs = id.size_schur;
	
	// reset right hand side
	const ID& dofID = theSOE->newDofID;
	const Vector& B = theSOE->getB();
	dV.resize(id.n);
	for(int i=0; i<dofID.Size(); i++) {
	    int id = dofID(i);
	    if(id >= 0) {
		dV(id) = B(i);
	    }
	}
	double* dV_ptr = &dV(0);
	id.rhs = dV_ptr;
	id.nrhs = 1;
    }
    ICNTL(26,1);
    id.job = JOB_SOLUTION;
    dmumps_c(&id);
    if(id.info[0] != 0) {
	opserr<<"WARNING: failed to get pressure rhs -- PFEMUnifiedSolver_Hybrid::solve\n";
	return -1;
    }

    // solve pressure using GPU
    if(myid == 0) {

	// cusolver handle
	cusolverStatus_t status;

	// host data
	double* h_S = id.schur;
	double* h_dP = &dP(0);

	// allocate on GPU
	double* d_S, *d_dP;
	cudaMalloc((void**) &d_S, id.size_schur*id.size_schur*sizeof(double));
	cudaMalloc((void**) &d_dP, id.size_schur*sizeof(double));

	// copy L and dP to GPU
	cudaMemcpy(d_S, h_S, id.size_schur*id.size_schur*sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(d_dP, h_dP, id.size_schur*sizeof(double), cudaMemcpyHostToDevice);

	// get workspace size
	int lda = id.size_schur, Lwork=0;
	status = cusolverDnDgetrf_bufferSize(handle,id.size_schur,id.size_schur,d_S,lda,&Lwork);
	if(status != CUSOLVER_STATUS_SUCCESS) {
	    opserr<<"WARNING: failed to get GPU workspace size--PFEMUnifiedSolver_Hybrid::solve\n";
	    return -1;
	}

	// set up workspace
	double* d_workspace;
	int* d_ipiv, *d_info;
	cudaMalloc((void**)&d_workspace, Lwork*sizeof(double));
	cudaMalloc((void**)&d_ipiv, id.size_schur*sizeof(int));
	cudaMalloc((void**)&d_info, sizeof(int));

	// LU decomposition
	status = cusolverDnDgetrf(handle,id.size_schur,id.size_schur,
				  d_S,lda,d_workspace,d_ipiv,d_info);
	if(status != CUSOLVER_STATUS_SUCCESS) {
	    opserr<<"WARNING: failed to do LU in GPU for pressure--PFEMUnifiedSolver_Hybrid::solve\n";
	    return -1;
	}

	// solve
	int ldb = lda, nrhs = 1;
	cublasOperation_t trans = CUBLAS_OP_N;
	status = cusolverDnDgetrs(handle,trans,lda,nrhs,d_S,lda,d_ipiv,d_dP,ldb,d_info);
	if(status != CUSOLVER_STATUS_SUCCESS) {
	    opserr<<"WARNING: failed to solve pressure in GPU--PFEMUnifiedSolver_Hybrid::solve\n";
	    return -1;
	}

	// check results
	cudaMemcpy(h_dP, d_dP, id.size_schur*sizeof(double), cudaMemcpyDeviceToHost);

	// clear
	cudaFree(d_S);
	cudaFree(d_dP);
	cudaFree(d_workspace);
	cudaFree(d_ipiv);
	cudaFree(d_info);
    }

    timer.pause();
    if(myid == 0) {
	opserr<<"dP = "<<dP.Norm()<<"\n";
	opserr<<"pressure time = "<<timer.getReal()<<"\n";
    }
    timer.start();

    // corrector
    ICNTL(26,2);
    id.job = JOB_SOLUTION;
    dmumps_c(&id);
    if(id.info[0] != 0) {
	opserr<<"WARNING: failed to solve corrector -- PFEMUnifiedSolver_Hybrid::solve\n";
	return -1;
    }

    timer.pause();
    if(myid==0) {
	opserr<<"corrector time = "<<timer.getReal()<<"\n";
    }
    
    // copy to X
    if(myid == 0) {
	const ID& dofID = theSOE->newDofID;
	Vector X = theSOE->getX();
	for(int i=0; i<dofID.Size(); i++) {
	    int id = dofID(i);
	    if(id >= 0) {
		X(i) = dV(id);
	    }
	}
	theSOE->setX(X);
    }

    return 0;
}

int PFEMUnifiedSolver_Hybrid::setSize()
{
    Timer timer;
    timer.start();

    // assembled format
    ICNTL(5,0);

    // input matrix is centralized on the host
    ICNTL(18,0);

    // non parallel ordering
    ICNTL(28,1);

    // workspace relaxation: 20%
    if(relax <= 0) {
	relax = 20;
    }
    ICNTL(14,relax);

    // Schur complement matrix
    ICNTL(19,2);

    // dense right hand side
    ICNTL(20,0);

    // centralized right hand side
    ICNTL(21,0);

    // No error messages
    if(err < 0) err = 0;
    ICNTL(1,err);
    ICNTL(2,err);
    ICNTL(3,err);
    ICNTL(4,err);

    // Pressure Schur Complement (PSC)
    Index psc, irn, jcn;
    if(myid == 0) {
	Index& rowInd = theSOE->rowInd;
	Index& colPtr = theSOE->colPtr;
	Value& A = theSOE->A;
	const ID& dofType = theSOE->getDofType();
	const ID& dofID = theSOE->newDofID;

	// matrix A
	id.n = theSOE->N;
	id.nz = A.size();
	if(id.n<=0 || id.nz<=0) return 0;
	id.a = &A[0];
	irn.resize(id.nz,1);
	jcn.resize(id.nz,1);
	for(int j=0; j<id.n; j++) {
	    for(int k=colPtr[j]; k<colPtr[j+1]; k++) {
		// increment row index for mumps fortran indexing
		irn[k] = rowInd[k]+1;
		jcn[k] = j+1;
	    }
	}
	id.irn = &irn[0];
	id.jcn = &jcn[0];

	// schur
	for(int i=0; i<dofType.Size(); i++) {
	    if(dofType(i) == 3) {
		psc.push_back(dofID(i)+1);
	    }
	}
	id.size_schur = psc.size();
	id.listvar_schur = &psc[0];
	id.nprow = 1;
	id.npcol = 1;
	id.mblock = 100;
	id.nblock = 100;
    }

    // call mumps
    id.job = JOB_ANALYSIS;
    dmumps_c(&id);

    if(id.info[0] != 0) {
	opserr<<"WARNING: failed to analyze -- PFEMUnifiedSolver_Hybrid::setSize\n";
	return -1;
    }

    timer.pause();
    if(myid == 0) {
	opserr<<"analysis time = "<<timer.getReal()<<"\n";
    }

    return 0;
}

int
PFEMUnifiedSolver_Hybrid::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMUnifiedSolver_Hybrid::recvSelf(int ctag,
				Channel &theChannel, 
				FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMUnifiedSolver_Hybrid::setLinearSOE(PFEMGeneralLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}


void
PFEMUnifiedSolver_Hybrid::ICNTL(int I, int val)
{
    id.icntl[I-1] = val;
}
