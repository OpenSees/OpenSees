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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h
//
// Written: Minjie 
// Created: Sep 17 2012
//

#include <PFEMSolver_Mumps.h>
#include <PFEMLinSOE.h>
#include <iostream>
#include <cmath>
#include <Timer.h>
#include <mpi.h>

PFEMSolver_Mumps::PFEMSolver_Mumps(int r, int e, int h, int s)
    :PFEMSolver(), theSOE(0), sid(), pid(), myid(0),
     relax(r), err(e), host(h), sym(s)
{
    // mumps id
    sid.job = JOB_INIT;
    sid.par = 1;
    sid.sym = sym;
    sid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&sid);

    sid.irn = 0;
    sid.jcn = 0;

    pid.job = JOB_INIT;
    pid.par = 1;
    pid.sym = sym;
    pid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&pid);

    // get process id
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

}

PFEMSolver_Mumps::~PFEMSolver_Mumps()
{
    sid.job = JOB_END;
    dmumps_c(&sid);

    if(sid.irn != 0) delete [] sid.irn;
    if(sid.jcn != 0) delete [] sid.jcn;

    pid.job = JOB_END;
    dmumps_c(&pid);
}

int
PFEMSolver_Mumps::solve()
{
    // Timer timer;
    // timer.start();
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
    cs* Qt = theSOE->Qt;
    Vector& Mhat = theSOE->Mhat;
    Vector& Mf = theSOE->Mf;
    Vector& X = theSOE->X;
    Vector& B = theSOE->B;
    ID& dofType = theSOE->dofType;
    ID& dofID = theSOE->dofID;

    int Msize = M->n;
    int Isize = Git->n;
    int Ssize = Msize-Isize;
    int Fsize = Mf.Size();
    int Psize = L->n;
    int Pisize = Mhat.Size();
    int size = X.Size();

    // numeric LU factorization of M
    // call mumps for factorization
    MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Isize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ssize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Fsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Psize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pisize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(Msize > 0) {
	sid.job = JOB_FACTORIZATION;
	dmumps_c(&sid);
    }
    
    if(sid.info[0] != 0) {
	opserr<<"WARNING: failed to factorize -- PFEMSolver_Mumps::solve\n";
	return -1;
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"factorization time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();
    
    // structure and interface predictor : deltaV1 = M^{-1} * rsi
    Vector deltaV1;
    if(Msize > 0) {
	if(myid == 0) {
	    deltaV1.resize(Msize);
	    deltaV1.Zero();
	
	    // rsi
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 2) {
		    deltaV1(rowid+Ssize) = B(i);   // rsi
		} else if(rowtype == 0) {
		    deltaV1(rowid) = B(i);         // rsi
		}
	    }
	    double* deltaV1_ptr = &deltaV1(0);
	    sid.rhs = deltaV1_ptr;
	    sid.nrhs = 1;
	}
	sid.job = JOB_SOLUTION;
	dmumps_c(&sid);
	if(sid.info[0] != 0) {
	    opserr<<"WARNING: failed to solve predictor -- PFEMSolver_Mumps::solve\n";
	    return -1;
	}
    }
    
    // fluid predictor: deltaVf1 = Mf^{-1} * rf
    Vector deltaVf1;
    if(Fsize > 0) {
	if(myid == 0) {
	    deltaVf1.resize(Fsize);
	    deltaVf1.Zero();
	    
	    // rf
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 1) {
		    if(Mf(rowid) == 0) {
			opserr<<"WANING: Zero Mf at location "<<rowid<<" ";
			opserr<<" - PFEMLinSOE::solve()\n";
			return -1;
		    }
		    deltaVf1(rowid) = B(i)/Mf(rowid);         // rf
		}
	    }
	}
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"dVs1 = "<<deltaV1.Norm()<<"\n";
    // 	opserr<<"dVf1 = "<<deltaVf1.Norm()<<"\n";
    // 	opserr<<"predictor time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();

    // Mi^{-1}, Msi^{-1}
    cs* invMi1 = cs_spalloc(Isize, Isize, 1, 1, 1);
    cs* invMsi1 = cs_spalloc(Ssize, Isize, 1, 1, 1);
    if (Isize > 0) {
	Vector eyes(Msize);
	double* eyes_ptr = &eyes(0);

	for (int j=0; j<Isize; j++) {
	    if (myid == 0) {
		// rhs
		eyes.Zero();
		eyes(j+Ssize) = 1.0;
		
		// M^{-1}*eyes
		sid.rhs = eyes_ptr;
		sid.nrhs = 1;		
	    }
	    sid.job = JOB_SOLUTION;
	    dmumps_c(&sid);
	    if(sid.info[0] != 0) {
		opserr<<"WARNING: failed to solve invMi -- PFEMSolver_Mumps::solve\n";
		return -1;
	    }
	    if (myid == 0) {
		// copy
		for (int i=0; i<Msize; i++) {
		    if (eyes(i) != 0.0) {
			if (i>=Ssize) {
			    cs_entry(invMi1,i-Ssize,j,eyes(i));
			} else {
			    cs_entry(invMsi1,i,j,eyes(i));
			}
		    }
		}

	    }
	}

    }
    cs* invMi = cs_compress(invMi1);
    cs* invMsi = cs_compress(invMsi1);
    cs_spfree(invMi1);
    cs_spfree(invMsi1);
    opserr<<"invMi = "<<cs_norm(invMi)<<"\n";

    // Gi, Mf^{-1}*Gf
    cs* Gi = 0;
    cs* Gf = 0;
    if(Fsize > 0) {
	if(myid == 0) {
	    Gi = cs_transpose(Git, 1);
	    Gf = cs_transpose(Gft, 1);
	    for(int j=0; j<Psize; j++) {
		for(int k=Gf->p[j]; k<Gf->p[j+1]; k++) {
		    int i = Gf->i[k];
		    double& x = Gf->x[k];
		    x /= Mf(i);
		}
	    }
	}
    }

    // solve for pressure
    Vector deltaP;
    if(Psize>0) {
	// assembled format
	ICNTL(pid,5,0);

	// input matrix is centralized on the host
	ICNTL(pid,18,0);

	// workspace relaxation: 20%
	ICNTL(pid,14,relax);

	// dense right hand side
	ICNTL(pid,20,0);

	// centralized right hand side
	ICNTL(pid,21,0);

	// error messages
	ICNTL(pid,1,err);
	ICNTL(pid,2,err);
	ICNTL(pid,3,err);
	ICNTL(pid,4,err);

	// get S
	cs* S = 0;
	if(myid == 0) {
	    deltaP.resize(Psize);
	    deltaP.Zero();
	    double* deltaP_ptr = &deltaP(0);

	    // Gft*deltaVf1
	    if(Fsize > 0) {
		double* deltaVf1_ptr = &deltaVf1(0);
		cs_gaxpy(Gft, deltaVf1_ptr, deltaP_ptr);
	    }

	    // Git*deltaVi1
	    if(Isize > 0) {
		double* deltaVi1_ptr = &deltaV1(0) + Ssize;
		cs_gaxpy(Git, deltaVi1_ptr, deltaP_ptr);
	    }

	    // rp-Git*deltaVi1-Gft*deltaVf1
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 3) {             // pressure
		    deltaP(rowid) = B(i)-deltaP(rowid);   // rp-Git*deltaVi1-Gft*deltaVf1
		}
	    }

	    // S = L + Git*Mi{-1}*Gi + Gft*Mf{-1}*Gf
	    if(Isize > 0) {
		cs* S1 = cs_multiply(Git, invMi);
		S = cs_multiply(S1, Gi);
		cs_spfree(S1);
	    }
	    if(Fsize > 0) {
		cs* S1 = cs_multiply(Gft, Gf);
		if(S == 0) {
		    S = S1;
		} else {
		    cs* S2 = cs_add(S, S1, 1.0, 1.0);
		    cs_spfree(S);
		    cs_spfree(S1);
		    S = S2;
		}
	    }

	    // S
	    if(S == 0) {
		S = L;
	    } else {
		cs* S1 = cs_add(S, L, 1.0, 1.0);
		cs_spfree(S);
		S = S1;
	    }

	    // set S
	    pid.n = S->n;
	    pid.nz = S->nzmax;
	    pid.a = S->x;
	    pid.irn = new int[pid.nz];
	    pid.jcn = new int[pid.nz];
	    for(int j=0; j<pid.n; j++) {
		for(int k=S->p[j]; k<S->p[j+1]; k++) {
		    pid.irn[k] = S->i[k]+1;
		    pid.jcn[k] = j+1;
		}
	    }

	    pid.rhs = deltaP_ptr;
	    pid.nrhs = 1;
	}
	// timer.pause();
	// if(myid == 0) {
	//     opserr<<"pressure setup time = "<<timer.getReal()<<"\n";
	// }
	
	// timer.start();
	
	// solve
	pid.job = JOB_SOLVE;
	dmumps_c(&pid);
	if(pid.info[0] != 0) {
	    opserr<<"WARNING: failed to solve pressure -- PFEMSolver_Mumps::solve\n";
	    return -1;
	}

	// release
	if (myid == 0) {
	    if(S != L) cs_spfree(S);
	    delete [] pid.irn;
	    delete [] pid.jcn;
	}
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"dP = "<<deltaP.Norm()<<"\n";
    // 	opserr<<"pressure  time = "<<timer.getReal()<<"\n";
    // }

    //timer.start();
    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    Vector deltaV;
    if(myid == 0) {
	deltaV.resize(Msize);
	deltaV.Zero();
    }
    if(Isize > 0) {
	if(myid == 0) {
	    // Gi*deltaP
	    Vector Gip(Isize);
	    double* Gip_ptr = &Gip(0);
	    if(Psize > 0) {
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Gi, deltaP_ptr, Gip_ptr);

		// Msi^{-1}*Gi*deltaP
		if(Ssize > 0) {
		    Vector vs(Ssize);
		    double* vs_ptr = &vs(0);
		    cs_gaxpy(invMsi, Gip_ptr, vs_ptr);
		    for(int i=0; i<Ssize; i++) {
			deltaV(i) += vs(i);
		    }
		}

		// Mi^{-1}*Gi*deltaP
		Vector vi(Isize);
		double* vi_ptr = &vi(0);
		cs_gaxpy(invMi, Gip_ptr, vi_ptr);
		for(int i=0; i<Isize; i++) {
		    deltaV(i+Ssize) = vi(i);
		}

	    }
	    cs_spfree(Gi);
	}
    }
    if(myid == 0) {
	deltaV += deltaV1;
    }
    cs_spfree(invMi);
    cs_spfree(invMsi);
    opserr<<"dVs = "<<deltaV.Norm()<<"\n";

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    Vector deltaVf(Fsize);
    if(Fsize > 0) {
        if(Psize > 0) {
	    if(myid == 0) {
		deltaVf.resize(Fsize);
		deltaVf.Zero();
		double* deltaVf_ptr = &deltaVf(0);
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);
		deltaVf += deltaVf1;
		cs_spfree(Gf);
	    }
        }
    }
    opserr<<"dVf = "<<deltaVf.Norm()<<"\n";

    // deltaPi = Mhatd^{-1}*(Qt*deltaP - rpi)
    Vector deltaPi;
    if(Pisize > 0) {
	if(myid == 0) {
	    deltaPi.resize(Pisize);
	    deltaPi.Zero();
	    
	    // Qt*deltaP
	    if(Psize > 0) {
		double* deltaPi_ptr = &deltaPi(0);
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Qt, deltaP_ptr, deltaPi_ptr);   // Qt*deltaP
	    }

	    // rpi
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype != 4) continue;     // rpi
		if(Mhat(rowid) == 0.0) {
		    opserr<<"Zero Mhat at location "<<rowid<<" ";
		    opserr<<" - PFEMLinSOE::solve()\n";
		    return -1;
		}
		deltaPi(rowid) = (B(i)-deltaPi(rowid))/Mhat(rowid); // rpi
	    }
	}
        
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"corrector  time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();
    // copy to X
    if(myid == 0) {
	X.Zero();

	for(int i=0; i<size; i++) {            // row
	    int rowtype = dofType(i);          // row type
	    int rowid = dofID(i);
	    if(rowtype == 0) {
		X(i) = deltaV(rowid);            
	    } else if(rowtype == 2) {
		X(i) = deltaV(rowid+Ssize);
	    } else if(rowtype == 1) {
		X(i) = deltaVf(rowid);
	    } else if(rowtype == 3) {
		X(i) = deltaP(rowid);
	    } else if(rowtype == 4) {
		X(i) = deltaPi(rowid);
	    }

	}
    }

    return 0;
}

int PFEMSolver_Mumps::setSize()
{
    // Timer timer;
    // timer.start();

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

    // check M size
    cs* M = theSOE->M;
    int Msize = M->n;
    MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(Msize <= 0) return 0;

    // host
    if(myid == 0) {
	sid.n = M->n;
	sid.nz = M->nzmax;
	sid.a = M->x;
	if(sid.irn != 0) delete [] sid.irn;
	if(sid.jcn != 0) delete [] sid.jcn;
	sid.irn = new int[sid.nz];
	sid.jcn = new int[sid.nz];

	for(int j=0; j<sid.n; j++) {
	    for(int k=M->p[j]; k<M->p[j+1]; k++) {
		sid.irn[k] = M->i[k]+1;
		sid.jcn[k] = j+1;
	    }
	}
    }

    // call mumps
    sid.job = JOB_ANALYSIS;
    dmumps_c(&sid);
    if(sid.info[0] != 0) {
	opserr<<"WARNING: failed to analyze -- PFEMSolver_Mumps::setSize\n";
	return -1;
    }

    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"analysis time = "<<timer.getReal()<<"\n";
    // }
    return 0;
}

int
PFEMSolver_Mumps::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_Mumps::recvSelf(int ctag,
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
{
    // nothing to do
    return 0;
}


int 
PFEMSolver_Mumps::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}

void
PFEMSolver_Mumps::ICNTL(DMUMPS_STRUC_C& id, int I, int val)
{
    id.icntl[I-1] = val;
}


