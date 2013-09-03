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

PFEMCompressibleSolver::PFEMCompressibleSolver()
    :LinearSOESolver(SOLVER_TAGS_PFEMCompressibleSolver), theSOE(0)
{
}

PFEMCompressibleSolver::~PFEMCompressibleSolver()
{

}

int
PFEMCompressibleSolver::solve()
{
    cs* M = theSOE->M;
    cs* Gt = theSOE->Gt;
    const Vector& B = theSOE->getB();
    Vector& Mp = theSOE->Mp;
    const ID& dofType = theSOE->getDofType();
    const ID& dofID = theSOE->getDofID();
    bool expl = theSOE->expl;
    
    int Vsize = M->n;
    int Psize = Mp.Size();
    int size = B.Size();

    if(Vsize<=0 || Psize<=0) {
        opserr<<"WARNING: Fsize or Psize or Pisize <= 0 -- ";
        opserr<<"PFEMCompressibleSolver::solve\n";
        return -1;
    }

    // solve velocity
    cs* GMp = cs_transpose(Gt,1);
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
    cs* H = cs_add(M, GMpGt, 1.0, 1.0);

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
            dV(rowid) += B(i);         // rm+GMp*rp
        }
    }

    cs_lusol(3, H, dV_ptr, 1e-16);
    cs_spfree(GMp);
    cs_spfree(GMpGt);
    cs_spfree(H);

    // solve pressure
    dP *= -1;                          // -rp
    cs_gaxpy(Gt, dV_ptr, dP_ptr);      // Gt*dV - rp
    for(int i=0; i<Psize; i++) {       // (rp - Gt*dV) / Mp
        dP(i) /= -Mp(i);
    }

    // pressure gradient
    // Vector dPi(Pisize);
    // double* dPi_ptr = &dPi(0);
    // cs_gaxpy(Qt, dP_ptr, dPi_ptr);
    // for(int i=0; i<size; i++) {        // row
    //     int rowtype = dofType(i);      // row type
    //     int rowid = dofID(i);          // row id
    //     if(rowtype == 4) {
    //         if(Mhat(rowid) == 0) {
    //             opserr<<"WARNING: Mhat is zero at "<<rowid<<"\n";
    //             return -1;
    //         }
    //         dPi(rowid) = (B(i) - dPi(rowid))/Mhat(rowid);         // (rpi-Qt*dP) / Mhat
    //     }
    // }

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
