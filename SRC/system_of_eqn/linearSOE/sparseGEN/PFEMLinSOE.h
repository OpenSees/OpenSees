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

#ifndef PFEMLinSOE_h
#define PFEMLinSOE_h


#include <LinearSOE.h>
#include <OPS_Stream.h>
#include <Vector.h>
#include <ID.h>
extern "C" {
#include <cs.h>
}
class PFEMSolver;

class PFEMLinSOE : public LinearSOE
{
  public:
    PFEMLinSOE(PFEMSolver &theSolver);
    PFEMLinSOE(int classTag);
    PFEMLinSOE();
    PFEMLinSOE(PFEMSolver& theSolver, int classTag);

    virtual ~PFEMLinSOE();

    virtual int solve(void);

    virtual int getNumEqn(void) const;
    virtual int setSize(Graph& theGraph);
    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual int addB(const Vector &, const ID &, double fact = 1.0);
    virtual int setB(const Vector &, double fact = 1.0);

    virtual void zeroA(void);
    virtual void zeroB(void);

    virtual const Vector &getX(void);
    virtual const Vector &getB(void);
    virtual double normRHS(void);

    virtual void setX(int loc, double value);
    virtual void setX(const Vector &x);
    int setPFEMSolver(PFEMSolver& newSolver);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    virtual const ID& getDofType()const {return dofType;}
    virtual const ID& getDofID()const {return dofID;}

    friend class PFEMSolver;
    friend class PFEMSolver_Mumps;
    friend class PFEMSolver_Umfpack;
    friend class PFEMSolver_Laplace;
    friend class PFEMSolver_LumpM;

    virtual bool isFluidID(const ID& id) const;
    void saveK(OPS_Stream& output);

private:

    virtual int setDofIDs(int size,int& Ssize, int&Fsize, int& Isize,int& Psize,int& Pisize);
    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize);

private:

    cs* M, *Gft, *Git, *L, *Qt;
    Vector X, B, Mhat, Mf;
    ID dofType, dofID;
};

#endif
