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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMCompressibleLinSOE.h,v $

// Written: Minjie Zhu
// Created: August 2012
//
// Description: This file contains the class definition for PFEMCompressibleLinSOE
// PFEMCompressibleLinSOE is a subclass of SparseGenColLinSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. It solves the equations using the Fractional Step Method in PFEM.
//
// What: "@(#) PFEMCompressibleLinSOE.h, revA"

#ifndef PFEMCompressibleLinSOE_h
#define PFEMCompressibleLinSOE_h


#include <PFEMLinSOE.h>

class PFEMCompressibleSolver;
class PFEMCompressibleSolver_Mumps;

class PFEMCompressibleLinSOE : public PFEMLinSOE
{
  public:
    PFEMCompressibleLinSOE(PFEMCompressibleSolver &theSolver, bool exp=false);        
    PFEMCompressibleLinSOE();        

    virtual ~PFEMCompressibleLinSOE();

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual void zeroA(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    virtual const ID& getDofID()const {return newDofID;}

    friend class PFEMCompressibleSolver;
    friend class PFEMCompressibleSolver_Mumps;

private:    

    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize);

private:

    cs* M, *Gt, *G;
    Vector Mp;
    ID newDofID;
};

#endif

