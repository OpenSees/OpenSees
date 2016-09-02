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
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMUnifiedLinSOE.h,v $

// Written: Minjie Zhu
// Created: September 2014
//
// Description: This file contains the class definition for PFEMUnifiedLinSOE
// PFEMUnifiedLinSOE is a subclass of SparseGenColLinSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the 
// matrix A. It solves the equations using the Fractional Step Method in PFEM.
//
// What: "@(#) PFEMUnifiedLinSOE.h, revA"

#ifndef PFEMUnifiedLinSOE_h
#define PFEMUnifiedLinSOE_h


#include <PFEMLinSOE.h>

class PFEMUnifiedSolver;

class PFEMUnifiedLinSOE : public PFEMLinSOE
{
  public:
    PFEMUnifiedLinSOE(PFEMUnifiedSolver &theSolver);
    PFEMUnifiedLinSOE();        

    virtual ~PFEMUnifiedLinSOE();

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual void zeroA(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    virtual const ID& getDofID()const {return newDofID;}

    friend class PFEMUnifiedSolver;

private:    

    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize);

private:

    cs* M, *G, *Gt, *Mp;
    ID newDofID;
};

#endif

