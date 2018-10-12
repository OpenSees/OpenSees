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

#ifndef PFEMDiaLinSOE_h
#define PFEMDiaLinSOE_h


#include <PFEMLinSOE.h>

class PFEMDiaSolver;

class PFEMDiaLinSOE : public PFEMLinSOE
{
  public:
    PFEMDiaLinSOE(PFEMDiaSolver &theSolver);        
    PFEMDiaLinSOE();        

    virtual ~PFEMDiaLinSOE();

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual void zeroA(void);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);    

    virtual const ID& getDofID()const {return newDofID;}

    friend class PFEMDiaSolver;

private:    

    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize,
			  int Isize, int Psize, int Pisize);

private:

    cs *Gt, *G;
    Vector M, Mp;
    ID newDofID;
};

#endif

