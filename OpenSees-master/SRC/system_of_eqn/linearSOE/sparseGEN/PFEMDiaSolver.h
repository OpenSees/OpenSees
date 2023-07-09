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
                                                                        
#ifndef PFEMDiaSolver_h
#define PFEMDiaSolver_h

#include <LinearSOESolver.h>
#include <Vector.h>


class PFEMDiaLinSOE;

class PFEMDiaSolver : public LinearSOESolver
{
public:
    PFEMDiaSolver();
    virtual ~PFEMDiaSolver();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMDiaLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:
    
    PFEMDiaLinSOE* theSOE;
};

#endif

