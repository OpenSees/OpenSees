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
// $Date: 2015-03-27 9:29:32 $

// Written: Minjie Zhu
// Created: March 2015
//
// Description: This file contains the class definition for PFEMGeneralLinSOE
// PFEMGeneralLinSOE is a subclass of LinearSOE. It stores the matrix equation
// Ax=b using the sparse column-compacted storage scheme for storing the
// matrix A.
//

#ifndef PFEMGeneralLinSOE_h
#define PFEMGeneralLinSOE_h


#include <PFEMLinSOE.h>
#include <Vector.h>
#include <ID.h>
#include <vector>

class PFEMUnifiedSolver_Hybrid;

class PFEMGeneralLinSOE : public PFEMLinSOE
{
public:
    PFEMGeneralLinSOE(PFEMUnifiedSolver_Hybrid &theSolver);
    PFEMGeneralLinSOE();

    virtual ~PFEMGeneralLinSOE();

    virtual int addA(const Matrix &, const ID &, double fact = 1.0);
    virtual void zeroA(void);

    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    friend class PFEMUnifiedSolver_Hybrid;

private:

    virtual int setMatIDs(Graph& theGraph, int Ssize, int Fsize, int Isize, int Psize, int Pisize);

private:

    typedef std::vector<int> Index;
    typedef std::vector<double> Value;

    int N;
    Index rowInd, colPtr;
    Value A;
    ID newDofID;
};

#endif
