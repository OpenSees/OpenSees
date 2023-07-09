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

// $Revision$
// $Date$
// $URL$

// Written: Minjie
// Created: May 22 2018

#ifndef PFEMSolver_LumpM_h
#define PFEMSolver_LumpM_h

#include <PFEMSolver.h>
extern "C" {
#include <cs.h>
}
#include "../../../../OTHER/UMFPACK/umfpack.h"

class PFEMLinSOE;

class PFEMSolver_LumpM : public PFEMSolver
{
public:
    PFEMSolver_LumpM(bool once);
    virtual ~PFEMSolver_LumpM();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private:

    void *MSym, *MNum, *LSym, *LNum;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    PFEMLinSOE* theSOE;
    bool numonce, factored;
};

#endif
