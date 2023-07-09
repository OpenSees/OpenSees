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
// Description: Solve PFEM with quasi-incompressible form, no global matrix operations
//
                                                                        
                                                                        
#ifndef PFEMQuasiSolver_h
#define PFEMQuasiSolver_h


#include <LinearSOESolver.h>
#include <Vector.h>
extern "C" {
#include <cs.h>
}
#include "../../../../OTHER/UMFPACK/umfpack.h"


class PFEMQuasiLinSOE;

class PFEMQuasiSolver : public LinearSOESolver
{
public:
    PFEMQuasiSolver();
    virtual ~PFEMQuasiSolver();

    int solve();
    int setSize();
    virtual int setLinearSOE(PFEMQuasiLinSOE& theSOE);

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);  

private:

    void reorder(cs* A);
    
    PFEMQuasiLinSOE* theSOE;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
};

#endif

