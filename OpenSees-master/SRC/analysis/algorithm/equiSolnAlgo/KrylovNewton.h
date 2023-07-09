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
                                                                        
// $Revision: 1.9 $
// $Date: 2007-04-02 23:41:13 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/KrylovNewton.h,v $
                                                                        
#ifndef KrylovNewton_h
#define KrylovNewton_h

// Written: MHS
// Created: June 2001
//
// Description: This file contains the class definition for 
// KrylovNewton.  KrylovNewton is a class which uses a Krylov
// subspace accelerator on the modified Newton method.
// The accelerator is described by Carlson and Miller in
// "Design and Application of a 1D GWMFE Code"
// from SIAM Journal of Scientific Computing (Vol. 19, No. 3,
// pp. 728-765, May 1998)

#include <EquiSolnAlgo.h>
#include <Vector.h>

class KrylovNewton: public EquiSolnAlgo
{
  public:
    KrylovNewton(int tangent = CURRENT_TANGENT, int maxDim = 3);    
    KrylovNewton(ConvergenceTest &theTest, int tangent = CURRENT_TANGENT, int maxDim = 3);
    ~KrylovNewton();

    int solveCurrentStep(void);    
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    int tangent;

    // Storage for update vectors
    Vector **v;

    // Storage for subspace vectors
    Vector **Av;

    // Array data sent to LAPACK subroutine
    double *AvData;
    double *rData;
    double *work;

    // Length of work array
    int lwork;

    // Size information
    int numEqns;
    int maxDimension;

    // Private lsq routine to do Krylov updates
    // dimension is the current dimension of the subspace
    int leastSquares(int dimension);
};

#endif


