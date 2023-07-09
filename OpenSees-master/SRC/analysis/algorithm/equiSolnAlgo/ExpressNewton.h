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
// $Date: 2019-01-28 17:53:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/algorithm/equiSolnAlgo/ExpressNewton.h,v $

// Written: Yuli Huang (yulee@berkeley.edu)
// Created: 02/2020
// Revision: A
//
// Description: This file contains the definition for the ExpressNewton class.
// ExpressNewton is a class which performs a ExpressNewton solution algorithm
// to solve the equations.
//
// Reference:
// Junjie Xu, Yuli Huang, Zhe Qu,
// An efficient and unconditionally stable numerical algorithm for nonlinear structural dynamics
// International Journal for Numerical Methods in Engineering,
// https://doi.org/10.1002/nme.6456.
// (https://onlinelibrary.wiley.com/doi/abs/10.1002/nme.6456)

// What: "@(#)ExpressNewton.h, revA"

#ifndef ExpressNewton_h
#define ExpressNewton_h

#include <EquiSolnAlgo.h>

class ExpressNewton: public EquiSolnAlgo
{
  public:
    ExpressNewton(int nIter, double kMultiplier, int tagent, int factorOnce);
    ~ExpressNewton();

    int solveCurrentStep(void);
    int setConvergenceTest(ConvergenceTest *theNewTest);
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag =0);    
    
  protected:
    
  private:
    int factorOnce;
    int nIter;
    double kMultiplier1, kMultiplier2;
};

#endif


