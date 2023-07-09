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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:30 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/SparseGenColLinSolver.cpp,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/SparseGenCol/SparseGenColLinSolver.C
//
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
// Revision: A
//
// Description: This file contains the implementation of SparseGenColLinSolver.
//
// What: "@(#) SparseGenColLinSolver.C, revA"

#include <SparseGenColLinSolver.h>
#include <SparseGenColLinSOE.h>

SparseGenColLinSolver::SparseGenColLinSolver(int theClassTag)    
:LinearSOESolver(theClassTag),
 theSOE(0)
{

}    

SparseGenColLinSolver::~SparseGenColLinSolver()    
{

}    

int 
SparseGenColLinSolver::setLinearSOE(SparseGenColLinSOE &theSparseGenColSOE)
{
    theSOE = &theSparseGenColSOE;
    return 0;
}








