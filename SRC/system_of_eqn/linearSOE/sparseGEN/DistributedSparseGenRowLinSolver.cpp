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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-12-06 22:20:18 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/DistributedSparseGenRowLinSolver.cpp,v $
                                                                        
// Written: fmk 
// Created: 04/05
//
// Description: This file contains the implementation for DistributedSparseGenRowSolver
//
// What: "@(#) DistributedSparseGenRowLinSolver.C, revA"

#include <DistributedSparseGenRowLinSolver.h>
#include <DistributedSparseGenRowLinSOE.h>

DistributedSparseGenRowLinSolver::DistributedSparseGenRowLinSolver(int theClassTag)    
:LinearSOESolver(theClassTag),
 theSOE(0)
{

}    

DistributedSparseGenRowLinSolver::~DistributedSparseGenRowLinSolver()    
{

}    

int 
DistributedSparseGenRowLinSolver::setLinearSOE(DistributedSparseGenRowLinSOE &theDistributedSparseGenRowSOE)
{
    theSOE = &theDistributedSparseGenRowSOE;
    return 0;
}








