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
// $Date: 2005-01-27 22:22:50 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/diagonal/DiagonalSolver.cpp,v $
                                                                        
// Written: fmk 
// Created: Jan 2005
// Revision: A
//
// Description: This file contains the class definition for DiagonalSolver.
// DiagonalSolver is an abstract base class and thus no objects of it's type
// can be instantiated. It has pure virtual functions which must be
// implemented in it's derived classes.  Instances of DiagonalSolver 
// are used to solve a system of equations of type DiagonalSOE.
//
// What: "@(#) DiagonalSolver.C, revA"

#include <DiagonalSolver.h>
#include <DiagonalSOE.h>

DiagonalSolver::DiagonalSolver(int theClassTag)    
:LinearSOESolver(theClassTag),
 theSOE(0)
{

}    

DiagonalSolver::~DiagonalSolver()    
{

}    

int 
DiagonalSolver::setLinearSOE(DiagonalSOE &theProfileSPDSOE)
{
    theSOE = &theProfileSPDSOE;
    return 0;
}

