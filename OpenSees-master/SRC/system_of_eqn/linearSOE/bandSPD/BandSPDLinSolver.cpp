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
                                                                        
// $Revision: 1.2 $
// $Date: 2009-05-11 20:55:31 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/bandSPD/BandSPDLinSolver.cpp,v $
                                                                        
// Written: fmk 
// Created: Tue Sep 26 16:27:47: 1996
//
// Description: This file contains the implementation of BandSPDLinSolver.
//
// What: "@(#) BandSPDLinSolver.C, revA"

#include <BandSPDLinSolver.h>
#include <BandSPDLinSOE.h>

BandSPDLinSolver::BandSPDLinSolver(int theClassTag)    
:LinearSOESolver(theClassTag),
 theSOE(0)
{

}    

BandSPDLinSolver::~BandSPDLinSolver()    
{

}    

int 
BandSPDLinSolver::setLinearSOE(BandSPDLinSOE &theBandSPDSOE)
{
    theSOE = &theBandSPDSOE;
    return 0;
}

