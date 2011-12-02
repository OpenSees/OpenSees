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
// $Date: 2000-09-15 08:23:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/analysis/analysis/TransientAnalysis.cpp,v $
                                                                        
                                                                        
// File: ~/analysis/analysis/TransientAnalysis.C
// 
// Written: fmk 
// Created: Sun Sept 15 11:47:47: 1996
// Revision: A
//
// Description: This file contains the implementation of TransientAnalysis.
// TransientAnalysis is an abstract base class, i.e. no objects of it's
// type can be created. 
//
// What: "@(#) TransientAnalysis.C, revA"

#include <TransientAnalysis.h>
#include <Domain.h>

TransientAnalysis::TransientAnalysis(Domain &theDom)
:Analysis(theDom)
{

}


TransientAnalysis::~TransientAnalysis()
{

}
















