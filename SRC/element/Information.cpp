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
// $Date: 2000-09-15 08:23:19 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Information.cpp,v $
                                                                        
                                                                        
// File: ~/element/Information.C
//
// Written: fmk 10/99
// Revised:
//
// Purpose: This file contains the class implementation for Information.

#include "Information.h"
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

Information::Information() 
  :theType(UnknownType), theID(0), theVector(0), theMatrix(0)
{
    // does nothing
}


Information::~Information() 
{
    if (theID != 0)
	delete theID;

    if (theVector != 0)
	delete theVector;
    
    if (theMatrix != 0)
	delete theMatrix;
}





