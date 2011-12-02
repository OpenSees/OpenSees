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
                                                                        
// $Revision: 1.4 $
// $Date: 2010-02-04 19:13:39 $
// $Source: /usr/local/cvs/OpenSees/SRC/G3Globals.h,v $
                                                                        
                                                                        
#ifndef G3Globals_h
#define G3Globals_h

// Written: fmk 
// Created: 11/99
//
// Description: This file contains global variables used in OpenSees files.
// if you change some of the variables, you must recompile ALL the code.

//#include <ErrorHandler.h>
class Domain;
class Element;

#define MAX_FILENAMELENGTH 50

//extern ErrorHandler *g3ErrorHandler;   // error handler for sending warning & fatal error messages
extern double   ops_Dt;                // current delta T for current domain doing an update
// extern double  *ops_Gravity;        // gravity factors for current domain undergoing an update
extern Domain  *ops_TheActiveDomain;   // current domain undergoing an update
extern Element *ops_TheActiveElement;  // current element undergoing an update

#endif
