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
                                                                        
// $Revision: 1.5 $
// $Date: 2008-04-14 21:56:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/ErrorHandler.h,v $
                                                                        
                                                                        
#ifndef ErrorHandler_h
#define ErrorHandler_h
 
// File: ~/handler/ErrorHandler.h
//
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class interface for ErrorHandler.
// An ErrorHandler is an object is used to process error information. The
// ErrorHandler class is an abstract class.
//
// What: "@(#) ErrorHandler.h, revA"

#include <iostream>
using std::ostream;

//class ostream;
#include <stdarg.h>

class ErrorHandler {
  public:
    ErrorHandler();
	virtual ~ErrorHandler();

	virtual void warning(const char *, ...) =0;
	virtual void fatal(const char *, ...) =0;

  protected:
	void outputMessage(ostream &theStream, const char *msg, va_list args);
    
  private:    

};

#endif


