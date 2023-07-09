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
// $Date: 2003-02-14 23:01:24 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/ConsoleErrorHandler.h,v $
                                                                        
                                                                        
#ifndef ConsoleErrorHandler_h
#define ConsoleErrorHandler_h

// File: ~/handler/ConsoleErrorHandler.h
//
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class interface for ConsoleErrorHandler.
// A ConsoleErrorHandler is an object which sends the error messages to the opserr stream.
//
// What: "@(#) ConsoleErrorHandler.h, revA"

#include "ErrorHandler.h"

class ConsoleErrorHandler : public ErrorHandler
{
  public:
    ConsoleErrorHandler();
    virtual ~ConsoleErrorHandler();

    void warning(const char *, ...);
    void fatal(const char *, ...);

  protected:
    
  private:    
};


#endif


