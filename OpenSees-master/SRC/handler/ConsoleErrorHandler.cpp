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
// $Source: /usr/local/cvs/OpenSees/SRC/handler/ConsoleErrorHandler.cpp,v $
                                                                        
                                                                        
// File: ~/handler/ConsoleErrorHandler.C
//
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for ConsoleErrorHandler.
//
// What: "@(#) ConsoleErrorHandler.C, revA"

#include "ConsoleErrorHandler.h"
#include <stdlib.h>

ConsoleErrorHandler::ConsoleErrorHandler()
  :ErrorHandler()
{

}

ConsoleErrorHandler::~ConsoleErrorHandler()
{
    // does nothing
}
 

void 
ConsoleErrorHandler::warning(const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  this->outputMessage(std::cerr, msg, args);
  va_end(args);
}  

void 
ConsoleErrorHandler::fatal(const char *msg, ...)
{
  va_list args;
  va_start(args, msg);
  this->outputMessage(std::cerr, msg, args);
  va_end(args);
  exit(-1);
}
