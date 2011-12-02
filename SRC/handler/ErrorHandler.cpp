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
// $Date: 2008-01-16 00:19:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/handler/ErrorHandler.cpp,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 11/99
// Revision: A
//
// Description: This file contains the class implementation for
// ErrorHandler. 
//
// What: "@(#) ErrorHandler.C, revA"

#include "ErrorHandler.h"

ErrorHandler::ErrorHandler()
{

}

ErrorHandler::~ErrorHandler()
{
    // does nothing
}
 
void
ErrorHandler::outputMessage(ostream &theStream, const char *msg, va_list args) 
{
  int dataInt;
  double dataDouble;
  char *dataString;
  char dataChar;

  int pos =0;
  int done =0;
  
  // parse the msg string until end of string '\0' is
  // encounterd, send to the output stream any character,
  // if we encounter a %d or a %f get the integer or double 
  // from the next arg on the va_list and send it to the stream

  while (done != 1) {

    // if reach string end then we are done
    if (msg[pos] == '\0')
      break;

    // otherwise parse string , looking for special %d and %f
    if (msg[pos] != '%') {
      dataChar = msg[pos];
      theStream << dataChar;
    } else {
      pos++;
      switch(msg[pos]) {
      case 'd':
	dataInt = va_arg(args,int);
	theStream << dataInt;
	break;
	
      case 'f':
	dataDouble = va_arg(args,double);
	theStream << dataDouble;
	break;	

      case 's':
	dataString = va_arg(args,char *);
	theStream << dataString;
	break;	

      default:
	;
      }
     }
    pos++;
  }
  theStream << std::endl;
}
