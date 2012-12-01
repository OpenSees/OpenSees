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
                                                                        
// Written: fmk 
// Created: Nov, 2012

// Description: This file contains the class definition for TclInterpreter
// TclInterpreter implements a DL_Interpreter for the tcl language
//

#ifndef TclInterpreter_h
#define TclInterpreter_h

#include "DL_Interpreter.h"
#include <tcl.h>


class TclInterpreter: public DL_Interpreter
{
  public:
    TclInterpreter(int argc, char **argv);
    virtual ~TclInterpreter();

    // method to run once the interpreter is set up
    virtual int run();

    // methods to add & remove additional commands
    virtual int addCommand(const char *, Command &);
    virtual int removeCommand(const char *);

    // methods for commands to parse the command line
    virtual int getNumRemainingInputArgs(void);
    virtual int getInt(int *, int numArgs);
    virtual int getDouble(double *, int numArgs);
    virtual int getString(char *cArray, int size);
    virtual int getStingCopy(char **stringPtr);

    // methods for interpreters to output results
    
  private:
    Tcl_Obj *resultPtr;
    Tcl_Obj *commandPtr;
    char buffer[1000], *args;
    int code, gotPartial, tty, length;
    int exitCode;
    Tcl_Channel inChannel, outChannel, errChannel;
    Tcl_DString argString;
    Tcl_Interp *interp;
};


#endif

