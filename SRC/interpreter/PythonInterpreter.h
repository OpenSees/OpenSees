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

// Description: This file contains the class definition for PythonInterpreter
// PythonInterpreter implements a DL_Interpreter for the tcl language
//

#ifndef PythonInterpreter_h
#define PythonInterpreter_h

#include "DL_Interpreter.h"
#include <Python.h>




class PythonInterpreter: public DL_Interpreter
{
  public:
    PythonInterpreter(int argc, char **argv);
    virtual ~PythonInterpreter();

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

};


#endif

