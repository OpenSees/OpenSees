/* *****************************************************************************
Copyright (c) 2012-2017, The Regents of the University of California (Regents).
All rights reserved.

Redistribution and use in source and binary forms, with or without 
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

The views and conclusions contained in the software and documentation are those
of the authors and should not be interpreted as representing official policies,
either expressed or implied, of the FreeBSD Project.

REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS 
PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, 
UPDATES, ENHANCEMENTS, OR MODIFICATIONS.

*************************************************************************** */

#ifndef TclInterpreter_h
#define TclInterpreter_h

#include "DL_Interpreter.h"
#include <tcl.h>
#include <OPS_Globals.h>
#include "TclWrapper.h"
#include "OpenSeesCommands.h"


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
    virtual const char* getString();
    virtual int getStringCopy(char **stringPtr);
    virtual void resetInput(int cArg);

    // methods for interpreters to output results
    virtual int setInt(int *, int numArgs, bool scalar);
    virtual int setDouble(double *, int numArgs, bool scalar);
    virtual int setString(const char*);
    
  private:
    Tcl_Obj *resultPtr;
    Tcl_Obj *commandPtr;
    char buffer[1000], *args;
    int code, gotPartial, tty, length;
    int exitCode;
    Tcl_Channel inChannel, outChannel, errChannel;
    Tcl_DString argString;
    Tcl_Interp *interp;

    TclWrapper wrapper;
    OpenSeesCommands cmds;
};


#endif

