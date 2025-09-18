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
                                                                        
// Written: Minjie

// Description: A tcl wrapper for OpenSees commands
//
#ifndef TclWrapper_h
#define TclWrapper_h

#include <tcl.h>
#include <OPS_Globals.h>

#include <vector>
#include <map>

class TclWrapper
{
public:

    TclWrapper();
    ~TclWrapper();

    // reset command line
    void resetCommandLine(int nArgs, int cArg, TCL_Char** argv);
    void resetCommandLine(int cArg);

    // wrapper commands
    void addOpenSeesCommands(Tcl_Interp* interp);
    void addCommand(Tcl_Interp* interp, const char* name, Tcl_CmdProc* proc);

    // get command line arguments
    TCL_Char** getCurrentArgv() {return currentArgv;}
    int getCurrentArg() const {return currentArg;}
    int getNumberArgs() const {return numberArgs;}
    void incrCurrentArg() {currentArg++;}

    // set outputs
    void setOutputs(Tcl_Interp* interp, int* data, int numArgs);
    void setOutputs(Tcl_Interp* interp, double* data, int numArgs);
    void setOutputs(Tcl_Interp* interp, const char* str);
    void setOutputs(Tcl_Interp* interp, std::vector<std::vector<int>>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, int>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, std::vector<int>>& data);
    void setOutputs(Tcl_Interp* interp, std::vector<std::vector<double>>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, double>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, std::vector<double>>& data);
    void setOutputs(Tcl_Interp* interp, std::vector<const char*>& data);
    void setOutputs(Tcl_Interp* interp, std::vector<std::vector<const char*>>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, const char*>& data);
    void setOutputs(Tcl_Interp* interp, std::map<const char*, std::vector<const char*>>& data);

private:

    // command line arguments
    TCL_Char** currentArgv;
    int currentArg;
    int numberArgs;

};

#endif
