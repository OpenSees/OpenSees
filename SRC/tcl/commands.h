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
                                                                        
// $Revision: 1.12 $
// $Date: 2004-07-12 21:21:58 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/commands.h,v $
                                                                        
                                                                        
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.
//
// What: "@(#) commands.C, revA"

#include <OPS_Globals.h>

int
g3AppInit(Tcl_Interp *interp);

int 
wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
wipeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
setTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
getTime(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
buildModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
printModel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
specifySOE(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
specifyNumberer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int
specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int
specifyCTest(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
addRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
addAlgoRecorder(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
addDatabase(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
playbackRecorders(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
playbackAlgorithmRecorders(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
int 
rigidLink(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
rigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
videoPlayer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
removeObject(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

// AddingSensitivity:BEGIN /////////////////////////////////////////////////
int 
nodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
computeGradients(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
sensNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
sensNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
sensitivityAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
sensitivityIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);
// AddingSensitivity:END ///////////////////////////////////////////////////


int 
startTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
stopTimer(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
rayleighDamping(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
addRegion(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);









