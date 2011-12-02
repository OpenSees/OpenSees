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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/tcl/commands.h,v $
                                                                        
                                                                        
// File: ~/tcl/commands.h
// 
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the functions that will be called by
// the interpreter when the appropriate command name is specified,
// see tkAppInit.C for command names.
//
// What: "@(#) commands.C, revA"

int
g3AppInit(Tcl_Interp *interp);

int 
wipeModel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
wipeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
resetModel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
initializeAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
setLoadConst(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
setTime(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
buildModel(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
analyzeModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	     char **argv);

int 
printModel(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv);
int 
specifyAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv);
int 
specifySOE(ClientData clientData, Tcl_Interp *interp, int argc, 
	   char **argv);

int 
specifyNumberer(ClientData clientData, Tcl_Interp *interp, int argc, 
		char **argv);
int 
specifyConstraintHandler(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv);
int
specifyAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc, 
		 char **argv);

int
specifyCTest(ClientData clientData, Tcl_Interp *interp, int argc, 
	     char **argv);

int 
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv);
int 
addRecorder(ClientData clientData, Tcl_Interp *interp, int argc, 
	    char **argv);
int 
addAlgoRecorder(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
addDatabase(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
playbackRecorders(ClientData clientData, Tcl_Interp *interp, int argc, 
		  char **argv);
int 
groundExcitation(ClientData clientData, Tcl_Interp *interp, int argc, 
		 char **argv);
int 
rigidLink(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
rigidDiaphragm(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
eigenAnalysis(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
videoPlayer(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
removeObject(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);












