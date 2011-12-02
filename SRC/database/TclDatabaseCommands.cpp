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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/TclDatabaseCommands.cpp,v $
                                                                        
                                                                        
// File: ~/recorders/TclRecordersCommand.C
// 
// Written: fmk 
// Created: 03/00
// Revision: A
//
// Description: This file contains the function that is invoked
// by the interpreter when the comand 'database' is invoked by the 
// user.
//
// What: "@(#) commands.C, revA"

#include <tcl.h>
#include <tk.h>

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <Domain.h>
#include <EquiSolnAlgo.h>

// databases
#include <FileDatastore.h>
#include <FEM_ObjectBroker.h>

static bool createdDatabaseCommands = false;
static FE_Datastore *theDatabase = 0;
static FEM_ObjectBroker *theBroker = 0;

int 
save(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);

int 
restore(ClientData clientData, Tcl_Interp *interp, int argc, char **argv);


int
TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc, 
	       char **argv, Domain &theDomain)
{
  if (createdDatabaseCommands == false) {

    // create the object broker
    theBroker = new FEM_ObjectBroker;

    // create the commands to commit and reset
    Tcl_CreateCommand(interp, "save", save,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "restore", restore,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        

    createdDatabaseCommands = true;
  }

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    interp->result = "WARNING need to specify a Database type ";
    return TCL_ERROR;
  }    

  //
  // check argv[1] for type of Database, parse in rest of arguments
  // needed for the type of Database, create the object and add to Domain
  //

  // a FileDatabase
  if (strcmp(argv[1],"File") == 0) {
    if (argc < 3) {
      cerr << "WARNING database File fileName? ";
      return TCL_ERROR;
    }    

    // delete the old database
    if (theDatabase != 0)
      delete theDatabase;

    theDatabase = new FileDatastore(argv[2], theDomain, *theBroker);
  }
  
  // no recorder type specified yet exists
  else {
    cerr << "WARNING No database type exists ";
    cerr << "for database of type:" << argv[1];
    return TCL_ERROR;
  }    

  // check we instantiated a recorder .. if not ran out of memory
  if (theDatabase == 0) {
    cerr << "WARNING ran out of memory - database " << argv[1]<< endl;
    return TCL_ERROR;
  } 

  // operation successfull
  return TCL_OK;
}




int 
save(ClientData clientData, Tcl_Interp *interp, int argc, 
     char **argv)
{

  if (theDatabase == 0) {
    cerr << "WARNING: save - no database has been constructed\n";
    return TCL_OK;
  }

     // make sure at least one other argument to contain type of system
    if (argc < 2) {
      cerr << "WARNING save no commit tag - want save commitTag?";
      return TCL_OK;
    }    

    // check argv[1] for commitTag
    int commitTag;
    if (Tcl_GetInt(interp, argv[1], &commitTag) != TCL_OK) {
      cerr << "WARNING - save could not read commitTag " << argv[1] << endl;
      return TCL_OK;	
    }	

    if (theDatabase->commitState(commitTag) < 0) {
      cerr << "WARNING - database failed to commitState \n";
      return TCL_ERROR;
    }
    
    return TCL_OK;
}


int 
restore(ClientData clientData, Tcl_Interp *interp, int argc, 
	char **argv)
{

  if (theDatabase == 0) {
    cerr << "WARNING: restore - no database has been constructed\n";
    return TCL_OK;
  }

     // make sure at least one other argument to contain type of system
    if (argc < 2) {
      cerr << "WARNING restore no commit tag - want restore commitTag?";
      return TCL_OK;
    }    

    // check argv[1] for commitTag
    int commitTag;
    if (Tcl_GetInt(interp, argv[1], &commitTag) != TCL_OK) {
      cerr << "WARNING - restore could not read commitTag " << argv[1] << endl;
      return TCL_OK;	
    }	

    if (theDatabase->restoreState(commitTag) < 0) {
      cerr << "WARNING - database failed to restoreState \n";
      return TCL_ERROR;
    }
    
    return TCL_OK;
}
