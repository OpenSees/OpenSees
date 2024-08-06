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

#include <OPS_Globals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Domain.h>
#include <EquiSolnAlgo.h>

// known databases
#include <FileDatastore.h>

// linked list of struct for other types of
// databases that can be added dynamically

#include <packages.h>

typedef struct databasePackageCommand {
  char *funcName;
  int (*funcPtr)(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv, Domain *, FEM_ObjectBroker *,
                 FE_Datastore **);
  struct databasePackageCommand *next;
} DatabasePackageCommand;

// static variables
static DatabasePackageCommand *theDatabasePackageCommands = NULL;
static bool createdDatabaseCommands = false;

int save(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);

int restore(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv);

extern FE_Datastore *theDatabase;

int
TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv,
               Domain &theDomain, FEM_ObjectBroker &theBroker)
{
  if (createdDatabaseCommands == false) {

    // create the commands to commit and reset
    Tcl_CreateCommand(interp, "save", save, (ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);
    Tcl_CreateCommand(interp, "restore", restore, (ClientData)NULL,
                      (Tcl_CmdDeleteProc *)NULL);

    createdDatabaseCommands = true;
  }

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify a Database type; valid type File, "
              "MySQL, BerkeleyDB \n";
    return TCL_ERROR;
  }

  //
  // check argv[1] for type of Database, parse in rest of arguments
  // needed for the type of Database, create the object and add to Domain
  //

  // a File Database
  if (strcmp(argv[1], "File") == 0) {
    if (argc < 3) {
      opserr << "WARNING database File fileName? ";
      return TCL_ERROR;
    }

    // delete the old database
    if (theDatabase != 0)
      delete theDatabase;

    theDatabase = new FileDatastore(argv[2], theDomain, theBroker);
    // check we instantiated a database .. if not ran out of memory
    if (theDatabase == nullptr) {
      opserr << "WARNING ran out of memory - database File " << argv[2]
             << endln;
      return TCL_ERROR;
    }

    return TCL_OK;
  } else {

    //
    // maybe a database package
    //

    // try existing loaded packages

    DatabasePackageCommand *dataCommands = theDatabasePackageCommands;
    bool found = false;
    while (dataCommands != NULL && found == false) {
      if (strcmp(argv[1], dataCommands->funcName) == 0) {
        int result = (*(dataCommands->funcPtr))(
            clientData, interp, argc, argv, &theDomain, &theBroker, &theDatabase);
        return result;
      } else
        dataCommands = dataCommands->next;
    }

    // load new package

    void *libHandle;
    int (*funcPtr)(ClientData, Tcl_Interp *, int ,
                   TCL_Char ** const argv, Domain *, FEM_ObjectBroker *,
                   FE_Datastore **);
    int databaseNameLength = strlen(argv[1]);
    char *tclFuncName = new char[databaseNameLength + 12];
    strcpy(tclFuncName, "TclCommand_");
    strcpy(&tclFuncName[11], argv[1]);

    int res =
        getLibraryFunction(argv[1], tclFuncName, &libHandle, (void **)&funcPtr);

    if (res == 0) {
      char *databaseName = new char[databaseNameLength + 1];
      strcpy(databaseName, argv[1]);
      DatabasePackageCommand *theDataCommand = new DatabasePackageCommand;
      theDataCommand->funcPtr = funcPtr;
      theDataCommand->funcName = databaseName;
      theDataCommand->next = theDatabasePackageCommands;
      theDatabasePackageCommands = theDataCommand;

      int result = (*funcPtr)(clientData, interp, argc, argv, &theDomain,
                              &theBroker, &theDatabase);
      return result;
    }
  }
  opserr << "WARNING No database type exists ";
  opserr << "for database of type:" << argv[1] << "valid database type File\n";

  return TCL_ERROR;
}

int
save(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  if (theDatabase == nullptr) {
    opserr << "WARNING: save - no database has been constructed\n";
    return TCL_OK;
  }

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING save no commit tag - want save commitTag?";
    return TCL_OK;
  }

  // check argv[1] for commitTag
  int commitTag;
  if (Tcl_GetInt(interp, argv[1], &commitTag) != TCL_OK) {
    opserr << "WARNING - save could not read commitTag " << argv[1] << endln;
    return TCL_OK;
  }

  if (theDatabase->commitState(commitTag) < 0) {
    opserr << "WARNING - database failed to commitState \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
restore(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  if (theDatabase == 0) {
    opserr << "WARNING: restore - no database has been constructed\n";
    return TCL_OK;
  }

  // make sure at least one other argument to contain type of system
  if (argc < 2) {
    opserr << "WARNING restore no commit tag - want restore commitTag?";
    return TCL_OK;
  }

  // check argv[1] for commitTag
  int commitTag;
  if (Tcl_GetInt(interp, argv[1], &commitTag) != TCL_OK) {
    opserr << "WARNING - restore could not read commitTag " << argv[1] << endln;
    return TCL_OK;
  }

  if (theDatabase->restoreState(commitTag) < 0) {
    opserr << "WARNING - database failed to restoreState \n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
