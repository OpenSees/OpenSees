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
// $Date: 2003-02-25 23:32:40 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/TclDatabaseCommands.cpp,v $
                                                                        
                                                                        
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

#include <OPS_Globals.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include <Domain.h>
#include <EquiSolnAlgo.h>

// databases
#include <FileDatastore.h>

#ifdef _MYSQL
#include <MySqlDatastore.h>
#endif
#ifdef _BERKELEYDB
#include <BerkeleyDbDatastore.h>
#endif
#include <FEM_ObjectBroker.h>

static bool createdDatabaseCommands = false;

int 
save(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

int 
restore(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv);

extern FE_Datastore *theDatabase;

int
TclAddDatabase(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, 
	       Domain &theDomain, 
	       FEM_ObjectBroker &theBroker)
{
  if (createdDatabaseCommands == false) {

    // create the commands to commit and reset
    Tcl_CreateCommand(interp, "save", save,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);    
    Tcl_CreateCommand(interp, "restore", restore,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);        

    createdDatabaseCommands = true;
  }

  // make sure at least one other argument to contain integrator
  if (argc < 2) {
    opserr << "WARNING need to specify a Database type; valid type File, MySQL, BerkeleyDB \n";
    return TCL_ERROR;
  }    

  //
  // check argv[1] for type of Database, parse in rest of arguments
  // needed for the type of Database, create the object and add to Domain
  //

  // a File Database
  if (strcmp(argv[1],"File") == 0) {
    if (argc < 3) {
      opserr << "WARNING database File fileName? ";
      return TCL_ERROR;
    }    

    // delete the old database
    if (theDatabase != 0)
      delete theDatabase;

    theDatabase = new FileDatastore(argv[2], theDomain, theBroker);
  }

#ifdef _MYSQL
  // a MySQL Database
  else if (strcmp(argv[1],"MySQL") == 0) {
    if (argc < 3) {
      opserr << "WARNING database MySql fileName? ";
      return TCL_ERROR;
    }    

    // delete the old database
    if (theDatabase != 0)
      delete theDatabase;

    theDatabase = new MySqlDatastore(argv[2], theDomain, theBroker);
  }
#endif

#ifdef _BERKELEYDB
  // a BerkeleyDB database
  else  if (strcmp(argv[1],"BerkeleyDB") == 0) {
    if (argc < 3) {
      opserr << "WARNING database BerkeleyDB fileName? ";
      return TCL_ERROR;
    }    


    // delete the old database 
    if (theDatabase != 0) 
      delete theDatabase;

    theDatabase = new BerkeleyDbDatastore(argv[2], theDomain, theBroker);
  }
#endif

  // no recorder type specified yet exists
  else {
    opserr << "WARNING No database type exists ";
    opserr << "for database of type:" << argv[1] << "valid database types File, BerkeleyDB and MySQL\n";
    return TCL_ERROR;
  }    

  // check we instantiated a database .. if not ran out of memory
  if (theDatabase == 0) {
    opserr << "WARNING ran out of memory - database " << argv[1]<< endln;
    return TCL_ERROR;
  } 

  // operation successfull
  return TCL_OK;
}




int 
save(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
{

  if (theDatabase == 0) {
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
restore(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv)
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
