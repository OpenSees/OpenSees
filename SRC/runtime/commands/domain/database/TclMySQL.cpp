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
//
// Description: This file contains the function invoked when the user invokes
// the MySQL command in the interpreter.
//
// What: "@(#) TclCommand_MySQL.C, revA"

#include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include <MySqlDatastore.h>

#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

extern "C" DllExport int
TclCommand_MySQL(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv, Domain *theDomain,
                 FEM_ObjectBroker *theBroker, FE_Datastore **theDatabase)
{

  if (argc < 3) {
    opserr << "WARNING database MySql dabaseName? ";
    return TCL_ERROR;
  }

  // delete the old database
  if (*theDatabase != 0)
    delete (*theDatabase);

  if (argc == 3)
    (*theDatabase) = new MySqlDatastore(argv[2], *theDomain, *theBroker);
  else {
    const char *database = argv[2];
    const char *host = NULL;
    const char *user = NULL;
    const char *passwd = NULL;
    const char *socket = NULL;
    int port = 0;
    int clientFlag = 0;

    int counter = 3;
    while (counter < argc) {
      if (strcmp(argv[counter], "-host") == 0) {
        host = argv[counter + 1];
        counter += 2;
      } else if (strcmp(argv[counter], "-user") == 0) {
        user = argv[counter + 1];
        counter += 2;
      } else if (strcmp(argv[counter], "-passwd") == 0) {
        passwd = argv[counter + 1];
        counter += 2;
      } else if (strcmp(argv[counter], "-socket") == 0) {
        socket = argv[counter + 1];
        counter += 2;
      } else if (strcmp(argv[counter], "-port") == 0) {
        if (Tcl_GetInt(interp, argv[counter + 1], &port) != TCL_OK)
          return TCL_ERROR;
        counter += 2;
      } else if (strcmp(argv[counter], "-clientFlag") == 0) {
        if (Tcl_GetInt(interp, argv[counter + 1], &clientFlag) != TCL_OK)
          return TCL_ERROR;
        counter += 2;
      } else {
        counter++;
      }
    }
    (*theDatabase) =
        new MySqlDatastore(database, host, user, passwd, port, socket,
                           clientFlag, *theDomain, *theBroker);
  }

  if (*theDatabase == 0) {
    opserr << "WARNING database MySql dabaseName? - out of memory\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
