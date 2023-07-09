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
                                                                        
// $Revision: 1.1 $
// $Date: 2005-07-25 18:06:09 $
// $Source: /usr/local/cvs/OpenSees/SRC/database/TclBerkeleyDB.cpp,v $
                                                                        
// Written: fmk

// Description: This file contains the function invoked when the user invokes
// the MySQL command in the interpreter. 
//
// What: "@(#) TclCommand_MySQL.C, revA"

#include <OPS_Globals.h>
#include <stdlib.h>
#include <string.h>
#include <tcl.h>

#include <BerkeleyDbDatastore.h>

#ifdef _USRDLL
#define DllExport _declspec(dllexport)
#else
#define DllExport
#endif

extern "C" DllExport int
TclCommand_BerkeleyDB(ClientData clientData, 
		      Tcl_Interp *interp,  
		      int argc, 
		      TCL_Char **argv, 
		      Domain *theDomain, 
		      FEM_ObjectBroker *theBroker,
		      FE_Datastore **theDatabase)
{
  
  // delete the old database
  if (*theDatabase != 0)
    delete (*theDatabase);
  
  (*theDatabase) = new BerkeleyDbDatastore(argv[2], *theDomain, *theBroker);

  if (*theDatabase == 0) {
    opserr << "WARNING database MySql dabaseName? - out of memory\n";
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

