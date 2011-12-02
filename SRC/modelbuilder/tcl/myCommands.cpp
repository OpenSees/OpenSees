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
// $Source: /usr/local/cvs/OpenSees/SRC/modelbuilder/tcl/myCommands.cpp,v $
                                                                        
                                                                        
// File: ~/modelbuilder/tcl/myCommands.C
// 
// Written: fmk 
// Created: 04/98
// Revision: A
//
// Description: This file contains the function myCommands().
// myCommands() is called in g3AppInit() - all new user commands
// are to be placed in here.
//
// What: "@(#) myCommands.C, revA"

#include <Domain.h>
#include "TclModelBuilder.h"

#include <tcl.h>
#include <tk.h>

#include <iostream.h>
#include <fstream.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

extern ModelBuilder *theBuilder;
extern Domain theDomain;

int
specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc, 
		    char **argv);

int myCommands(Tcl_Interp *interp) {
    Tcl_CreateCommand(interp, "model", specifyModelBuilder,
		      (ClientData)NULL, (Tcl_CmdDeleteProc *)NULL);
    return 0;
}

int
specifyModelBuilder(ClientData clientData, Tcl_Interp *interp, int argc, 
			 char **argv)
{
    // make sure at least one other argument to contain model builder type given
    if (argc < 2) {
      cerr << "WARNING need to specify a model type, valid types:\n";
      cerr << "\tBasicBuilder\n";
      return TCL_ERROR;
    }    

    // check argv[1] for type of ModelBuilder and create the object 
    if (strcmp(argv[1],"basic") == 0 || strcmp(argv[1],"BasicBuilder") == 0) {
      int ndm =0;
      int ndf = 0;
      
      if (argc < 4) {
	cerr << "WARNING incorrect number of command arguments\n";
	cerr << "model modelBuilderType -ndm ndm? <-ndf ndf?> \n";
	return TCL_ERROR;
      }

      int argPos = 2;
      while (argPos < argc) {
	if (strcmp(argv[argPos],"-ndm") == 0 ||
	                strcmp(argv[argPos],"-NDM") == 0) {	
	  argPos++;
	  if (argPos < argc)
	    if (Tcl_GetInt(interp, argv[argPos], &ndm) != TCL_OK) {
	      cerr << "WARNING error reading ndm: " << argv[argPos];
	      cerr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
	      return TCL_ERROR;
	    }	  
	  argPos++;
	}

	else if (strcmp(argv[argPos],"-ndf") == 0 ||
		        strcmp(argv[argPos],"-NDF") == 0) {	
	  argPos++;
	  if (argPos < argc)
	    if (Tcl_GetInt(interp, argv[argPos], &ndf) != TCL_OK) {
	      cerr << "WARNING error reading ndf: " << argv[argPos];
	      cerr << "\nmodel modelBuilderType -ndm ndm? <-ndf ndf?>\n";
	      return TCL_ERROR;
	    }	  
	  argPos++;
	}

	else // Advance to next input argument if there are no matches -- MHS
	  argPos++;
      }

      // check that ndm was specified
      if (ndm == 0) {
	cerr << "WARNING need to specify ndm\n";
	cerr << "model modelBuilderType -ndm ndm? <-ndf ndf?>\n";
	return TCL_ERROR;
      }

      // check for ndf, if not assume one
      if (ndf == 0) {
	if (ndm == 1) 
	  ndf = 1;
	else if (ndm == 2)
	  ndf = 3;
	else if (ndm == 3)
	  ndf = 6;
	else {
	  cerr << "WARNING specified ndm, " << ndm << ", will not work\n";
	  cerr << "with any elements in BasicBuilder\n";
	  return TCL_ERROR;
	}
      }

      // create the model builder
      theBuilder = new TclModelBuilder(theDomain, interp, ndm, ndf);
      if (theBuilder == 0) {
	cerr << "WARNING ran out of memory in creating BasicBuilder model\n";
	return TCL_ERROR;
      }
    }
    else
    {
      cerr << "WARNING model builder type " << argv[1]
	   << " not supported\n";
      return TCL_ERROR;
    }
    
    return TCL_OK;
}





