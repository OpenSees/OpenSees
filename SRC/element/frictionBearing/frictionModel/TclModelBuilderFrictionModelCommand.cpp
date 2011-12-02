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
// $Date: 2009-04-17 23:02:41 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/frictionBearing/frictionModel/TclModelBuilderFrictionModelCommand.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the frictionModel command in the interpreter. 
//
// What: "@(#) TclModelBuilderFrictionModelCommand.cpp, revA"

#include <TclModelBuilder.h>
#include <FrictionModel.h>

#include <CoulombFriction.h>
#include <VDependentFriction.h>
#include <VPDependentFriction.h>

#include <ID.h>
#include <Vector.h>
#include <string.h>


static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
		opserr << argv[i] << " ";
    opserr << endln;
} 


int TclModelBuilderFrictionModelCommand(ClientData clientData, Tcl_Interp *interp, int argc,
    TCL_Char **argv, TclModelBuilder *theTclBuilder, Domain *theDomain)
{
    // make sure there is a minimum number of arguments
    if (argc < 3)  {
		opserr << "WARNING insufficient number of friction model arguments\n";
		opserr << "Want: frictionModel type tag <specific friction model args>\n";
		return TCL_ERROR;
    }
    
    // pointer to a friction model that will be added to the model builder
    FrictionModel *theFrnMdl = 0;
	
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"Coulomb") == 0)  {
		if (argc != 4)  {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: frictionModel Coulomb tag mu\n";
			return TCL_ERROR;
		}    
		
		int tag;
        double mu;
		
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
			opserr << "WARNING invalid Coulomb friction model tag\n";
			return TCL_ERROR;		
		}
		if (Tcl_GetDouble(interp, argv[3], &mu) != TCL_OK)  {
			opserr << "WARNING invalid mu\n";
			opserr << "Coulomb friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		
		// parsing was successful, allocate the friction model
		theFrnMdl = new CoulombFriction(tag, mu);
    }
	
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VDependent") == 0)  {
		if (argc != 6)  {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: frictionModel VDependent tag muSlow muFast transRate\n";
			return TCL_ERROR;
		}    
		
		int tag;
        double muSlow, muFast, transRate;
		
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
			opserr << "WARNING invalid VDependent friction model tag\n";
			return TCL_ERROR;		
		}
		if (Tcl_GetDouble(interp, argv[3], &muSlow) != TCL_OK)  {
			opserr << "WARNING invalid muSlow\n";
			opserr << "VDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[4], &muFast) != TCL_OK)  {
			opserr << "WARNING invalid muFast\n";
			opserr << "VDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[5], &transRate) != TCL_OK)  {
			opserr << "WARNING invalid transRate\n";
			opserr << "VDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		
		// parsing was successful, allocate the friction model
		theFrnMdl = new VDependentFriction(tag, muSlow, muFast, transRate);
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VPDependent") == 0)  {
		if (argc != 9)  {
			opserr << "WARNING invalid number of arguments\n";
			printCommand(argc,argv);
			opserr << "Want: frictionModel VPDependent tag muSlow muFast0 A deltaMu alpha transRate\n";
			return TCL_ERROR;
		}    
		
		int tag;
        double muSlow, muFast0, A, deltaMu, alpha, transRate;
		
		if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
			opserr << "WARNING invalid VPDependent friction model tag\n";
			return TCL_ERROR;		
		}
		if (Tcl_GetDouble(interp, argv[3], &muSlow) != TCL_OK)  {
			opserr << "WARNING invalid muSlow\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[4], &muFast0) != TCL_OK)  {
			opserr << "WARNING invalid muFast0\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)  {
			opserr << "WARNING invalid A\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[6], &deltaMu) != TCL_OK)  {
			opserr << "WARNING invalid deltaMu\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[7], &alpha) != TCL_OK)  {
			opserr << "WARNING invalid alpha\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		if (Tcl_GetDouble(interp, argv[8], &transRate) != TCL_OK)  {
			opserr << "WARNING invalid transRate\n";
			opserr << "VPDependent friction model: " << tag << endln;
			return TCL_ERROR;	
		}
		
		// parsing was successful, allocate the friction model
		theFrnMdl = new VPDependentFriction(tag, muSlow, muFast0, A, deltaMu, alpha, transRate);
    }
	
    // ----------------------------------------------------------------------------	
	if (theFrnMdl == 0)  {
		opserr << "WARNING could not create friction model " << argv[1] << endln;
		return TCL_ERROR;
	}
	
	// now add the friction model to the modelBuilder
	if (theTclBuilder->addFrictionModel(*theFrnMdl) < 0)  {
		opserr << "WARNING could not add friction model to the domain\n";
		opserr << *theFrnMdl << endln;
		delete theFrnMdl; // invoke the destructor, otherwise mem leak
		return TCL_ERROR;
	}
	
	return TCL_OK;
}
