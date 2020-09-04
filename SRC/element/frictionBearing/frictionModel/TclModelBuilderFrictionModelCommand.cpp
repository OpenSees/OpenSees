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

// $Revision$
// $Date$
// $URL$

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function invoked when the user invokes
// the frictionModel command in the interpreter.

#include <FrictionModel.h>
#include <tcl.h>
#include <elementAPI.h>

extern void *OPS_Coulomb();
extern void *OPS_VelDependent();
extern void *OPS_VelDepMultiLinear();
extern void *OPS_VelNormalFrcDep();
extern void *OPS_VelPressureDep();


static void printCommand(int argc, TCL_Char **argv)
{
    opserr << "Input command: ";
    for (int i=0; i<argc; i++)
        opserr << argv[i] << " ";
    opserr << endln;
} 


int TclModelBuilderFrictionModelCommand(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theDomain)
{
    // make sure there is a minimum number of arguments
    if (argc < 3)  {
        opserr << "WARNING insufficient number of friction model arguments\n";
        opserr << "Want: frictionModel type tag <specific friction model args>\n";
        return TCL_ERROR;
    }
    
    OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, theDomain);
    
    // pointer to a friction model that will be added to the model builder
    FrictionModel *theFrnMdl = 0;
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"Coulomb") == 0 || strcmp(argv[1],"Constant") == 0)  {
        void *theFrn = OPS_Coulomb();
        if (theFrn != 0)
            theFrnMdl = (FrictionModel *)theFrn;
        else
            return TCL_ERROR;
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelDependent") == 0 || strcmp(argv[1],"VDependent") == 0)  {
        void *theFrn = OPS_VelDependent();
        if (theFrn != 0)
            theFrnMdl = (FrictionModel *)theFrn;
        else
            return TCL_ERROR;
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1], "VelDepMultiLinear") == 0 || strcmp(argv[1], "VDependentMultiLinear") == 0) {
        void *theFrn = OPS_VelDepMultiLinear();
        if (theFrn != 0)
            theFrnMdl = (FrictionModel *)theFrn;
        else
            return TCL_ERROR;
    }

    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1], "VelNormalFrcDep") == 0 || strcmp(argv[1], "VNDependent") == 0) {
        void *theFrn = OPS_VelNormalFrcDep();
        if (theFrn != 0)
            theFrnMdl = (FrictionModel *)theFrn;
        else
            return TCL_ERROR;
    }

    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelPressureDep") == 0 || strcmp(argv[1],"VPDependent") == 0)  {
        void *theFrn = OPS_VelPressureDep();
        if (theFrn != 0)
            theFrnMdl = (FrictionModel *)theFrn;
        else
            return TCL_ERROR;
    }
    
    // ----------------------------------------------------------------------------	
    if (theFrnMdl == 0)  {
        opserr << "WARNING could not create friction model " << argv[1] << endln;
        return TCL_ERROR;
    }
    
    // now add the friction model to the modelBuilder
    if (OPS_addFrictionModel(theFrnMdl) == false)  {
        opserr << "WARNING could not add friction model to the domain\n";
        opserr << *theFrnMdl << endln;
        delete theFrnMdl; // invoke the destructor, otherwise mem leak
        return TCL_ERROR;
    }
    
    return TCL_OK;
}
