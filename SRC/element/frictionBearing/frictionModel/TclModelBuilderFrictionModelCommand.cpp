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

#include <Coulomb.h>
#include <VelDependent.h>
#include <VelPressureDep.h>
#include <VelDepMultiLinear.h>
#include <VelNormalFrcDep.h>

#include <ID.h>
#include <Vector.h>
#include <string.h>
#include <FrictionModel.h>


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
    
    // pointer to a friction model that will be added to the model builder
    FrictionModel *theFrnMdl = 0;
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"Coulomb") == 0 || strcmp(argv[1],"Constant") == 0)  {
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
        theFrnMdl = new Coulomb(tag, mu);
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelDependent") == 0 || strcmp(argv[1],"VDependent") == 0)  {
        if (argc != 6)  {
            opserr << "WARNING invalid number of arguments\n";
            printCommand(argc,argv);
            opserr << "Want: frictionModel VelDependent tag muSlow muFast transRate\n";
            return TCL_ERROR;
        }    
        
        int tag;
        double muSlow, muFast, transRate;
        
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
            opserr << "WARNING invalid VelDependent friction model tag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[3], &muSlow) != TCL_OK)  {
            opserr << "WARNING invalid muSlow\n";
            opserr << "VelDependent friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4], &muFast) != TCL_OK)  {
            opserr << "WARNING invalid muFast\n";
            opserr << "VelDependent friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5], &transRate) != TCL_OK)  {
            opserr << "WARNING invalid transRate\n";
            opserr << "VelDependent friction model: " << tag << endln;
            return TCL_ERROR;
        }
        
        // parsing was successful, allocate the friction model
        theFrnMdl = new VelDependent(tag, muSlow, muFast, transRate);
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelPressureDep") == 0 || strcmp(argv[1],"VPDependent") == 0)  {
        if (argc != 9)  {
            opserr << "WARNING invalid number of arguments\n";
            printCommand(argc,argv);
            opserr << "Want: frictionModel VelPressureDep tag muSlow muFast0 A deltaMu alpha transRate\n";
            return TCL_ERROR;
        }    
        
        int tag;
        double muSlow, muFast0, A, deltaMu, alpha, transRate;
        
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
            opserr << "WARNING invalid VelPressureDep friction model tag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[3], &muSlow) != TCL_OK)  {
            opserr << "WARNING invalid muSlow\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4], &muFast0) != TCL_OK)  {
            opserr << "WARNING invalid muFast0\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5], &A) != TCL_OK)  {
            opserr << "WARNING invalid A\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6], &deltaMu) != TCL_OK)  {
            opserr << "WARNING invalid deltaMu\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[7], &alpha) != TCL_OK)  {
            opserr << "WARNING invalid alpha\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[8], &transRate) != TCL_OK)  {
            opserr << "WARNING invalid transRate\n";
            opserr << "VelPressureDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        
        // parsing was successful, allocate the friction model
        theFrnMdl = new VelPressureDep(tag, muSlow, muFast0, A, deltaMu, alpha, transRate);
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelDepMultiLinear") == 0 || strcmp(argv[1],"VDependentMultiLinear") == 0)  {
        if (argc < 9)  {
            opserr << "WARNING invalid number of arguments\n";
            printCommand(argc,argv);
            opserr << "Want: frictionModel VelDepMultiLinear tag ";
            opserr << "-vel velocityPoints -frn frictionPoints  ";
            opserr << "(with at least two friction-velocity points)";
            return TCL_ERROR;
        }    
        
        int tag, numVelPts, numFrnPts, i;
        double velData[64];
        double frnData[64];
        
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
            opserr << "WARNING invalid VelDepMultiLinear friction model tag\n";
            return TCL_ERROR;
        }
        
        // get velocity data points
        i = 3;
        if (strcmp(argv[i],"-vel") == 0)  {
            i++;
            numVelPts = 0;
            while (i < argc && strcmp(argv[i],"-frn") != 0)  {
                if (Tcl_GetDouble(interp, argv[i], (velData+numVelPts)) != TCL_OK)  {
                    opserr << "WARNING invalid velocity value\n";
                    opserr << "VelDepMultiLinear friction model: " << tag << endln;
                    return TCL_ERROR;
                }
                numVelPts++;
                i++;
            }
        } else  {
            opserr << "WARNING expecting -vel but got " << argv[i] << endln;
            opserr << "VelDepMultiLinear friction model: " << tag << endln;
            return TCL_ERROR;
        }
        Vector velocityPts(velData,numVelPts);
        
        // get friction data points
        if (strcmp(argv[i],"-frn") == 0)  {
            i++;
            numFrnPts = 0;
            while (i < argc)  {
                if (Tcl_GetDouble(interp, argv[i], (frnData+numFrnPts)) != TCL_OK)  {
                    opserr << "WARNING invalid friction value\n";
                    opserr << "VelDepMultiLinear friction model: " << tag << endln;
                    return TCL_ERROR;
                }
                numFrnPts++;
                i++;
            }
        } else  {
            opserr << "WARNING expecting -frn but got " << argv[i] << endln;
            opserr << "VelDepMultiLinear friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (numVelPts != numFrnPts)  {
            opserr << "WARNING velocity and friction arrays have different length\n";
            opserr << "VelDepMultiLinear friction model: " << tag << endln;
            return TCL_ERROR;
        }
        Vector frictionPts(frnData,numFrnPts);
        
        // parsing was successful, allocate the friction model
        theFrnMdl = new VelDepMultiLinear(tag, velocityPts, frictionPts);
    }
    
    // ----------------------------------------------------------------------------	
    if (strcmp(argv[1],"VelNormalFrcDep") == 0 || strcmp(argv[1],"VNDependent") == 0)  {
        if (argc != 11)  {
            opserr << "WARNING invalid number of arguments\n";
            printCommand(argc,argv);
            opserr << "Want: frictionModel VelNormalFrcDep tag aSlow nSlow aFast nFast alpha0 alpha1 alpha2 maxMuFact\n";
            return TCL_ERROR;
        }    
        
        int tag;
        double aSlow, nSlow, aFast, nFast;
        double alpha0, alpha1, alpha2, maxMuFact;
        
        if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK)  {
            opserr << "WARNING invalid VelNormalFrcDep friction model tag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[3], &aSlow) != TCL_OK)  {
            opserr << "WARNING invalid aSlow\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4], &nSlow) != TCL_OK)  {
            opserr << "WARNING invalid nSlow\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5], &aFast) != TCL_OK)  {
            opserr << "WARNING invalid aFast\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6], &nFast) != TCL_OK)  {
            opserr << "WARNING invalid nFast\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[7], &alpha0) != TCL_OK)  {
            opserr << "WARNING invalid alpha0\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[8], &alpha1) != TCL_OK)  {
            opserr << "WARNING invalid alpha1\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[9], &alpha2) != TCL_OK)  {
            opserr << "WARNING invalid alpha2\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[10], &maxMuFact) != TCL_OK)  {
            opserr << "WARNING invalid maxMuFact\n";
            opserr << "VelNormalFrcDep friction model: " << tag << endln;
            return TCL_ERROR;
        }
        
        // parsing was successful, allocate the friction model
        theFrnMdl = new VelNormalFrcDep(tag, aSlow, nSlow, aFast, nFast,
            alpha0, alpha1, alpha2, maxMuFact);
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
