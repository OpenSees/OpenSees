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
// Created: 12/13
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the RJWatsonEqsBearing element.


#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <RJWatsonEQS2d.h>
#include <RJWatsonEQS3d.h>
#include <FrictionModel.h>
#include <UniaxialMaterial.h>

extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addRJWatsonEqsBearing(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - RJWatsonEqsBearing\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    int ndf = theTclBuilder->getNDF();
    int tag;
    
    if (ndm == 2)  {
        // check plane frame problem has 3 dof per node
        if (ndf != 3)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for plane problem need 3 - RJWatsonEqsBearing\n";
            return TCL_ERROR;
        }
        
        // check the number of arguments is correct
        if ((argc-eleArgStart) < 12)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: RJWatsonEqsBearing eleTag iNode jNode frnMdlTag kInit -P matTag -Vy matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
            return TCL_ERROR;
        }
        
        // get the id and end nodes
        int iNode, jNode, frnMdlTag, matTag, argi, i, j;
        int recvMat = 0;
        double kInit;
        double shearDistI = 1.0;
        int doRayleigh = 0;
        double mass = 0.0;
        int maxIter = 25;
        double tol = 1E-12;
        double kFactUplift = 1E-12;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid RJWatsonEqsBearing eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[4+eleArgStart], &frnMdlTag) != TCL_OK)  {
            opserr << "WARNING invalid frnMdlTag\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        FrictionModel *theFrnMdl = OPS_getFrictionModel(frnMdlTag);
        if (theFrnMdl == 0)  {
            opserr << "WARNING friction model not found\n";
            opserr << "frictionModel: " << frnMdlTag << endln;
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &kInit) != TCL_OK)  {
            opserr << "WARNING invalid kInit\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[3];
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                theMaterials[0] = 0;
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid axial matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++) {
            if (i+1 < argc && strcmp(argv[i], "-Vy") == 0) {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK) {
                    opserr << "WARNING invalid shear y matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0) {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment z matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[2] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[2] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 3)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 3 materials\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // check for optional arguments
        Vector x = 0;
        Vector y = 0;
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-shearDist") != 0 &&
                    strcmp(argv[j],"-doRayleigh") != 0 &&
                    strcmp(argv[j],"-mass") != 0 &&
                    strcmp(argv[j],"-iter") != 0 &&
                    strcmp(argv[j],"-kFactUplift") != 0)  {
                    numOrient++;
                    j++;
                }
                if (numOrient == 6)  {
                    argi = i+1;
                    x.resize(3);
                    y.resize(3);
                    double value;
                    // read the x values
                    for (j=0; j<3; j++)  {
                        if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK)  {
                            opserr << "WARNING invalid -orient value\n";
                            opserr << "RJWatsonEqsBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            x(j) = value;
                        }
                    }
                    // read the y values
                    for (j=0; j<3; j++)  {
                        if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK)  {
                            opserr << "WARNING invalid -orient value\n";
                            opserr << "RJWatsonEqsBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-shearDist") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &shearDistI) != TCL_OK)  {
                    opserr << "WARNING invalid -shearDist value\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i], "-doRayleigh") == 0)
                doRayleigh = 1;
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid -mass value\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+2 < argc && strcmp(argv[i], "-iter") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &maxIter) != TCL_OK)  {
                    opserr << "WARNING invalid maxIter\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                if (Tcl_GetDouble(interp, argv[i+2], &tol) != TCL_OK)  {
                    opserr << "WARNING invalid tol\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-kFactUplift") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &kFactUplift) != TCL_OK)  {
                    opserr << "WARNING invalid kFactUplift\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the RJWatsonEqsBearing
        theElement = new RJWatsonEQS2d(tag, iNode, jNode, *theFrnMdl, kInit,
            theMaterials, y, x, shearDistI, doRayleigh, mass, maxIter, tol,
            kFactUplift);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the RJWatsonEqsBearing to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }
    }
    
    else if (ndm == 3)  {
        // check space frame problem has 6 dof per node
        if (ndf != 6)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for space problem need 6 - RJWatsonEqsBearing \n";
            return TCL_ERROR;
        }
        
        // check the number of arguments is correct
        if ((argc-eleArgStart) < 18)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: RJWatsonEqsBearing eleTag iNode jNode frnMdlTag kInit -P matTag -Vy matTag -Vz matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
            return TCL_ERROR;
        }
        
        // get the id and end nodes 
        int iNode, jNode, frnMdlTag, matTag, argi, i, j;
        int recvMat = 0;
        double kInit;
        double shearDistI = 1.0;
        int doRayleigh = 0;
        double mass = 0.0;
        int maxIter = 25;
        double tol = 1E-12;
        double kFactUplift = 1E-12;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid RJWatsonEqsBearing eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[4+eleArgStart], &frnMdlTag) != TCL_OK)  {
            opserr << "WARNING invalid frnMdlTag\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        FrictionModel *theFrnMdl = OPS_getFrictionModel(frnMdlTag);
        if (theFrnMdl == 0)  {
            opserr << "WARNING friction model not found\n";
            opserr << "frictionModel: " << frnMdlTag << endln;
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &kInit) != TCL_OK)  {
            opserr << "WARNING invalid kInit\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[6];
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid axial matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Vy") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid shear y matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Vz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid shear z matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[2] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[2] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-T") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid torsional matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[3] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[3] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-My") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment y matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[4] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[4] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment z matTag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[5] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[5] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 6)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 6 materials\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // check for optional arguments
        Vector x(0);
        Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-shearDist") != 0 &&
                    strcmp(argv[j],"-doRayleigh") != 0 &&
                    strcmp(argv[j],"-mass") != 0 &&
                    strcmp(argv[j],"-iter") != 0 &&
                    strcmp(argv[j],"-kFactUplift") != 0)  {
                    numOrient++;
                    j++;
                }
                if (numOrient == 3)  {
                    argi = i+1;
                    double value;
                    // read the y values
                    for (j=0; j<3; j++)  {
                        if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK)  {
                            opserr << "WARNING invalid -orient value\n";
                            opserr << "RJWatsonEqsBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;
                        }
                    }
                }
                else if (numOrient == 6)  {
                    argi = i+1;
                    x.resize(3);
                    double value;
                    // read the x values
                    for (j=0; j<3; j++)  {
                        if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK)  {
                            opserr << "WARNING invalid -orient value\n";
                            opserr << "RJWatsonEqsBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            x(j) = value;
                        }
                    }
                    // read the y values
                    for (j=0; j<3; j++)  {
                        if (Tcl_GetDouble(interp, argv[argi], &value) != TCL_OK)  {
                            opserr << "WARNING invalid -orient value\n";
                            opserr << "RJWatsonEqsBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-shearDist") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &shearDistI) != TCL_OK)  {
                    opserr << "WARNING invalid -shearDist value\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-doRayleigh") == 0)
                doRayleigh = 1;
        }
        for (i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid -mass value\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+2 < argc && strcmp(argv[i], "-iter") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &maxIter) != TCL_OK)  {
                    opserr << "WARNING invalid maxIter\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                if (Tcl_GetDouble(interp, argv[i+2], &tol) != TCL_OK)  {
                    opserr << "WARNING invalid tol\n";
                    opserr << "RJWatsonEqsBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 6+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-kFactUplift") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &kFactUplift) != TCL_OK)  {
                    opserr << "WARNING invalid kFactUplift\n";
                    opserr << "singleFPBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the RJWatsonEqsBearing
        theElement = new RJWatsonEQS3d(tag, iNode, jNode, *theFrnMdl, kInit,
            theMaterials, y, x, shearDistI, doRayleigh, mass, maxIter, tol,
            kFactUplift);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the RJWatsonEqsBearing to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "RJWatsonEqsBearing element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }
    }
    
    else  {
        opserr << "WARNING RJWatsonEqsBearing command only works when ndm is 2 or 3, ndm: ";
        opserr << ndm << endln;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the RJWatsonEqsBearing and added it to the domain
    return TCL_OK;
}
