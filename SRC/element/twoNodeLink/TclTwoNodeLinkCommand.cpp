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
// Created: 08/08
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the TwoNodeLink element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>
#include <UniaxialMaterial.h>

#include <TwoNodeLink.h>


extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addTwoNodeLink(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - twoNodeLink\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    int ndf = theTclBuilder->getNDF();
    
    // check the number of arguments is correct
    if ((argc-eleArgStart) < 8)  {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: twoNodeLink eleTag iNode jNode -mat matTags -dir dirs <-orient <x1 x2 x3> y1 y2 y3> <-pDelta Mratios> <-shearDist sDratios> <-doRayleigh> <-mass m>\n";
        return TCL_ERROR;
    }
    
    // get the id and end nodes 
    int tag, iNode, jNode, numMat, matTag, numDir, dirID, i, j;
    int argi = eleArgStart+1;
    Vector Mratio(0), shearDistI(0);
    int doRayleigh = 0;
    double mass = 0.0;
    
    if (Tcl_GetInt(interp, argv[argi], &tag) != TCL_OK)  {
        opserr << "WARNING invalid twoNodeLink eleTag\n";
        return TCL_ERROR;
    }
    argi++;
    if (Tcl_GetInt(interp, argv[argi], &iNode) != TCL_OK)  {
        opserr << "WARNING invalid iNode\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    if (Tcl_GetInt(interp, argv[argi], &jNode) != TCL_OK)  {
        opserr << "WARNING invalid jNode\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    // read the number of materials
    numMat = 0;
    if (strcmp(argv[argi],"-mat") != 0)  {
        opserr << "WARNING expecting -mat flag\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;	
    }
    argi++;
    i = argi;
    while (i < argc &&
        strcmp(argv[i],"-dir") != 0)  {
            numMat++;
            i++;
    }
    if (numMat == 0)  {
        opserr << "WARNING no directions specified\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    // create array of uniaxial materials
    UniaxialMaterial **theMaterials = new UniaxialMaterial* [numMat];
    for (i=0; i<numMat; i++)  {
        theMaterials[i] = 0;
        if (Tcl_GetInt(interp, argv[argi], &matTag) != TCL_OK)  {
            opserr << "WARNING invalid matTag\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
        }
        theMaterials[i] = OPS_getUniaxialMaterial(matTag);
        if (theMaterials[i] == 0)  {
            opserr << "WARNING material model not found\n";
            opserr << "uniaxialMaterial " << matTag << endln;
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
        }
        argi++;
    }
    // read the number of directions
    numDir = 0;
    if (strcmp(argv[argi],"-dir") != 0)  {
        opserr << "WARNING expecting -dir flag\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    i = argi;
    while (i < argc &&
        strcmp(argv[i],"-orient") != 0 &&
        strcmp(argv[i],"-pDelta") != 0 &&
        strcmp(argv[i],"-shearDist") != 0 &&
        strcmp(argv[i],"-doRayleigh") != 0 &&
        strcmp(argv[i],"-mass") != 0)  {
            numDir++;
            i++;
    }
    if (numDir != numMat)  {
        opserr << "WARNING wrong number of directions specified\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    // create the ID array to hold the direction IDs
    ID theDirIDs(numDir);
    // fill in the directions
    for (i=0; i<numDir; i++)  {
        if (Tcl_GetInt(interp, argv[argi], &dirID) != TCL_OK)  {
            opserr << "WARNING invalid direction ID\n";
            opserr << "twoNodeLink element: " << tag << endln;
            return TCL_ERROR;
        }
        if (dirID < 1 || dirID > ndf)  {
            opserr << "WARNING invalid direction ID: ";
            opserr << "dir = " << dirID << " > ndf = " << ndf;
            opserr << "\ntwoNodeLink element: " << tag << endln;
            return TCL_ERROR;
        }
        theDirIDs(i) = dirID-1;
        argi++;
    }
    
    // check for optional arguments
    Vector x(0), y(0);
    for (i = argi; i < argc; i++)  {
        if (strcmp(argv[i],"-orient") == 0)  {
            j = i+1;
            int numOrient = 0;
            while (j < argc &&
                strcmp(argv[j],"-pDelta") != 0 &&
                strcmp(argv[j],"-shearDist") != 0 &&
                strcmp(argv[j],"-doRayleigh") != 0 &&
                strcmp(argv[j],"-mass") != 0)  {
                    numOrient++;
                    j++;
            }
            if (numOrient == 3)  {
                y.resize(3);
                double value;
                // read the y values
                for (j=0; j<3; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+1+j], &value) != TCL_OK)  {
                        opserr << "WARNING invalid -orient value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    } else  {
                        y(j) = value;
                    }
                }
            }
            else if (numOrient == 6)  {
                x.resize(3);
                y.resize(3);
                double value;
                // read the x values
                for (j=0; j<3; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+1+j], &value) != TCL_OK)  {
                        opserr << "WARNING invalid -orient value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    } else  {
                        x(j) = value;
                    }
                }
                // read the y values
                for (j=0; j<3; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+4+j], &value) != TCL_OK)  {
                        opserr << "WARNING invalid -orient value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    } else  {
                        y(j) = value;
                    }
                }
            }
            else  {
                opserr << "WARNING insufficient arguments after -orient flag\n";
                opserr << "twoNodeLink element: " << tag << endln;
                return TCL_ERROR;
            }
        }
    }
    for (i=argi; i<argc; i++)  {
        if (i+1 < argc && strcmp(argv[i], "-pDelta") == 0)  {
            double Mr;
            Mratio.resize(4);
            if (ndm == 2)  {
                Mratio.Zero();
                for (j=0; j<2; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+1+j], &Mr) != TCL_OK)  {
                        opserr << "WARNING invalid -pDelta value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    }
                    Mratio(2+j) = Mr;
                }
            } else if (ndm == 3)  {
                for (j=0; j<4; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+1+j], &Mr) != TCL_OK)  {
                        opserr << "WARNING invalid -pDelta value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    }
                    Mratio(j) = Mr;
                }
            }
        }
    }
    for (i=argi; i<argc; i++)  {
        if (i+1 < argc && strcmp(argv[i], "-shearDist") == 0)  {
            double sDI;
            shearDistI.resize(2);
            if (ndm == 2)  {
                if (Tcl_GetDouble(interp, argv[i+1], &sDI) != TCL_OK)  {
                    opserr << "WARNING invalid -shearDist value\n";
                    opserr << "twoNodeLink element: " << tag << endln;
                    return TCL_ERROR;
                }
                shearDistI(0) = sDI;
                shearDistI(1) = 0.5;
            } else if (ndm == 3)  {
                for (j=0; j<2; j++)  {
                    if (Tcl_GetDouble(interp, argv[i+1+j], &sDI) != TCL_OK)  {
                        opserr << "WARNING invalid -shearDist value\n";
                        opserr << "twoNodeLink element: " << tag << endln;
                        return TCL_ERROR;
                    }
                    shearDistI(j) = sDI;
                }
            }
        }
    }
    for (i=argi; i<argc; i++)  {
        if (strcmp(argv[i], "-doRayleigh") == 0)
            doRayleigh = 1;
    }
    for (i=argi; i<argc; i++)  {
        if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
            if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                opserr << "WARNING invalid -mass value\n";
                opserr << "twoNodeLink element: " << tag << endln;
                return TCL_ERROR;
            }
        }
    }
    
    // now create the twoNodeLink
    theElement = new TwoNodeLink(tag, ndm, iNode, jNode, theDirIDs,
        theMaterials, y, x, Mratio, shearDistI, doRayleigh, mass);
    
    // cleanup dynamic memory
    if (theMaterials != 0)
        delete [] theMaterials;
    
    if (theElement == 0)  {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "twoNodeLink element: " << tag << endln;
        return TCL_ERROR;
    }
    
    // then add the twoNodeLink to the domain
    if (theTclDomain->addElement(theElement) == false)  {
        opserr << "WARNING could not add element to the domain\n";
        opserr << "twoNodeLink element: " << tag << endln;
        delete theElement;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the twoNodeLink and added it to the domain
    return TCL_OK;
}
