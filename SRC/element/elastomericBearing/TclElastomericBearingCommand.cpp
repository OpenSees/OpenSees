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

// $Revision: 1.3 $
// $Date: 2009-04-17 23:00:48 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/elastomericBearing/TclElastomericBearingCommand.cpp,v $

// Written: Andreas Schellenberg (andreas.schellenberg@gmx.net)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the elastomericBearing element.
//
// What: "@(#) TclElastomericBearingCommand.cpp, revA"


#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <ElastomericBearing2d.h>
#include <ElastomericBearing3d.h>


extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addElastomericBearing(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - elastomericBearing\n";    
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
            opserr << ", for plane problem need 3 - elastomericBearing\n";    
            return TCL_ERROR;
        } 
        
        // check the number of arguments is correct
        if ((argc-eleArgStart) < 11)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: elastomericBearing eleTag iNode jNode ke fy alpha -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-mass m>\n";
            return TCL_ERROR;
        }    
        
        // get the id and end nodes
        int iNode, jNode, matTag, argi, i, j;
        int recvMat = 0;
        double ke, fy, alpha;
        double mass = 0.0;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid elastomericBearing eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4+eleArgStart], &ke) != TCL_OK)  {
            opserr << "WARNING invalid ke\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &fy) != TCL_OK)  {
            opserr << "WARNING invalid fy\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6+eleArgStart], &alpha) != TCL_OK)  {
            opserr << "WARNING invalid alpha\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[2];
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                theMaterials[0] = 0;
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 2)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 2 materials\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
       
        // check for optional arguments
        Vector x = 0;
        Vector y = 0;
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-mass") != 0)  {
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
                            opserr << "elastomericBearing element: " << tag << endln;
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
                            opserr << "elastomericBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;		
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid mass\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the elastomericBearing
        theElement = new ElastomericBearing2d(tag, iNode, jNode, ke, fy, alpha, theMaterials, y, x, mass);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the elastomericBearing to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "elastomericBearing element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }       
    }

    else if (ndm == 3)  {
        // check space frame problem has 6 dof per node
        if (ndf != 6)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for space problem need 6 - elastomericBearing \n";    
            return TCL_ERROR;
        } 

        // check the number of arguments is correct
        if ((argc-eleArgStart) < 15)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: elastomericBearing eleTag iNode jNode ke fy alpha -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-mass m>\n";
            return TCL_ERROR;
        }
        
        // get the id and end nodes
        int iNode, jNode, matTag, argi, i, j;
        int recvMat = 0;
        double ke, fy, alpha;
        double mass = 0.0;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid elastomericBearing eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4+eleArgStart], &ke) != TCL_OK)  {
            opserr << "WARNING invalid ke\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &fy) != TCL_OK)  {
            opserr << "WARNING invalid fy\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6+eleArgStart], &alpha) != TCL_OK)  {
            opserr << "WARNING invalid alpha\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[4];
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid axial matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-T") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid torsional matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-My") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment y matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[2] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[2] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment z matTag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[3] = theTclBuilder->getUniaxialMaterial(matTag);
                if (theMaterials[3] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 4)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 4 materials\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // check for optional arguments
        Vector x(0);
        Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-mass") != 0)  {
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
                            opserr << "elastomericBearing element: " << tag << endln;
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
                            opserr << "elastomericBearing element: " << tag << endln;
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
                            opserr << "elastomericBearing element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;		
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 7+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid mass\n";
                    opserr << "elastomericBearing element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the elastomericBearing
        theElement = new ElastomericBearing3d(tag, iNode, jNode, ke, fy, alpha, theMaterials, y, x, mass);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "elastomericBearing element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the elastomericBearing to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "elastomericBearing element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }       
    }
    
    else  {
        opserr << "WARNING elastomericBearing command only works when ndm is 2 or 3, ndm: ";
        opserr << ndm << endln;
        return TCL_ERROR;
    }
    
    // if get here we have sucessfully created the elastomericBearing and added it to the domain
    return TCL_OK;
}
