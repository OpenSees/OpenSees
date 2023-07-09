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

// $Revision: 4945 $
// $Date: 2012-07-27 16:05:47 -0700 (Fri, 27 Jul 2012) $
// $URL: svn://opensees.berkeley.edu/usr/local/svn/OpenSees/trunk/SRC/element/elastomericBearing/TclElastomericBearingBoucWenCommand.cpp $

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 02/06
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the elastomericBearingBoucWen element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <ID.h>
#include <Vector.h>

#include <ElastomericBearingBoucWen2d.h>
#include <ElastomericBearingBoucWen3d.h>
#include <UniaxialMaterial.h>

extern void printCommand(int argc, TCL_Char **argv);


int TclModelBuilder_addElastomericBearingBoucWen(ClientData clientData,
    Tcl_Interp *interp, int argc, TCL_Char **argv, Domain *theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - elastomericBearingBoucWen\n";
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
            opserr << ", for plane problem need 3 - elastomericBearingBoucWen\n";
            return TCL_ERROR;
        } 
        
        // check the number of arguments is correct
        if ((argc-eleArgStart) < 16)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: elastomericBearingBoucWen eleTag iNode jNode kInit qd alpha1 alpha2 mu eta beta gamma -P matTag -Mz matTag <-orient x1 x2 x3 y1 y2 y3> <-shearDist sDratio> <-doRayleigh> <-mass m> <-iter maxIter tol>\n";
            return TCL_ERROR;
        }    
        
        // get the id and end nodes
        int iNode, jNode, matTag, argi, i, j;
        int recvMat = 0;
        double kInit, qd, alpha1;
        double alpha2 = 0.0;
        double mu = 2.0;
        double eta = 1.0;
        double beta = 0.5;
        double gamma = 0.5;
        double shearDistI = 0.5;
        int doRayleigh = 0;
        double mass = 0.0;
        int maxIter = 25;
        double tol = 1E-12;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid elastomericBearingBoucWen eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4+eleArgStart], &kInit) != TCL_OK)  {
            opserr << "WARNING invalid kInit\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &qd) != TCL_OK)  {
            opserr << "WARNING invalid qd\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6+eleArgStart], &alpha1) != TCL_OK)  {
            opserr << "WARNING invalid alpha1\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[7+eleArgStart], &alpha2) != TCL_OK)  {
            opserr << "WARNING invalid alpha2\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[8+eleArgStart], &mu) != TCL_OK)  {
            opserr << "WARNING invalid mu\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[9+eleArgStart], &eta) != TCL_OK)  {
            opserr << "WARNING invalid eta\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[10+eleArgStart], &beta) != TCL_OK)  {
            opserr << "WARNING invalid beta\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[11+eleArgStart], &gamma) != TCL_OK)  {
            opserr << "WARNING invalid gamma\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[2];
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                theMaterials[0] = 0;
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 2)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 2 materials\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // check for optional arguments
        Vector x = 0;
        Vector y = 0;
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-shearDist") != 0 &&
                    strcmp(argv[j],"-doRayleigh") != 0 &&
                    strcmp(argv[j],"-mass") != 0 &&
                    strcmp(argv[j],"-iter") != 0)  {
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
                            opserr << "elastomericBearingBoucWen element: " << tag << endln;
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
                            opserr << "elastomericBearingBoucWen element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-shearDist") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &shearDistI) != TCL_OK)  {
                    opserr << "WARNING invalid -shearDist value\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 12+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i], "-doRayleigh") == 0)
                doRayleigh = 1;
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid -mass value\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 12+eleArgStart; i < argc; i++)  {
            if (i+2 < argc && strcmp(argv[i], "-iter") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &maxIter) != TCL_OK)  {
                    opserr << "WARNING invalid maxIter\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                if (Tcl_GetDouble(interp, argv[i+2], &tol) != TCL_OK)  {
                    opserr << "WARNING invalid tol\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the elastomericBearingBoucWen
        theElement = new ElastomericBearingBoucWen2d(tag, iNode, jNode, kInit, qd,
            alpha1, theMaterials, y, x, alpha2, mu, eta, beta, gamma, shearDistI,
            doRayleigh, mass, maxIter, tol);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the elastomericBearingBoucWen to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }       
    }
    
    else if (ndm == 3)  {
        // check space frame problem has 6 dof per node
        if (ndf != 6)  {
            opserr << "WARNING invalid ndf: " << ndf;
            opserr << ", for space problem need 6 - elastomericBearingBoucWen \n";
            return TCL_ERROR;
        } 
        
        // check the number of arguments is correct
        if ((argc-eleArgStart) < 18)  {
            opserr << "WARNING insufficient arguments\n";
            printCommand(argc, argv);
            opserr << "Want: elastomericBearingBoucWen eleTag iNode jNode kInit qd alpha1 alpha2 mu eta beta gamma -P matTag -T matTag -My matTag -Mz matTag <-orient <x1 x2 x3> y1 y2 y3> <-shearDist sDratio> <-mass m> <-iter maxIter tol>\n";
            return TCL_ERROR;
        }
        
        // get the id and end nodes
        int iNode, jNode, matTag, argi, i, j;
        int recvMat = 0;
        double kInit, qd, alpha1;
        double alpha2 = 0.0;
        double mu = 2.0;
        double eta = 1.0;
        double beta = 0.5;
        double gamma = 0.5;
        double shearDistI = 0.5;
        int doRayleigh = 0;
        double mass = 0.0;
        int maxIter = 25;
        double tol = 1E-12;
        
        if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
            opserr << "WARNING invalid elastomericBearingBoucWen eleTag\n";
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[2+eleArgStart], &iNode) != TCL_OK)  {
            opserr << "WARNING invalid iNode\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetInt(interp, argv[3+eleArgStart], &jNode) != TCL_OK)  {
            opserr << "WARNING invalid jNode\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[4+eleArgStart], &kInit) != TCL_OK)  {
            opserr << "WARNING invalid kInit\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[5+eleArgStart], &qd) != TCL_OK)  {
            opserr << "WARNING invalid qd\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[6+eleArgStart], &alpha1) != TCL_OK)  {
            opserr << "WARNING invalid alpha1\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[7+eleArgStart], &alpha2) != TCL_OK)  {
            opserr << "WARNING invalid alpha2\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[8+eleArgStart], &mu) != TCL_OK)  {
            opserr << "WARNING invalid mu\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[9+eleArgStart], &eta) != TCL_OK)  {
            opserr << "WARNING invalid eta\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[10+eleArgStart], &beta) != TCL_OK)  {
            opserr << "WARNING invalid beta\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[11+eleArgStart], &gamma) != TCL_OK)  {
            opserr << "WARNING invalid gamma\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        UniaxialMaterial *theMaterials[4];
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-P") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid axial matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[0] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[0] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-T") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid torsional matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[1] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[1] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-My") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment y matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[2] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[2] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-Mz") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &matTag) != TCL_OK)  {
                    opserr << "WARNING invalid moment z matTag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                theMaterials[3] = OPS_getUniaxialMaterial(matTag);
                if (theMaterials[3] == 0)  {
                    opserr << "WARNING material model not found\n";
                    opserr << "uniaxialMaterial: " << matTag << endln;
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                recvMat++;
            }
        }
        if (recvMat != 4)  {
            opserr << "WARNING wrong number of materials\n";
            opserr << "got " << recvMat << " materials, but want 4 materials\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // check for optional arguments
        Vector x(0);
        Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (strcmp(argv[i],"-orient") == 0)  {
                j = i+1;
                int numOrient = 0;
                while (j < argc &&
                    strcmp(argv[j],"-shearDist") != 0 &&
                    strcmp(argv[j],"-doRayleigh") != 0 &&
                    strcmp(argv[j],"-mass") != 0 &&
                    strcmp(argv[j],"-iter") != 0)  {
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
                            opserr << "elastomericBearingBoucWen element: " << tag << endln;
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
                            opserr << "elastomericBearingBoucWen element: " << tag << endln;
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
                            opserr << "elastomericBearingBoucWen element: " << tag << endln;
                            return TCL_ERROR;
                        } else  {
                            argi++;
                            y(j) = value;		
                        }
                    }
                }
                else  {
                    opserr << "WARNING insufficient arguments after -orient flag\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-shearDist") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &shearDistI) != TCL_OK)  {
                    opserr << "WARNING invalid -shearDist value\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-doRayleigh") == 0)
                doRayleigh = 1;
        }
        for (i = 12+eleArgStart; i < argc; i++)  {
            if (i+1 < argc && strcmp(argv[i], "-mass") == 0)  {
                if (Tcl_GetDouble(interp, argv[i+1], &mass) != TCL_OK)  {
                    opserr << "WARNING invalid -mass value\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        for (int i = 12+eleArgStart; i < argc; i++)  {
            if (i+2 < argc && strcmp(argv[i], "-iter") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &maxIter) != TCL_OK)  {
                    opserr << "WARNING invalid maxIter\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
                if (Tcl_GetDouble(interp, argv[i+2], &tol) != TCL_OK)  {
                    opserr << "WARNING invalid tol\n";
                    opserr << "elastomericBearingBoucWen element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
        
        // now create the elastomericBearingBoucWen
        theElement = new ElastomericBearingBoucWen3d(tag, iNode, jNode, kInit, qd,
            alpha1, theMaterials, y, x, alpha2, mu, eta, beta, gamma, shearDistI,
            doRayleigh, mass, maxIter, tol);
        
        if (theElement == 0)  {
            opserr << "WARNING ran out of memory creating element\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            return TCL_ERROR;
        }
        
        // then add the elastomericBearingBoucWen to the domain
        if (theTclDomain->addElement(theElement) == false)  {
            opserr << "WARNING could not add element to the domain\n";
            opserr << "elastomericBearingBoucWen element: " << tag << endln;
            delete theElement;
            return TCL_ERROR;
        }       
    }
    
    else  {
        opserr << "WARNING elastomericBearingBoucWen command only works when ndm is 2 or 3, ndm: ";
        opserr << ndm << endln;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the elastomericBearingBoucWen and added it to the domain
    return TCL_OK;
}
