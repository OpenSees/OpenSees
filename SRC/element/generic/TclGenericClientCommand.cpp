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
// Created: 11/06
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the genericClient element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <GenericClient.h>


extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addGenericClient(ClientData clientData, Tcl_Interp *interp,  int argc, 
    TCL_Char **argv, Domain*theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - genericClient\n";
        return TCL_ERROR;
    }
    
    // check the number of arguments is correct
    if ((argc-eleArgStart) < 8)  {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: element genericClient eleTag -node Ndi Ndj ... -dof dofNdi -dof dofNdj ... -server ipPort <ipAddr> <-ssl> <-udp> <-dataSize size> <-noRayleigh>\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    
    // get the id and end nodes
    int tag, node, dof, ipPort, argi, i, j;
    int numNodes = 0, numDOFj = 0, numDOF = 0;
    char *ipAddr = 0;
    int ssl = 0, udp = 0;
    int dataSize = 256;
    int doRayleigh = 1;
    
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
        opserr << "WARNING invalid genericClient eleTag\n";
        return TCL_ERROR;
    }
    // read the number of nodes
    if (strcmp(argv[2+eleArgStart], "-node") != 0)  {
        opserr << "WARNING expecting -node flag\n";
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
    }
    argi = 3+eleArgStart;
    i = argi;
    while (strcmp(argv[i], "-dof") != 0  && i < argc)  {
        numNodes++;
        i++;
    }
    if (numNodes == 0)  {
        opserr << "WARNING no nodes specified\n";
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
    }
    // create the ID arrays to hold the nodes and dofs
    ID nodes(numNodes);
    ID *dofs = new ID [numNodes];
    if (dofs == 0)  {
        opserr << "WARNING out of memory\n";
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
    }
    // fill in the nodes ID
    for (i=0; i<numNodes; i++)  {
        if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK)  {
            opserr << "WARNING invalid node\n";
            opserr << "genericClient element: " << tag << endln;
            return TCL_ERROR;
        }
        nodes(i) = node;
        argi++; 
    }
    for (j=0; j<numNodes; j++)  {
        // read the number of dofs per node j
        numDOFj = 0;
        if (strcmp(argv[argi], "-dof") != 0)  {
            opserr << "WARNING expect -dof\n";
            opserr << "genericClient element: " << tag << endln;
            return TCL_ERROR;
        }
        argi++;
        i = argi;
        while (strcmp(argv[i], "-dof") != 0 &&
            strcmp(argv[i], "-server") != 0 &&
            strcmp(argv[i], "-doRayleigh") != 0 &&
            strcmp(argv[i], "-noRayleigh") != 0 &&
            i < argc)  {
                numDOFj++;
                numDOF++;
                i++;
        }
        // fill in the dofs ID array
        ID dofsj(numDOFj);
        for (i=0; i<numDOFj; i++)  {
            if (Tcl_GetInt(interp, argv[argi], &dof) != TCL_OK)  {
                opserr << "WARNING invalid dof\n";
                opserr << "genericClient element: " << tag << endln;
                return TCL_ERROR;
            }
            dofsj(i) = dof-1;
            argi++;
        }
        dofs[j] = dofsj;
    }
    if (strcmp(argv[argi], "-server") == 0)  {
        argi++;
        if (Tcl_GetInt(interp, argv[argi], &ipPort) != TCL_OK)  {
            opserr << "WARNING invalid ipPort\n";
            opserr << "genericClient element: " << tag << endln;
            return TCL_ERROR;
        }
        argi++;
        if (argi < argc &&
            strcmp(argv[argi], "-doRayleigh") != 0 &&
            strcmp(argv[argi], "-noRayleigh") != 0 &&
            strcmp(argv[argi], "-dataSize") != 0 &&
            strcmp(argv[argi], "-ssl") != 0 &&
            strcmp(argv[argi], "-udp") != 0)  {
                ipAddr = new char [strlen(argv[argi])+1];
                strcpy(ipAddr,argv[argi]);
                argi++;
        }
        else  {
            ipAddr = new char [9+1];
            strcpy(ipAddr,"127.0.0.1");
        }
        for (i = argi; i < argc; i++)  {
            if (strcmp(argv[i], "-ssl") == 0)  {
                ssl = 1; udp = 0;
            }
            else if (strcmp(argv[i], "-udp") == 0)  {
                udp = 1; ssl = 0;
            }
            else if (strcmp(argv[i], "-dataSize") == 0)  {
                if (Tcl_GetInt(interp, argv[i+1], &dataSize) != TCL_OK)  {
                    opserr << "WARNING invalid dataSize\n";
                    opserr << "genericClient element: " << tag << endln;
                    return TCL_ERROR;
                }
            }
        }
    }
    else  {
        opserr << "WARNING expecting -server string but got ";
        opserr << argv[argi] << endln;
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
    }
    for (i = argi; i < argc; i++)  {
        if (strcmp(argv[i], "-doRayleigh") == 0)  {
            doRayleigh = 1;
        } else if (strcmp(argv[i], "-noRayleigh") == 0)  {
            doRayleigh = 0;
        }
    }
    
    // now create the GenericClient
    theElement = new GenericClient(tag, nodes, dofs, ipPort, ipAddr,
        ssl, udp, dataSize, doRayleigh);
    
    // cleanup dynamic memory
    if (dofs != 0)
        delete [] dofs;
    
    if (theElement == 0)  {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "genericClient element: " << tag << endln;
        return TCL_ERROR;
    }
    
    // then add the GenericClient to the domain
    if (theTclDomain->addElement(theElement) == false)  {
        opserr << "WARNING could not add element to the domain\n";
        opserr << "genericClient element: " << tag << endln;
        delete theElement;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the genericClient and added it to the domain
    return TCL_OK;
}
