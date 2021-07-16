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
// Created: 09/07
// Revision: A
//
// Description: This file contains the function to parse the TCL input
// for the adapter element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Adapter.h>


extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addAdapter(ClientData clientData, Tcl_Interp *interp,  int argc, 
    TCL_Char **argv, Domain*theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0) {
        opserr << "WARNING builder has been destroyed - adapter\n";
        return TCL_ERROR;
    }
    
    // check the number of arguments is correct
    if ((argc-eleArgStart) < 8) {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: element adapter eleTag -node Ndi Ndj ... -dof dofNdi -dof dofNdj ... -stif Kij ipPort <-doRayleigh> <-mass Mij>\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    
    // get the id and end nodes 
    int tag, node, dof, ipPort, argi, i, j, k;
    int numNodes = 0, numDOFj = 0, numDOF = 0;
    int doRayleigh = 0;
    Matrix *mass = 0;
    
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK) {
        opserr << "WARNING invalid adapter eleTag" << endln;
        return TCL_ERROR;
    }
    // read the number of nodes
    if (strcmp(argv[2+eleArgStart], "-node") != 0)  {
        opserr << "WARNING expecting -node flag\n";
        opserr << "adapter element: " << tag << endln;
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
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;
    }
    // create the ID arrays to hold the nodes and dofs
    ID nodes(numNodes);
    ID *dofs = new ID [numNodes];
    if (dofs == 0)  {
        opserr << "WARNING out of memory\n";
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;
    }
    // fill in the nodes ID
    for (i=0; i<numNodes; i++)  {
        if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK)  {
            opserr << "WARNING invalid node\n";
            opserr << "adapter element: " << tag << endln;
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
            opserr << "adapter element: " << tag << endln;
            return TCL_ERROR;
        }
        argi++;
        i = argi;
        while (strcmp(argv[i], "-dof") != 0 && 
            strcmp(argv[i], "-stif") != 0 && 
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
                opserr << "adapter element: " << tag << endln;
                return TCL_ERROR;
            }
            dofsj(i) = dof-1;
            argi++; 
        }
        dofs[j] = dofsj;
    }
    // get stiffness matrix
    Matrix kb(numDOF,numDOF);
    if (strcmp(argv[argi], "-stif") != 0)  {
        opserr << "WARNING expecting -stif flag\n";
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    if (argc-1 < argi+numDOF*numDOF)  {
        opserr << "WARNING incorrect number of stiffness terms\n";
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;      
    }
    double stif;
    for (j=0; j<numDOF; j++)  {
        for (k=0; k<numDOF; k++)  {
            if (Tcl_GetDouble(interp, argv[argi], &stif) != TCL_OK)  {
                opserr << "WARNING invalid stiffness term\n";
                opserr << "adapter element: " << tag << endln;
                return TCL_ERROR;
            }
            kb(j,k) = stif;
            argi++;
        }
    }
    // get ip-port
    if (Tcl_GetInt(interp, argv[argi], &ipPort) != TCL_OK) {
        opserr << "WARNING invalid ipPort\n";
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    // get optional rayleigh flag
    for (int i = argi; i < argc; i++)  {
        if (strcmp(argv[i], "-doRayleigh") == 0)
            doRayleigh = 1;
    }
    // get optional mass matrix
    for (int i = argi; i < argc; i++) {
        if (strcmp(argv[i], "-mass") == 0) {
            if (argc-1 < i+numDOF*numDOF)  {
                opserr << "WARNING incorrect number of mass terms\n";
                opserr << "adapter element: " << tag << endln;
                return TCL_ERROR;      
            }
            mass = new Matrix(numDOF,numDOF);
            double m;
            for (j=0; j<numDOF; j++)  {
                for (k=0; k<numDOF; k++)  {
                    if (Tcl_GetDouble(interp, argv[i+1 + numDOF*j+k], &m) != TCL_OK)  {
                        opserr << "WARNING invalid mass term\n";
                        opserr << "adapter element: " << tag << endln;
                        return TCL_ERROR;
                    }
                    (*mass)(j,k) = m;
                }
            }
        }
    }
    
    // now create the adapter and add it to the Domain
    if (mass == 0)
        theElement = new Adapter(tag, nodes, dofs, kb, ipPort, doRayleigh);
    else
        theElement = new Adapter(tag, nodes, dofs, kb, ipPort, doRayleigh, mass);
    
    // cleanup dynamic memory
    if (dofs != 0)
        delete [] dofs;
    
    if (theElement == 0) {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "adapter element: " << tag << endln;
        return TCL_ERROR;
    }
    
    if (theTclDomain->addElement(theElement) == false) {
        opserr << "WARNING could not add element to the domain\n";
        opserr << "adapter element: " << tag << endln;
        delete theElement;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the adapter and added it to the domain
    return TCL_OK;
}
