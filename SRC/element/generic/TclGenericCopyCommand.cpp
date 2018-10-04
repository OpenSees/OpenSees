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
// for the GenericCopy element.

#include <TclModelBuilder.h>

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <GenericCopy.h>


extern void printCommand(int argc, TCL_Char **argv);

int TclModelBuilder_addGenericCopy(ClientData clientData, Tcl_Interp *interp,  int argc, 
    TCL_Char **argv, Domain*theTclDomain,
    TclModelBuilder *theTclBuilder, int eleArgStart)
{
    // ensure the destructor has not been called
    if (theTclBuilder == 0)  {
        opserr << "WARNING builder has been destroyed - expElement genericCopy\n";
        return TCL_ERROR;
    }
    
    // check the number of arguments is correct
    if ((argc-eleArgStart) < 6)  {
        opserr << "WARNING insufficient arguments\n";
        printCommand(argc, argv);
        opserr << "Want: expElement genericCopy eleTag -node Ndi ... -src srcTag\n";
        return TCL_ERROR;
    }
    
    Element *theElement = 0;
    int ndm = theTclBuilder->getNDM();
    
    // get the id and end nodes
    int tag, node, srcTag, argi, i;
    int numNodes = 0;
    
    if (Tcl_GetInt(interp, argv[1+eleArgStart], &tag) != TCL_OK)  {
        opserr << "WARNING invalid genericCopy eleTag\n";
        return TCL_ERROR;
    }
    // read the number of nodes
    if (strcmp(argv[2+eleArgStart], "-node") != 0)  {
        opserr << "WARNING expecting -node flag\n";
        opserr << "genericCopy element: " << tag << endln;
        return TCL_ERROR;
    }
    argi = 3+eleArgStart;
    i = argi;
    while (strcmp(argv[i], "-src") != 0  && i < argc)  {
        numNodes++;
        i++;
    }
    if (numNodes == 0)  {
        opserr << "WARNING no nodes specified\n";
        opserr << "genericCopy element: " << tag << endln;
        return TCL_ERROR;
    }
    // create and fill in the ID array to hold the nodes
    ID nodes(numNodes);
    for (i=0; i<numNodes; i++)  {
        if (Tcl_GetInt(interp, argv[argi], &node) != TCL_OK)  {
            opserr << "WARNING invalid node\n";
            opserr << "genericCopy element: " << tag << endln;
            return TCL_ERROR;
        }
        nodes(i) = node;
        argi++; 
    }
    if (strcmp(argv[argi], "-src") != 0)  {
        opserr << "WARNING expect -src\n";
        opserr << "genericCopy element: " << tag << endln;
        return TCL_ERROR;
    }
    argi++;
    if (Tcl_GetInt(interp, argv[argi], &srcTag) != TCL_OK)  {
        opserr << "WARNING invalid srcTag\n";
        opserr << "genericCopy element: " << tag << endln;
        return TCL_ERROR;
    }
    
    // now create the GenericCopy
    theElement = new GenericCopy(tag, nodes, srcTag);
    
    if (theElement == 0)  {
        opserr << "WARNING ran out of memory creating element\n";
        opserr << "genericCopy element: " << tag << endln;
        return TCL_ERROR;
    }
    
    // then add the GenericCopy to the domain
    if (theTclDomain->addElement(theElement) == false)  {
        opserr << "WARNING could not add element to the domain\n";
        opserr << "genericCopy element: " << tag << endln;
        delete theElement;
        return TCL_ERROR;
    }
    
    // if get here we have successfully created the GenericCopy and added it to the domain
    return TCL_OK;
}
