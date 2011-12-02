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
                                                                        

// $Revision: 1.7 $
// $Date: 2003-04-29 18:43:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/TclJoint2dCommand.cpp,v $

// Written: Arash	Created: 04/01
// Revision: 
//				AAA		02/03
//
// Description: This file contains the implementation of the TclModelBuilder_addJoint2D()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Joint2D.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addJoint2D(ClientData clientData, Tcl_Interp *interp,  
				int argc, 
				TCL_Char **argv, 
				Domain *theTclDomain,
				TclModelBuilder *theTclBuilder)
{
	// ensure the destructor has not been called
	if (theTclBuilder == 0) {
		opserr << "WARNING builder has been destroyed\n";
		return TCL_ERROR;
	}
	
	if (theTclBuilder->getNDM() != 2 || theTclBuilder->getNDF() != 3) {
		opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with BCconnect element\n";
		return TCL_ERROR;
	}
	
	// check the number of arguments is correct
	int argStart = 2;
	
	if ( (argc-argStart) != 8 && (argc-argStart) != 10 && (argc-argStart) != 12 && (argc-argStart) != 18) {
		opserr << "WARNING incorrect number of arguments\n";
		printCommand(argc, argv);
		opserr << "Want:\n";
		opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? LrgDsp?\n";
		opserr << "or:\n";
		opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatC? LrgDsp? -damage DmgTag?\n";
		opserr << "or:\n";
		opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? MatK? MatL? MatC? LrgDsp?\n";
		opserr << "or:\n";
		opserr << "element Joint2D Tag? NodI? NodJ? NodK? NodL? NodC? MatI? MatJ? MatK? MatL? MatC? LrgDsp? -damage DmgI DmgJ DmgK DmgL DmgC\n";
		return TCL_ERROR;
	}
	
	// get the id and end nodes
	int Joint2DId, iNode, jNode, kNode, lNode;
	if (Tcl_GetInt(interp, argv[argStart], &Joint2DId) != TCL_OK) {
		opserr << "WARNING invalid Joint2D eleTag" << endln;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
		opserr << "WARNING invalid iNode\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
		opserr << "WARNING invalid jNode\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
		opserr << "WARNING invalid kNode\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}
	
	if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
		opserr << "WARNING invalid lNode\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}

	// Get the center node
	int CenterNodeTag;
	if (Tcl_GetInt(interp, argv[5+argStart], &CenterNodeTag) != TCL_OK) {
		opserr << "WARNING invalid tag for center node\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}
	
	// check domain for existence of internal node tag
	Node *CenterNode = theTclDomain->getNode(CenterNodeTag);
	if (CenterNode != 0) {
		opserr << "WARNING node tag specified for the center node already exists.\n";
		opserr << "Use a new node tag.\n";
		opserr << "Joint2D element: " << Joint2DId << endln;
		return TCL_ERROR;
	}
	
	UniaxialMaterial *MatI = NULL;
	UniaxialMaterial *MatJ = NULL;
	UniaxialMaterial *MatK = NULL;
	UniaxialMaterial *MatL = NULL;
	UniaxialMaterial *PanelMaterial = NULL;
	Joint2D *theJoint2D;
	int LargeDisp;


	// Decide to use which constructor, based on the number of arguments
	if ( (argc-argStart) == 8 || (argc-argStart) == 12 ) {

		// Using Joint2D constructor without damage 
		
		if ( (argc-argStart) == 8  )
		{
			int PanelMatId;
			if (Tcl_GetInt(interp, argv[6+argStart], &PanelMatId) != TCL_OK) {
				opserr << "WARNING invalid matID\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
			
			if (Tcl_GetInt(interp, argv[7+argStart], &LargeDisp) != TCL_OK) {
				// use 0 as default
				LargeDisp = 0;
			}

			PanelMaterial = theTclBuilder->getUniaxialMaterial(PanelMatId);
			
			if ( PanelMaterial == 0 ) {
				opserr << "WARNING material not found\n";
				opserr << "Material: " << PanelMatId;
				opserr << "\nJoint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
		}
		
		else			// if ( (argc-argStart) == 12  )
		{
			int MatIid;
			if (Tcl_GetInt(interp, argv[6+argStart], &MatIid) != TCL_OK) {
				opserr << "WARNING invalid material ID for spring I\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
					
			if ( MatIid != 0 ) {
				MatI = theTclBuilder->getUniaxialMaterial(MatIid);
				
				if ( MatI == NULL )
				{
					opserr << "WARNING material not found\n";
					opserr << "Material: " << MatIid;
					opserr << "\nJoint2D element: " << Joint2DId << endln;
					return TCL_ERROR;
				}
			} else MatI = NULL;
			
			int MatJid;
			if (Tcl_GetInt(interp, argv[7+argStart], &MatJid) != TCL_OK) {
				opserr << "WARNING invalid material ID for spring J\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
			
			if ( MatJid != 0 ) {
				MatJ = theTclBuilder->getUniaxialMaterial(MatJid);
				
				if ( MatJ == NULL )
				{
					opserr << "WARNING material not found\n";
					opserr << "Material: " << MatJid;
					opserr << "\nJoint2D element: " << Joint2DId << endln;
					return TCL_ERROR;
				}
			} else MatJ = NULL;

			
			int MatKid;
			if (Tcl_GetInt(interp, argv[8+argStart], &MatKid) != TCL_OK) {
				opserr << "WARNING invalid material ID for spring K\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				
				return TCL_ERROR;
			}
			if ( MatKid != 0 ) {
				MatK = theTclBuilder->getUniaxialMaterial(MatKid);
				
				if ( MatK == NULL )
				{
					opserr << "WARNING material not found\n";
					opserr << "Material: " << MatKid;
					opserr << "\nJoint2D element: " << Joint2DId << endln;
					return TCL_ERROR;
				}
			} else MatK = NULL;
			
			int MatLid;
			if (Tcl_GetInt(interp, argv[9+argStart], &MatLid) != TCL_OK) {
				opserr << "WARNING invalid material ID for spring L\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
			if ( MatLid != 0 ) {
				MatL = theTclBuilder->getUniaxialMaterial(MatLid);
				
				if ( MatL == NULL )
				{
					opserr << "WARNING material not found\n";
					opserr << "Material: " << MatLid;
					opserr << "\nJoint2D element: " << Joint2DId << endln;
					return TCL_ERROR;
				}
			} else MatL = NULL;
			
			int PanelMatId;
			if (Tcl_GetInt(interp, argv[10+argStart], &PanelMatId) != TCL_OK) {
				opserr << "WARNING invalid matID\n";
				opserr << "Joint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}
			PanelMaterial = theTclBuilder->getUniaxialMaterial(PanelMatId);
			
			if ( PanelMaterial == 0 ) {
				opserr << "WARNING material not found\n";
				opserr << "Material: " << PanelMatId;
				opserr << "\nJoint2D element: " << Joint2DId << endln;
				return TCL_ERROR;
			}

			if (Tcl_GetInt(interp, argv[11+argStart], &LargeDisp) != TCL_OK) {
				// use 0 as default
				LargeDisp = 0;
			}
		}

		theJoint2D = new Joint2D( Joint2DId,
			iNode,jNode,kNode,lNode,CenterNodeTag,
			*MatI,*MatJ,*MatK,*MatL,*PanelMaterial, 
			theTclDomain, 
			LargeDisp);
		
		if (theJoint2D == 0) {
			opserr << "WARNING ran out of memory creating element\n";
			opserr << "Joint2D element: " << Joint2DId << endln;
			return TCL_ERROR;
		}
		
		if (theTclDomain->addElement(theJoint2D) == false) {
			opserr << "WARNING could not add element to the domain\n";
			opserr << "Joint2D element: " << Joint2DId << endln;
			delete theJoint2D;
			return TCL_ERROR;
		}
		
		// if get here we have sucessfully created the element and added it to the domain
		return TCL_OK;
	}
	
	else if ( (argc-argStart) == 10 || (argc-argStart) == 18 )
	{ 
		// Using Joint2D constructor with damage 
		// not implemented in this version
		return TCL_ERROR;
	}
	else 
	{
		return TCL_ERROR;
	}
}

