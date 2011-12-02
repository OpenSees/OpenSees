/* ****************************************************************** **
**    Opensees - Open System for Earthquake Engineering Simulation    **
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
                                                                        
// $Revision: 1.2 $
// $Date: 2003-06-11 18:33:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/TclBeamColumnJointCommand.cpp,v $

// Written: NM (nmitra@u.washington.edu)
// Created: April 2002
//
// Description: This file contains the implementation of the commands used 
// to add beam column joint to a model.


#include <stdlib.h>
#include <string.h>

#include <BeamColumnJoint2d.h>
#include <BeamColumnJoint3d.h>
#include <Information.h>
#include <ElementResponse.h>

#include <Domain.h>
#include <Node.h>
#include <Channel.h>  // not yet implemented 
#include <FEM_ObjectBroker.h>  
#include <Renderer.h>
#include <UniaxialMaterial.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);
static Domain *theTclModelBuilderDomain = 0;       
static TclModelBuilder *theTclModelBuilder =0;     

int TclModelBuilder_addBeamColumnJoint(ClientData clientData, Tcl_Interp *interp, int argc,
									TCL_Char **argv, Domain* theTclDomain,
									TclModelBuilder *theTclBuilder, int eleArgStart)
{
	// ensure destructor not called
	if (theTclBuilder == 0) {
		opserr << "Warning builder has been destroyed \n";
		return TCL_ERROR;
	}

    theTclModelBuilderDomain = theTclDomain;
    theTclModelBuilder = theTclBuilder;
    
    int NDM, NDF;
     
    NDM = theTclModelBuilder->getNDM();   // dimension of the structure (1d, 2d, or 3d)
    NDF = theTclModelBuilder->getNDF();   // number of degrees of freedom per node

	
	if ((NDM == 2 && NDF == 3) || (NDM == 3 && NDF == 6)) {
	
	// check no of arguments
	if ((argc-eleArgStart) < 18)
	{
		opserr << "WARNING insufficient arguments\n";
		printCommand(argc,argv);
		opserr << "Want: element beamColumnJoint eleTag? node1? node2? node3? node4? matTag1? matTag2? matTag3?\n";
		opserr<< "matTag4? matTag5? matTag6? matTag7? matTag8? matTag9? matTag10? matTag11? matTag12? matTag13?\n";

		return TCL_ERROR;
	}
	
	int id, nd1, nd2, nd3, nd4, matId1, matId2, matId3, matId4, matId5, matId6, matId7, matId8, matId9, matId10;
	int matId11, matId12, matId13;


	if (Tcl_GetInt(interp, argv[1+eleArgStart], &id) != TCL_OK)
	{
		opserr << "WARNING invalid beamColumnJoint eleTag" << endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[2+eleArgStart], &nd1) != TCL_OK)
	{
		opserr << "WARNING invalid Node 1\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[3+eleArgStart], &nd2) != TCL_OK)
	{
		opserr << "WARNING invalid Node 2\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[4+eleArgStart], &nd3) != TCL_OK)
	{
		opserr << "WARNING invalid Node 3\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[5+eleArgStart], &nd4) != TCL_OK)
	{
		opserr << "WARNING invalid Node 4\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[6+eleArgStart], &matId1) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 1\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[7+eleArgStart], &matId2) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 2\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[8+eleArgStart], &matId3) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 3\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[9+eleArgStart], &matId4) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 4\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[10+eleArgStart], &matId5) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 5\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[11+eleArgStart], &matId6) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 6\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[12+eleArgStart], &matId7) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 7\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[13+eleArgStart], &matId8) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 8\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[14+eleArgStart], &matId9) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 9\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[15+eleArgStart], &matId10) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 10\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[16+eleArgStart], &matId11) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 11\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[17+eleArgStart], &matId12) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 12\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
	if (Tcl_GetInt(interp, argv[18+eleArgStart], &matId13) != TCL_OK)
	{
		opserr << "WARNING invalid Material Tag 13\n";
		opserr << "beamColumnJoint Element: " << id <<endln;
		return TCL_ERROR;
	}
 

	UniaxialMaterial *theMaterial1 = theTclBuilder->getUniaxialMaterial(matId1);
      
    if (theMaterial1 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId1;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial2 = theTclBuilder->getUniaxialMaterial(matId2);
      
    if (theMaterial2 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId2;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial3 = theTclBuilder->getUniaxialMaterial(matId3);
      
    if (theMaterial3 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId3;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial4 = theTclBuilder->getUniaxialMaterial(matId4);
      
    if (theMaterial4 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId4;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }
	UniaxialMaterial *theMaterial5 = theTclBuilder->getUniaxialMaterial(matId5);
      
    if (theMaterial5 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId5;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial6 = theTclBuilder->getUniaxialMaterial(matId6);
      
    if (theMaterial6 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId6;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }
	UniaxialMaterial *theMaterial7 = theTclBuilder->getUniaxialMaterial(matId7);
      
    if (theMaterial7 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId7;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial8 = theTclBuilder->getUniaxialMaterial(matId8);
      
    if (theMaterial8 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId8;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }
	UniaxialMaterial *theMaterial9 = theTclBuilder->getUniaxialMaterial(matId9);
      
    if (theMaterial9 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId9;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial10 = theTclBuilder->getUniaxialMaterial(matId10);
      
    if (theMaterial10 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId10;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial11 = theTclBuilder->getUniaxialMaterial(matId11);
      
    if (theMaterial11 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId11;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }
	UniaxialMaterial *theMaterial12 = theTclBuilder->getUniaxialMaterial(matId12);
      
    if (theMaterial12 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId12;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	UniaxialMaterial *theMaterial13 = theTclBuilder->getUniaxialMaterial(matId13);
      
    if (theMaterial13 == 0) {
	 opserr << "WARNING material not found\n";
	 opserr << "Material: " << matId13;
	 opserr << "\nbeamColumnJoint element: " << id << endln;
	 return TCL_ERROR;
    }

	// create Beam Column Joint Element and add it to Domain
	Element *theBeamColumnJoint = 0;

	if (NDM == 2)
	{
		theBeamColumnJoint = new BeamColumnJoint2d(id,nd1,nd2,nd3,nd4,*theMaterial1,*theMaterial2,
			*theMaterial3,*theMaterial4,*theMaterial5,*theMaterial6,*theMaterial7,*theMaterial8,
			*theMaterial9,*theMaterial10,*theMaterial11,*theMaterial12,*theMaterial13);
	}
	else if (NDM == 3)
	{
		theBeamColumnJoint = new BeamColumnJoint3d(id,nd1,nd2,nd3,nd4,*theMaterial1,*theMaterial2,
			*theMaterial3,*theMaterial4,*theMaterial5,*theMaterial6,*theMaterial7,*theMaterial8,
			*theMaterial9,*theMaterial10,*theMaterial11,*theMaterial12,*theMaterial13);
	}

	if (theBeamColumnJoint == 0)
	{
		opserr << "WARNING ran out of memory creating elements\n";
		opserr << "beamColumnJoint element: " << id <<endln;
		return TCL_ERROR;
	}

	if (theTclDomain->addElement(theBeamColumnJoint) == false)
	{
		opserr << "WARNING could not add element to the domain\n";
		opserr << "beamColumnJoint element: " << id << endln;
		delete theBeamColumnJoint;
		return TCL_ERROR;
	}

	}
	else {
		opserr << "WARNING NDM = " << NDM << " and NDF = " << NDF << "is imcompatible with available frame elements";
		return TCL_ERROR;
	}      

	//the element succesfully created and added to the domain
	return TCL_OK;
}
