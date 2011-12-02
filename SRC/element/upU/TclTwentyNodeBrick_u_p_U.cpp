///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:             
// CLASS:            
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class for coupled system
// RETURN:
// VERSION:
// LANGUAGE:          C++.ver >= 3.0
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Xiaoyan Wu
// PROGRAMMER:        Boris Jeremic, Xiaoyan Wu
// DATE:              Sept. 2001
//
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
//
///////////////////////////////////////////////////////////////////////////////

//
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <TwentyNodeBrick.h>
#include <TwentyNodeBrick_u_p_U.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addTwentyNodeBrick_u_p_U(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  TCL_Char **argv, Domain *theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element Brick20N_u_p_U - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  // should have 36 arguments.  
  //  element  Brick20N_u_p_U eleTag?  node1? node2? .. node20?  matTag?   bforce1?  bforce2?  bforce3?\n
  //  argv[0]	argv[1]	      argv[2] argv[3]			 argv[23]  argv[24]  argv[25]  argv[26]
  //  porosity?  alpha?  solidDensity? fluidDensity? x_permeability? y_permeability? z_permeability?
  //  argv[27]  argv[28]   argv[29]	 argv[30]       argv[31]          argv[32]      argv[33]
  //  solid_bulk_modulus? fluid_bulk_modulus?   pressure?
  //       argv[34]            argv[35]     	 argv[36] 
  // Xiaoyan added this comments. 01/07/2002
  if ((argc-eleArgStart) < 35) {
    opserr << "command: element Brick20N_u_p_U - insufficient args - want " << 
      "element Brick20N_u_p_U eleTag? node1? node2? .. node20?  matTag? bforce1? bforce2? bforce3?\n" <<
      "porosity? alpha?  solidDensity? fluidDensity? \n" <<
      "permeability_in_x_dir? permeability_in_y_dir? permeability_in_z_dir?" << 
      "solid_bulk_modulus, fluid_bulk_modulus? pressure?\n";
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[20];
  double bodyforces[3], porosity, alpha, solidDensity, fluidDensity;
  double perm_x, perm_y, perm_z,kks, kkf;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    opserr << "command: element Brick20N_u_p_U - invalid integer tag " << argv[1+eleArgStart] << "\n";
			      

      return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<20; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid integer tag " << argv[2+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[22+eleArgStart], &matID) != TCL_OK) {	     // argv[23]  Xiaoyan 01/07/2002
    opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid matID tag " << argv[22+eleArgStart] << "\n";      
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - no NDMaterial with tag " << argv[22+eleArgStart] << " exists\n";
    return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];

  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[23+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid bodyforces tag " << argv[23+i+eleArgStart] << "\n";   
	return TCL_ERROR;
      }
  }

  // now get the porosity
  if (Tcl_GetDouble(interp, argv[26+eleArgStart], &porosity) != TCL_OK) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid porosity " << argv[26+eleArgStart] << endln;
    return TCL_ERROR;
  } 
 
  // now get the alpha for solid alpha=1.0
   if (Tcl_GetDouble(interp, argv[27+eleArgStart], &alpha) != TCL_OK) {
     opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid alpha " << argv[27+eleArgStart] << endln;
     return TCL_ERROR;
  }  
 
   // now get the solidDensity
   if (Tcl_GetDouble(interp, argv[28+eleArgStart], &solidDensity) != TCL_OK) {
     opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid solidDensity " << argv[28+eleArgStart] << endln;
     return TCL_ERROR;
  }  
 
   // now get the fluidDensity
   if (Tcl_GetDouble(interp, argv[29+eleArgStart], &fluidDensity) != TCL_OK) {
     opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid fluidDensity " << argv[29+eleArgStart] << endln;
     return TCL_ERROR;
  }  
  // permeability in x direction
  if (Tcl_GetDouble(interp, argv[30+eleArgStart], &perm_x) != TCL_OK) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid permeability in x direction " << argv[30+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in y direction
  if (Tcl_GetDouble(interp, argv[31+eleArgStart], &perm_y) != TCL_OK) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid permeability in y direction " << argv[31+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in z direction
  if (Tcl_GetDouble(interp, argv[32+eleArgStart], &perm_z) != TCL_OK) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - invalid permeability in z direction " << argv[32+eleArgStart] << endln;
    return TCL_ERROR;
  }  
  // now get the bulk modulus of solid
  if (Tcl_GetDouble(interp, argv[33+eleArgStart], &kks) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << "  - invalid bulk modulus of solid " << argv[33+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // now get the bulk modulus of fluid
  if (Tcl_GetDouble(interp, argv[34+eleArgStart], &kkf) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << "  - invalid bulk modulus of fluid " << argv[34+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  // now create the EightNodeBrick and add it to the Domain
  TwentyNodeBrick_u_p_U *theEle = new TwentyNodeBrick_u_p_U(eleID, nodes[ 0], nodes [1], nodes[ 2], nodes[ 3], nodes[ 4],
							    nodes[ 5], nodes [6], nodes[ 7], nodes[ 8], nodes[ 9],
							    nodes[10], nodes[11], nodes[12], nodes[13], nodes[14],
							    nodes[15], nodes[16], nodes[17], nodes[18], nodes[19],
							    theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], 
							    porosity, alpha, solidDensity, fluidDensity,
							    perm_x, perm_y, perm_z, kks, kkf, 0.0);
  
  if (theEle == 0) {
    opserr << "command: element Brick20N_u_p_U " << eleID << "  - out of memory\n";      
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element Brick20N_u_p_U " << eleID << " - could not add ele to domain\n"; 
    delete theEle;
    return TCL_ERROR;
  }
  
  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



