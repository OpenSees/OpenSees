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

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <EightNodeBrick_u_p_U.h>
#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addEightNodeBrick_u_p_U(ClientData clientData, Tcl_Interp *interp,  int argc, 
					TCL_Char **argv, Domain *theTclDomain,
					TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element Brick8_u_p_U - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
 
  //  element  Brick20N_u_p_U eleTag?  node1? node2? .. node8?  matTag?   bforce1?  bforce2?  bforce3?\n
  //  argv[0]	argv[1]	      argv[2] argv[3]		       argv[11]  argv[12]  argv[13]  argv[14]
  //  porosity?  alpha?  solidDensity? fluidDensity?  x_permeability? y_permeability? z_permeability? 
  //  argv[15]  argv[16]   argv[17]	 argv[18]        argv[19]          argv[20]      argv[21]     
  //  solid_bulk_modulus? fluid_bulk_modulus?   pressure?
  //       argv[22]            argv[23]     	 argv[24] 
  // Xiaoyan added this comments. 01/07/2002

  if ((argc-eleArgStart) < 23) {
      opserr << "command: element Brick8_u_p_U - insufficient args - want %s",
          "element Brick8_u_p_U eleTag? node1? node2? ... node8? matTag? bforce1? bforce2? bforce3? \n"
	  " porosity? alpha? solidDensity? fluidDensity? \n"
 	  "permeability_in_x_dir? permeability_in_y_dir? permeability_in_z_dir?"
	  "solid_bulk_modulus? fluid_bulk_modulus? pressure?\n";
      return TCL_ERROR;
  }    

  // get the id and end nodes 
  int eleID, matID;
  int nodes[8];					  
  double bodyforces[3], porosity, alpha, solidDensity, fluidDensity;
  double perm_x, perm_y, perm_z, kks, kkf;
  //char *type;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
      opserr << "command: element Brick8_u_p_U - invalid integer tag " <<      
	argv[1+eleArgStart] << endln;

      return TCL_ERROR;
  }
  
  // read the 8 node tags
  int i;
  for (i=0; i<8; i++) {
    if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
      opserr << "command: element Brick8_u_p_U " << eleID << " - invalid integer tag " <<      
	argv[2+i+eleArgStart] << endln;
      return TCL_ERROR;
    }
  }
  
  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid matID tag " <<      
      argv[10+eleArgStart] << endln;      
    return TCL_ERROR;
  }
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element Brick8_u_p_U " << eleID << " - no NDMaterial with tag %s exists" <<
      argv[10+eleArgStart] << endln;      
    return TCL_ERROR;      
  }
  
  //type = argv[11+eleArgStart];
  
  // read the 3 bodyforce accel's 
  for (i=0; i<3; i++) {
    if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
      opserr << "command: element Brick8_u_p_U " << eleID << " - invalid bodyforces tag " <<      
	argv[11+i+eleArgStart] << endln;
      return TCL_ERROR;
    }
  }
  
  // now get the porosity
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &porosity) != TCL_OK) {	       // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid porosity " <<      
      argv[14+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  // now get the alpha
  if (Tcl_GetDouble(interp, argv[15+eleArgStart], &alpha) != TCL_OK) {	      // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid alpha " <<      
      argv[15+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  // now get the soldDensity
  
  if (Tcl_GetDouble(interp, argv[16+eleArgStart], &solidDensity) != TCL_OK) {	      // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid solidDensity " <<      
      argv[16+eleArgStart] << endln;      
    return TCL_ERROR;
  }    
  
  // now get the fludiDensity
  if (Tcl_GetDouble(interp, argv[17+eleArgStart], &fluidDensity) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid fluidDensity " <<      
      argv[17+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  
  // now get the permeability in x direction
  if (Tcl_GetDouble(interp, argv[18+eleArgStart], &perm_x) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid permeability in x direction " <<      
      argv[18+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // now get the permeability in y direction
  if (Tcl_GetDouble(interp, argv[19+eleArgStart], &perm_y) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid permeability in y direction " <<      
      argv[19+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // now get the permeability in z direction
  if (Tcl_GetDouble(interp, argv[20+eleArgStart], &perm_z) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid permeability in z direction " <<      
      argv[20+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  // now get the bulk modulus of solid
  if (Tcl_GetDouble(interp, argv[21+eleArgStart], &kks) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid bulk modulus of solid " <<      
      argv[21+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  
  // now get the bulk modulus of fluid
  if (Tcl_GetDouble(interp, argv[22+eleArgStart], &kkf) != TCL_OK) {        // wxy added 01/07/2002
    opserr << "command: element Brick8_u_p_U " << eleID << " - invalid bulk modulus of fluid " <<      
      argv[21+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // now create the EightNodeBrick and add it to the Domain
  EightNodeBrick_u_p_U *theEle = new EightNodeBrick_u_p_U(eleID,nodes[0], nodes[1], nodes[2], nodes[3], nodes[4],
							  nodes[5],nodes[6], nodes[7], theMaterial, 
							  bodyforces[0], bodyforces[1], bodyforces[2], 
							  porosity, alpha, solidDensity, fluidDensity, 
							  perm_x, perm_y, perm_z, kks, kkf,0.0);
  
  if (theEle == 0) {
    opserr << "command: element Brick8_u_p_U " << eleID << " - out of memory\n";
    return TCL_ERROR;
  }
  
  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element Brick8_u_p_U " << eleID << " - could not add ele to domain\n";
    delete theEle;
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



