///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              EightNodeBrick_u_p.cpp
// CLASS:             EightNodeBrick_u_p
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, Boris Jeremic
// DATE:              Aug. 2006
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////

//  "Coupled system": Z.Cheng & B.Jeremic @ UCDavis
//                    u-- Solid displacement
//                    p-- Pore pressure

//
#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <ErrorHandler.h>
//#include <TwentyNodeBrick.h>

#include <EightNode_LDBrick_u_p.h>  //large deformation
#include <EightNode_Brick_u_p.h>    //small deformation

#include <TclModelBuilder.h>

extern void printCommand(int argc, TCL_Char **argv);


int
TclModelBuilder_addEightNode_LDBrick_u_p(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  TCL_Char **argv, Domain *theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element EightNode_LDBrick_u_p - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  // should have 35 arguments.  
  if ((argc-eleArgStart) < 21) {
    opserr << "command: element EightNode_LDBrick_u_p - insufficient args ";
      return TCL_ERROR;
  }    
  // get the id and end nodes
  int eleID = 0;
  int matID = 0;
  int nodes[8] = {0,0,0,0,0,0,0,0};
  double bodyforces[3] = {0.0, 0.0, 0.0};
  double fluidfraction = 0.0;
  double solidDensity = 0.0;
  double fluidDensity = 0.0;
  double perm_x = 0.0; 
  double perm_y = 0.0;
  double perm_z = 0.0;
  double kkf = 0.0;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p - invalid integer tag " << argv[1+eleArgStart] << "\n";			     
    return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<8; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	opserr << "command: EightNode_LDBrick_u_p " << eleID << " - invalid integer tag " << argv[2+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid matID tag " << argv[10+eleArgStart] << "\n";      
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - no NDMaterial with tag " << argv[10+eleArgStart] << " exists\n";
    return TCL_ERROR;      
  }

  // read the 3 bodyforces 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid body forces tag " << argv[11+i+eleArgStart] << "\n";   
	return TCL_ERROR;
      }
  }

  // now get the void fraction
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &fluidfraction) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid void fraction " << argv[14+eleArgStart] << endln;
    return TCL_ERROR;
  }  
 
  // now get the solid density
  if (Tcl_GetDouble(interp, argv[15+eleArgStart], &solidDensity) != TCL_OK) {
     opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid solidDensity " << argv[15+eleArgStart] << endln;
     return TCL_ERROR;
  }  
 
   // now get the fluid density
   if (Tcl_GetDouble(interp, argv[16+eleArgStart], &fluidDensity) != TCL_OK) {
     opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid fluidDensity " << argv[16+eleArgStart] << endln;
     return TCL_ERROR;
  }  

  // permeability in x direction
  if (Tcl_GetDouble(interp, argv[17+eleArgStart], &perm_x) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid permeability in x direction " << argv[17+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in y direction
  if (Tcl_GetDouble(interp, argv[18+eleArgStart], &perm_y) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid permeability in y direction " << argv[18+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in z direction
  if (Tcl_GetDouble(interp, argv[19+eleArgStart], &perm_z) != TCL_OK) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - invalid permeability in z direction " << argv[19+eleArgStart] << endln;
    return TCL_ERROR;
  }  

 
  // now get the bulk modulus of fluid
  if (Tcl_GetDouble(interp, argv[20+eleArgStart], &kkf) != TCL_OK) { 
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << "  - invalid bulk modulus of fluid " << argv[20+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  

  // now create the Twenty_EightNode_LDBrick_u_p and add it to the Domain
  EightNode_LDBrick_u_p *theEle = new EightNode_LDBrick_u_p(eleID, 
                                nodes[0], nodes[1], nodes[2], nodes[3],
				nodes[4], nodes[5], nodes[6], nodes[7],
				theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], 
				fluidfraction, solidDensity, fluidDensity,
				perm_x, perm_y, perm_z, kkf); 
  if (theEle == 0) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << "  - out of memory\n";      
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element EightNode_LDBrick_u_p " << eleID << " - could not add ele to domain\n"; 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}



int
TclModelBuilder_addEightNode_Brick_u_p(ClientData clientData, Tcl_Interp *interp,  int argc, 
				  TCL_Char **argv, Domain *theTclDomain,
				  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
      opserr << "command: element EightNode_Brick_u_p - no modelbuilder\n";
      return TCL_ERROR;
  }

  // check the number of arguments is correct
  // should have 35 arguments.  
  if ((argc-eleArgStart) < 23) {
    opserr << "command: element EightNode_Brick_u_p - insufficient args ";
      return TCL_ERROR;
  }    
  // get the id and end nodes
  int eleID = 0;
  int matID = 0;
  int nodes[8] = {0,0,0,0,0,0,0,0};
  double bodyforces[3] = {0.0, 0.0, 0.0};
  double fluidfraction = 0.0;
  double solidDensity = 0.0;
  double fluidDensity = 0.0;
  double perm_x = 0.0; 
  double perm_y = 0.0;
  double perm_z = 0.0;
  double kks = 0.0;
  double kkf = 0.0;
  double alphaf = 1.0;
  
  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p - invalid integer tag " << argv[1+eleArgStart] << "\n";			     
    return TCL_ERROR;
  }
  
  // read the 20 node tags
  int i;
  for (i=0; i<8; i++) {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK) {
	opserr << "command: EightNode_Brick_u_p " << eleID << " - invalid integer tag " << argv[2+i+eleArgStart] << endln;
	return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[10+eleArgStart], &matID) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid matID tag " << argv[10+eleArgStart] << "\n";      
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);
  
  if (theMaterial == 0) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - no NDMaterial with tag " << argv[10+eleArgStart] << " exists\n";
    return TCL_ERROR;      
  }

  // read the 3 bodyforces 
  for (i=0; i<3; i++) {
      if (Tcl_GetDouble(interp, argv[11+i+eleArgStart], &bodyforces[i]) != TCL_OK) {
	opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid body forces tag " << argv[11+i+eleArgStart] << "\n";   
	return TCL_ERROR;
      }
  }

  // now get the void fraction
  if (Tcl_GetDouble(interp, argv[14+eleArgStart], &fluidfraction) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid void fraction " << argv[14+eleArgStart] << endln;
    return TCL_ERROR;
  }
    
  // now get the alpha
  if (Tcl_GetDouble(interp, argv[15+eleArgStart], &alphaf) != TCL_OK) { 
    opserr << "command: element EightNode_Brick_u_p " << eleID << "  - invalid alpha " << argv[15+eleArgStart] << endln;      
    return TCL_ERROR;
  }
   
  // now get the solid density
  if (Tcl_GetDouble(interp, argv[16+eleArgStart], &solidDensity) != TCL_OK) {
     opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid solidDensity " << argv[16+eleArgStart] << endln;
     return TCL_ERROR;
  }  
 
   // now get the fluid density
   if (Tcl_GetDouble(interp, argv[17+eleArgStart], &fluidDensity) != TCL_OK) {
     opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid fluidDensity " << argv[17+eleArgStart] << endln;
     return TCL_ERROR;
  }  

  // permeability in x direction
  if (Tcl_GetDouble(interp, argv[18+eleArgStart], &perm_x) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid permeability in x direction " << argv[18+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in y direction
  if (Tcl_GetDouble(interp, argv[19+eleArgStart], &perm_y) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid permeability in y direction " << argv[19+eleArgStart] << endln;      
    return TCL_ERROR;
  }  
  // permeability in z direction
  if (Tcl_GetDouble(interp, argv[20+eleArgStart], &perm_z) != TCL_OK) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - invalid permeability in z direction " << argv[20+eleArgStart] << endln;
    return TCL_ERROR;
  }  

  // now get the bulk modulus of solid
  if (Tcl_GetDouble(interp, argv[21+eleArgStart], &kks) != TCL_OK) { 
    opserr << "command: element EightNode_Brick_u_p " << eleID << "  - invalid bulk modulus of solid " << argv[21+eleArgStart] << endln;      
    return TCL_ERROR;
  } 
 
  // now get the bulk modulus of fluid
  if (Tcl_GetDouble(interp, argv[22+eleArgStart], &kkf) != TCL_OK) { 
    opserr << "command: element EightNode_Brick_u_p " << eleID << "  - invalid bulk modulus of fluid " << argv[22+eleArgStart] << endln;      
    return TCL_ERROR;
  }     

  // now create the Twenty_EightNode_LDBrick_u_p and add it to the Domain
  EightNode_Brick_u_p *theEle = new EightNode_Brick_u_p(eleID, 
                                nodes[0], nodes[1], nodes[2], nodes[3],
				nodes[4], nodes[5], nodes[6], nodes[7],
				theMaterial, bodyforces[0], bodyforces[1], bodyforces[2], 
				fluidfraction, alphaf, solidDensity, fluidDensity,
				perm_x, perm_y, perm_z, kks, kkf); 
  if (theEle == 0) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << "  - out of memory\n";      
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false) {
    opserr << "command: element EightNode_Brick_u_p " << eleID << " - could not add ele to domain\n"; 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

