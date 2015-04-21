
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// TclAcousticCommand.cpp
//------------------------------------------


#include <stdlib.h>
#include <string.h>
#include <Domain.h>
#include <AC3D8HexWithSensitivity.h>
#include <TclModelBuilder.h>
  
#include <AC3D8Hex.h>
#include <ASI3D8Quad.h>  
#include <AV3D4Quad.h>
#include <ASI3D8QuadWithSensitivity.h>
#include <AV3D4QuadWithSensitivity.h>

extern void printCommand(int argc, TCL_Char **argv);




int 
TclModelBuilder_addAC3D8Hex(ClientData clientData, Tcl_Interp *interp, 
  int argc, TCL_Char **argv, Domain *theTclDomain, 
  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    printf("WARNING builder has been destroyed\n");    
    return TCL_ERROR;
  }
  
  int eleID, numNodes, matTag;
  int nodes[8];
  int i;
  
  if(argc - eleArgStart < 11) {
    printf("WARNING insufficient arguments!\n");
    printCommand(argc, argv);
    printf("WANT: element AC3D8 eleTag? Node1? Node2? ... NodeX matTag?\n");
    return TCL_ERROR;
  }
  
  
  if(Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    printf("WARNING invalid AC3D8 eleTag\n");
    return TCL_ERROR;
  }

  // printf("Element id is %d!\n", eleID);
  
  numNodes = 8;
  for(i = 0; i < numNodes; i++) {
    if(Tcl_GetInt(interp, argv[i+2+eleArgStart], &nodes[i]) != TCL_OK) {
      printf("Command: element AC3D8 %d - invalid integer tag %s\n", 
        eleID, argv[i+2+eleArgStart]);
      return TCL_ERROR;
    }
    // printf("Nodes[%d] is %d!\n", i+1, nodes[i]);
  }
  
  if(Tcl_GetInt(interp, argv[11], &matTag) != TCL_OK) {
    printf("WARNING invalid AC3D8 matTag\n");
    return TCL_ERROR;
  }

  // printf("Material tag is %d!\n", matTag);
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matTag);
  if (theMaterial == 0) {
    printf("command: element AC3D8 %d - no NDMaterial with tag %s exists!\n", 
      eleID, argv[11+eleArgStart]);
    return TCL_ERROR;      
  }
  
  // printf("Get input data for AC3D8 ok...\n");
  //AC3D8Hex *theEle = new AC3D8Hex(eleID, nodes[0], nodes[1], nodes[2],   
  //  nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], theMaterial);

  Element *theEle = 0;
  
  if (strcmp(argv[1],"AC3D8") == 0) {
  theEle = new AC3D8Hex(eleID, nodes[0], nodes[1], nodes[2],     
    nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], theMaterial);
  }
  else if (strcmp(argv[1],"AC3D8WithSensitivity") == 0){
    theEle = new AC3D8HexWithSensitivity(eleID, nodes[0], nodes[1], nodes[2],      //
    nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], theMaterial);
  }


  else {
    opserr << "WARNING element " << argv[1] << " type not recognized\n";
    return TCL_ERROR;
  }

  //  AC3D8Hex *theEle = new AC3D8Hex(eleID, nodes[0], nodes[1], nodes[2],   
   // nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], theMaterial);
  
  if(theTclDomain->addElement(theEle) == false) {
    printf("command: element AC3D8 %d - could not add element to domain!\n", 
      eleID); 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


int 
TclModelBuilder_addASI3D8Quad(ClientData clientData, Tcl_Interp *interp, 
  int argc, TCL_Char **argv, Domain *theTclDomain, 
  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    printf("WARNING builder has been destroyed\n");    
    return TCL_ERROR;
  }
  
  int eleID, numNodes;
  int nodes[8];
  int i;
  
  if(argc - eleArgStart < 10) {
    printf("WARNING insufficient arguments!\n");
    printCommand(argc, argv);
    printf("WANT: element ASI3D8 eleTag? Node1? Node2? ... NodeX?\n");
    return TCL_ERROR;
  }
  
  if(Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    printf("WARNING invalid ASI3D8 eleTag\n");
    return TCL_ERROR;
  }
  
  numNodes = 8;
  for(i = 0; i < numNodes; i++) {
    if(Tcl_GetInt(interp, argv[i+2+eleArgStart], &nodes[i]) != TCL_OK) {
      printf("Command: element ASI3D8 %d - invalid integer tag %s\n", 
        eleID, argv[i+2+eleArgStart]);
      return TCL_ERROR;
    }
  }
  
  // if(Tcl_GetInt(interp, argv[11+eleArgStart], &matTag) != TCL_OK) {
  //   printf("WARNING invalid AC3D8 matTag\n");
  //   return TCL_ERROR;
  // }

  Element *theEle = 0;
  if (strcmp(argv[1],"ASI3D8") == 0){
    theEle = new ASI3D8Quad(eleID, nodes[0], nodes[1], nodes[2], 
    nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
  }
      else if (strcmp(argv[1],"ASI3D8WithSensitivity") == 0){
    theEle = new ASI3D8QuadWithSensitivity(eleID, nodes[0], nodes[1], nodes[2],      //
    nodes[3], nodes[4], nodes[5], nodes[6], nodes[7]);
  }

  else {
    opserr << "WARNING element " << argv[1] << " type not recognized\n";
    return TCL_ERROR;
  }

  
  if(theTclDomain->addElement(theEle) == false) {
    printf("command: element ASI3D8 %d - could not add element to domain!\n", 
      eleID); 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


int 
TclModelBuilder_addAV3D4Quad(ClientData clientData, Tcl_Interp *interp, 
  int argc, TCL_Char **argv, Domain *theTclDomain, 
  TclModelBuilder *theTclBuilder, int eleArgStart)
{
  // ensure the destructor has not been called - 
  if (theTclBuilder == 0) {
    printf("WARNING builder has been destroyed\n");    
    return TCL_ERROR;
  }
  
  int eleID, numNodes, matTag;
  int nodes[4];
  int i;
  
  if(argc - eleArgStart < 7) {
    printf("WARNING insufficient arguments!\n");
    printCommand(argc, argv);
    printf("WANT: element AV3D4 eleTag? Node1? Node2? ... NodeX matTag?\n");
    return TCL_ERROR;
  }
  
  
  if(Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    printf("WARNING invalid AV3D4 eleTag\n");
    return TCL_ERROR;
  }

  // printf("Element id is %d!\n", eleID);
  
  numNodes = 4;
  for(i = 0; i < numNodes; i++) {
    if(Tcl_GetInt(interp, argv[i+2+eleArgStart], &nodes[i]) != TCL_OK) {
      printf("Command: element AV3D4 %d - invalid integer tag %s\n", 
        eleID, argv[i+2+eleArgStart]);
      return TCL_ERROR;
    }
    // printf("Nodes[%d] is %d!\n", i+1, nodes[i]);
  }
  
  if(Tcl_GetInt(interp, argv[7], &matTag) != TCL_OK) {
    printf("WARNING invalid AV3D4 matTag\n");
    return TCL_ERROR;
  }

  // printf("Material tag is %d!\n", matTag);
  
  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matTag);
  if (theMaterial == 0) {
    printf("command: element AV3D4 %d - no NDMaterial with tag %s exists!\n", 
      eleID, argv[7+eleArgStart]);
    return TCL_ERROR;
  }

  // printf("Get input data for AV3D4 ok...\n");

  Element *theEle = 0;
  
  if (strcmp(argv[1],"AV3D4") == 0) {

  theEle = new AV3D4Quad(eleID, nodes[0], nodes[1], nodes[2], nodes[3], theMaterial);
  
  }

  else if (strcmp(argv[1],"AV3D4WithSensitivity") == 0){

  theEle = new AV3D4QuadWithSensitivity(eleID, nodes[0], nodes[1], nodes[2], nodes[3], theMaterial);
  
  }

  else {
    opserr << "WARNING element " << argv[1] << " type not recognized\n";
    return TCL_ERROR;
  }
  
  if(theTclDomain->addElement(theEle) == false) {
    printf("command: element AV3D4 %d - could not add element to domain!\n", 
      eleID); 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

