
// Written: Quan Gu, Yichao Gao and Zhijian Qiu  
// Created: 2015/01/25 
// TclAcousticCommand.cpp
//------------------------------------------

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <TclModelBuilder.h>


#include <VS3D4Quad.h>
#include <VS3D4QuadWithSensitivity.h>

extern void printCommand(int argc, TCL_Char **argv);

int 
TclModelBuilder_addVS3D4Quad(ClientData clientData, Tcl_Interp *interp, 
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
  double _E = 0.0;
  double _G = 0.0;
  double _rho = 1.0;
  double _R = 1.0;
  double _alphaN = 1.33;
  double _alphaT = 0.67;
  
  
  if(argc - eleArgStart < 9) {
    printf("WARNING insufficient arguments!\n");
    printCommand(argc, argv);
    printf("WANT: element VS3D4 eleTag? Node1? Node2? ... NodeX G rho R <alphaN> <alphaT>?\n");
    return TCL_ERROR;
  }
  
  
  if(Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK) {
    printf("WARNING invalid VS3D4 eleTag\n");
    return TCL_ERROR;
  }

  // printf("Element id is %d!\n", eleID);
  
  numNodes = 4;
  for(i = 0; i < numNodes; i++) {
    if(Tcl_GetInt(interp, argv[i+2+eleArgStart], &nodes[i]) != TCL_OK) {
      printf("Command: element VS3D4 %d - invalid integer tag %s\n", 
        eleID, argv[i+2+eleArgStart]);
      return TCL_ERROR;
    }
    // printf("Nodes[%d] is %d!\n", i+1, nodes[i]);
  }
  
  if(Tcl_GetDouble(interp, argv[7], &_E) != TCL_OK) {
    printf("WARNING invalid VS3D4 E\n");
    return TCL_ERROR;
  }
  
  if(Tcl_GetDouble(interp, argv[8], &_G) != TCL_OK) {
    printf("WARNING invalid VS3D4 G\n");
    return TCL_ERROR;
  }
  
  if(Tcl_GetDouble(interp, argv[9], &_rho) != TCL_OK) {
    printf("WARNING invalid VS3D4 rho\n");
    return TCL_ERROR;
  }
  
  if(Tcl_GetDouble(interp, argv[10], &_R) != TCL_OK) {
    printf("WARNING invalid VS3D4 R\n");
    return TCL_ERROR;
  }
  
  if(argc > 11) {
    if(Tcl_GetDouble(interp, argv[11], &_alphaN) != TCL_OK) {
       printf("WARNING invalid VS3D4 alphaN\n");
       return TCL_ERROR;
    }
  }
  
  if(argc > 12) {
    if(Tcl_GetDouble(interp, argv[12], &_alphaT) != TCL_OK) {
       printf("WARNING invalid VS3D4 alphaT\n");
       return TCL_ERROR;
    }
  }
  

  
Element *theEle = 0;
  
  if (strcmp(argv[1],"VS3D4") == 0) {

  theEle = new VS3D4Quad(eleID, nodes[0], nodes[1], nodes[2], 
    nodes[3], _E, _G, _rho, _R, _alphaN, _alphaT);
  
  }

  else if (strcmp(argv[1],"VS3D4WithSensitivity") == 0){

  theEle = new VS3D4QuadWithSensitivity(eleID, nodes[0], nodes[1], nodes[2], 
    nodes[3], _E, _G, _rho, _R, _alphaN, _alphaT);
  
  }
  
  if(theTclDomain->addElement(theEle) == false) {
    printf("command: element VS3D4 %d - could not add element to domain!\n", 
      eleID); 
    delete theEle;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}

  






