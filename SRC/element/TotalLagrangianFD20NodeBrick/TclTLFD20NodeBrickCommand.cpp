//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              Sept2003
//# UPDATE HISTORY:
//#
//#
//===============================================================================
#include <stdlib.h>
#include <string.h>
#include <OPS_Stream.h>
#include <Domain.h>

#include <ErrorHandler.h>
#include <TotalLagrangianFD20NodeBrick.h>
#include <TclModelBuilder.h>

#define NumNodes 20
#define NumDof 3

extern void printCommand(int argc, TCL_Char **argv);

int
TclModelBuilder_addTLFD20nBrick(ClientData clientData,
                                Tcl_Interp *interp,
                                int argc,
                                TCL_Char **argv,
                                Domain*theTclDomain,
                                TclModelBuilder *theTclBuilder,
                                int eleArgStart)
{
  // ensure the destructor has not been called -
  if (theTclBuilder == 0)
  {
    opserr << "command: element TotalLagrangianFD20nbrick - no modelbuilder\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  if ((argc-eleArgStart) < 26)
  {
    opserr << "command: element TotalLagrangianFD20nbrick - insufficient args - want " <<
      "element TotalLagrangianFD20nbrick eleTag? node1? node2? .. node20? matTag? bforce1? bforce2? bforce3? \n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int eleID, matID;
  int nodes[NumNodes];
  double bodyforces[NumDof];

  // read the eleTag
  if (Tcl_GetInt(interp, argv[1+eleArgStart], &eleID) != TCL_OK)
  {
    opserr << "command: element TotalLagrangianFD20nbrick - invalid integer tag " << argv[1+eleArgStart] << endln;
    return TCL_ERROR;
  }

  // read the 20 node tags
  for (int i=0; i<NumNodes; i++)
  {
      if (Tcl_GetInt(interp, argv[2+i+eleArgStart], &nodes[i]) != TCL_OK)
      {
  opserr << "command: element TotalLagrangianFD20nbrick " << eleID << " - invalid integer tag " <<
    argv[2+i+eleArgStart] << endln;
  return TCL_ERROR;
      }
  }

  // read in material tag & check the material exists in the model builder
  if (Tcl_GetInt(interp, argv[22+eleArgStart], &matID) != TCL_OK)
  {
    opserr << "command: element TotalLagrangianFD20nbrick " << eleID << " - invalid matID tag " <<
      argv[22+eleArgStart] << endln;
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = theTclBuilder->getNDMaterial(matID);

  if (theMaterial == 0)
  {
    opserr << "command: element TotalLagrangianFD20nbrick " << eleID <<
      " - no NDMaterial with tag " << argv[22+eleArgStart] << "exists\n";
    return TCL_ERROR;
  }

  // read the 3 bodyforce accel's
  for (int j=0; j<NumDof; j++)
  {
      if (Tcl_GetDouble(interp, argv[23+j+eleArgStart], &bodyforces[j]) != TCL_OK)
      {
  opserr << "command: element TotalLagrangianFD20nbrick " << eleID << " - invalid bodyforces tag " <<
    argv[23+j+eleArgStart] << endln;
  return TCL_ERROR;
      }
  }

  // now create the TwentyNodeBrick and add it to the Domain
  TotalLagrangianFD20NodeBrick *theEle = new TotalLagrangianFD20NodeBrick(eleID,
                                                                          nodes[ 0],
                                                                          nodes[ 1],
                                                                          nodes[ 2],
                                                                          nodes[ 3],
                                                                          nodes[ 4],
                                                                          nodes[ 5],
                                                                          nodes[ 6],
                                                                          nodes[ 7],
                                                                          nodes[ 8],
                                                                          nodes[ 9],
                                                                          nodes[10],
                                                                          nodes[11],
                                                                          nodes[12],
                                                                          nodes[13],
                                                                          nodes[14],
                                                                          nodes[15],
                                                                          nodes[16],
                                                                          nodes[17],
                                                                          nodes[18],
                                                                          nodes[19],
                                                                         *theMaterial,
                                                                          bodyforces[0],
                                                                          bodyforces[1],
                                                                          bodyforces[2]);

  if (theEle == 0)
  {
    opserr << "command: element TotalLagrangianFD20nbrick " << eleID << " - out of memory\n";
    return TCL_ERROR;
  }

  if (theTclDomain->addElement(theEle) == false)
  {
    opserr << "command: element TotalLagrangianFD20nbrick  - could not add ele: " << eleID << " to domain\n";
    delete theEle;
    return TCL_ERROR;
  }

   // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}




