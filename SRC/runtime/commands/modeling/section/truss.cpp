//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements a command to add a uniaxial section.
// Written: cmp
//
#include <set>
#include <tcl.h>
#include <string.h>
#include <Parsing.h>
#include <Logging.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>
#include <UniaxialMaterial.h>
#include <SectionAggregator.h>


int 
TclCommand_addUniaxialSection(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // section Uniaxial tag? material? code?
  if (argc < 5) {
    opserr << "WARNING insufficient arguments\n";
    return TCL_ERROR;
  }
  
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING invalid section tag\n";
    return TCL_ERROR;
  } 

  int mat;
  if (Tcl_GetInt(interp, argv[3], &mat) != TCL_OK) {
    opserr << "WARNING invalid uniaxial material tag\n";
    return TCL_ERROR;
  } 

  
  UniaxialMaterial* material = builder->getTypedObject<UniaxialMaterial>(mat);
  if (material == nullptr)
    return TCL_ERROR;
  
  UniaxialMaterial *theMats[1] = {material};

  //
  int code;
  const char* type = argv[4];
  if (strcmp(type,"Mz") == 0)
    code = SECTION_RESPONSE_MZ;
  else if (strcmp(type,"P") == 0)
    code = SECTION_RESPONSE_P;
  else if (strcmp(type,"Vy") == 0)
    code = SECTION_RESPONSE_VY;
  else if (strcmp(type,"My") == 0)
    code = SECTION_RESPONSE_MY;
  else if (strcmp(type,"Vz") == 0)
    code = SECTION_RESPONSE_VZ;
  else if (strcmp(type,"T") == 0)
    code = SECTION_RESPONSE_T;
  else {
    opserr << "WARNING invalid code " << type << "\n";
    return TCL_ERROR;
  }

  ID codeID(1);
  codeID(0) = code;
  return builder->addTaggedObject<FrameSection>(*new SectionAggregator(tag, 1, theMats, codeID));
}

