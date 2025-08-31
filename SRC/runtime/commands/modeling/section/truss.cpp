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


#if 0

int
TclCommand_addTrussSection(ClientData clientData, Tcl_Interp *interp,
                              int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  enum class Positions : int {
    Material, Area, End
  };
  ArgumentTracker<Positions> tracker;
  std::set<int> positional;


  // section Type? $tag $material $area
  if (argc < 5) {
    opserr << OpenSees::PromptValueError
           << "incorrect number of arguments\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError
           << "failed to read integer tag\n";
    return TCL_ERROR;
  }

  int mtag;
  double area;
  for (int i=3; i<argc; ++i) {
    if (strcmp(argv[i], "-material") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &mtag) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read integer material tag\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Material);

    }
    else if (strcmp(argv[i], "-area") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read area\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Area);
    }
    else {
      positional.insert(i);
    }
  }
  for (int i: positional) {
    switch (tracker.current()) {
      case Positions::Material:
        if (Tcl_GetInt(interp, argv[i], &mtag) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read integer material tag\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Positions::Area:
        if (Tcl_GetDouble(interp, argv[i], &area) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read area\n";
          return TCL_ERROR;
        }
        tracker.increment();
        break;

      case Positions::End:
      default:
        opserr << OpenSees::PromptValueError
               << "unexpected argument\n";
        return TCL_ERROR;
    }
  }

  // if (Tcl_GetInt(interp, argv[3], &mtag) != TCL_OK) {
  //   opserr << OpenSees::PromptValueError
  //          << "failed to read integer material tag\n";
  //   return TCL_ERROR;
  // }

  // if (Tcl_GetDouble(interp, argv[4], &area) != TCL_OK) {
  //   opserr << OpenSees::PromptValueError
  //          << "failed to read area\n";
  //   return TCL_ERROR;
  // }

  UniaxialMaterial* mptr = builder->getTypedObject<UniaxialMaterial>(mtag);
  if (mptr == nullptr)
    return TCL_ERROR;

  UniaxialMaterial* pptr = nullptr;
  if (strcmp(argv[1], "Uniaxial") == 0) {
    area = 1.0;
  }
  else if (strcmp(argv[1], "PlaneStress") == 0) {
    if (!(pptr = mptr->getCopy("PlaneStress"))) {
      pptr = new PlaneStressMaterial(tag, *mptr);
    }
  }
  else {
    opserr << OpenSees::PromptValueError
           << "unknown plane section\n";
    return TCL_ERROR;
  }

  builder->addTaggedObject<TrussSection<UniaxialMaterial>>(*(new TrussSection<UniaxialMaterial>(tag, *pptr, area)));

  return TCL_OK;
}
#endif