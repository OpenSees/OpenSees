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
// Written: cmp
//
#include <tcl.h>
#include <string.h>
#include <Parsing.h>
#include <Logging.h>
#include <BasicModelBuilder.h>
#include "PlaneSection.h"
#include <ArgumentTracker.h>
#include <PlaneStrainMaterial.h>
#include <PlaneStressMaterial.h>

int
TclCommand_addPlaneSection(ClientData clientData, Tcl_Interp *interp,
                           int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  enum class Positions : int {
    Material, Thickness, End
  };
  ArgumentTracker<Positions> tracker;
  std::set<int> positional;

  // section Plane[Strain|Stress] $tag $material $thickness
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
  double thickness;
  for (int i=3; i<argc; ++i) {
    if (strcmp(argv[i], "-material") == 0) {
      if (argc == ++i || Tcl_GetInt(interp, argv[i], &mtag) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read integer material tag\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Material);

    }
    else if (strcmp(argv[i], "-thickness") == 0) {
      if (argc == ++i || Tcl_GetDouble(interp, argv[i], &thickness) != TCL_OK) {
        opserr << OpenSees::PromptValueError
               << "failed to read thickness\n";
        return TCL_ERROR;
      }
      tracker.consume(Positions::Thickness);
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

      case Positions::Thickness:
        if (Tcl_GetDouble(interp, argv[i], &thickness) != TCL_OK) {
          opserr << OpenSees::PromptValueError
                 << "failed to read thickness\n";
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

  // if (Tcl_GetDouble(interp, argv[4], &thickness) != TCL_OK) {
  //   opserr << OpenSees::PromptValueError
  //          << "failed to read thickness\n";
  //   return TCL_ERROR;
  // }

  NDMaterial* mptr = builder->getTypedObject<NDMaterial>(mtag);
  if (mptr == nullptr)
    return TCL_ERROR;

  NDMaterial* pptr = nullptr;
  if (strcmp(argv[1], "PlaneStrain") == 0) {
    if (!(pptr = mptr->getCopy("PlaneStrain"))) {
      pptr = new PlaneStrainMaterial(tag, *mptr);
    }
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

  builder->addTaggedObject<PlaneSection<NDMaterial>>(*(new PlaneSection<NDMaterial>(tag, *pptr, thickness)));

  return TCL_OK;
}
