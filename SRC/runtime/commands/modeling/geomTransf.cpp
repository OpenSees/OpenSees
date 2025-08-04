//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
//===----------------------------------------------------------------------===//

/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
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

//
// Description: Geometric transformation command
//
// cmp
//
#include <tcl.h>
#include <string>
#include <string.h>
#include <assert.h>
#include <Logging.h>
#include <BasicModelBuilder.h>

#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>
#include <LinearCrdTransf2dInt.h>
#include <CorotCrdTransfWarping2d.h>


#include <BasicFrameTransf.h>
#include <transform/FrameTransformBuilder.hpp>

using namespace OpenSees;

int 
TclCommand_addTransformBuilder(ClientData clientData, Tcl_Interp *interp, int argc,
                         const char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << OpenSees::PromptValueError << "insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int ndm = builder->getNDM();
  
  if (ndm != 2 && ndm != 3)
    return TCL_ERROR;

  int tag;
  const char *name = argv[1];
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid tag\n";
    return TCL_ERROR;
  }

  FrameTransformBuilder& transform = *new FrameTransformBuilder(ndm, tag, name);

  // Parse orientation vector
  int argi = 3;
  bool parsed_xz = (ndm == 2);
  int argxz = 0;
  while (argi != argc) {
    if (strcmp(argv[argi], "-jntOffset") == 0) {
      argi++;
      transform.offsets[1].zero();
      for (int i = 0; i < ndm; ++i) {
        if (argi == argc ||
            Tcl_GetDouble(interp, argv[argi++], &transform.offsets[1][i]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid offset at end I\n";
          return TCL_ERROR;
        }
      }

      transform.offsets[2].zero();
      for (int i = 0; i < ndm; ++i) {
        if (argi == argc ||
            Tcl_GetDouble(interp, argv[argi++], &transform.offsets[2][i]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid offset at end J\n";
          return TCL_ERROR;
        }
      }
    }

    else if (strcmp(argv[argi], "-orient") == 0) {
      // -orient {x y z}
      argi++;
      if (parsed_xz) {
        opserr << OpenSees::PromptValueError 
               << "orientation already provided\n";
        return TCL_ERROR;
      }

      if (argi == argc) {
        opserr << OpenSees::PromptValueError 
               << "missing orientation vector\n";
        return TCL_ERROR;
      }

      const char ** xzarg;
      int xznum;
      Tcl_SplitList(interp, argv[argi], &xznum, &xzarg);
      if (xznum != 3) {
        Tcl_Free((char *)xzarg);
        opserr << OpenSees::PromptValueError 
               << "invalid orientation vector\n";
        return TCL_ERROR;
      }
      for (int i=0; i<3; ++i)
        if (Tcl_GetDouble(interp, xzarg[i], &transform.vz[i]) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "failed  to parse vecxz\n";
          return TCL_ERROR;
        }
      
      Tcl_Free((char *)xzarg);
      argi++;
      parsed_xz = true;
    }

    else if (strcmp(argv[argi], "-pull") == 0) {
      // -pull iter|incr|init|default
      argi++;
      if (argi == argc) {
        opserr << OpenSees::PromptValueError 
               << "missing pull type\n";
        return TCL_ERROR;
      }

      if (strcmp(argv[argi], "iter") == 0) {
        transform.offset_flags |= LogIter;
      } else if (strcmp(argv[argi], "incr") == 0) {
        transform.offset_flags |= LogIncr;
      } else if (strcmp(argv[argi], "init") == 0) {
        transform.offset_flags |= LogInit;
      } else if (strcmp(argv[argi], "default") == 0) {
        transform.offset_flags |= LogDefault;
      } else {
        opserr << OpenSees::PromptValueError 
               << "invalid pull type\n";
        return TCL_ERROR;
      }
    }

    else if (strcmp(argv[argi], "-offset-local") == 0) {
      transform.offset_flags |= OffsetLocal;
      argi++;
    }
    else if (strcmp(argv[argi], "-offset-length") == 0) {
      transform.offset_flags |= OffsetNormalized;
      argi++;
    }

    else if (strcmp(argv[argi], "-offset") == 0) {
      // -offset {1 {x y z}; 2 {x y z}}
      argi++;
      if (argi == argc) {
        opserr << OpenSees::PromptValueError << "missing offset block.\n";
        return TCL_ERROR;
      }
      constexpr static const char * const offset_nodes[] = {"1", "2", "3", "4", "5", "6", "7", "8", "9"};
      constexpr static int MAX_OFFSETS = sizeof(offset_nodes)/sizeof(offset_nodes[0]);
      auto TclCommand_addOffset  = [](ClientData cd, Tcl_Interp* interp, int oargc, const char** const oargv) ->int {
        
        FrameTransformBuilder* transform = static_cast<FrameTransformBuilder*>(cd);
        assert(transform != nullptr);
        if (oargc < 2) {
          opserr << OpenSees::PromptValueError 
                 << "insufficient number of offset arguments\n";
          return TCL_ERROR;
        }
        int offset_tag;
        if (Tcl_GetInt(interp, oargv[0], &offset_tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid offset tag\n";
          return TCL_ERROR;
        }
        if (offset_tag < 1 || offset_tag > MAX_OFFSETS) {
          opserr << OpenSees::PromptValueError 
                 << "invalid offset tag\n";
          return TCL_ERROR;
        }
        if (oargc < 2) {
          opserr << OpenSees::PromptValueError 
                 << "missing offset vector\n";
          return TCL_ERROR;
        }

        const char ** xzarg;
        int xznum;
        if (oargc == 2)
          Tcl_SplitList(interp, oargv[1], &xznum, &xzarg);
        else {
          xznum = oargc - 1;
          xzarg = oargv + 1;
        }

        for (int i=0; i<xznum; ++i)
          if (Tcl_GetDouble(interp, xzarg[i], &transform->offsets[offset_tag][i]) != TCL_OK) {
            if (oargc == 2)
              Tcl_Free((char *)xzarg);
            opserr << OpenSees::PromptValueError 
                   << "failed to parse offset vector\n";
            return TCL_ERROR;
          }
        if (oargc == 2)
          Tcl_Free((char *)xzarg);
        return TCL_OK;
      };

      for (int i=0; i<MAX_OFFSETS; i++) {
        Tcl_CreateCommand(interp, offset_nodes[i], 
                          TclCommand_addOffset, &transform, nullptr);
      }

      if (Tcl_Eval(interp, argv[argi]) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "failed to parse offset block\n";
        return TCL_ERROR;
      } else {
        for (int i=0; i<MAX_OFFSETS; i++) {
          Tcl_DeleteCommand(interp, offset_nodes[i]);
        }
      }
      argi++;
    }
    else {
      if (argxz == 0) {
        argxz = argi;
      }
      argi++;
    }
  } // keyword arguments


  // Ensure orientation was provided
  if (!parsed_xz) {
    if (argxz == 0 || argxz >= argc) {
      opserr << OpenSees::PromptValueError << "missing orientation vector\n";
      return TCL_ERROR;
    }
    const char ** xzarg;
    int xznum;
    Tcl_SplitList(interp, argv[argxz], &xznum, &xzarg);
    if (xznum == 3) {
      for (int i=0; i<3; ++i)
          if (Tcl_GetDouble(interp, xzarg[i], &transform.vz[i]) != TCL_OK) {
            Tcl_Free((char *)xzarg);
            opserr << OpenSees::PromptValueError 
                   << "Failed to parse vectxz\n";
            return TCL_ERROR;
          }
      argi++;
      parsed_xz = true;
    }
    Tcl_Free((char *)xzarg);
  } 

  if (!parsed_xz) {
    if (argxz+3 > argc) {
      opserr << OpenSees::PromptValueError 
             << "missing orientation vector\n";
      return TCL_ERROR;
    }
    for (int i=0; i<3; i++)
      if (Tcl_GetDouble(interp, argv[argxz++], &transform.vz[i]) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid vecxz component\n";
        return TCL_ERROR;
      }
  }

  if (builder->addTaggedObject<FrameTransformBuilder>(transform) != TCL_OK)
    return TCL_ERROR;
  
  return TCL_OK;
}


int
TclCommand_addGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,
                         const char ** const argv)

{
  if (TclCommand_addTransformBuilder(clientData, interp, argc, argv) != TCL_OK)
    return TCL_ERROR;
  
  if (strcmp(argv[0], "transform") == 0)
    return TCL_OK;


  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);


  // Make sure there is a minimum number of arguments
  if (argc < 3) {
    opserr << OpenSees::PromptValueError 
           << "insufficient number of arguments\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid tag"
           << "\n";
    return TCL_ERROR;
  }

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm == 3 && !getenv("CRD")) {
    auto tb = builder->getTypedObject<FrameTransformBuilder>(tag);
    if (tb == nullptr) {
      opserr << OpenSees::PromptValueError 
             << "transformation not found with tag " << tag 
             << "\n";
      return TCL_ERROR;
    }
    CrdTransf* t = new BasicFrameTransf3d(tb->template create<2,6>());
    return builder->addTaggedObject<CrdTransf>(*t);
  }

  //
  // 2D Case
  //
  if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4)) {
    Vector jntOffsetI(2),
           jntOffsetJ(2);

    int argi = 3;

    // Additional options at end of command

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "invalid jntOffset value"
                   << "\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << OpenSees::PromptValueError 
                   << "invalid jntOffset value"
                   << "\n";
            return TCL_ERROR;
          }
        }
      }

      else {
        opserr << OpenSees::PromptValueError 
               << "unexpected argument " << argv[argi] 
               << "\n";
        return TCL_ERROR;
      }
    }

    //
    // construct the transformation
    //

    CrdTransf *crdTransf2d = nullptr;

    if (strcmp(argv[1], "Linear") == 0)
        crdTransf2d = new LinearCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "LinearInt") == 0)
      crdTransf2d =
          new LinearCrdTransf2dInt(tag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
      crdTransf2d = new PDeltaCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if ((strcmp(argv[1], "Corotational") == 0 || strcmp(argv[1], "Corotational02") == 0) && ndf == 3)
      crdTransf2d = new CorotCrdTransf2d(tag, jntOffsetI, jntOffsetJ);

    else if ((strcmp(argv[1], "Corotational") == 0 || strcmp(argv[1], "Corotational02") == 0) && ndf == 4)
      crdTransf2d =
          new CorotCrdTransfWarping2d(tag, jntOffsetI, jntOffsetJ);

    else {
      opserr << OpenSees::PromptValueError 
             << "invalid Type: " << argv[1]
             << "\n";
      return TCL_ERROR;
    }

    //
    if (builder->addTaggedObject<CrdTransf>(*crdTransf2d) != TCL_OK)
      return TCL_ERROR;
  }

  else if (ndm == 3 && ndf >= 6) {
    // vector that defines local xz plane
    Vector vecxzPlane(3);
    // joint offsets in global coordinates
    Vector jntOffsetI(3), jntOffsetJ(3);

    if (argc < 6) {
      opserr << OpenSees::PromptValueError 
             << "insufficient arguments"
             << "\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &tag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid tag\n";
      return TCL_ERROR;
    }

    // parse orientation vector
    bool parsed_xz = false;
    if (!parsed_xz) {
      const char ** xzarg;
      int xznum;
      Tcl_SplitList(interp, argv[argi], &xznum, &xzarg);
      if (xznum == 3) {
        for (int i=0; i<3; ++i)
           if (Tcl_GetDouble(interp, xzarg[i], &vecxzPlane(i)) != TCL_OK) {
             opserr << OpenSees::PromptValueError 
                    << "Failed  to parse vectxz\n";
             return TCL_ERROR;
           }

        argi++;
        parsed_xz = true;
      }
    } 

    if (!parsed_xz) {
      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(0)) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid vecxzPlaneX\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(1)) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid vecxzPlaneY\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(2)) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
               << "invalid vecxzPlaneZ\n";
        return TCL_ERROR;
      }
    }

    //
    // Additional keyword options
    //

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid jntOffset\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid jntOffset\n";
            return TCL_ERROR;
          }
        }
      } else {
        opserr << OpenSees::PromptValueError << "unexpected argument: " << argv[argi] << "\n";
        return TCL_ERROR;
      }
    }

    //
    // construct the transformation
    //
    CrdTransf *crdTransf3d=nullptr;

    if (strcmp(argv[1], "Linear") == 0)
        crdTransf3d = new LinearCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
        crdTransf3d = new PDeltaCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0)
      crdTransf3d = new CorotCrdTransf3d(tag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else {
      opserr << OpenSees::PromptValueError << "invalid Type\n";
      return TCL_ERROR;
    }

    if (crdTransf3d == nullptr) {
      opserr << OpenSees::PromptValueError << "Failed to create transform\n";
      return TCL_ERROR;
    }

    // add the transformation to the modelBuilder
    if (builder->addTaggedObject<CrdTransf>(*crdTransf3d) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "Failed to add transformation to model\n";
      return TCL_ERROR;
    }

  }

  else {
    opserr << OpenSees::PromptValueError 
           << "ndm = " << ndm << " and ndf = " << ndf
           << " is incompatible with available frame elements\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
