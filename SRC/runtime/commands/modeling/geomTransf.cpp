//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: Geometric transformation command
//
// cmp
//
#include <tcl.h>
#include <string.h>
#include <assert.h>
#include <BasicModelBuilder.h>
#include <G3_Logging.h>

#include <LinearCrdTransf2d.h>
#include <LinearCrdTransf3d.h>
#include <PDeltaCrdTransf2d.h>
#include <PDeltaCrdTransf3d.h>
#include <CorotCrdTransf2d.h>
#include <CorotCrdTransf3d.h>
//#include <LinearCrdTransf2d02.h>
#include <LinearCrdTransf2dInt.h>
#include <CorotCrdTransfWarping2d.h>

#include <LinearFrameTransf3d.h>
#include <PDeltaFrameTransf3d.h>
#include <CorotFrameTransf3d.h>

//
// Create a coordinate transformation
//
int
TclCommand_addGeomTransf(ClientData clientData, Tcl_Interp *interp, int argc,
                         const char ** const argv)

{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  // Make sure there is a minimum number of arguments
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "insufficient number of geomTransf arguments\n";
    opserr << "    Expected: geomTransf type? tag? <specific transf args>\n";
    return TCL_ERROR;
  }

  int ndm = builder->getNDM();
  int ndf = builder->getNDF(); // number of degrees of freedom per node

  // create 2d coordinate transformation
  if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4)) {

    int crdTransfTag;
    Vector jntOffsetI(2),
           jntOffsetJ(2);

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "insufficient arguments\n"
        "    Expected: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tag\n"
        "    Expected: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>\n";
      return TCL_ERROR;
    }

    // Additional options at end of command

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n"
              "    Expected: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 2; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n"
              "    Expected: geomTransf type? tag? <-jntOffset dXi? dYi? dXj? dYj?>\n";
            return TCL_ERROR;
          }
        }
      }

      else {
        opserr << G3_ERROR_PROMPT << "bad command\n    Expected: geomTransf type? tag? "
                  "<-jntOffset dXi? dYi? dXj? dYj?>\n";
        opserr << "invalid: " << argv[argi] << "\n";
        return TCL_ERROR;
      }
    }

    // construct the transformation

    FrameTransform2d *crdTransf2d = nullptr;

    if (strcmp(argv[1], "Linear") == 0)
//    if (getenv("CRD")) {
//      crdTransf2d = new LinearCrdTransf2d02(crdTransfTag, jntOffsetI, jntOffsetJ);
//    } else
        crdTransf2d = new LinearCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "LinearInt") == 0)
      crdTransf2d =
          new LinearCrdTransf2dInt(crdTransfTag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
      crdTransf2d = new PDeltaCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0 && ndf == 3)
      crdTransf2d = new CorotCrdTransf2d(crdTransfTag, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0 && ndf == 4)
      crdTransf2d =
          new CorotCrdTransfWarping2d(crdTransfTag, jntOffsetI, jntOffsetJ);

    else {
      opserr << G3_ERROR_PROMPT << "invalid Type\n";
      opserr << argv[1] << "\n";
      return TCL_ERROR;
    }

    // add the transformation to the modelBuilder
    if (builder->addTaggedObject<FrameTransform2d>(*crdTransf2d) != TCL_OK)
      return TCL_ERROR;

  } else if (ndm == 3 && ndf == 6) {
    int crdTransfTag;
    Vector vecxzPlane(3);                // vector that defines local xz plane
    Vector jntOffsetI(3), jntOffsetJ(3); // joint offsets in global coordinates

    if (argc < 6) {
      opserr << G3_ERROR_PROMPT 
             << "insufficient arguments\n"
             << "    Expected: geomTransf type? tag? "
                "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? "
                "dZi? dXj? dYj? dZj? >\n";
      return TCL_ERROR;
    }

    int argi = 2;
    if (Tcl_GetInt(interp, argv[argi++], &crdTransfTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid tag\n    Expected: geomTransf type? tag? "
                "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? "
                "dZi? dXj? dYj? dZj? >\n";
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
             opserr << G3_ERROR_PROMPT << "Failed  to parse vectxz\n";
             return TCL_ERROR;
           }

        argi++;
        parsed_xz = true;
      }
    } 

    if (!parsed_xz) {
      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(0)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneX\n    Expected: geomTransf type? tag? "
                  "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? "
                  "dZi? dXj? dYj? dZj? >\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(1)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneY\n    Expected: geomTransf type? tag? "
                  "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? "
                  "dZi? dXj? dYj? dZj? >\n";
        return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[argi++], &vecxzPlane(2)) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid vecxzPlaneZ\n    Expected: geomTransf type? tag? "
                  "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? dYi? "
                  "dZi? dXj? dYj? dZj? >\n";
        return TCL_ERROR;
      }
    }

    // additional keyword options at end of command

    while (argi != argc) {
      if (strcmp(argv[argi], "-jntOffset") == 0) {
        argi++;
        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetI(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n    Expected: geomTransf "
                      "type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  "
                      "<-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >\n";
            return TCL_ERROR;
          }
        }

        for (int i = 0; i < 3; ++i) {
          if (argi == argc ||
              Tcl_GetDouble(interp, argv[argi++], &jntOffsetJ(i)) != TCL_OK) {
            opserr << G3_ERROR_PROMPT << "invalid jntOffset value\n    Expected: geomTransf "
                      "type? tag? vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  "
                      "<-jntOffset dXi? dYi? dZi? dXj? dYj? dZj? >\n";
            return TCL_ERROR;
          }
        }
      } else {
        opserr << G3_ERROR_PROMPT << "bad command\n    Expected: geomTransf type? tag? "
                  "vecxzPlaneX? vecxzPlaneY? vecxzPlaneZ?  <-jntOffset dXi? "
                  "dYi? dZi? dXj? dYj? dZj? > \n";
        opserr << "    invalid: " << argv[argi] << "\n";
        return TCL_ERROR;
      }
    }

    // construct the transformation object
    FrameTransform3d *crdTransf3d=nullptr;

    if (strcmp(argv[1], "Linear") == 0)
      if (!getenv("CRD"))
        crdTransf3d = new LinearFrameTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
      else
        crdTransf3d = new LinearCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "PDelta") == 0 ||
             strcmp(argv[1], "LinearWithPDelta") == 0)
      if (!getenv("CRD"))
        crdTransf3d = new PDeltaFrameTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
      else
        crdTransf3d = new PDeltaCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else if (strcmp(argv[1], "Corotational") == 0)
      // By default use new faster version
      if (!getenv("CRD"))
        crdTransf3d = new CorotFrameTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);
      else
        crdTransf3d = new CorotCrdTransf3d(crdTransfTag, vecxzPlane, jntOffsetI, jntOffsetJ);

    else {
      opserr << G3_ERROR_PROMPT << "invalid Type\n";
      return TCL_ERROR;
    }

    if (crdTransf3d == nullptr) {
      opserr << G3_ERROR_PROMPT << "Failed to create transform\n";
      return TCL_ERROR;
    }

    // add the transformation to the modelBuilder
    if (builder->addTaggedObject<FrameTransform3d>(*crdTransf3d) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "could not add "
                "geometric transformation to model Builder\n";
      return TCL_ERROR;
    }

  } else {
    opserr << G3_ERROR_PROMPT << "ndm = " << ndm << " and ndf = " << ndf
           << " is incompatible with available frame elements\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
