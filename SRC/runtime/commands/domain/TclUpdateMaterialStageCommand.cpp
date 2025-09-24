//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
// 
// Description: This command is used to update a PressureDependMultiYield,
// PressureDependMultiYield02, PressureIndependMultiYield, or FluidSolidPorous
// material. To conduct a seismic analysis, two stages should be followed.
// First, during the application of gravity load (and static loads if any), set
// material stage to 0, and material behavior is linear elastic (with Gr and Br
// as elastic moduli). A FluidSolidPorous material does not contribute to the
// material response if its stage is set to 0. After the application of gravity
// load, set material stage to 1 or 2. In case of stage 2, all the elastic
// material properties are then internally determined at the current effective
// confinement, and remain constant thereafter. In the subsequent dynamic
// (fast) loading phase(s), the deviatoric stress-strain response is
// elastic-plastic (stage 1) or linear-elastic (stage 2), and the volumetric
// response remains linear-elastic. Please visit
// http://cyclic.ucsd.edu/opensees for examples.
//
// Written: ZHY
//
// $Revision: 1.16 $
// $Date: 2007-10-16 00:15:07 $
//
#include <ID.h>
#include <BasicModelBuilder.h>
#include <PressureIndependMultiYield.h>
#include <PressureDependMultiYield.h>
#include <FluidSolidPorousMaterial.h>
#include <Information.h>
#include <UniaxialMaterial.h>
#include <MaterialStageParameter.h>
#include <tcl.h>
#include <PyLiq1.h>
#include <TzLiq1.h>
#include <QzLiq1.h>

#include <string.h>


int
TclCommand_updateMaterialStage(ClientData clientData,
                               Tcl_Interp *interp, int argc,
                               TCL_Char ** const argv)
{

  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;
  Domain* domain = builder->getDomain();

  if (argc < 5) {
    opserr << "WARNING insufficient number of UpdateMaterialStage arguments\n";
    opserr << "Want: UpdateMaterialStage material matTag? stage value?"
           << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "-material") != 0) {
    opserr << "WARNING UpdateMaterialStage: Only accept parameter '-material' "
              "for now"
           << endln;
    return TCL_ERROR;
  }

  int materialTag, value;

  if (Tcl_GetInt(interp, argv[2], &materialTag) != TCL_OK) {
    opserr << "WARNING MYSstage: invalid material tag" << endln;
    return TCL_ERROR;
  }

  int parTag = domain->getNumParameters();
  parTag++;

  if (argc > 6) {
    if (strcmp(argv[5], "-parameter") == 0) {
      if (Tcl_GetInt(interp, argv[6], &parTag) != TCL_OK) {
        opserr << "WARNING UpdateMaterialStage: invalid parameter tag used"
               << endln;
        return TCL_ERROR;
      }
    }
  }

  MaterialStageParameter *theParameter =
      new MaterialStageParameter(parTag, materialTag);
  if (domain->addParameter(theParameter) == false) {
    opserr << "WARNING could not add updateMaterialStage - "
              "MaterialStageParameter to domain"
           << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[3], "-stage") != 0) {
    opserr
        << "WARNING UpdateMaterialStage: Only accept parameter '-stage' for now"
        << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4], &value) != TCL_OK) {
    opserr << "WARNING UpdateMaterialStage: invalid parameter value" << endln;
    return TCL_ERROR;
  }

  domain->updateParameter(parTag, value);

  domain->removeParameter(parTag);

  delete theParameter;

  return TCL_OK;
}


int
TclBasicBuilderUpdateParameterCommand(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{

  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 5) {
    opserr << "WARNING insufficient number of updateParameter arguments\n";
    opserr << "Want: updateParameter -material matNum? -param? newValue?"
           << endln;
    return TCL_ERROR;
  }

  if (strcmp(argv[1], "-material") != 0) {
    opserr
        << "WARNING UpdateParameter: Only accept parameter '-material' for now"
        << endln;
    return TCL_ERROR;
  }

  int tag, id;
  double value;

  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << "WARNING UpdateParameter: invalid material tag" << endln;
    return TCL_ERROR;
  }

  // TODO: This will print an error message if not found; maybe
  // getTypedObject should accept flag that tells it not to print
  NDMaterial *a = builder->getTypedObject<NDMaterial>(tag);

  if (a == 0) {
    // opserr << "WARNING UpdateParameter: couldn't get NDmaterial tagged: " <<
    // tag << endln; return TCL_ERROR;
    UniaxialMaterial *a = builder->getTypedObject<UniaxialMaterial>(tag);
    if (a == 0) {
      opserr
          << "WARNING UpdateParameter: couldn't get Uniaxialmaterial tagged: "
          << tag << endln;
      return TCL_ERROR;
    }
    if (strcmp(argv[3], "-E") == 0) {
      if (Tcl_GetDouble(interp, argv[4], &value) != TCL_OK) {
        opserr << "WARNING UpdateParameter: invalid parameter value" << endln;
        return TCL_ERROR;
      }
      Information info;
      info.setDouble(value);
      a->updateParameter(0, info);
    } else if (strcmp(argv[3], "-fy") == 0) {
      if (Tcl_GetDouble(interp, argv[4], &value) != TCL_OK) {
        opserr << "WARNING UpdateParameter: invalid parameter value" << endln;
        return TCL_ERROR;
      }
      Information info;
      info.setDouble(value);
      a->updateParameter(1, info);
    } else {
      opserr << "WARNING UpdateParameter: Only accept parameter '-E' or '-fy' "
                "for now"
             << endln;
      return TCL_ERROR;
    }
    return TCL_OK;
  }

  if (strcmp(argv[3], "-refG") == 0)
    id = 10;
  else if (strcmp(argv[3], "-refB") == 0)
    id = 11;
  else {
    opserr << "WARNING UpdateParameter: Only accept parameter '-refG' or "
              "'-refB' for now"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[4], &value) != TCL_OK) {
    opserr << "WARNING UpdateParameter: invalid parameter value" << endln;
    return TCL_ERROR;
  }

  const char *c = a->getType();

  if (strcmp(c, "PlaneStrain") == 0 || strcmp(c, "ThreeDimensional") == 0) {
    Information info;
    info.setDouble(value);
    a->updateParameter(id, info);
  } else {
    opserr << "WARNING UpdateParameter: The tagged is not a " << endln;
    opserr << "PressureDependMultiYield/PressureIndependMultiYield/"
              "FluidSolidPorous material. "
           << endln;
    return TCL_ERROR;
  }

  return TCL_OK;
}
