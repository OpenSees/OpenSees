//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Written: cmp
//

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <tcl.h>
#include <Vector.h>
#include <DummyStream.h>
#include <G3_Logging.h>
#include <Response.h>
#include <FrameSection.h>
#include <BasicModelBuilder.h>

static Tcl_CmdProc SectionTest_setStrainSection;
static Tcl_CmdProc SectionTest_getStressSection;
static Tcl_CmdProc SectionTest_getTangSection;
static Tcl_CmdProc SectionTest_getResponseSection;


static int count;
static int countsTillCommit;

// invoke Section $tag $commands
int
TclCommand_useCrossSection(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  assert(clientData != nullptr);
  // TODO: Parse tag properly
  SectionForceDeformation *theSection = 
    ((BasicModelBuilder*)clientData)->getTypedObject<FrameSection>(std::atoi(argv[2]));

  if (theSection == nullptr) {
    opserr << G3_ERROR_PROMPT << "no section found with tag '" << argv[2] << "'\n";
    return TCL_ERROR;
  } else {
    // theSection = theSection->getCopy();
  }

  //
  //
  //
  Tcl_CreateCommand(interp, "update",
                    SectionTest_setStrainSection, (ClientData)theSection, NULL);

  Tcl_CreateCommand(interp, "stress",
                    SectionTest_getStressSection, (ClientData)theSection, NULL);

  Tcl_CreateCommand(interp, "tangent", SectionTest_getTangSection,
                    (ClientData)theSection, NULL);

  Tcl_CreateCommand(interp, "responseSectionTest",
                    SectionTest_getResponseSection, (ClientData)theSection, NULL);

  Tcl_CreateCommand(interp, "response",
                    SectionTest_getResponseSection, (ClientData)theSection, NULL);
  //
  //
  Tcl_Eval(interp, argv[3]);
  //
  //

  Tcl_DeleteCommand(interp, "strain");
  Tcl_DeleteCommand(interp, "stress");
  Tcl_DeleteCommand(interp, "tangent");
  Tcl_DeleteCommand(interp, "responseSectionTest");

  return TCL_OK;
}

static int
SectionTest_setStrainSection(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  SectionForceDeformation *theSection = (SectionForceDeformation*)clientData;

  // check number of arguments in command line
  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "bad command - want: strainSectionTest strain?\n";
    return TCL_ERROR;
  }

  // get the sectionID form command line
  // Need to set the data based on argc, otherwise it crashes when setting
  // "data(i-1) = strain"
  static Vector data(argc - 1);
  double strain;
  for (int i = 1; i < argc; ++i) {
    if (Tcl_GetDouble(interp, argv[i], &strain) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "could not read strain: strainSectionTest strain1? "
                "strain2? ... strainN?\n";
      return TCL_ERROR;
    }
    data(i - 1) = strain;
  }

  theSection->setTrialSectionDeformation(data);

  if (count == countsTillCommit) {
    theSection->commitState();
    count = 1;
  } else
    count++;

  return TCL_OK;
}

static int
SectionTest_getStressSection(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  SectionForceDeformation *theSection = (SectionForceDeformation*)clientData;
  const Vector &stress = theSection->getStressResultant();
  for (int i = 0; i < stress.Size(); ++i) {
    char buffer[40];
    sprintf(buffer, "%.10e ", stress(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }
  return TCL_OK;
}

static int
SectionTest_getTangSection(ClientData clientData, Tcl_Interp *interp,
                                int argc, TCL_Char ** const argv)
{
  SectionForceDeformation *theSection = (SectionForceDeformation*)clientData;

  const Matrix &tangent = theSection->getSectionTangent();
  for (int i = 0; i < tangent.noRows(); ++i)
    for (int j = 0; j < tangent.noCols(); j++) {
      char buffer[40];
      sprintf(buffer, "%.10e ", tangent(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  return TCL_OK;
}

static int
SectionTest_getResponseSection(ClientData clientData, Tcl_Interp *interp,
                                    int argc, TCL_Char ** const argv)
{
  SectionForceDeformation *theSection = (SectionForceDeformation*)clientData;
  DummyStream dummy;
  Response *theResponse =
      theSection->setResponse(argv + 1, argc - 1, dummy);

  if (theResponse == nullptr) {
    opserr << G3_ERROR_PROMPT << "Response returned a null pointer\n";
    return TCL_ERROR;
  }

  if (theResponse->getResponse() < 0) {
    delete theResponse;
    opserr << G3_ERROR_PROMPT << "Failed to get response\n";
    return TCL_ERROR;
  }

  Information &eleInfo = theResponse->getInformation();
  const Vector &data = eleInfo.getData();

  for (int i = 0; i < data.Size(); ++i) {
    char buffer[40];
    sprintf(buffer, "%.10e ", data(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;
  return TCL_OK;
}
