//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Matrix.h>
#include <Domain.h>
#include <Logging.h>
#include <Response.h>

#include <DummyStream.h>
#include <Element.h>

int
basicDeformation(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - basicDeformation eleTag? \n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicDeformation eleTag? dofNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicDeformation eleTag? dofNum? - could not read dofNum?
  \n"; return TCL_ERROR;
  }
  */

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "basicDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicDeformation";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; ++i) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int
basicForce(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - basicForce eleTag? \n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicForce eleTag? dofNum? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicDeformation eleTag? dofNum? - could not read dofNum \n";
    return TCL_ERROR;
  }
  */

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "basicDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicForce";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  int nbf = theVec.Size();

  char buffer[200];
  for (int i = 0; i < nbf; ++i) {
    sprintf(buffer, "%12.8f ", theVec(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  delete theResponse;

  return TCL_OK;
}

int
basicStiffness(ClientData clientData, Tcl_Interp *interp, int argc,
               TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - basicStiffness eleTag? \n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicStiffness eleTag? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  /*
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "basicDeformation eleTag? dofNum? - could not read dofNum?\n";
    return TCL_ERROR;
  }
  */

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "basicStiffness element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "basicStiffness";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMatrix = *(info.theMatrix);
  int nbf = theMatrix.noCols();

  char buffer[200];
  for (int i = 0; i < nbf; ++i) {
    for (int j = 0; j < nbf; j++) {
      sprintf(buffer, "%12.8f ", theMatrix(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

int
sectionForce(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "want - sectionForce eleTag? <secNum?> dof? \n";
    return TCL_ERROR;
  }

  int tag, dof;
  int secNum = 0;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionForce eleTag? secNum? dof? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }

  // Make this work for zeroLengthSection too
  int currentArg = 2;
  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[currentArg++], &secNum) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "sectionForce eleTag? secNum? dof? - could not read "
                "secNum? \n";
      return TCL_ERROR;
    }
  }
  if (Tcl_GetInt(interp, argv[currentArg++], &dof) != TCL_OK) {
    opserr
        << G3_ERROR_PROMPT << "sectionForce eleTag? secNum? dof? - could not read dof? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionForce element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "force";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;
  if (argc < 4) { // For zeroLengthSection
    argcc = 2;
    argvv[1] = c;
  }

  DummyStream dummy;
  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theVec(dof-1)));

  delete theResponse;

  return TCL_OK;
}

int
sectionDeformation(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 4) {
    opserr << G3_ERROR_PROMPT << "want - sectionDeformation eleTag? secNum? dof? \n";
    return TCL_ERROR;
  }

  int tag, secNum, dof;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionDeformation eleTag? secNum? dof? - could not "
              "read eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionDeformation eleTag? secNum? dof? - could not "
              "read secNum? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionDeformation eleTag? secNum? dof? - could not "
              "read dof? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionDeformation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "deformation";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;
  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theVec(dof-1)));

  delete theResponse;

  return TCL_OK;
}

int
sectionLocation(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "want - sectionLocation eleTag? secNum? \n";
    return TCL_ERROR;
  }

  int tag, secNum;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionLocation eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionLocation eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionLocation element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "integrationPoints";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theVec(secNum - 1)));

  delete theResponse;

  return TCL_OK;
}

int
sectionWeight(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "want - sectionWeight eleTag? secNum? \n";
    return TCL_ERROR;
  }

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << G3_ERROR_PROMPT << "sectionWeight eleTag? secNum? - could not read eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr
        << G3_ERROR_PROMPT << "sectionWeight eleTag? secNum? - could not read secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionWeight element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 1;
  char a[80] = "integrationWeights";
  const char *argvv[1];
  argvv[0] = a;

  DummyStream dummy;
  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Vector &theVec = *(info.theVector);

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theVec(secNum - 1)));

  delete theResponse;

  return TCL_OK;
}

int
sectionStiffness(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "want - sectionStiffness eleTag? secNum? \n";
    return TCL_ERROR;
  }

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionStiffness eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionStiffness eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionStiffness element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "stiffness";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; ++i) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer, "%12.8g ", theMat(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

int
sectionFlexibility(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "want - sectionFlexibility eleTag? secNum? \n";
    return TCL_ERROR;
  }

  int tag, secNum;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionFlexibility eleTag? secNum? - could not read "
              "eleTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "sectionFlexibility eleTag? secNum? - could not read "
              "secNum? \n";
    return TCL_ERROR;
  }

  Element *theElement = the_domain->getElement(tag);
  if (theElement == nullptr) {
    opserr << G3_ERROR_PROMPT << "sectionFlexibility element with tag " << tag
           << " not found in domain \n";
    return TCL_ERROR;
  }

  int argcc = 3;
  char a[80] = "section";
  char b[80];
  sprintf(b, "%d", secNum);
  char c[80] = "flexibility";
  const char *argvv[3];
  argvv[0] = a;
  argvv[1] = b;
  argvv[2] = c;

  DummyStream dummy;

  Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
  if (theResponse == nullptr) {
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(0.0));
    return TCL_OK;
  }

  theResponse->getResponse();
  Information &info = theResponse->getInformation();

  const Matrix &theMat = *(info.theMatrix);
  int nsdof = theMat.noCols();

  char buffer[200];
  for (int i = 0; i < nsdof; ++i) {
    for (int j = 0; j < nsdof; j++) {
      sprintf(buffer, "%12.8g ", theMat(i, j));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  delete theResponse;

  return TCL_OK;
}

int
sectionTag(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
    assert(clientData != nullptr);
    Domain *the_domain = (Domain*)clientData;

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "want - sectionTag eleTag? secNum? \n";
      return TCL_ERROR;
    }

    int tag, secNum;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "sectionTag eleTag? secNum? - could not read "
                "eleTag? \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "sectionTag eleTag? secNum? - could not read "
                "secNum? \n";
      return TCL_ERROR;
    }

    Element *theElement = the_domain->getElement(tag);
    if (theElement == nullptr) {
      opserr << G3_ERROR_PROMPT << "sectionFlexibility element with tag " << tag
             << " not found in domain \n";
      return TCL_ERROR;
    }

    int argcc = 1;
    char a[80] = "sectionTags";
    const char *argvv[1];
    argvv[0] = a;

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == nullptr) {
        return TCL_ERROR;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const ID &theID = *(info.theID);
    int Np = theID.Size();

    if (secNum > 0 && secNum <= Np) { // One IP
      int value = theID(secNum-1);
      Tcl_SetObjResult(interp, Tcl_NewIntObj(value));

    } else {
      // All IPs in a list
      Tcl_Obj* list = Tcl_NewListObj(Np, nullptr);

      for (int i = 0; i < Np; i++)
        Tcl_ListObjAppendElement(interp, list, Tcl_NewIntObj(theID(i)));


      Tcl_SetObjResult(interp, list);
    }

    delete theResponse;

    return TCL_OK;
}

int
sectionDisplacement(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
    assert(clientData != nullptr);
    Domain *theDomain = (Domain*)clientData;

    if (argc < 3) {
      opserr << G3_ERROR_PROMPT << "want - sectionLocation eleTag? secNum? \n";
      return TCL_ERROR;
    }

    int tag, secNum;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "sectionLocation eleTag? secNum? - could not read "
                "eleTag? \n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[2], &secNum) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "sectionLocation eleTag? secNum? - could not read "
                "secNum? \n";
      return TCL_ERROR;
    }


    bool local = false;
    if (argc > 4) {
      if (strstr(argv[4], "local") != 0)
        local = true;
    }


    Element *theElement = theDomain->getElement(tag);
    if (theElement == nullptr) {
        opserr << "WARNING sectionDisplacement element with tag " << tag << " not found in domain \n";
        return TCL_ERROR;
    }

    int argcc = 2;
    char a[80] = "sectionDisplacements";
    const char *argvv[2];
    argvv[0] = a;
    if (local)
      argvv[1] = "local";
    else
      argvv[1] = "global";

    DummyStream dummy;

    Response *theResponse = theElement->setResponse(argvv, argcc, dummy);
    if (theResponse == nullptr) {
        return TCL_ERROR;
    }

    theResponse->getResponse();
    Information &info = theResponse->getInformation();

    const Matrix &theMatrix = *(info.theMatrix);
    if (secNum <= 0 || secNum > theMatrix.noRows()) {
        opserr << "WARNING invalid secNum\n";
        delete theResponse;
        return TCL_ERROR;
    }


    const int nc = theMatrix.noCols();
    Tcl_Obj* list = Tcl_NewListObj(nc, nullptr);
    for (int i=0; i<nc; i++)
      Tcl_ListObjAppendElement(interp, list, Tcl_NewDoubleObj(theMatrix(secNum-1, i)));

    Tcl_SetObjResult(interp, list);

    delete theResponse;

    return TCL_OK;
}
