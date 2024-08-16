//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements commands for interacting with nodes
// in the domain.
//
// All commands assume a Domain* is passed as clientData.
//
#include <math.h>
#include <assert.h>
#include <set>
#include <vector>
#include <algorithm>
#include <tcl.h>
#include <Logging.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Domain.h>
#include <DOF_Group.h>
#include <Node.h>
#include <NodeIter.h>
#include <Pressure_Constraint.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>

// TODO(cmp): Remove global vars
static char *resDataPtr  = nullptr;
static int   resDataSize = 0;


int
getNodeTags(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  NodeIter &nodeIter = the_domain->getNodes();

  Node *node;
  char buffer[20];
  while ((node = nodeIter()) != nullptr) {
    sprintf(buffer, "%d ", node->getTag());
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
findID(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *theDomain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - findNodesWithID ?id\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING eleForce eleTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  NodeIter &theNodes = theDomain->getNodes();
  Node *theNode;
  char buffer[20] = {0};

  while ((theNode = theNodes()) != nullptr) {
    DOF_Group *theGroup = theNode->getDOF_GroupPtr();
    if (theGroup != nullptr) {
      const ID &nodeID = theGroup->getID();
      for (int i = 0; i < nodeID.Size(); ++i) {
        if (nodeID(i) == tag) {
          sprintf(buffer, "%d ", theNode->getTag());
          Tcl_AppendResult(interp, buffer, NULL);
          break;
        }
      }
    }
  }

  return TCL_OK;
}

int
setNodeCoord(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 4) {
    opserr << "WARNING want - setNodeCoord nodeTag? dim? value?\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  int dim;
  double value;

  if (Tcl_GetInt(interp, argv[2], &dim) != TCL_OK) {
    opserr
        << "WARNING setNodeCoord nodeTag? dim? value? - could not read dim? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING setNodeCoord nodeTag? dim? value? - could not read "
              "value? \n";
    return TCL_ERROR;
  }

  Node *theNode = domain->getNode(tag);

  if (theNode == nullptr) {
    // TODO: add error message
    return TCL_ERROR;
  }

  //
  // TODO: Check dimensions

  Vector coords(theNode->getCrds());
  coords(dim - 1) = value;
  theNode->setCrds(coords);

  return TCL_OK;
}

int
nodeDisp(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeDisp nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeDisp nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeDisp nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, NodeData::Disp);

  if (nodalResponse == nullptr)
    // TODO: add error message
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeDisp nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  } else {
    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
nodeMass(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << "WARNING want - nodeMass nodeTag? nodeDOF?\n";
    return TCL_ERROR;
  }

  int tag, dof;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeMass nodeTag? nodeDOF? \n";
    return TCL_ERROR;
  }

  char buffer[40];

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "WARNING nodeMass node " << tag << " not found" << endln;
    return TCL_ERROR;
  }
  int numDOF = theNode->getNumberDOF();
  if (dof < 1 || dof > numDOF) {
    opserr << "WARNING nodeMass dof " << dof << " not in range" << endln;
    return TCL_ERROR;
  } else {
    const Matrix &mass = theNode->getMass();
    sprintf(buffer, "%35.20f", mass(dof - 1, dof - 1));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

int
nodePressure(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING: want - nodePressure nodeTag?\n";
    return TCL_ERROR;
  }
  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING: nodePressure " << argv[1] << "\n";
    return TCL_ERROR;
  }
  double pressure = 0.0;
  Pressure_Constraint *thePC = the_domain->getPressure_Constraint(tag);
  if (thePC != nullptr) {
    pressure = thePC->getPressure();
  }
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(pressure));
  return TCL_OK;
}

int
nodeBounds(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  const int requiredDataSize = 20*6;
  if (requiredDataSize > resDataSize) {
    if (resDataPtr != 0) {
      delete[] resDataPtr;
    }
    resDataPtr = new char[requiredDataSize];
    resDataSize = requiredDataSize;
  }

  for (int i = 0; i < requiredDataSize; ++i)
    resDataPtr[i] = '\n';

  const Vector &bounds = the_domain->getPhysicalBounds();

  int cnt = 0;
  for (int j = 0; j < 6; j++) {
    cnt += sprintf(&resDataPtr[cnt], "%.6e  ", bounds(j));
  }

  Tcl_SetResult(interp, resDataPtr, TCL_STATIC);

  return TCL_OK;
}

int
nodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeVel nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeVel nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeVel nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = the_domain->getNodeResponse(tag, NodeData::Vel);

  if (nodalResponse == nullptr)
    // TODO: add error message
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {
    if (size < dof)
    // TODO: add error message
      return TCL_ERROR;

    double value = (*nodalResponse)(dof);

    // now we copy the value to the tcl string that is returned
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  } else {

    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
setNodeVel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 4) {
    opserr << "WARNING want - setNodeVel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "WARNING setNodeVel -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING setNodeVel nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr
        << "WARNING setNodeVel nodeTag? dof? value?- could not read value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getVel();
    vel(dof) = value;
    theNode->setTrialVel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

int
setNodeDisp(ClientData clientData, Tcl_Interp *interp, int argc,
            TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 4) {
    opserr << "WARNING want - setNodeDisp nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeDisp nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "WARNING setNodeDisp -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getDisp();
    vel(dof) = value;
    theNode->setTrialDisp(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}


int
setNodeAccel(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  Domain *the_domain = (Domain *)clientData;

  if (argc < 4) {
    opserr << "WARNING want - setNodeAccel nodeTag? dof? value? <-commit>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;
  double value = 0.0;
  bool commit = false;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read "
              "nodeTag? \n";
    return TCL_ERROR;
  }

  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << "WARNING setNodeAccel -- node with tag " << tag << " not found"
           << endln;
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr
        << "WARNING setNodeDisp nodeTag? dof? value?- could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &value) != TCL_OK) {
    opserr << "WARNING setNodeAccel nodeTag? dof? value?- could not read "
              "value? \n";
    return TCL_ERROR;
  }
  if (argc > 4 && strcmp(argv[4], "-commit") == 0)
    commit = true;

  dof--;

  int numDOF = theNode->getNumberDOF();

  if (dof >= 0 && dof < numDOF) {
    Vector vel(numDOF);
    vel = theNode->getAccel();
    vel(dof) = value;
    theNode->setTrialAccel(vel);
  }
  if (commit)
    theNode->commitState();

  return TCL_OK;
}

int
nodeAccel(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  Domain *the_domain = (Domain *)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeAccel nodeTag? dof?\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeAccel nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeAccel nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = the_domain->getNodeResponse(tag, NodeData::Accel);
  if (nodalResponse == nullptr)
    // TODO: add error message
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {
    if (size < dof)
      return TCL_ERROR;

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  } else {
    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
nodeUnbalance(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeUnbalance nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeUnbalance nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeUnbalance nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, NodeData::UnbalancedLoad);

  if (nodalResponse == nullptr)
    // TODO: add error message
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeUnbalance nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  } else {
    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
nodeResponse(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain *)clientData;

  if (argc < 4) {
    opserr << "WARNING want - nodeResponse nodeTag? dof? responseID?\n";
    return TCL_ERROR;
  }

  int tag, dof, responseID;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? - could not read dof? \n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &responseID) != TCL_OK) {
    opserr << "WARNING nodeResponse nodeTag? dof? responseID? - could not read "
              "responseID? \n";
    return TCL_ERROR;
  }

  dof--;

  const Vector *nodalResponse =
      the_domain->getNodeResponse(tag, (NodeData)responseID);

  if (nodalResponse == 0 || nodalResponse->Size() < dof || dof < 0)
    // TODO: add error message
    return TCL_ERROR;

  Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  return TCL_OK;
}

int
nodeEigenvector(ClientData clientData, Tcl_Interp *interp, int argc,
                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 3) {
    opserr << "WARNING want - nodeEigenVector nodeTag? eigenVector? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int eigenvector = 0;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr
        << "WARNING nodeEigenvector nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2], &eigenvector) != TCL_OK) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
    return TCL_ERROR;
  }

  if (argc > 3) {
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK) {
      opserr
          << "WARNING nodeEigenvector nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;
  eigenvector--;
  Node *theNode = domain->getNode(tag);
  const Matrix &theEigenvectors = theNode->getEigenvectors();

  int size = theEigenvectors.noRows();
  int numEigen = theEigenvectors.noCols();

  if (eigenvector < 0 || eigenvector >= numEigen) {
    opserr << "WARNING nodeEigenvector nodeTag? dof? - eigenvecor too large\n";
    return TCL_ERROR;
  }

  if (dof >= 0) {
    if (dof >= size) {
      opserr << "WARNING nodeEigenvector nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    double value = theEigenvectors(dof, eigenvector);

    // now we copy the value to the Tcl string that is returned
    char buffer[40];
    sprintf(buffer, "%35.20f", value);
    Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  } else {

    char buffer[40];
    for (int i = 0; i < size; ++i) {
      double value = theEigenvectors(i, eigenvector);
      sprintf(buffer, "%35.20f", value);
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
calculateNodalReactions(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain *)clientData;

  int incInertia = 0;

  if (argc == 2) {
    if ((strcmp(argv[1], "-incInertia") == 0) ||
        (strcmp(argv[1], "-dynamical") == 0) ||
        (strcmp(argv[1], "-Dynamic") == 0) ||
        (strcmp(argv[1], "-dynamic") == 0))

      incInertia = 1;

    else if ((strcmp(argv[1], "-rayleigh") == 0))

      incInertia = 2;
  }

  the_domain->calculateNodalReactions(incInertia);

  return TCL_OK;
}

int
nodeReaction(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << "WARNING want - nodeReaction nodeTag? <dof?>\n";
    return TCL_ERROR;
  }

  int tag;
  int dof = -1;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << "WARNING nodeReaction nodeTag? dof? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  if (argc > 2) {
    if (Tcl_GetInt(interp, argv[2], &dof) != TCL_OK) {
      opserr << "WARNING nodeReaction nodeTag? dof? - could not read dof? \n";
      return TCL_ERROR;
    }
  }

  dof--;

  const Vector *nodalResponse = domain->getNodeResponse(tag, NodeData::Reaction);

  if (nodalResponse == nullptr)
    // TODO: add error message
    return TCL_ERROR;

  int size = nodalResponse->Size();

  if (dof >= 0) {

    if (dof >= size) {
      opserr << "WARNING nodeReaction nodeTag? dof? - dofTag? too large\n";
      return TCL_ERROR;
    }

    Tcl_SetObjResult(interp, Tcl_NewDoubleObj((*nodalResponse)(dof)));

  } else {
    char buffer[40];
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", (*nodalResponse)(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
  }

  return TCL_OK;
}

int
nodeCoord(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT << "want - nodeCoord nodeTag? <dim?>\n";
    return TCL_ERROR;
  }

  int tag;

  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "nodeCoord nodeTag? dim? - could not read nodeTag? \n";
    return TCL_ERROR;
  }

  int dim = -1;

  if (argc > 2) {
    if (strcmp(argv[2], "X") == 0 || strcmp(argv[2], "x") == 0 ||
        strcmp(argv[2], "1") == 0)
      dim = 0;
    else if (strcmp(argv[2], "Y") == 0 || strcmp(argv[2], "y") == 0 ||
             strcmp(argv[2], "2") == 0)
      dim = 1;
    else if (strcmp(argv[2], "Z") == 0 || strcmp(argv[2], "z") == 0 ||
             strcmp(argv[2], "3") == 0)
      dim = 2;
    else {
      opserr << G3_ERROR_PROMPT << "" << "nodeCoord nodeTag? dim? - could not read dim? \n";
      return TCL_ERROR;
    }
  }

  Node *theNode = the_domain->getNode(tag);

  if (theNode == nullptr) {
    opserr << G3_ERROR_PROMPT << "Unable to retrieve node with tag '" << tag << "'\n";
    return TCL_ERROR;
  }

  const Vector &coords = theNode->getCrds();

  char buffer[40];
  int size = coords.Size();
  if (dim == -1) {
    for (int i = 0; i < size; ++i) {
      sprintf(buffer, "%35.20f", coords(i));
      Tcl_AppendResult(interp, buffer, NULL);
    }
    return TCL_OK;

  } else if (dim < size) {
    double value = coords(dim);
    Tcl_SetObjResult(interp, Tcl_NewDoubleObj(value));
    return TCL_OK;
  }

  return TCL_ERROR;
}

int
retainedNodes(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *domain = (Domain*)clientData;
  bool all = 1;
  int cNode;
  if (argc > 1) {
    if (Tcl_GetInt(interp, argv[1], &cNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "retainedNodes <cNode?> - could not read cNode? \n";
      return TCL_ERROR;
    }
    all = 0;
  }

  MP_Constraint *theMP;
  MP_ConstraintIter &mpIter = domain->getMPs();

  // get unique constrained nodes with set
  std::set<int> tags;
  int tag;
  while ((theMP = mpIter()) != nullptr) {
    tag = theMP->getNodeRetained();
    if (all || cNode == theMP->getNodeConstrained()) {
      tags.insert(tag);
    }
  }

  // assign set to vector and sort
  std::vector<int> tagv;
  tagv.assign(tags.begin(), tags.end());
  std::sort(tagv.begin(), tagv.end());

  // loop through unique, sorted tags, adding to output
  char buffer[20];
  for (int tag : tagv) {
    sprintf(buffer, "%d ", tag);
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}



int
nodeDOFs(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  if (argc != 2) {
    opserr << G3_ERROR_PROMPT << "expected - nodeDOFs nodeTag?\n";
    return TCL_ERROR;
  }

  int tag;
  if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "nodeDOFs nodeTag?\n";
    return TCL_ERROR;
  }


  Node *theNode = the_domain->getNode(tag);
  if (theNode == nullptr) {
    opserr << G3_ERROR_PROMPT << "nodeDOFs node " << tag << " not found" << endln;
    return TCL_ERROR;
  }

  int numDOF = theNode->getNumberDOF();
  DOF_Group *theDOFgroup = theNode->getDOF_GroupPtr();
  if (theDOFgroup == nullptr) {
    opserr << G3_ERROR_PROMPT << "nodeDOFs DOF group null" << endln;
    return -1;
  }

  char buffer[40];
  const ID &eqnNumbers = theDOFgroup->getID();
  for (int i = 0; i < numDOF; ++i) {
    sprintf(buffer, "%d ", eqnNumbers(i));
    Tcl_AppendResult(interp, buffer, NULL);
  }

  return TCL_OK;
}

