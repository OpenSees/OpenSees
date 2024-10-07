//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
// 
// Description: This file implements commands that configure Node objects
// for an analysis.
//
// Author: cmp
//
#include <assert.h>
#include <string.h>
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Node.h>
#include <NodeND.h>
#include <Matrix.h>
#include <Domain.h>
#include <BasicModelBuilder.h>

#define G3_MAX_NUM_DOFS 1000000000000
#define G3_NUM_DOF_BUFFER 20

int
TclCommand_addNode(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);

  Domain *theTclDomain = builder->getDomain();

  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  // make sure corect number of arguments on command line
  if (argc < 2 + ndm) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments, expected:\n";
    opserr << "      node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  Node *theNode = 0;

  // read the node id
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid nodeTag\n";
    opserr << "        Want: node nodeTag? [ndm coordinates?] <-mass [ndf values?]>\n";
    return TCL_ERROR;
  }

  // read in the coordinates and create the node
  double xLoc, yLoc, zLoc;
  if (ndm == 1) {
    // create a node in 1d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid coordinate\n";
      return TCL_ERROR;
    }
  }

  else if (ndm == 2) {
    // create a node in 2d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid 1st coordinate\n";
      opserr << "node: " << nodeId << "\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid 2nd coordinate\n";
      opserr << "node: " << nodeId << "\n";
      return TCL_ERROR;
    }
  }

  else if (ndm == 3) {
    // create a node in 3d space
    if (Tcl_GetDouble(interp, argv[2], &xLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid 1st coordinate\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[3], &yLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid 2nd coordinate\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[4], &zLoc) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid 3rd coordinate\n";
      return TCL_ERROR;
    }

  } else {
    opserr << G3_ERROR_PROMPT << "unsupported model dimension\n";
    return TCL_ERROR;
  }

  // check for -ndf override option
  int currentArg = 2 + ndm;
  if (currentArg < argc && strcmp(argv[currentArg], "-ndf") == 0) {
    if (Tcl_GetInt(interp, argv[currentArg + 1], &ndf) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nodal ndf given for node " << nodeId << "\n";
      return TCL_ERROR;
    }
    currentArg += 2;
  }

  //
  // create the node
  //
  switch (ndm) {
  case 1:
    theNode = new Node(nodeId, ndf, xLoc);
    break;
  case 2:
    theNode = new Node(nodeId, ndf, xLoc, yLoc);
    break;
  case 3:
    if (getenv("NODE")) {
      switch (ndf) {
        case 3:
          theNode = new NodeND<3, 3>(nodeId, xLoc, yLoc, zLoc);
          break;
        case 6:
          theNode = new NodeND<3, 6>(nodeId, xLoc, yLoc, zLoc);
          break;
        default:
          theNode = new Node(nodeId, ndf, xLoc, yLoc, zLoc);
          break;
      }
    } else
      theNode = new Node(nodeId, ndf, xLoc, yLoc, zLoc);
    break;
  }

  while (currentArg < argc) {
    if (strcmp(argv[currentArg], "-mass") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << G3_ERROR_PROMPT << "incorrect number of nodal mass terms\n";
        opserr << "node: " << nodeId << "\n";
        return TCL_ERROR;
      }

      double theMass;
      Matrix mass(ndf, ndf);
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theMass) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid nodal mass term";
          opserr << " at dof " << i + 1 << "\n";
          return TCL_ERROR;
        }
        mass(i, i) = theMass;
      }
      theNode->setMass(mass);

    } else if (strcmp(argv[currentArg], "-dispLoc") == 0) {
      currentArg++;
      if (argc < currentArg + ndm) {
        opserr << G3_ERROR_PROMPT << "incorrect number of nodal display location terms, "
                  "need ndm\n";
        return TCL_ERROR;
      }
      Vector displayLoc(ndm);
      double theCrd;
      for (int i = 0; i < ndm; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theCrd) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << "\n";
          return TCL_ERROR;
        }
        displayLoc(i) = theCrd;
      }
      theNode->setDisplayCrds(displayLoc);

    } else if (strcmp(argv[currentArg], "-disp") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << G3_ERROR_PROMPT << "incorrect number of nodal disp terms\n";
        opserr << "node: " << nodeId << "\n";
        return TCL_ERROR;
      }
      Vector disp(ndf);
      double theDisp;
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid nodal disp term\n";
          opserr << "node: " << nodeId << ", dof: " << i + 1 << "\n";
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialDisp(disp);
      theNode->commitState();

    } else if (strcmp(argv[currentArg], "-vel") == 0) {
      currentArg++;
      if (argc < currentArg + ndf) {
        opserr << G3_ERROR_PROMPT << "incorrect number of nodal vel terms, ";
        opserr << "expected " << ndf << "\n";
        return TCL_ERROR;
      }

      double theDisp;
      Vector disp(ndf);
      for (int i = 0; i < ndf; ++i) {
        if (Tcl_GetDouble(interp, argv[currentArg++], &theDisp) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid nodal vel term at ";
          opserr << " dof " << i + 1 << "\n";
          return TCL_ERROR;
        }
        disp(i) = theDisp;
      }
      theNode->setTrialVel(disp);
      theNode->commitState();

    } else
      currentArg++;
  }

  //
  // add the node to the domain
  //
  if (theTclDomain->addNode(theNode) == false) {
    opserr << G3_ERROR_PROMPT << "failed to add node to the domain\n";
    delete theNode;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclCommand_addNodalMass(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theTclDomain = builder->getDomain();

  int ndf = argc - 2;

  // make sure at least one other argument
  if (argc < (1 + ndf)) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments, expected:\n"
              "      mass nodeId <" << ndf << " mass values>\n"; 
    return TCL_ERROR;
  }

  // get the id of the node
  int nodeId;
  if (Tcl_GetInt(interp, argv[1], &nodeId) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid nodeId: " << argv[1];
    opserr << " - mass nodeId " << ndf << " forces\n";
    return TCL_ERROR;
  }

  // check for mass terms
  Matrix mass(ndf,ndf);
  for (int i=0; i<ndf; ++i) {
     double theMass;
     if (Tcl_GetDouble(interp, argv[i+2], &theMass) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid nodal mass term\n";
          opserr << "node: " << nodeId << ", dof: " << i+1 << "\n";
          return TCL_ERROR;
      }
      mass(i,i) = theMass;
  }

  if (theTclDomain->setMass(mass, nodeId) != 0) {
    opserr << G3_ERROR_PROMPT << "failed to set mass at node " << nodeId << "\n";
    return TCL_ERROR;
  }

  // if get here we have sucessfully created the node and added it to the domain
  return TCL_OK;
}



int
TclCommand_getNDM(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *the_domain = builder->getDomain();

  int ndm;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "ndm nodeTag? \n";
      return TCL_ERROR;
    }

    Node *theNode = the_domain->getNode(tag);
    if (theNode == nullptr) {
      opserr << G3_ERROR_PROMPT << "nodeTag " << tag << " does not exist \n";
      return TCL_ERROR;
    }
    const Vector &coords = theNode->getCrds();
    ndm = coords.Size();

  } else {
      ndm = builder->getNDM();
  }

  Tcl_SetObjResult(interp, Tcl_NewIntObj(ndm));
  return TCL_OK;
}

int
TclCommand_getNDF(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *the_domain = builder->getDomain();
  int ndf;

  if (argc > 1) {
    int tag;
    if (Tcl_GetInt(interp, argv[1], &tag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "ndf nodeTag? \n";
      return TCL_ERROR;
    }
    Node *theNode = the_domain->getNode(tag);
    if (theNode == nullptr) {
      opserr << G3_ERROR_PROMPT << "nodeTag " << tag << " does not exist \n";
      return TCL_ERROR;
    }
    ndf = theNode->getNumberDOF();

  } else {
    ndf = builder->getNDF();
  }

  char buffer[G3_NUM_DOF_BUFFER];
  sprintf(buffer, "%d", ndf);

  Tcl_AppendResult(interp, buffer, NULL);

  return TCL_OK;
}
