//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// NOTE: This doesnt really need access to the model builder, it would
// work with just the domain
//
#ifdef _TCL85
#  define TCL_Char const char
#elif _TCL84
#  define TCL_Char const char
#else
#  define TCL_Char char
#endif
#include <tcl.h>
#include <assert.h>
#include <Domain.h>
#include <Matrix.h>
#include <Node.h>
#include <NodeND.h>
#include <ID.h>
#include <G3_Logging.h>
#include <BasicModelBuilder.h>
#include "Block2D.h"
#include "Block3D.h"


int
TclCommand_doBlock2D(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{ 
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theTclDomain = builder->getDomain();
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm < 2) {
    opserr << G3_ERROR_PROMPT
           << "model dimension (ndm) must be at leat 2 " << "\n";
    return TCL_ERROR;
  }

  if (argc < 8) {
    opserr << G3_ERROR_PROMPT << "incorrect number of args, expected:"
      "\n\tblock2D numX? numY? startNode? startEle? eleType? eleArgs?"; 
    return TCL_ERROR;
  }
  int numX, numY, startNodeNum, startEleNum;
  if (Tcl_GetInt(interp, argv[1], &numX) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid numX (" << argv[1] << "), expected:\n\t"
      << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &numY) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &startNodeNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[3] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &startEleNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[4] << "\n";
    return TCL_ERROR;
  }


  static Matrix Coordinates(9,3);
  Coordinates.Zero();

  static ID     haveNode(9);
  for (int k=0; k<9; k++)
    haveNode(k) = -1;

  int numNodes = 4;
  if (argc == 10) {
    if (strcmp(argv[7],"-numEleNodes") == 0)
      if (Tcl_GetInt(interp, argv[8], &numNodes) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType eleArgs?";
        opserr << " -numEleNodes numNodes?: invalid numNodes: " << argv[8] << "\n"; 
        return TCL_ERROR;
      }
    if (numNodes != 4 && numNodes != 9) {
      opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
      opserr << "-numEleNodes numNodes?: invalid numNodes: " << argv[8] << " 4 or 9 only\n";
      return TCL_ERROR;
    }

    if (numNodes == 9) {
      if (((numX % 2) != 0) || ((numY % 2) != 0)) {
        opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs? ";
        opserr << "numX and numY must both be even when using -numEleNodes 9\n";
        return TCL_ERROR;
      }
    }
  }

  int nodalInfo = 9;
  if (numNodes == 4)
    nodalInfo = 7;

  TCL_Char ** argvNodes;
  int  argcNodes;

  Tcl_SplitList(interp, argv[nodalInfo], &argcNodes, &argvNodes);


  int count = 0;
  while (count < argcNodes) {
    int nodeTag;
    double value;
    if ((count + ndm + 1) >  argcNodes) {
      opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid number of node args: " << argv[7] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "block2D numX? numY? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node tag: " << argvNodes[count] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (nodeTag < 1 || nodeTag > 9) {
      opserr << G3_ERROR_PROMPT
             << " : invalid node tag out of bounds [1,9]: " << argvNodes[count] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; ++i) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
        opserr << G3_ERROR_PROMPT
               << " : invalid node coordinate for node: " << argvNodes[count] << "\n";
        Tcl_Free((char *)argvNodes); return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block2D  theBlock(numX, numY, haveNode, Coordinates, numNodes);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int jj;
  for (jj=0; jj<=numY; jj++) {
    for (int ii=0; ii<=numX; ii++) {

      Vector3D nodeCoords = theBlock.getNodalCoords(ii,jj);

      Node *theNode = nullptr;

      if (ndm == 2) {
        theNode = new Node(nodeID, ndf, nodeCoords(0), nodeCoords(1));

      } else if (ndm == 3) {
        if (getenv("NODE")) {
          switch (ndf) {
            case 3:
              theNode = new NodeND<3, 3>(nodeID, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
            case 6:
              theNode = new NodeND<3, 6>(nodeID, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
            default:
              theNode = new Node(nodeID, ndf, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
          }
        } else
          theNode = new Node(nodeID, ndf, nodeCoords(0), nodeCoords(1), nodeCoords(2));
        // theNode = new Node(nodeID,ndf, nodeCoords(0), nodeCoords(1), nodeCoords(2));
      }

      if (theTclDomain->addNode(theNode) == false) {
        opserr << G3_ERROR_PROMPT << "failed to add node to the domain\n";
        opserr << "node: " << nodeID << "\n";
        delete theNode;
        return TCL_ERROR;
      }
      nodeID++;
    }
  }

  // create the elements: numX*numY elements to be created if 4 node elements
  //                      numX/2 * numY /2 nodes to be v=created if 9 node elements 
  TCL_Char *eleType = argv[5]; 
  TCL_Char *additionalEleArgs = argv[6];
  //  const ID &nodeTags = theBlock.getElementNodes(0,0);
  //  int numNodes = nodeTags.Size();

  // assumes 15 is largest string for individual nodeTags
  count = int(10 + strlen(eleType) + strlen(additionalEleArgs) + 15 *(numNodes+1));
  char *eleCommand = new char[count];
  int initialCount = int(8 + strlen(eleType));

  int  eleID = startEleNum;
  if (numNodes == 9) {
    numX /= 2;
    numY /= 2;
  }

  for (int jj=0; jj<numY; jj++) {
    for (int ii=0; ii<numX; ii++) {
      count = initialCount;
      const ID &nodeTags = theBlock.getElementNodes(ii,jj);

      // create the string to be evaluated
      strcpy(eleCommand, "element ");
      strcpy(&eleCommand[8], eleType);
      count += sprintf(&eleCommand[count], " %d ", eleID);
      for (int i=0; i<numNodes; ++i) {
        int nodeTag = nodeTags(i)+startNodeNum;
        count += sprintf(&eleCommand[count], " %d ", nodeTag);
      }
      strcat(eleCommand, additionalEleArgs);

      // now to create the element we get the string evaluated
      if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
          delete [] eleCommand;
        return TCL_ERROR;
      }
      eleID++;
    }
  }
  delete [] eleCommand;
  return TCL_OK;
}


int
TclCommand_doBlock3D(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *theTclDomain = builder->getDomain();
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  if (ndm < 3) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? startNode? startEle? eleType? eleArgs?";
    opserr << " : model dimension (ndm) must be at leat 2 " << "\n";
    return TCL_ERROR;
  }

  int numX, numY, numZ, startNodeNum, startEleNum;
  if (Tcl_GetInt(interp, argv[1], &numX) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numX: " << argv[1] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &numY) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numY: " << argv[2] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[3], &numZ) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid numZ: " << argv[3] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[4], &startNodeNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startNode: " << argv[4] << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[5], &startEleNum) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
    opserr << " : invalid startEle: " << argv[5] << "\n";
    return TCL_ERROR;
  }

  static Matrix Coordinates(27,3);
  static ID     haveNode(27);
  Coordinates.Zero();
  for (int k=0; k<27; k++) haveNode(k) = -1;

  TCL_Char *nodalInfo = argv[8];
  TCL_Char ** argvNodes;
  int  argcNodes;

  Tcl_SplitList(interp, nodalInfo, &argcNodes, &argvNodes);


  int count = 0;
  while (count < argcNodes) {
    if ((count + ndm + 1) > argcNodes) {
      opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType eleArgs?";
      opserr << " : invalid number of node args: " << argv[8] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    int nodeTag;
    double value;
    if (Tcl_GetInt(interp, argvNodes[count], &nodeTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : invalid node id in node args: " << argvNodes[count] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    if (nodeTag < 1 || nodeTag > 27) {
      opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
      opserr << " : node tag out of bounds [1, 27]: " << argvNodes[count] << "\n";
      Tcl_Free((char *)argvNodes);
      return TCL_ERROR;
    }
    for (int i=0; i<ndm; ++i) {
      if (Tcl_GetDouble(interp, argvNodes[count+1+i], &value) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "block3D numX? numY? numZ? startNode? startEle? eleType? eleArgs?";
        opserr << " : invalid coordinate in node args: " << argvNodes[count] << "\n";
        Tcl_Free((char *)argvNodes);
        return TCL_ERROR;
      }
      Coordinates(nodeTag-1,i) = value;
      haveNode(nodeTag-1) = nodeTag;
    }
    count += 1 + ndm;
  }

  Tcl_Free((char *)argvNodes);

  Block3D  theBlock(numX, numY, numZ, haveNode, Coordinates);

  // create the nodes: (numX+1)*(numY+1) nodes to be created
  int nodeID = startNodeNum;
  int kk;
  for (kk=0; kk<=numZ; kk++) {
    for (int jj=0; jj<=numY; jj++) {
      for (int ii=0; ii<=numX; ii++) {

        const Vector &nodeCoords = theBlock.getNodalCoords(ii,jj,kk);
        //
        Node* theNode = nullptr;
        if (getenv("NODE")) {
          switch (ndf) {
            case 3:
              theNode = new NodeND<3, 3>(nodeID, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
            case 6:
              theNode = new NodeND<3, 6>(nodeID, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
            default:
              theNode = new Node(nodeID, ndf, nodeCoords(0), nodeCoords(1), nodeCoords(2));
              break;
          }
        } else
          theNode = new Node(nodeID, ndf, nodeCoords(0), nodeCoords(1), nodeCoords(2));

        // theNode = new Node(nodeID,ndf,nodeCoords(0),nodeCoords(1),nodeCoords(2));

        if (theTclDomain->addNode(theNode) == false) {
          opserr << G3_ERROR_PROMPT << "failed to add node to the domain\n";
          opserr << "node: " << nodeID << "\n";
          delete theNode;
          return TCL_ERROR;
        }

        nodeID++;
      }
    }
  }

  // create the elements: numX*numY elements to be created
  TCL_Char *eleType = argv[6];
  TCL_Char *additionalEleArgs = argv[7];
  const ID &nodeTags = theBlock.getElementNodes(0,0,0);
  int numNodes = nodeTags.Size();

  // TODO: assumes 15 is largest string for individual nodeTags
  count = int(10 + strlen(eleType) + strlen(additionalEleArgs) + 15 * (numNodes+1));
  char *eleCommand = new char[count];
  int initialCount = int(8 + strlen(eleType));

  int  eleID = startEleNum;
  for (kk=0; kk<numZ; kk++) {
    for (int jj=0; jj<numY; jj++) {
      for (int ii=0; ii<numX; ii++) {
        count = initialCount;

        const ID &nodeTags = theBlock.getElementNodes(ii,jj,kk);

        // create the string to be evaluated
        strcpy(eleCommand, "element ");
        strcpy(&eleCommand[8], eleType);
        count += sprintf(&eleCommand[count], " %d ", eleID);
        for (int i=0; i<numNodes; ++i) {
          int nodeTag = nodeTags(i)+startNodeNum;
          count += sprintf(&eleCommand[count], " %d ", nodeTag);
        }
        strcat(eleCommand, additionalEleArgs);

        // now to create the element we get the string evaluated
        if (Tcl_Eval(interp, eleCommand) != TCL_OK) {
          delete [] eleCommand;
          return TCL_ERROR;
        }
        eleID++;
      }
    }
  }

  delete [] eleCommand;
  return TCL_OK;
}

int
TclCommand_doBlock(ClientData clientData, Tcl_Interp *interp, int argc,
                          TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  int ndm = builder->getNDM();

  if (argc < 1) {
    opserr << G3_ERROR_PROMPT << "block <type> {args}\n";
    return TCL_ERROR;
  }

  if (ndm == 2)
    return TclCommand_doBlock2D(clientData, interp, argc, argv);

  else if (strcmp(argv[1], "2d") == 0)
    return TclCommand_doBlock2D(clientData, interp, argc-1, argv+1);

  else
    return TclCommand_doBlock3D(clientData, interp, argc, argv);

  return TCL_OK;
}
