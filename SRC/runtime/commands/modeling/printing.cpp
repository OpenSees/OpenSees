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
// Description: Commands that are used to print out the domain
//
// Author: cmp
//
#ifdef _WIN32
#  include <io.h>
#  define isatty _isatty
#  define STDOUT_FILENO _fileno(stdout)
#else
#  include <unistd.h>               
#endif
#include <assert.h>
#include <tcl.h>
#include <Logging.h>
#include <FileStream.h>
#include <DummyStream.h>

#include <BasicModelBuilder.h>

#include <ID.h>
#include <Domain.h>
#include <LoadPattern.h>
#include <Parameter.h>
#include <SP_Constraint.h>
#include <SP_ConstraintIter.h>
#include <MP_Constraint.h>
#include <MP_ConstraintIter.h>

#include <Parameter.h>
#include <ParameterIter.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <SectionForceDeformation.h>
#include <FrameSection.h>

#include <Pressure_Constraint.h>
#include <Element.h>
#ifdef OPS_USE_DAMPING
#include <damping/Damping.h>
#endif
#include <ElementIter.h>

#include <Node.h>
#include <NodeIter.h>

#include <LoadPattern.h>
#include <LoadPatternIter.h>

#include <CrdTransf.h>

int printElement(ClientData clientData, Tcl_Interp *interp, int argc,
                 TCL_Char ** const argv, OPS_Stream &output);

int printNode(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv, OPS_Stream &output);

int printIntegrator(ClientData clientData, Tcl_Interp *interp, int argc,
                    TCL_Char ** const argv, OPS_Stream &output);

int printAlgorithm(ClientData clientData, Tcl_Interp *interp, int argc,
                   TCL_Char ** const argv, OPS_Stream &output);


int
TclCommand_classType(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char** const argv)
{

  assert(clientData != nullptr);
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  if (argc < 3) {
    opserr << "ERROR want - classType objectType tag?\n";
    return TCL_ERROR;
  }

  std::string type = argv[1];

  MovableObject* theObject = nullptr;
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) < 0) {
    opserr << OpenSees::PromptValueError << "classType objectType tag? - unable to read tag" << "\n";
    return TCL_ERROR;
  }

  if (type == "uniaxialMaterial")
    theObject = builder->getTypedObject<UniaxialMaterial>(tag);

  else if (type == "section")
    theObject = builder->getTypedObject<SectionForceDeformation>(tag);
#ifdef OPS_USE_DAMPING
  else if (type == "damping")
    theObject = builder->getTypedObject<Damping>(tag);
#endif
  else {
    opserr << OpenSees::PromptValueError << "classType - " << type.c_str() << " not yet supported" << "\n";
    return TCL_ERROR;
  }

  std::string classType = theObject->getClassType();
  
  Tcl_SetObjResult(interp, Tcl_NewStringObj(classType.c_str(), strlen(classType.c_str())));

  return TCL_OK;
}

template <typename T>
static int
printRegistryObject(const BasicModelBuilder& builder, int tag, int flag, OPS_Stream *output)
{
  TaggedObject* object = builder.getTypedObject<T>(tag);
  object->Print(*output, flag);
  return TCL_OK;
}

static int
printRegistry(const BasicModelBuilder& builder, TCL_Char* type, int flag, OPS_Stream *output)
{
  if (type == nullptr)
    builder.printRegistry<BasicModelBuilder>(*output, flag);
  return TCL_OK;
}


static void
printDomain(OPS_Stream &s, BasicModelBuilder* builder, int flag) 
{

  Domain* theDomain = builder->getDomain();

  const char* tab = "  ";
  // TODO: maybe add a method called countRegistry<>
  // to BasicModelBuilder

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "{\n";
    s << "\"StructuralAnalysisModel\": {\n";

    s << tab << "\"properties\": {\n";
    //
    {
      s << tab << tab << "\"sections\": [\n";        
      int n = builder->printRegistry<SectionForceDeformation>(s, flag);

      DummyStream dummy;
      if (builder->printRegistry<FrameSection>(dummy, flag) > 0) {
        if (n > 0)
          s << ",\n";
        builder->printRegistry<FrameSection>(s, flag);
      }
      s << "\n" << tab << tab << "]";
    }
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"nDMaterials\": [\n";        
      builder->printRegistry<NDMaterial>(s, flag);
      s << "\n" << tab << tab << "]";
    }
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"uniaxialMaterials\": [\n";        
      builder->printRegistry<UniaxialMaterial>(s, flag);
      s << "\n" << tab << tab << "]";
    }
    s << ",\n";
    //
    s << tab << tab << "\"crdTransformations\": [\n";
    {
      int n = builder->printRegistry<CrdTransf>(s, flag);

      DummyStream dummy;
      if (builder->printRegistry<CrdTransf>(dummy, flag) > 0) {
        if (n > 0)
          s << ",\n";
        builder->printRegistry<CrdTransf>(s, flag);
      }
      s << "\n" << tab << tab << "]";
    }
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"patterns\": [\n";
      LoadPatternIter &patterns = theDomain->getLoadPatterns();
      LoadPattern *p;
      bool first_mp = true;
      while ((p = patterns()) != nullptr) {
        if (!first_mp)
          s << ",\n";

        p->Print(s, flag);
        first_mp = false;
      }
      s << "\n" << tab << tab << "]";
    }
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"parameters\": [\n";
      ParameterIter &params = theDomain->getParameters();
      Parameter *param;
      bool first_mp = true;
      while ((param = params()) != nullptr) {
        if (!first_mp)
          s << ",\n";

        param->Print(s, flag);
        first_mp = false;
      }
      s << "\n" << tab << tab << "]\n";
    }
    //
    s << "\n";
    //
    //
    s << tab << "},\n";
    //
    //
    s << tab << "\"geometry\": {\n";
    int numPrinted = 0;
    int numToPrint = theDomain->getNumNodes();
    NodeIter &theNodess = theDomain->getNodes();
    Node *theNode;
    s << tab << tab << "\"nodes\": [\n";
    while ((theNode = theNodess()) != nullptr) {
      theNode->Print(s, flag);
      numPrinted += 1;
      if (numPrinted < numToPrint)
        s << ",\n";
    }
    s << "\n" << tab << tab << "]";
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"elements\": [\n";
      Element *theEle;
      ElementIter &theElementss = theDomain->getElements();
      numToPrint = theDomain->getNumElements();
      numPrinted = 0;
      while ((theEle = theElementss()) != nullptr) {
        theEle->Print(s, flag);
        numPrinted += 1;
        if (numPrinted < numToPrint)
          s << ",\n";
      }
      s << "\n" << tab << tab << "]";
    }
    //
    s << ",\n";
    //
    {
      s << tab << tab << "\"constraints\": [\n";
      MP_ConstraintIter &theMPs = theDomain->getMPs();
      MP_Constraint *theMP;
      bool first_mp = true;
      while ((theMP = theMPs()) != nullptr) {
        if (!first_mp)
          s << ",\n";
        theMP->Print(s, flag);
        first_mp = false;
      }

      SP_ConstraintIter &theSPs = theDomain->getSPs();
      SP_Constraint *theSP;
      bool first_sp = true;
      while ((theSP = theSPs()) != nullptr) {
        if (!first_sp || !first_mp)
          s << ",\n";
        theSP->Print(s, flag);
        first_sp = false;
      }
      s << "\n" << tab << tab << "]";
    }

    // END
    s << "\n";

    s << tab << "}\n";
    s << "}\n";
    s << "}\n";

    return;
  }
}

int
TclCommand_print(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);

  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  Domain * domain = builder->getDomain();

  int currentArg = 1;
  int res = TCL_OK;

  int flag = OPS_PRINT_CURRENTSTATE;

  FileStream outputFile;
  OPS_Stream *output = nullptr;
  bool close_file = false;

  // if called interactively, print to console
  if (isatty(STDOUT_FILENO)) {
    output = &opserr;
  } else {
#ifdef _WIN32
    outputFile.setFile("CON", openMode::APPEND);
#else
    outputFile.setFile("/dev/stdout", openMode::APPEND);
#endif
    output = &outputFile;
  }

  bool done = false;

  while (currentArg < argc) {

    // if 'print ele i j k..' print out some elements
    if ((strcmp(argv[currentArg], "-ele") == 0) ||
        (strcmp(argv[currentArg], "-element") == 0)  ||
        (strcmp(argv[currentArg], "ele") == 0)) {
      currentArg++;
      res = printElement((ClientData)domain, interp, argc - currentArg, argv + currentArg, *output);
      done = true;
    }

    // if 'print node i j k ..' print out some nodes
    else if ((strcmp(argv[currentArg], "-node") == 0) ||
             (strcmp(argv[currentArg], "node") == 0)) {
      currentArg++;
      res = printNode((ClientData)domain, interp, argc - currentArg, argv + currentArg,
                      *output);
      done = true;
    }

    // if 'print material i j k ..' print out some nodes
    else if ((strcmp(argv[currentArg], "-material") == 0)) {
      currentArg++;
      if (currentArg == argc) {
        opserr << OpenSees::PromptValueError << "print -material <tag> .. - no tag specified\n";
        return TCL_ERROR;
      }
      for (int i = currentArg; i < argc; i++) {
        int tag;
        if (Tcl_GetInt(interp, argv[i], &tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "print -material failed to get integer tag: " << argv[i]
                 << "\n";
          return TCL_ERROR;
        }
        res += printRegistryObject<NDMaterial>(*((BasicModelBuilder*)clientData), tag, OPS_PRINT_PRINTMODEL_JSON, output);
      }
      done = true;
    }


    else if ((strcmp(argv[currentArg], "-registry") == 0)) {
      currentArg++;
      if (currentArg == argc)
        res = printRegistry(*((BasicModelBuilder*)clientData), nullptr, flag, output);
      else
        res = printRegistry(*((BasicModelBuilder*)clientData), argv[currentArg++], flag, output);
      done = true;
    }

    // if 'print integrator flag' print out the integrator
    else if ((strcmp(argv[currentArg], "integrator") == 0) ||
             (strcmp(argv[currentArg], "-integrator") == 0)) {
      currentArg++;
      res = printIntegrator((ClientData)domain, interp, argc - currentArg,
                            argv + currentArg, *output);
      done = true;
    }

    // if 'print algorithm flag' print out the algorithm
    else if ((strcmp(argv[currentArg], "algorithm") == 0) ||
             (strcmp(argv[currentArg], "-algorithm") == 0)) {
      currentArg++;

      Tcl_CmdInfo info;
      if (Tcl_GetCommandInfo(interp, "analyze", &info)==1) {
        res = printAlgorithm(info.clientData, interp, argc - currentArg,
                             argv + currentArg, *output);
      } else {
        opserr << OpenSees::PromptValueError << "Cannot print algorithm\n";
      }
      done = true;
    }

    else if ((strcmp(argv[currentArg], "-JSON") == 0) ||
             (strcmp(argv[currentArg], "-json") == 0)) {
      currentArg++;
      flag = OPS_PRINT_PRINTMODEL_JSON;
    }

    else {
      if ((strcmp(argv[currentArg], "file") == 0) ||
          (strcmp(argv[currentArg], "-file") == 0))
        currentArg++;

      openMode mode = openMode::APPEND;
      if (flag == OPS_PRINT_PRINTMODEL_JSON)
        mode = openMode::OVERWRITE;

      if (currentArg < argc) {
        if (outputFile.setFile(argv[currentArg], mode) != 0) {
          opserr << "print <filename> .. - failed to open file: "
                 << argv[currentArg] << "\n";
          return TCL_ERROR;
        } else {
          close_file = true;
          output = &outputFile;
        }
      }

      currentArg++;
    }
  }

  // if just 'print <filename>' then print out the entire domain to eof
  if (!done) {
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
      // simulationInfo.Print(*output, flag);
      printDomain(*output, builder, flag);
    } else {
      domain->Print(*output, flag);
      // Domain doesnt leave a new line
      *output << "\n";
    }

    res = TCL_OK;
    done = true;
  }

  // close the output file if one has been opened
  if (close_file)
    outputFile.close();

  return res;
}

int
printElement(ClientData clientData, Tcl_Interp *interp, int argc,
             TCL_Char ** const argv, OPS_Stream &output)
{
  assert(clientData != nullptr);
  Domain * the_domain = (Domain*)clientData;

  int flag   = 0; // default flag sent to a nodes Print() method
  int eleArg = 0;

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> Element flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) ||
      (strcmp(argv[0], "-flag")) == 0) { // get the specified flag
    if (argc < 2) {
      opserr << OpenSees::PromptValueError << "print <filename> ele <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "print ele failed to get integer flag: \n";
      opserr << argv[eleArg] << "\n";
      return TCL_ERROR;
    }
    eleArg += 2;
  }

  // now print the Elements with the specified flag, 0 by default
  if (argc == eleArg) {
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0)
      theElement->Print(output, flag);
    return TCL_OK;

  } else {

    // otherwise print out the specified elements i j k .. with flag
    int numEle = argc - eleArg;
    ID *theEle = new ID(numEle);
    for (int i = 0; i < numEle; ++i) {
      int eleTag;
      if (Tcl_GetInt(interp, argv[i + eleArg], &eleTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "print -ele failed to get integer: " << argv[i]
               << "\n";
        return TCL_ERROR;
      }
      (*theEle)(i) = eleTag;
    }

    the_domain->Print(output, 0, theEle, flag);
    delete theEle;
  }

  return TCL_OK;
}

// function to print out the nodal information conatined in line
//     print <filename> node <flag int> <int int int>
// Parameters
//   nodeArg: integer equal to arg count to node plus 1
//   output:  output stream to which the results are sent
//
int
printNode(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv,
          OPS_Stream &output)
{
  int flag = 0; // default flag sent to a nodes Print() method
  int nodeArg = 0;

  assert(clientData != nullptr);
  Domain* domain = (Domain*)clientData; 

  // if just 'print <filename> node' print all the nodes - no flag
  if (argc == 0) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != nullptr)
      theNode->Print(output);
    return TCL_OK;
  }

  // if 'print <filename> node flag int <int int ..>' get the flag
  if ((strcmp(argv[0], "flag") == 0) || (strcmp(argv[0], "-flag") == 0)) {
    // get the specified flag
    if (argc <= nodeArg) {
      opserr << OpenSees::PromptValueError << "print <filename> node <flag int> no int specified \n";
      return TCL_ERROR;
    }
    if (Tcl_GetInt(interp, argv[1], &flag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "print node failed to get integer flag: \n";
      opserr << argv[nodeArg] << "\n";
      return TCL_ERROR;
    }
    nodeArg += 2;
  }

  // now print the nodes with the specified flag, 0 by default

  // if 'print <filename> node flag'
  //     print out all the nodes in the domain with flag
  if (nodeArg == argc) {
    NodeIter &theNodes = domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != nullptr)
      theNode->Print(output, flag);
    return TCL_OK;
  } else {
    // otherwise print out the specified nodes i j k .. with flag
    int numNodes = argc - nodeArg;
    ID *theNodes = new ID(numNodes);
    for (int i = 0; i < numNodes; ++i) {
      int nodeTag;
      if (Tcl_GetInt(interp, argv[nodeArg], &nodeTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "print node failed to get integer: " << argv[nodeArg]
               << "\n";
        return TCL_ERROR;
      }
      (*theNodes)(i) = nodeTag;
      nodeArg++;
    }

    domain->Print(output, theNodes, 0, flag);
    delete theNodes;
  }

  return TCL_OK;
}

// Print domain in GiD format
// Author: Talledo
int
printModelGID(ClientData clientData, Tcl_Interp *interp, int argc,
              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  Domain *the_domain = (Domain*)clientData;

  enum ObjectTypes {
    Linear = 1<<0,
    Tri3   = 1<<1,
    Quad4  = 1<<2,
    Quad8  = 1<<3,
    Quad9  = 1<<4,
    Brick  = 1<<5
  };

  // This function print's a file with node and elements in a format useful for
  // GID
  int res = 0;
  bool hasLinear = false;
  bool hasTri3  = false;
  bool hasQuad4 = false;
//bool hasQuad8 = false;
  bool hasQuad9 = false;
  bool hasBrick = false;
  int startEle = 1;
  int endEle = 1;
  int eleRange = 0;
  int i = 2;

  FileStream outputFile;

  if (argc < 2) {
    opserr << OpenSees::PromptValueError << "printGID fileName? - no filename supplied\n";
    return TCL_ERROR;
  }
  openMode mode = openMode::OVERWRITE;
  if (argc >= 3) {
    if (strcmp(argv[i], "-append") == 0) {
      mode = openMode::APPEND;
      i++;
    }
    if (strcmp(argv[i], "-eleRange") == 0) {

      eleRange = 1;
      if (Tcl_GetInt(interp, argv[i + 1], &startEle) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "print node failed to get integer: " << argv[i + 1]
               << "\n";
        return TCL_ERROR;
      }
      if (Tcl_GetInt(interp, argv[i + 2], &endEle) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "print node failed to get integer: " << argv[i + 2]
               << "\n";
        return TCL_ERROR;
      }

    }
  }

  if (outputFile.setFile(argv[1], mode) < 0) {
    opserr << OpenSees::PromptValueError << "printGID " << argv[1] << " failed to set the file\n";
    return TCL_ERROR;
  }

  // Cycle over Elements to understand what type of elements are there
  ElementIter &theElements = the_domain->getElements();
  Element *theElement;
  while ((theElement = theElements()) != 0) {

    // Check type of Element with Number of Nodes
    // if 2 Nodes print the Element
    switch (theElement->getNumExternalNodes()) {
    case (2):
      hasLinear = true;
      break;
    case (4):
      hasQuad4 = true;
      break;
    case (3):
      hasTri3 = true;
      break;
    case (9):
      hasQuad9 = true;
      break;
    case (8):
      if (strcmp(theElement->getClassType(), "Brick") == 0) {
        hasBrick = true;
      } else {
        ;// hasQuad8 = true;
      }
    }
  }

  //
  // **** Linear Elements - 2 Nodes
  //
  if (hasLinear == 1) {
    // Print HEADER
    outputFile << "MESH \"2NMESH\" dimension 3 ElemType Linear Nnode 2"
               << "\n";
    outputFile << "#color 0 0 255" << "\n\n";

    // Print node coordinates
    outputFile << "Coordinates" << "\n";
    NodeIter &theNodes = the_domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << "\n";
    }
    outputFile << "End coordinates" << endln << "\n";

    // Print elements connectivity
    outputFile << "Elements" << "\n";
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {
        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 2) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; ++i) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; ++i) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << "\n";
        }
      }
    }
    outputFile << "End elements" << "\n";
  }
  //
  // **** Quadrilateral Elements - 4 Nodes
  //
  if (hasQuad4 == 1) {
    // Print HEADER
    outputFile << "MESH \"4NMESH\" dimension 3 ElemType Quadrilateral Nnode 4"
               << "\n";
    outputFile << "#color 0 255 0" << endln << "\n";

    // Print node coordinates
    outputFile << "Coordinates" << "\n";
    NodeIter &theNodes = the_domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << "\n";
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << "\n";
    }
    outputFile << "End coordinates" << endln << "\n";

    // Print elements connectivity
    outputFile << "Elements" << "\n";
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 4) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; ++i) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; ++i) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << "\n";
        }
      }
    }
    outputFile << "End elements" << "\n";
  }
  //
  // **** Triangular Elements - 3 Nodes
  //
  if (hasTri3 == 1) {
    // Print HEADER
    outputFile << "MESH \"3NMESH\" dimension 3 ElemType Triangle Nnode 3"
               << "\n";
    outputFile << "#color 0 255 0" << endln << "\n";

    // Print node coordinates
    outputFile << "Coordinates" << "\n";
    NodeIter &theNodes = the_domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << "\n";
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << "\n";
    }
    outputFile << "End coordinates" << endln << "\n";

    // Print elements connectivity
    outputFile << "Elements" << "\n";
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 3 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 3) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; ++i) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; ++i) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << "\n";
        }
      }
    }
    outputFile << "End elements" << "\n";
  }
  //
  // **** Quadrilateral Elements - 9 Nodes
  //
  if (hasQuad9 == 1) {
    // Print HEADER
    outputFile << "MESH \"9NMESH\" dimension 3 ElemType Linear Nnode 9"
               << "\n";
    outputFile << "#color 0 255 0" << endln << "\n";

    // Print node coordinates
    outputFile << "Coordinates" << "\n";
    NodeIter &theNodes = the_domain->getNodes();
    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();

      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << "\n";
    }
    outputFile << "End coordinates" << endln << "\n";

    // Print elements connectivity
    outputFile << "Elements" << "\n";
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 9) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; ++i) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; ++i) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << "\n";
        }
      }
    }
    outputFile << "End elements" << "\n";
  }
  //
  // **** Hexahedra Elements - 8 Nodes
  //
  if (hasBrick == 1) {
    // Print HEADER
    outputFile << "MESH \"8NMESH\" dimension 3 ElemType Hexahedra Nnode 8"
               << "\n";
    outputFile << "#color 255 0 0" << endln << "\n";

    // Print node coordinates
    outputFile << "Coordinates" << "\n";
    NodeIter &theNodes = the_domain->getNodes();

    Node *theNode;
    while ((theNode = theNodes()) != 0) {
      int tag = theNode->getTag();
      const Vector &crds = theNode->getCrds();
      // outputFile << tag << "\t\t" << crds(0) << "\t" << crds(1) << "\t" <<
      // crds(2) << "\n";
      int l_tmp = crds.Size();
      outputFile << tag << "\t\t";
      for (int ii = 0; ii < l_tmp; ii++) {
        outputFile << crds(ii) << "\t";
      }
      for (int ii = l_tmp; ii < 3; ii++) {
        outputFile << 0.0 << "\t";
      }
      outputFile << "\n";
    }
    outputFile << "End coordinates" << endln << "\n";

    // Print elements connectivity
    outputFile << "Elements" << "\n";
    ElementIter &theElements = the_domain->getElements();
    Element *theElement;
    while ((theElement = theElements()) != 0) {
      int tag = theElement->getTag();
      // Check if element tag is inside theRange
      if (((tag <= endEle) & (tag >= startEle)) || (eleRange == 0)) {

        // Check type of Element with Number of Nodes
        // if 2 Nodes print the Element
        int nNode = theElement->getNumExternalNodes();
        if (nNode == 8) {
          Node **NodePtrs;
          NodePtrs = theElement->getNodePtrs();
          ID tagNodes(nNode);
          for (int i = 0; i < nNode; ++i) {
            tagNodes(i) = NodePtrs[i]->getTag();
          }
          outputFile << tag << "\t\t";
          for (int i = 0; i < nNode; ++i) {
            outputFile << tagNodes(i) << "\t";
          }
          outputFile << "\n";
        }
      }
    }
    outputFile << "End elements" << "\n";
  }

  outputFile.close();
  return res;
}
