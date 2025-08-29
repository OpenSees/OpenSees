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
#include <tcl.h>
#include <set>
#include <assert.h>
#include <unordered_map>
#include <Logging.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <BasicModelBuilder.h>
#include <element/Shell/ASDShellQ4.h>
#include <element/Shell/ShellANDeS.h>
#include <element/Shell/ShellDKGQ.h>
#include <element/Shell/ShellDKGT.h>
#include <element/Shell/ShellMITC4.h>
#include <element/Shell/ShellMITC9.h>
#include <element/Shell/ShellNLDKGQ.h>
#include <element/Shell/ShellMITC4Thermal.h>
#include <element/Shell/ShellNLDKGQThermal.h>
#include <element/Shell/ShellNLDKGT.h>
using namespace OpenSees;

#include <algorithm>
#include <string>
  #ifdef _MSC_VER 
  #  include <string.h>
  #  define strcasecmp _stricmp
  #else
  #  include <strings.h>
  #endif

static
std::string
toLower( const std::string & s )
{
  std::string copy = s;
  transform( copy.begin( ), copy.end( ), copy.begin( ), 
      [](unsigned char c) { return std::tolower(c); });
  return copy;
}

static bool 
equalsIgnoreCase( const std::string & lhs, const std::string & rhs )
{
  return toLower( lhs ) == toLower( rhs );
}

class CaseInsensitive
{
  public:
    size_t operator( ) ( const std::string & s ) const
    {  
        static std::hash<std::string> hf;
        return hf( toLower( s ) );
    }
    
    bool operator( ) ( const std::string & lhs, const std::string & rhs ) const
    {
        return equalsIgnoreCase( lhs, rhs );
    }
};

using namespace OpenSees;

static std::unordered_map<std::string, int, CaseInsensitive, CaseInsensitive> 
NodeCounts = {
  {"ASDShellQ4",         4},
  {"ShellMITC4",         4},
  {"ShellMITC9",         9},
  {"ShellDKGQ",          4},
  {"ShellNLDKGQ",        4},
  {"ShellDKGT",          3},
  // {"ASDShellT3",         3}, // TODO
  {"ShellNLDKGT",        3},
  {"ShellANDeS",         4},
  {"ShellMITC4Thermal",  4},
  {"ShellNLDKGQThermal", 4},
};

int
TclBasicBuilder_addShell(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  int tag;
  if (argc < 4) {
    opserr << OpenSees::PromptValueError 
           << "insufficient arguments for element " << argv[1] 
           << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError 
           << "invalid element tag " << argv[2] 
           << "\n";
    return TCL_ERROR;
  }

  // Determine the number of nodes
  int nen = -1;
  auto it = NodeCounts.find(argv[1]);
  if (it != NodeCounts.end())
    nen = it->second;

  int argi = 3;
  // Parse node tags
  std::vector<int> multi_nodes;
  {
    int list_argc;
    TCL_Char **list_argv;
    if (Tcl_SplitList(interp, argv[argi], &list_argc, &list_argv) == TCL_OK && list_argc >= 2) {
      for (int i = 0; i < list_argc; ++i) {
        int node;
        if (Tcl_GetInt(interp, list_argv[i], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid node " << list_argv[i]
                 << "\n";
          return TCL_ERROR;
        }
        multi_nodes.push_back(node); 
      }
      nen = multi_nodes.size();
      Tcl_Free((char *)list_argv);
      argi += 1;

    } else {
      if (nen == -1) {
        opserr << OpenSees::PromptValueError 
               << "Nodes must be supplied in a list for element type " << argv[1] 
               << "\n";
        return TCL_ERROR;
      }
      if (argi + nen > argc) {
        opserr << OpenSees::PromptValueError 
               << "expected " << nen << " nodes for element type " << argv[1] 
               << "\n";
        return TCL_ERROR;
      }
      for (int i=0; i<nen; i++) {
        int node;
        if (Tcl_GetInt(interp, argv[argi++], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError 
                 << "invalid node tag " << argv[argi-1] 
                 << "\n";
          return TCL_ERROR;
        }
        multi_nodes.push_back(node);
      }
    }
  }


  // Section/Material
  int mat_tag;
  SectionForceDeformation *section = nullptr;
  double b[3] = {0.0, 0.0, 0.0}; // body forces

  Vector3D local_x{};
  bool updateBasis = false;
  bool corotational = false;
  bool use_eas = true;
  bool use_drill_stab = false;
  double drilling_stab = 0.01;
  ASDShellQ4::DrillingDOFMode drill_mode = ASDShellQ4::DrillingDOF_Elastic;

  enum class Position : int {
    Section, B1, B2, B3, End
  };
  ArgumentTracker<Position> tracker;
  std::set<int> positional;

  //
  // Keywords
  //
  for (int i=argi; i<argc; i++) {
    if (strcmp(argv[i], "-section") == 0) {
      i++;
      if (i== argc) {
        opserr << OpenSees::PromptValueError 
                << "-section requires argument"
                << "\n";
        return TCL_ERROR;
      }

      int stag;
      if (Tcl_GetInt(interp, argv[i], &stag) != TCL_OK) {
        opserr << OpenSees::PromptValueError 
                << "failed to read section tag"
                << "\n";
        return TCL_ERROR;
      }
      section = builder->getTypedObject<SectionForceDeformation>(stag);
      if (section == nullptr)
        return TCL_ERROR;

      tracker.consume(Position::Section);
    }
    else if (strcmp(argv[i], "-updateBasis") == 0) {
      updateBasis = true;
    }

    else if ((strcasecmp(argv[i], "-corotational") == 0))
      corotational = true;

    else if (strcmp(argv[i], "-noeas") == 0) {
      use_eas = false;
    }

    else if (strcmp(argv[i], "-drillingStab") == 0) {
      if (drill_mode != ASDShellQ4::DrillingDOF_Elastic) {
          opserr << "Error: element ASDShellQ4: -drillingStab and -drillingNL options are mutually exclusive\n";
          return 0;
      }
      if (argc < i + 2) {
          opserr << "Error: drilling stabilization parameter not provided with -drillingStab option\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i+1], &drilling_stab) != TCL_OK) {
          opserr << "Error: cannot get drilling stabilization parameter with -drillingStab option\n";
          return TCL_ERROR;
      }
      drilling_stab = std::max(0.0, std::min(1.0, drilling_stab));
      drill_mode = ASDShellQ4::DrillingDOF_Elastic;
      use_drill_stab = true;
      i++;
    }
    else if (strcmp(argv[i], "-drillingNL") == 0) {
      if (use_drill_stab) {
          opserr << "Error: -drillingStab and -drillingNL options are mutually exclusive\n";
          return 0;
      }
      drill_mode = ASDShellQ4::DrillingDOF_NonLinear;
      drilling_stab = 1.0;
    }
    // else if (strcmp(argv[i], "-local") == 0) {
    //     if (OPS_GetNumRemainingInputArgs() < 3) {
    //         opserr << "Error: element ASDShellQ4: not enough arguments for -local options (3 components are required)\n";
    //         return 0;
    //     }
    //     for (int i = 0; i < 3; ++i) {
    //         double local_x_com;
    //         if (OPS_GetDoubleInput(&numData, &local_x_com) == 0) {
    //             local_x(i) = local_x_com;
    //         }
    //         else {
    //             opserr << "Error: element ASDShellQ4: cannot get the component " << i + 1 << " for the local X axis\n";
    //             return 0;
    //         }
    //     }
    // }

    else
      positional.insert(i);
  }

  //
  // Positional arguments
  //
  for (int i : positional) {
    switch (tracker.current()) {
      case Position::Section:
        if (Tcl_GetInt(interp, argv[i], &mat_tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid material tag " << argv[i] << "\n";
          return TCL_ERROR;
        } else {
          section = builder->getTypedObject<SectionForceDeformation>(mat_tag);
          if (section == nullptr)
            return TCL_ERROR;
        }
        tracker.increment();
        break;
  
      case Position::End:
      default:
        opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << "\n";
        return TCL_ERROR;
    }
    //
    // Check required positional arguments
    //
    while (tracker.current() != Position::End) {
      switch (tracker.current()) {
        case Position::Section:
          opserr << OpenSees::PromptValueError 
                 << "missing required positional argument section\n";
          return TCL_ERROR;
        default:
          break;
      }
      tracker.increment();
    }
  }

  //
  // Create the element
  //
  Element* theElement = nullptr;
  if (nen == 3) {

    std::array<int, 3> nodes;
    for (int i=0; i<3; i++)
      nodes[i] = multi_nodes[i];

    if (strcasecmp(argv[1], "ShellDKGT") == 0) {
      theElement = new ShellDKGT(tag, nodes[0], nodes[1], nodes[2], *section, b[0], b[1], b[2]);

    } else if (strcasecmp(argv[1], "ShellNLDKGT") == 0) {
      theElement = new ShellNLDKGT(tag, nodes[0], nodes[1], nodes[2], *section);

    }
  }
  else if (nen == 4) {
    std::array<int, 4> nodes;
    for (int i=0; i<4; i++)
      nodes[i] = multi_nodes[i];

    if ((strcasecmp(argv[1], "ShellMITC4") == 0) || 
        (strcasecmp(argv[1], "Shell") == 0)) {
      theElement = new ShellMITC4(tag, nodes[0], nodes[1], nodes[2], nodes[3], *section, updateBasis);

    } else if (strcasecmp(argv[1], "ShellDKGQ") == 0) {
      theElement = new ShellDKGQ(tag, nodes[0], nodes[1], nodes[2], nodes[3], *section);

    } else if (strcasecmp(argv[1], "ShellNLDKGQ") == 0) {
      theElement = new ShellNLDKGQ(tag, nodes[0], nodes[1], nodes[2], nodes[3], *section);

    } else if (strcasecmp(argv[1], "ASDShellQ4") == 0) {

      theElement = new ASDShellQ4(tag, 
                                  nodes[0], nodes[1], nodes[2], nodes[3],
                                  section, local_x, corotational, use_eas, drill_mode, drilling_stab);

    } else if (strcasecmp(argv[1], "ShellNLDKGQThermal") == 0) {
      theElement = new ShellNLDKGQThermal(tag, nodes[0], nodes[1], nodes[2], nodes[3], *section);

    } else if (strcasecmp(argv[1], "ShellMITC4Thermal") == 0) {
      theElement = new ShellMITC4Thermal(tag, nodes[0], nodes[1], nodes[2], nodes[3], *section);
    }
  }
  else if (nen == 9) {
    std::array<int, 9> nodes;
    for (int i=0; i<9; i++)
      nodes[i] = multi_nodes[i];

    if (strcasecmp(argv[1], "ShellMITC9") == 0) {
      theElement = new ShellMITC9(tag, nodes[0], nodes[1], nodes[2], nodes[3],
        nodes[4], nodes[5], nodes[6], nodes[7], nodes[8], *section);
    }
  }

  if (theElement == nullptr) {
    opserr << OpenSees::PromptValueError << "failed to create element\n";
    return TCL_ERROR;
  }

  //
  //
  //
  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }
  
  return TCL_OK;
}


#include <elementAPI.h>
Element*
TclDispatch_newShellANDeS(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{

  if (argc < 6) {
    opserr << "Want: element ShellANDeS $tag $iNode $jNode $kNode $thick $E $nu $rho";
    return nullptr;
  }

  int numArgs = OPS_GetNumRemainingInputArgs();

  int iData[4];
  int numData = 4;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  double dData[11];
  numArgs = OPS_GetNumRemainingInputArgs();
  if (OPS_GetDoubleInput(&numArgs, dData) != 0) {
    opserr << "WARNING invalid double thickness: element ShellANDeS \n";
    return nullptr;
  }

  Element *theElement = nullptr;

  if (numArgs == 4) {
    theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3],
                                dData[0], dData[1], dData[2], dData[3]);
  } else if (numArgs == 11) {
    theElement =
        new ShellANDeS(iData[0], iData[1], iData[2], iData[3], dData[0],
                       dData[1], dData[2], dData[3], dData[4], dData[5],
                       dData[6], dData[7], dData[8], dData[9], dData[10]);
  }

  return theElement;
}

#if 0
Element*
TclDispatch_newShellDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  if (argc < 6) {
    opserr << "Want: element ShellDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4], *theSection);
}

int
TclDispatch_newShellMITC4(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  bool updateBasis = false;
  Element *theElement = nullptr;

  if (argc < 6) {
    opserr << "Want: element ShellMITC4 $tag $iNode $jNode $kNode $lNode $secTag <-updateBasis>";
    return TCL_ERROR;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return TCL_ERROR;
  }

  if (argc == 7) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-updateBasis") == 0)
      updateBasis = true;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return TCL_ERROR;

  theElement = new ShellMITC4(iData[0], iData[1], iData[2], iData[3], iData[4],
                              *theSection, updateBasis);

  if (builder->getDomain()->addElement(theElement) == false)
    return TCL_ERROR;
  return TCL_OK;
}



Element*
TclDispatch_newShellMITC9(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 11) {
    opserr << "Want: element ShellMITC9 $tag $node1 $node2 .... $node9 $secTag";
    return nullptr;
  }

  int iData[11];
  int numData = 11;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[10]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellMITC9(iData[0], iData[1], iData[2], iData[3], iData[4], iData[5],
                     iData[6], iData[7], iData[8], iData[9], *theSection);
}


Element*
TclDispatch_newShellNLDKGQ(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr
        << "Want: element ShellNLDKGQ $tag $iNode $jNoe $kNode $lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQ(iData[0], iData[1], iData[2], iData[3], iData[4],
                               *theSection);

}

Element*
TclDispatch_newShellDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;


  double b_data[3] = {0, 0, 0};

  int num_remaining_args = OPS_GetNumRemainingInputArgs();

  if (num_remaining_args > 3) {
    num_remaining_args = 3;
  }
  if (num_remaining_args > 0) {
    if (OPS_GetDoubleInput(&num_remaining_args, b_data) < 0) {
      opserr << "WARNING: invalid double b_data\n";
      return nullptr;
    }
  }

  return new ShellDKGT(iData[0], iData[1], iData[2], iData[3],
                       *theSection, b_data[0], b_data[1], b_data[2]);
}

Element*
TclDispatch_newShellMITC4Thermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellMITC4Thermal $tag $iNode $jNoe $kNode $lNode "
              "$secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellMITC4Thermal(iData[0], iData[1], iData[2], iData[3], iData[4], *theSection);
}



Element*
TclDispatch_newShellNLDKGQThermal(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 6) {
    opserr << "Want: element ShellNLDKGQThermal $tag $iNode $jNoe $kNode "
              "$lNode $secTag";
    return nullptr;
  }

  int iData[6];
  int numData = 6;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[5]);
  if (theSection == nullptr)
    return nullptr;

  return new ShellNLDKGQThermal(iData[0], iData[1], iData[2], iData[3],
                                      iData[4], *theSection);

}


Element*
TclDispatch_newShellNLDKGT(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder* builder = (BasicModelBuilder*)clientData;

  int numArgs = OPS_GetNumRemainingInputArgs();

  if (numArgs < 5) {
    opserr << "Want: element ShellNLDKGT $tag $iNode $jNoe $kNode $secTag";
    return nullptr;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag\n";
    return nullptr;
  }

  SectionForceDeformation *theSection = builder->getTypedObject<SectionForceDeformation>(iData[4]);
  if (theSection == nullptr)
    return nullptr;

  return
      new ShellNLDKGT(iData[0], iData[1], iData[2], iData[3], *theSection);

}
#endif
