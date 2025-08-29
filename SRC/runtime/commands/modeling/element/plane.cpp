//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the implementation of the
//              TclBasicBuilder_addFourNodeQuad() command.
//
// Written: fmk
// Created: 07/99
//
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <set>
#include <unordered_map>
#ifdef _MSC_VER 
#  include <string.h>
#  define strcasecmp _stricmp
#else
#  include <strings.h>
#endif
#include <tcl.h>
#include <Logging.h>
#include <Parsing.h>
#include <Domain.h>
#include <ArgumentTracker.h>
#include <section/PlaneSection.h>
#include <FourNodeQuad.h>
#include <FourNodeQuad3d.h>
#include <FourNodeQuadWithSensitivity.h>
#include <ConstantPressureVolumeQuad.h>
#include <EnhancedQuad.h>
#include <NineNodeMixedQuad.h>
#include <NineNodeQuad.h>
#include <EightNodeQuad.h>
#include <LagrangeQuad.h>
//
#include <Tri31.h>
#include <SixNodeTri.h>
#include <BasicModelBuilder.h>

#include <algorithm>
#include <string>

namespace {
static
std::string toLower( const std::string & s )
{
    std::string copy = s;
    transform( copy.begin( ), copy.end( ), copy.begin( ), 
        [](unsigned char c) { return std::tolower(c); });
    return copy;
}

static bool 
equalsIgnoreCase(const std::string & lhs, const std::string & rhs )
{
    return toLower(lhs) == toLower( rhs );
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
} // namespace

using namespace OpenSees;

static std::unordered_map<std::string, int, CaseInsensitive, CaseInsensitive> 
NodeCounts = {
  {"Quad",                        4},
  {"FourNodeQuad",                4},
  {"quad8n",                      8},
  {"EightNodeQuad",               8},
  {"FourNodeQuad3d",              4},
  {"FourNodeQuadWithSensitivity", 4},
  {"ConstantPressureVolumeQuad",  4},
  {"EnhancedQuad",                4},
  {"NineNodeQuad",                9},
  {"quad9n",                      9},
  {"NineNodeMixedQuad",           9},
  {"LagrangeQuad",                4},
  {"Tri31",                       3},
  {"CST",                         3},
  {"T3",                          3},
  {"SixNodeTri",                  6}
};

int
TclBasicBuilder_addFourNodeQuad(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  //         0     1     2       3      4      5      6     7        8      9        10      11  12  13
  //  10  element name eleTag? iNode? jNode? kNode? lNode? thk?    type?  matTag?
  //  14  element name eleTag? iNode? jNode? kNode? lNode? thk?    type?  matTag? <pressure? rho? b1? b2?>
  //
  //      element name eleTag? iNode? jNode? kNode? lNode? sec? <pressure? rho?      b1?     b2?>

  if (builder->getNDM() != 2 || (builder->getNDF() != 2 && builder->getNDF() != 3)) {
    opserr << OpenSees::PromptValueError 
           << "model dimensions and/or nodal DOF not compatible with quad element\n";
    return TCL_ERROR;
  }


  int tag;
  if (argc < 6) {
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
  int nen = -1;

  auto it = NodeCounts.find(argv[1]);
  if (it != NodeCounts.end())
    nen = it->second;

  int argi = 3;
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



  int mat_tag;
  NDMaterial *nd_mat = nullptr;
  double thickness = 1.0;
  double p = 0.0;   // uniform normal traction (pressure)
  double rho = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;
  TCL_Char *type = nullptr;
  if (true) {
    enum class Position : int {
      Thickness, Type, Material, Pressure, Density, B1, B2, End
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
        PlaneSection<NDMaterial>* section = builder->getTypedObject<PlaneSection<NDMaterial>>(stag);
        if (section == nullptr)
          return TCL_ERROR;

        thickness = section->getThickness();
        nd_mat  = section->getMaterial();
        mat_tag = nd_mat->getTag();
        tracker.consume(Position::Material);
        tracker.consume(Position::Type);
        tracker.consume(Position::Thickness);
      }
      else
        // continue;
        positional.insert(i);
    }

    //
    // Positional arguments
    //
    for (int i : positional) {
      switch (tracker.current()) {
        case Position::Material:
          if (Tcl_GetInt(interp, argv[i], &mat_tag) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid material tag " << argv[i] << "\n";
            return TCL_ERROR;
          } else {
            nd_mat = builder->getTypedObject<NDMaterial>(mat_tag);
            if (nd_mat == nullptr)
              return TCL_ERROR;
      
            nd_mat = nd_mat->getCopy(type);
            if (nd_mat == nullptr) {
              opserr << OpenSees::PromptValueError << "invalid material\n";
              return TCL_ERROR;
            }
          }
          tracker.increment();
          break;
    
        case Position::Thickness:
          if (Tcl_GetDouble(interp, argv[i], &thickness) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid thickness\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::Type:
          type = argv[i];

          if (strcmp(type,"PlaneStrain") != 0   && 
              strcmp(type,"PlaneStress") != 0   &&
              strcmp(type,"PlaneStrain2D") != 0 && 
              strcmp(type,"PlaneStress2D") != 0) {
            opserr << OpenSees::PromptValueError 
                   << "improper material type: " << type << "\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;
        
        case Position::Pressure:
          if (Tcl_GetDouble(interp, argv[i], &p) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid pressure\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::Density:
          if (Tcl_GetDouble(interp, argv[i], &rho) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid density\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::B1:
          if (Tcl_GetDouble(interp, argv[i], &b1) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid b1\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::B2:
          if (Tcl_GetDouble(interp, argv[i], &b2) != TCL_OK) {
            opserr << OpenSees::PromptValueError << "invalid b2\n";
            return TCL_ERROR;
          }
          tracker.increment();
          break;

        case Position::End:
        default:
          opserr << OpenSees::PromptParseError << "unexpected argument " << argv[i] << "\n";
          return TCL_ERROR;
      }
    }
    //
    // Check required positional arguments
    //
    while (tracker.current() != Position::End) {
      switch (tracker.current()) {
        case Position::Thickness:
          opserr << OpenSees::PromptValueError 
                 << "missing required positional argument thickness\n";
          return TCL_ERROR;
        case Position::Type:
          opserr << OpenSees::PromptValueError 
                 << "missing required positional argument type\n";
          return TCL_ERROR;
        case Position::Material:
          opserr << OpenSees::PromptValueError 
                 << "missing required positional argument material\n";
          return TCL_ERROR;
        default:
          break;
      }
      tracker.increment();
    }
  }


  //
  //
  //
  Element* theElement = nullptr;
  if (nen == 3) {
    if (nd_mat == nullptr) {
      opserr << OpenSees::PromptValueError << "invalid material\n";
      return TCL_ERROR;
    }

    std::array<int, 3> nodes;
    for (int i=0; i<3; i++)
      nodes[i] = multi_nodes[i];

    if (strcasecmp(argv[1], "Tri31") == 0) {
      theElement = 
          new Tri31(tag, nodes, *nd_mat, thickness, p, rho, b1, b2);
    }
  }

  else if (nen == 4) {
    std::array<int, 4> nodes;
    for (int i=0; i<4; i++)
      nodes[i] = multi_nodes[i];

    if (strcasecmp(argv[1], "LagrangeQuad") == 0) {
      Mate<2> *mat_2d = builder->getTypedObject<Mate<2>>(mat_tag);
      if (mat_2d == nullptr)
        return TCL_ERROR;

      theElement =
          new LagrangeQuad<4,4>(tag, nodes, *mat_2d,
                                thickness, p, rho, b1, b2);

    } else {
      if (nd_mat == nullptr) {
        opserr << OpenSees::PromptValueError 
               << "invalid material" 
               << "\n";
        return TCL_ERROR;
      }
      if (strcasecmp(argv[1], "EnhancedQuad") == 0) {
        theElement =
            new EnhancedQuad(tag, nodes, *nd_mat, thickness);
      }
      else if (strcasecmp(argv[1], "bbarQuad") == 0 || 
               strcasecmp(argv[1], "mixedQuad") == 0) {
        theElement = new ConstantPressureVolumeQuad(tag, nodes[0], nodes[1], nodes[2], nodes[3], *nd_mat, thickness);

      }
      else 
        theElement =
            new FourNodeQuad(tag, nodes, *nd_mat, thickness, p, rho, b1, b2);
    }
  }

  else if (nen == 8) {
    std::array<int, 8> nodes;
    for (int i=0; i<8; i++)
      nodes[i] = multi_nodes[i];
    if ((strcasecmp(argv[1], "eightnodequad") == 0) ||
        (strcasecmp(argv[1], "quad8n") == 0) || 
        (strcasecmp(argv[1], "quad") == 0)) {
      if (nd_mat == nullptr) {
        opserr << OpenSees::PromptValueError 
               << "invalid material" 
               << "\n";
        return TCL_ERROR;
      }
      theElement =
          new EightNodeQuad(tag, nodes, *nd_mat, thickness, p, rho, b1, b2);
    }
  }
  else if (nen == 9) {
    std::array<int, 9> nodes;
    for (int i=0; i<9; i++)
      nodes[i] = multi_nodes[i];
    if ((strcasecmp(argv[1], "ninenodequad") == 0) ||
        (strcasecmp(argv[1], "quad9n") == 0) || 
        (strcasecmp(argv[1], "quad") == 0)) {

      if (nd_mat == nullptr) {
        opserr << OpenSees::PromptValueError 
                << "invalid material" 
                << "\n";
        return TCL_ERROR;
      }
      theElement =
          new NineNodeQuad(tag, nodes, *nd_mat, thickness, p, rho, b1, b2);
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

  if (type != nullptr)
    delete nd_mat;
  return TCL_OK;
}


int
TclBasicBuilder_addConstantPressureVolumeQuad(ClientData clientData,
                                              Tcl_Interp *interp, int argc,
                                              TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError 
           << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  int argStart = 2;

  // check the number of arguments is correct
  if ((argc - argStart) < 7) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element ConstantPressureVolumeQuad eleTag? iNode? jNode? kNode? lNode? thk? matTag?\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, iNode, jNode, kNode, lNode, matID;
  double thickness = 1.0;

  if (Tcl_GetInt(interp, argv[argStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the ConstantPressureVolumeQuad and add it to the Domain
  ConstantPressureVolumeQuad *theConstantPressureVolumeQuad =
      new ConstantPressureVolumeQuad(tag, 
                                     iNode, jNode, kNode, lNode, *theMaterial, thickness);


  if (builder->getDomain()->addElement(theConstantPressureVolumeQuad) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theConstantPressureVolumeQuad;
    return TCL_ERROR;
  }

  return TCL_OK;
}


/*  *****************************************************************************

    N I N E   N O D E   M I X E D  Q U A D

    *****************************************************************************
 */

int
TclBasicBuilder_addNineNodeMixedQuad(ClientData clientData, Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 11) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element NineNodeMixedQuad  eleTag?"
           << " iNode? jNode? kNode? lNode? mNode, nNode, pNode, qNode, centerNode  matTag?\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, iNode, jNode, kNode, lNode;
  int mNode, nNode, pNode, qNode;
  int centerNode;
  int matID;

  if (Tcl_GetInt(interp, argv[argStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid NineNodeMixedQuad eleTag" << "\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[5 + argStart], &mNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid mNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &nNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid nNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[7 + argStart], &pNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid pNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8 + argStart], &qNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid qNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[9 + argStart], &centerNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid centerNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[10 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the NineNodeMixedQuad and add it to the Domain
  NineNodeMixedQuad *theNineNodeMixed = new NineNodeMixedQuad(
      tag, iNode, jNode, kNode, lNode, mNode, nNode, pNode,
      qNode, centerNode, *theMaterial);

  if (builder->getDomain()->addElement(theNineNodeMixed) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theNineNodeMixed;
    return TCL_ERROR;
  }

  return TCL_OK;
}

int
TclBasicBuilder_addFourNodeQuadWithSensitivity(ClientData clientData,
                                               Tcl_Interp *interp, int argc,
                                               TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder == 0 || clientData == 0) {
    opserr << OpenSees::PromptValueError << "builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 8) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element FourNodeQuad eleTag? iNode? jNode? kNode? lNode? "
              "thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, iNode, jNode, kNode, lNode, matID;
  double thickness = 1.0;
  double p = 0.0; // uniform normal traction (pressure)
  double r = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element tag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  TCL_Char *type = argv[6 + argStart];

  if (Tcl_GetInt(interp, argv[7 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) > 11) {
    if (Tcl_GetDouble(interp, argv[8 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[9 + argStart], &r) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid rho\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the FourNodeQuad and add it to the Domain
  FourNodeQuadWithSensitivity *theFourNodeQuadWithSensitivity =
      new FourNodeQuadWithSensitivity(tag, iNode, jNode, kNode,
                                      lNode, *theMaterial, type, thickness, p,
                                      r, b1, b2);

  if (builder->getDomain()->addElement(theFourNodeQuadWithSensitivity) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theFourNodeQuadWithSensitivity;
    return TCL_ERROR;
  }

  return TCL_OK;
}

//
// Description: This file contains the implementation of
// TclBasicBuilder_addFourNodeQuadUP() ,
// TclBasicBuilder_addNineFourNodeQuadUP() ,
// TclBasicBuilder_addBBarFourNodeQuadUP(),
//
// Zhaohui Yang and Jinchi Lu (September 2009)
//
#include <stdlib.h>
#include <Domain.h>

#include <FourNodeQuadUP.h>
#include <Nine_Four_Node_QuadUP.h>
#include <BBarFourNodeQuadUP.h>


/*  *****************************************************************************

    Q U A D  U_P

    *****************************************************************************
 */

int
TclBasicBuilder_addFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 3) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 11) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element FourNodeQuadUP eleTag? iNode? jNode? kNode? "
              "lNode? thk? matTag? bulk? rho? perm_x? perm_y? <b1? b2? "
              "pressure? dM? dK?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int tag, iNode, jNode, kNode, lNode, matID;
  double thickness, bk, r, perm1, perm2;
  double p = 0.0; // uniform normal traction (pressure)
  double b1 = 0.0;
  double b2 = 0.0;

  // TCL_Char *type;
  if (Tcl_GetInt(interp, argv[argStart], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid FourNodeQuadUP eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7 + argStart], &bk) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8 + argStart], &r) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9 + argStart], &perm1) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lateral permeability\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &perm2) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid vertical permeability\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 12) {
    if (Tcl_GetDouble(interp, argv[11 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 13) {
    if (Tcl_GetDouble(interp, argv[12 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 14) {
    if (Tcl_GetDouble(interp, argv[13 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the FourNodeQuadUP and add it to the Domain
  FourNodeQuadUP *theFourNodeQuadUP = new FourNodeQuadUP(
      tag, iNode, jNode, kNode, lNode, *theMaterial, "PlaneStrain",
      thickness, bk, r, perm1, perm2, b1, b2, p);

  if (builder->getDomain()->addElement(theFourNodeQuadUP) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theFourNodeQuadUP;
    return TCL_ERROR;
  }

  return TCL_OK;
}


/*  *****************************************************************************

    9-4-N O D E  Q U A D  U_P

    *****************************************************************************
 */

int
TclBasicBuilder_addNineFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder == 0 || clientData == 0) {
    opserr << OpenSees::PromptValueError << "builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions not compatible with 9-4-NodeQuadUP "
              "element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 16) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr
        << "Want: element FourNodeQuadUP eleTag? Node1? ... Node9? thk? "
           "matTag? bulk? rho? perm_x? perm_y? <b1? b2? pressure? dM? dK?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int Ninetag, Node[9], matID;
  double thickness, bk, r, perm1, perm2;
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &Ninetag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid FourNodeQuadUP eleTag" << "\n";
    return TCL_ERROR;
  }
  for (int i = 1; i <= 9; i++) {
    if (Tcl_GetInt(interp, argv[i + argStart], &Node[i - 1]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid Node\n";
      return TCL_ERROR;
    }
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[11 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[12 + argStart], &bk) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[13 + argStart], &r) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[14 + argStart], &perm1) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lateral permeability\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[15 + argStart], &perm2) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid vertical permeability\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 17) {
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 18) {
    if (Tcl_GetDouble(interp, argv[17 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the FourNodeQuadUP and add it to the Domain
  NineFourNodeQuadUP *theNineFourNodeQuadUP = new NineFourNodeQuadUP(
      Ninetag, Node[0], Node[1], Node[2], Node[3], Node[4],
      Node[5], Node[6], Node[7], Node[8], *theMaterial, "PlaneStrain",
      thickness, bk, r, perm1, perm2, b1, b2);
  if (theNineFourNodeQuadUP == 0) {
    opserr << OpenSees::PromptValueError << "ran out of memory creating element\n";
    return TCL_ERROR;
  }

  if (builder->getDomain()->addElement(theNineFourNodeQuadUP) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theNineFourNodeQuadUP;
    return TCL_ERROR;
  }

  return TCL_OK;
}


/*  *****************************************************************************

    B B A R  Q U A D  U_P

    *****************************************************************************
 */

int
TclBasicBuilder_addBBarFourNodeQuadUP(ClientData clientData, Tcl_Interp *interp,
                                      int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder == 0 || clientData == 0) {
    opserr << OpenSees::PromptValueError << "builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 2 || builder->getNDF() != 3) {
    opserr << OpenSees::PromptValueError 
           << "model dimensions and/or nodal DOF not compatible "
              "with QuadUP element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 11) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element bbarQuadUP eleTag? iNode? jNode? kNode? lNode? "
              "thk? matTag? bulk? rho? perm_x? perm_y? <b1? b2? "
              "pressure? dM? dK?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int BBartag, iNode, jNode, kNode, lNode, matID;
  double thickness, bk, r, perm1, perm2;
  double p = 0.0; // uniform normal traction (pressure)
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &BBartag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid BBarFourNodeQuadUP eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[5 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[7 + argStart], &bk) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid bulk modulus\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[8 + argStart], &r) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid fluid mass density\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9 + argStart], &perm1) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lateral permeability\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &perm2) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid vertical permeability\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) >= 12) {
    if (Tcl_GetDouble(interp, argv[11 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 13) {
    if (Tcl_GetDouble(interp, argv[12 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }
  if ((argc - argStart) >= 14) {
    if (Tcl_GetDouble(interp, argv[13 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the BBarFourNodeQuadUP and add it to the Domain
  BBarFourNodeQuadUP *theBBarFourNodeQuadUP = new BBarFourNodeQuadUP(
      BBartag, iNode, jNode, kNode, lNode, *theMaterial,
      "PlaneStrain", thickness, bk, r, perm1, perm2, b1, b2, p);

  if (theBBarFourNodeQuadUP == nullptr) {
    opserr << OpenSees::PromptValueError << "ran out of memory creating element\n";
    return TCL_ERROR;
  }

  if (builder->getDomain()->addElement(theBBarFourNodeQuadUP) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theBBarFourNodeQuadUP;
    return TCL_ERROR;
  }

  return TCL_OK;
}

#if 0
int
TclDispatch_newTri31(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char **const argv)
{
  BasicModelBuilder* builder = static_cast<BasicModelBuilder*>(clientData);
  
  if (argc < 9) {
    opserr << "Invalid #args, want: "
    //            0      1      2       3      4      5    6     7      8        9       10  11  12
           << "element Tri31 eleTag? iNode? jNode? kNode? thk? type? matTag? <pressure? rho? b1? b2?>\n";
    return TCL_ERROR;
  }

  int tag, matID;
  std::array<int,3> nodes;
  char *type;
  double thickness,
         pressure=0, 
         density=0,
         b1 = 0,
         b2 = 0;
  
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element tag\n";
    return TCL_ERROR;
  }
  for (int i=0; i<3; i++) {
    if (Tcl_GetInt(interp, argv[i+3], &nodes[i]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid node tag\n";
      return TCL_ERROR;
    }
  }

  if (Tcl_GetDouble(interp, argv[6], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element thickness\n";
    return TCL_ERROR;
  }


  type = strdup(argv[7]);
  if (   strcmp(type,"PlaneStrain") != 0 
      && strcmp(type,"PlaneStress") != 0
      && strcmp(type,"PlaneStrain2D") != 0 
      && strcmp(type,"PlaneStress2D") != 0) {
        opserr << OpenSees::PromptValueError 
               << "improper material type: " << type << "for Tri31"
               << OpenSees::SignalMessageEnd;
        return TCL_ERROR;
  }


  if (Tcl_GetInt(interp, argv[8], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid material tag\n";
    return TCL_ERROR;
  }
  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr) {
    return TCL_ERROR;
  }
  
  if (argc > 9  && Tcl_GetDouble(interp, argv[ 9], &pressure) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element pressure\n";
    return TCL_ERROR;
  }
  if (argc > 10 && Tcl_GetDouble(interp, argv[10], &density) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element density\n";
    return TCL_ERROR;
  }
  if (argc > 11 && Tcl_GetDouble(interp, argv[11], &b1) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element load b1\n";
    return TCL_ERROR;
  }
  if (argc > 12 && Tcl_GetDouble(interp, argv[12], &b2) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid element load b2\n";
    return TCL_ERROR;
  }

  // parsing was successful, create the element

  Element *theElement = new Tri31(tag, 
                         nodes,
                         *theMaterial, 
                         type,
                         thickness, 
                         pressure,
                         density, 
                         b1, b2);

  if (builder->getDomain()->addElement(theElement) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theElement;
    return TCL_ERROR;
  }

  free(type);

  return TCL_OK;
}


// Regular nine node quad

int
TclBasicBuilder_addNineNodeQuad(ClientData clientData, Tcl_Interp *interp, int argc,
                                TCL_Char ** const argv)
{
  // TODO: assertions, clean up
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder == 0 || clientData == 0) {
    opserr << OpenSees::PromptValueError << "builder has been destroyed\n";
    return TCL_ERROR;
  }

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 13) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element NineNodeQuad eleTag? iNode? jNode? kNode? lNode? "
              "nNode? mNode? pNode? qNode? cNode? thk? type? matTag? "
              "<pressure? rho? b1? b2?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int NineNodeQuadId;
  std::array<int,9> nodes{};
  int matID;
  double thickness = 1.0;
  double p = 0.0;   // uniform normal traction (pressure)
  double rho = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &NineNodeQuadId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid NineNodeQuad eleTag" << "\n";
    return TCL_ERROR;
  }

  for (int i=0; i<9; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart+i], &nodes[i]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid node\n";
      return TCL_ERROR;
    }

  if (Tcl_GetDouble(interp, argv[10 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  TCL_Char *type = argv[11 + argStart];

  if (Tcl_GetInt(interp, argv[12 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) > 16) {
    if (Tcl_GetDouble(interp, argv[13 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[14 + argStart], &rho) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[16 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the NineNodeQuad and add it to the Domain
  NineNodeQuad *theNineNodeQuad = new NineNodeQuad(
      NineNodeQuadId, nodes, *theMaterial, type, thickness, p, rho, b1, b2);


  if (builder->getDomain()->addElement(theNineNodeQuad) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theNineNodeQuad;
    return TCL_ERROR;
  }

  return TCL_OK;
}


//
// Regular eight node quad
//
int
TclBasicBuilder_addEightNodeQuad(ClientData clientData, Tcl_Interp *interp,
                                 int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 12) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element EightNodeQuad eleTag? iNode? jNode? kNode? lNode? "
              "nNode? mNode? pNode? qNode? thk? type? matTag? <pressure? rho? "
              "b1? b2?>\n";
    return TCL_ERROR;
  }

  // get the id and end nodes
  int EightNodeQuadId;
  int iNode, jNode, kNode, lNode, nNode, mNode, pNode, qNode;
  int matID;
  double thickness = 1.0;
  double p = 0.0;   // uniform normal traction (pressure)
  double rho = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &EightNodeQuadId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid EightNodeQuad eleTag" << "\n";
    return TCL_ERROR;
  }
  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid iNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid jNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid kNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid lNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[5 + argStart], &nNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid nNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[6 + argStart], &mNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid mNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[7 + argStart], &pNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid pNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetInt(interp, argv[8 + argStart], &qNode) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid qNode\n";
    return TCL_ERROR;
  }

  if (Tcl_GetDouble(interp, argv[9 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    return TCL_ERROR;
  }

  TCL_Char *type = argv[10 + argStart];

  if (Tcl_GetInt(interp, argv[11 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) > 15) {
    if (Tcl_GetDouble(interp, argv[12 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13 + argStart], &rho) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[14 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[15 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the EightNodeQuad and add it to the Domain
  EightNodeQuad *theEightNodeQuad = new EightNodeQuad(
      EightNodeQuadId, iNode, jNode, kNode, lNode, nNode, mNode, pNode, qNode,
      *theMaterial, type, thickness, p, rho, b1, b2);

  if (builder->getDomain()->addElement(theEightNodeQuad) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theEightNodeQuad;
    return TCL_ERROR;
  }

  return TCL_OK;
}

#endif


int
TclBasicBuilder_addSixNodeTri(ClientData clientData, Tcl_Interp *interp, int argc,
                              TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  if (builder->getNDM() != 2 || builder->getNDF() != 2) {
    opserr << OpenSees::PromptValueError << "-- model dimensions and/or nodal DOF not compatible "
              "with quad element\n";
    return TCL_ERROR;
  }

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) < 10) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    opserr << "Want: element SixNodeTri eleTag? iNode? jNode? kNode? lNode? "
              "nNode? mNode? pNode? qNode? thk? type? matTag? <pressure? rho? "
              "b1? b2?>\n";
    return TCL_ERROR;
  }

  int SixNodeTriId;
  std::array<int,6> nodes;
  int matID;
  double thickness = 1.0;
  double p = 0.0;   // uniform normal traction (pressure)
  double rho = 0.0; // mass density
  double b1 = 0.0;
  double b2 = 0.0;

  if (Tcl_GetInt(interp, argv[argStart], &SixNodeTriId) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid SixNodeTri eleTag" << "\n";
    return TCL_ERROR;
  }
  
  // Nodes
  for (int i=0; i<6; i++)
    if (Tcl_GetInt(interp, argv[1 + argStart + i], &nodes[i]) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid node\n";
      return TCL_ERROR;
    }

  if (Tcl_GetDouble(interp, argv[7 + argStart], &thickness) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid thickness\n";
    opserr << "SixNodeTri element: " << SixNodeTriId << "\n";
    return TCL_ERROR;
  }

  TCL_Char *type = argv[8 + argStart];

  if (Tcl_GetInt(interp, argv[9 + argStart], &matID) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid matID\n";
    return TCL_ERROR;
  }

  if ((argc - argStart) > 13) {
    if (Tcl_GetDouble(interp, argv[10 + argStart], &p) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid pressure\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[11 + argStart], &rho) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[12 + argStart], &b1) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b1\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[13 + argStart], &b2) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid b2\n";
      return TCL_ERROR;
    }
  }

  NDMaterial *theMaterial = builder->getTypedObject<NDMaterial>(matID);
  if (theMaterial == nullptr)
    return TCL_ERROR;


  // now create the SixNodeTri and add it to the Domain
  SixNodeTri *theSixNodeTri =
      new SixNodeTri(SixNodeTriId, nodes,
                     *theMaterial, type, thickness, p, rho, b1, b2);


  if (builder->getDomain()->addElement(theSixNodeTri) == false) {
    opserr << OpenSees::PromptValueError << "could not add element to the domain\n";
    delete theSixNodeTri;
    return TCL_ERROR;
  }

  return TCL_OK;
}
