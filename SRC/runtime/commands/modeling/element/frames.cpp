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
//
// Written: cmp, mhs, rms, fmk
//
// Created: Feb 2023
//
// Standard library
  #include <string>
  #include <array>
  #include <algorithm>
  #include <vector>
  #include <utility>
  #include <stdlib.h>
  #include <string.h>
  #include <assert.h>
  #include <math.h>
  #ifdef _MSC_VER 
  #  include <string.h>
  #  define strcasecmp _stricmp
  #else
  #  include <strings.h>
  #endif
  #define strcmp strcasecmp
  
  // Parsing
  #include <tcl.h>
  #include <Logging.h>
  #include <Parsing.h>
  #include <ArgumentTracker.h>
  
  // Model
  #include <Node.h>
  #include <Domain.h>
  #include <BasicModelBuilder.h>
  
  // Sections
  #include <FrameSection.h>
  #include <ElasticSection2d.h>
  #include <ElasticSection3d.h>
  
  // Elements
  #include "ElasticBeam2d.h"
  #include "ElasticBeam2d.h"
  #include "ElasticBeam3d.h"
  #include "ElasticBeam3d.h"
  #include "PrismFrame2d.h"
  #include "PrismFrame2d.h"
  #include "PrismFrame3d.h"
  #include "PrismFrame3d.h"
  
  #include <CubicFrame3d.h>
  #include <ForceFrame3d.h>
  #include <ForceDeltaFrame3d.h>
  #include <EulerFrame3d.h>
  #include <EulerDeltaFrame3d.h>
  #include <ExactFrame3d.h>
  
  #include <DispBeamColumn2d.h>
  #include <DispBeamColumn2dThermal.h>
  #include <DispBeamColumn3d.h>
  #include <DispBeamColumn3dThermal.h>
  #include <DispBeamColumnNL2d.h>
  
  #include <ElasticForceBeamColumn2d.h>
  #include <ElasticForceBeamColumn3d.h>
  #include <ElasticForceBeamColumnWarping2d.h>
  
  #include <ForceBeamColumn2d.h>
  #include <ForceBeamColumn2d.h>
  #include <ForceBeamColumn2dThermal.h>
  #include <ForceBeamColumn3d.h>
  #include <ForceBeamColumnCBDI2d.h>
  #include <ForceBeamColumnCBDI3d.h>
  #include <ForceBeamColumnWarping2d.h>
  #include <TimoshenkoBeamColumn2d.h>
  
  // Quadrature
  #include <BeamIntegration.h>
  #include <LobattoBeamIntegration.h>
  #include <LegendreBeamIntegration.h>
  #include <HingeEndpointBeamIntegration.h>
  #include <HingeMidpointBeamIntegration.h>
  #include <HingeRadauBeamIntegration.h>
  #include <HingeRadauTwoBeamIntegration.h>

  #include <transform/FrameTransformBuilder.hpp>

using namespace OpenSees;

struct Options {
  int mass_flag;
  int use_mass;
  int shear_flag;
  int geom_flag;
};


extern BeamIntegration*     GetBeamIntegration(TCL_Char* type, int);
extern BeamIntegrationRule* GetHingeStencil(int argc, TCL_Char ** const argv);

static inline int
CheckTransformation(Domain& domain, int iNode, int jNode, CrdTransf& transform)
{
  Node* ni = domain.getNode(iNode);
  Node* nj = domain.getNode(jNode);
  if (ni == nullptr || nj == nullptr) {
    opserr << OpenSees::PromptValueError << "nodes not found with tags "
           << iNode << " and " << jNode
           << OpenSees::SignalMessageEnd;
  }

  if (transform.initialize(ni, nj) != 0) {
    if (transform.getInitialLength() <= 0.0) {
      opserr << OpenSees::PromptValueError 
            << "element has zero or negative initial length "
            << transform.getInitialLength()
            << "; check for duplicate nodes"
            << OpenSees::SignalMessageEnd;
    }
    else {
      opserr << OpenSees::PromptValueError 
            << "transformation with tag " << transform.getTag()
            << " could not be initialized with nodes "
            << iNode << " and " << jNode
            << "; check orientation"
            << OpenSees::SignalMessageEnd;
    }
    return TCL_ERROR;
  }
  return TCL_OK;
}


template <int ndm, typename Transform, typename Section>
static Element*
CreateFrame(BasicModelBuilder& builder, 
            const char* name,
            int tag,
            std::vector<int>& nodev,
            int transfTag,
            const std::vector<int>& section_tags,
            BeamIntegration& beamIntegr,
            double mass, int max_iter, double tol,
            Options& options) 
{

  std::vector<Section*> sections;

  // Finalize sections
  assert(section_tags.size() != 0);
  for (int tag : section_tags) {
    Section *section = builder.getTypedObject<Section>(tag);
    if (section == nullptr)
      return nullptr;
    sections.push_back(section);
  }
  int nIP = sections.size();

  SectionForceDeformation** secptrs = (SectionForceDeformation**)(sections.data());

  if (options.shear_flag == -1) {
    options.shear_flag = 0;
    const ID& resultants = sections[0]->getType();
    for (int i=0; i< sections[0]->getOrder(); i++)
      if (resultants(i) == FrameStress::Vy)
        options.shear_flag = 1;
  }

  // Finalize the coordinate transform
  CrdTransf* theTransf = builder.getTypedObject<CrdTransf>(transfTag);
  if (theTransf == nullptr) {
    opserr << OpenSees::PromptValueError << "transformation not found with tag " << transfTag << "\n";
    return nullptr;
  }

  //
  // Create the element
  //
  Element  *theElement   = nullptr;

  int iNode = nodev[0],
      jNode = nodev[1];
  bool use_mass = options.use_mass;

  if constexpr (ndm == 2) {

    if (strcmp(name, "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);
    else if (strcmp(name, "timoshenkoBeamColumn") == 0)
      theElement =
          new TimoshenkoBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);
    else if (strcmp(name, "dispBeamColumn") == 0)
      theElement =
          new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass, options.mass_flag);
    else if (strcmp(name, "dispBeamColumnNL") == 0)
      theElement =
          new DispBeamColumnNL2d(tag, iNode, jNode, nIP, secptrs,
                                 beamIntegr, *theTransf, mass);

    else if (strcmp(name, "dispBeamColumnThermal") == 0)
      theElement = new DispBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);


    else if (strcmp(name, "dispBeamColumnWithSensitivity") == 0)
      theElement = new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);


    // Force formulations
    else if (strcmp(name, "forceBeamColumnCBDI") == 0)
      theElement = new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                             beamIntegr, *theTransf, mass, false);

    else if (strcmp(name, "forceBeamColumnCSBDI") == 0)
      return  new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                              beamIntegr, *theTransf, mass, true);

    else if (strcmp(name, "forceBeamColumnWarping") == 0)
      return
          new ForceBeamColumnWarping2d(tag, iNode, jNode, nIP,
                                       secptrs, beamIntegr, *theTransf);

    else if (strcmp(name, "forceBeamColumnThermal") == 0)
      return new ForceBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcmp(name, "elasticForceBeamColumnWarping") == 0)
      return new ElasticForceBeamColumnWarping2d(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf);
    else
      theElement =
          new ForceBeamColumn2d(tag, iNode, jNode, nIP, secptrs,
                                beamIntegr, *theTransf, mass, max_iter, tol);


  } else {
    //
    // ndm == 3
    //

    if (strstr(name, "Frame") != nullptr) {
      if (strstr(name, "Exact") == nullptr) {

        if (CheckTransformation(*builder.getDomain(), nodev[0], nodev[nodev.size()-1], *theTransf) != TCL_OK)
          return nullptr;
        std::array<int, 2> nodes {nodev[0], nodev[1]};

        FrameTransformBuilder* tb = builder.getTypedObject<FrameTransformBuilder>(transfTag);

        if (!tb) {
          opserr << OpenSees::PromptValueError << "invalid transform\n";
          return nullptr;
        }


        if (strcmp(name, "EulerFrame") == 0) {
            theElement = new EulerFrame3d(tag, nodes, nIP, 
                                          sections.data(),
                                          beamIntegr, 
                                          *tb, 
                                          mass, options.mass_flag);
        }

        if (strcmp(name, "CubicFrame") == 0) {
          if (options.shear_flag)
            theElement = new CubicFrame3d<true,0>(tag, nodes, 
                                          sections,
                                          beamIntegr, 
                                          *theTransf, // TODO: Use FrameTransformBuilder
                                          mass);
          else
            theElement = new CubicFrame3d<false,0>(tag, nodes, 
                                          sections,
                                          beamIntegr, 
                                          *theTransf, // TODO: Use FrameTransformBuilder
                                          mass);
        } 

        else if (strcmp(name, "DisplFrame") == 0) {
          theElement =  new EulerDeltaFrame3d(tag, nodes, sections,
                                              beamIntegr, *theTransf, 
                                              mass, 
                                              options.mass_flag, 
                                              use_mass);
        }

        else if ((strstr(name, "Force") != 0) ||
                (strcmp(name, "MixedFrame") == 0)) {
          if (strcmp(name, "ForceDeltaFrame") == 0 || options.geom_flag) {
            if (!options.shear_flag)
              static_loop<2,6>([&](auto nip) constexpr {
                if (nip.value == sections.size())
                  theElement = new ForceDeltaFrame3d<nip.value, 4>(tag, nodes, sections,
                                                beamIntegr, *tb, 
                                                mass, 
                                                options.mass_flag, 
                                                use_mass,
                                                max_iter, tol,
                                                options.shear_flag
                                                );
              });
            else
              static_loop<2,6>([&](auto nip) constexpr {
                if (nip.value == sections.size())
                  theElement = new ForceDeltaFrame3d<nip.value, 6>(tag, nodes, sections,
                                                beamIntegr, *tb, 
                                                mass, 
                                                options.mass_flag, 
                                                use_mass,
                                                max_iter, tol,
                                                options.shear_flag
                                                );
              });
          } else {
            int ndf = builder.getNDF();

            static_loop<0, 3>([&](auto nwm) constexpr {
              if (nwm.value + 6 == ndf) {
                if (!options.shear_flag) {
                  static_loop<2,30>([&](auto nip) constexpr {
                    if (nip.value == sections.size())
                      theElement = new ForceFrame3d<nip.value, 4+nwm.value*2, nwm.value>(tag, 
                                                    nodes, sections,
                                                    beamIntegr, *tb,
                                                    mass, options.mass_flag, use_mass,
                                                    max_iter, tol
                                                    );
                    });
                }
                else
                  theElement = new ForceFrame3d<20, 6+nwm.value*2, nwm.value>(tag, 
                                                nodes, sections,
                                                beamIntegr, *tb,
                                                mass, options.mass_flag, use_mass,
                                                max_iter, tol
                                                );
              }
            });
          }
        }
      }

      else if (strcmp(name, "ExactFrame") == 0) {
        if (!options.shear_flag) {
          opserr << OpenSees::PromptValueError 
                 << "ExactFrame3d requires shear formulation"
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
        int ndf = builder.getNDF();
        if (sections.size() < nodev.size()-1)
          for (unsigned i = 0; i < nodev.size()-1; ++i)
            sections.push_back(sections[0]);
        
        unsigned nen = nodev.size();
        static_loop<2,6>([&](auto nn) constexpr {
          if (nn.value == nen) {
            std::array<int, nn.value> nodes;
            std::copy_n(nodev.begin(), nn.value, nodes.begin());
            static_loop<0,4>([&](auto nwm) constexpr {
              if (nwm.value+6 == ndf)
                theElement = new ExactFrame3d<nn.value, nwm.value>(tag, nodes, sections.data(), *theTransf);
            });
          }
        });
        if (theElement == nullptr) {
          opserr << OpenSees::PromptValueError 
                 << "invalid number of dofs for ExactFrame; got " << ndf 
                 << OpenSees::SignalMessageEnd;
          return nullptr;
        }
      }
    }

    else if (strcmp(name, "elasticForceBeamColumn") == 0)
      theElement = new ElasticForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs, 
                                                beamIntegr, *theTransf, mass);

    else if (strcasecmp(name, "dispBeamColumn") == 0)
      theElement = new DispBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                        beamIntegr, *theTransf, 
                                        mass, options.mass_flag);

    else if (strcmp(name, "dispBeamColumnWithSensitivity") == 0)
      theElement = new DispBeamColumn3d(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcmp(name, "dispBeamColumnThermal") == 0)
      theElement = new DispBeamColumn3dThermal(
          tag, iNode, jNode, nIP, secptrs, beamIntegr, *theTransf, mass);

    else if (strcmp(name, "forceBeamColumnCBDI") == 0)
      theElement = new ForceBeamColumnCBDI3d(tag, iNode, jNode, nIP, secptrs,
                                             beamIntegr, *theTransf, 
                                             mass, false, max_iter, tol);
    else
      theElement = new ForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                         beamIntegr, *theTransf, mass, max_iter, tol);
  }
  return theElement;
}


// 0       1    2 3  4
// element beam 1 $i $j 0 1 2
//
//  a)
//     0       1     2    3    4    5    6
//                                  0    1
//     element $type $tag $ndi $ndj $trn "Gauss arg1 arg2 ..." 
//             <-mass $mass> <-iter $iter $tol>
//
//  b)
//                                  0    1
//     element(type, tag, ndi, ndj, trn, itag, 
//              iter=(10, 1e-12), mass=0.0)
//
//  c) "Original/Obsolete"
//                                  0    1    2
//     element $type $tag $ndi $ndj $nip $sec $trn 
//             <-mass $mass> <-iter $iter $tol> <-integration $ityp>
//
//  d) 
//                                  0    1         ... (2 + nIP)
//     element $type $tag $ndi $ndj $nip -sections ... $trn 
//             <-mass $massDens> <-cMass> <-integration $ityp>
//
//  e)
//                                   0
//     element $type $tag $ndi $ndj  $trn 
//             -sections {...}
//             <-mass $massDens> <-cMass> <-integration $ityp>
//
//
// Integration may be specitied as either 
//   i  ) a single name, 
//   ii ) a pattern spec, or 
//   iii) a tag for a pattern
// 
// if a list of sections is given with -sections, or nIP is provided, then 
// we must have an integration with the form (i)
//
// 1) Parse common keyword args: 
//      "-mass" $mass, 
//      "-cMass"/"-lMass"
//      "-mass-form" $form
//
//      "-iter" $iter $tol, 
//      
//      "-integration" $Integration
//        - first try parsing $Integration as integer ($itag, form (iii))
//          if successfull, populate section_tags and continue
//        - next try parsing $Integration as basic quadrature ($ityp, form (i))
//        - finally, try parsing as full pattern spec
//
//      "-section" $Tag
//
//      "-transform" $Tag
//      "-vertical" {}
//      "-horizontal" {}
//
//      "-sections" $SectionTags
//        - if cannot split $SectionTags as list, then
//          mark "-sections" as positional and continue
//          with keyword loop
//        - Check if "-integration" was provided already; if so, it must have been in form (i);
//          otherwise throw an error.
//        - Parse $SectionTags
//          -sections {...} may occur anywhere
//          -sections ...   must occur after nIP is obtained
//
//
//  2) 
//     If pos[1] == "-sections" then command is Form (d):
//        nIP = pos[0]
//        trn = pos[2+nIP]
//
//     else
//
//     switch (pos.size())
//     case 1: // Form (e)
//        trn = pos[0]
//        if (section_tags.size() == 0)
//          ERROR
//
//     case 2: 
//        // Form (a) or (b)
//        trn = pos[0]
//        if GetInt(interp, pos[1]):
//           itag = pos[1]
//        else
//           ParseHingeScheme(pos[1])
//
//      case 3: 
//         // Form (c)
//         nip = int(pos[0])
//         sec = int(pos[1])
//         trn = int(pos[2])
//
int
TclBasicBuilder_addForceBeamColumn(ClientData clientData, Tcl_Interp *interp,
                                   int argc, TCL_Char **const argv)
{
  assert(clientData != nullptr);
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;
  Domain *domain = builder->getDomain();
  assert(domain != nullptr);

  int status = TCL_OK;

//enum class GaussMethod {
//  None,
//  Rule, // Quadrature name + section locations; Hinge methods
//  Quad, // Quadrature name
//} gauss_method = GaussMethod::None;

  // collect positional arguments
  std::vector<int> positions;


  // 
  // Preliminary checks
  //
  int ndm = builder->getNDM();
  int ndf = builder->getNDF();

  { // Check dimension and DOFs of problem
    int ok = 0;
    if ((ndm == 2 && ndf == 3) || (ndm == 2 && ndf == 4))
      ok = 1;
    if (ndm == 3 && ndf >= 6)
      ok = 1;

    if (ok == 0) {
      opserr << OpenSees::PromptValueError << "ndm = " << ndm << " and ndf = " << ndf
             << " not compatible with Frame element" << "\n";
      return TCL_ERROR;
    }
  }

  if (argc < 6) {
    opserr << OpenSees::PromptValueError << "insufficient arguments\n";
    return TCL_ERROR;
  }

  //
  // Required positional arguments
  //
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << OpenSees::PromptValueError << "invalid " << " tag " << tag << "\n";
    return TCL_ERROR;
  }

  int argi = 5;
  // int iNode=0, jNode=0;
  // bool multi_node = false;
  std::vector<int> multi_nodes;
  {
    int list_argc;
    TCL_Char **list_argv;
    if (Tcl_SplitList(interp, argv[3], &list_argc, &list_argv) == TCL_OK && list_argc >= 2) {
      argi -= 1;
      
      for (int i = 0; i < list_argc; ++i) {
        int node;
        if (Tcl_GetInt(interp, list_argv[i], &node) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid node\n";
          return TCL_ERROR;
        }
        multi_nodes.push_back(node); 
      }
      Tcl_Free((char *)list_argv);
    }
    else {
      int iNode, jNode;
      if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid iNode\n";
        return TCL_ERROR;
      }
  
      if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid jNode\n";
        return TCL_ERROR;
      }
      multi_nodes.push_back(iNode);
      multi_nodes.push_back(jNode);
    }
  }

  int max_iter = 10;
  double tol  = 1.0e-12;
  double mass = 0.0;
  bool use_mass = false;
  int transfTag;
  std::vector<int> section_tags;
  const char* integration_type = nullptr;
  BeamIntegration   *beamIntegr   = nullptr;
  BeamIntegrationRule  *theRule   = nullptr;
  int itg_tag;

  // If we get a BeamIntegration from a BeamIntegrationRule
  // then we dont own it and can't delete it
  bool deleteBeamIntegr = true;
  bool removeHingeIntegr = false;

  //
  // Defaults
  //
  struct Options options;
  options.mass_flag  =  0;
  options.shear_flag = -1;
  options.geom_flag  =  0;
  if (strcasecmp(argv[1], "elasticBeamColumn") == 0) {
    options.shear_flag = 0;
  }
  if (strcasecmp(argv[1], "dispBeamColumn") == 0 || 
      strcasecmp(argv[1], "nonlinearBeamColumn") == 0) {
    options.shear_flag = 0;
  }
  else if (strcasecmp(argv[1], "timoshenkoBeamColumn") == 0) {
    options.shear_flag = 1;
  }


  //
  // Parse positions
  //
  {
    while (argi < argc) {
      // Shear
      if (strcmp(argv[argi], "-shear") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -shear $flag\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.shear_flag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid shear_flag, expected integer\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 2;
      }

      // Geometry
      else if (strcmp(argv[argi], "-order") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -order $flag\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.geom_flag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid geom_flag, expected integer\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 2;
      }

      // -iter $max_iter $tol 
      else if (strcmp(argv[argi], "-iter") == 0) {
        if (argc < argi + 3) {
          opserr << OpenSees::PromptValueError << "not enough -iter args need -iter max_iter? tol?\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &max_iter) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid max_iter\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetDouble(interp, argv[argi + 2], &tol) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid tol\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 3;
      }
      // mass
      else if (strcmp(argv[argi], "-mass") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -mass $mass\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid mass\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 2;
        use_mass = true;

      // mass type
      } else if ((strcmp(argv[argi], "-lMass") == 0) ||
                 (strcmp(argv[argi], "lMass") == 0)) {
        options.mass_flag = 0;
        argi++;
      }
      else if ((strcmp(argv[argi], "-cMass") == 0) ||
                 (strcmp(argv[argi], "cMass") == 0)) {
        options.mass_flag = 1;
        argi++;
      }

      // Quadrature
      else if (strcmp(argv[argi], "-integration") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -integration $integration\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        integration_type = argv[argi];
        // beamIntegr = GetBeamIntegration(argv[argi]);

        // if (beamIntegr == nullptr) {
        //   opserr << OpenSees::PromptValueError << "invalid integration type\n";
        //   status = TCL_ERROR;
        //   goto clean_up;
        // }
        argi++;
      }

      // Transform
      else if (strcmp(argv[argi], "-transform") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -transform $transform\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        if (Tcl_GetInt(interp, argv[argi], &transfTag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid transform\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi++;
      }

      // Section
      else if (strcmp(argv[argi], "-section") == 0) {
        if (argc < argi + 2) {
          opserr << OpenSees::PromptValueError << "not enough arguments, expected -section $section\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        int sec_tag;
        if (Tcl_GetInt(interp, argv[argi], &sec_tag) != TCL_OK) {
          opserr << OpenSees::PromptValueError << "invalid sec_tag\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        section_tags.push_back(sec_tag);

        argi++;
      }


    //else if (strcmp(argv[argi], "-sections") == 0) {
    // split possible lists present in argv
    //char *List = Tcl_Merge(inArgc, inArgv);
    //if (List == nullptr) {
    //  opserr << OpenSees::PromptValueError << "problem merging list\n";
    //  return TCL_ERROR;
    //}
    //  int secc;
    //  TCL_Char ** secv;
    //  if (Tcl_SplitList(interp, argv[positions[2]], &secc, &secv) != TCL_OK) {
    //    opserr << OpenSees::PromptValueError << "problem splitting list\n";
    //    return TCL_ERROR;
    //  }
    //Tcl_Free((char *)List);

    //}
      else {
        positions.push_back(argi);
        argi++;
      }
    }
  }

  //
  // II Parse Positional Arguments
  //

  // Version d)
  // positional arguments are:
  //   0: nIP
  //   1: -sections
  //   2: secTag1
  //   3: secTag2...
  if (positions.size() > 1 && strcmp(argv[positions[1]], "-sections") == 0) {
    
    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nIP\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    // TODO: Make sure 2+nIP < positions.size()
  
    // Get section tags
    for (int i = 0; i < nIP; i++) {
      int secTag;
      if (Tcl_GetInt(interp, argv[positions[2+i]], &secTag) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "invalid section\n";
        status = TCL_ERROR;
        goto clean_up;
      }
      section_tags.push_back(secTag);
    }

    if (Tcl_GetInt(interp, argv[positions[2+nIP]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid transform\n";
      status = TCL_ERROR;
      goto clean_up;
    }

  }

  // Version e) ?
  else if (positions.size() == 1) {
    if (section_tags.empty()) {
      status = TCL_ERROR;
      goto clean_up;
    }
  }
  
  // Version a or b
  else if (positions.size() == 2 || positions.size() > 3) {
    // Here we create a BeamIntegrationRule (theRule) which is a pair of
    // section tags and a BeamIntegration. In this case we do not
    // delete the BeamIntegration because it is owned by theRule.

    // Geometric transformation
    if (Tcl_GetInt(interp, argv[positions[0]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError 
             << "invalid transform " << argv[positions[0]] 
             << OpenSees::SignalMessageEnd;
      status = TCL_ERROR;
      goto clean_up;
    }
    
    // Version b)
    if (Tcl_GetInt(interp, argv[positions[1]], &itg_tag) == TCL_OK) {
      deleteBeamIntegr = false;
      removeHingeIntegr = false;
    }

    // Version a)
    else {
      // If we fail to parse an integer tag, treat it like an inline definition
      builder->findFreeTag<BeamIntegrationRule>(itg_tag);
      std::string integrCommand{argv[positions[1]]};
      integrCommand.insert(integrCommand.find(" "), " "+std::to_string(itg_tag)+" ");
      integrCommand.insert(0, "beamIntegration ");
      if (Tcl_Eval(interp, integrCommand.c_str()) != TCL_OK) {
        opserr << OpenSees::PromptValueError << "failed to parse integration\n";
        status = TCL_ERROR;
        goto clean_up;
      }

      deleteBeamIntegr = false;
      removeHingeIntegr = true;
    }

    theRule = builder->getTypedObject<BeamIntegrationRule>(itg_tag);
    if (theRule == nullptr) {
      status = TCL_ERROR;
      goto clean_up;
    }

    beamIntegr = theRule->getBeamIntegration();
    const ID& secTags = theRule->getSectionTags();

    for (int i=0; i < secTags.Size(); i++)
      section_tags.push_back(secTags(i));
  }

  // Version c)
  //
  // .. nip section transf
  else if (positions.size() == 3) {

    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid nIP\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    if (nIP <= 0) {
      opserr << OpenSees::PromptValueError << "invalid nIP, must be > 0\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    //
    int secTag;
    if (Tcl_GetInt(interp, argv[positions[1]], &secTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid secTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    for (int i=0; i < nIP; i++)
      section_tags.push_back(secTag);

    // Transform
    if (Tcl_GetInt(interp, argv[positions[2]], &transfTag) != TCL_OK) {
      opserr << OpenSees::PromptValueError << "invalid transform\n";
      status = TCL_ERROR;
      goto clean_up;
    }
  }

  //
  // Finalize the quadrature
  //
  // TODO
  if (section_tags.size() == 1 && theRule == nullptr) {
    if (strstr(argv[1], "isp") == 0) {
      section_tags.resize(5, section_tags[0]);
    } else {
      section_tags.resize(5, section_tags[0]);
    }
  }

  if (beamIntegr == nullptr) {
    if (integration_type == nullptr) {
      if (strstr(argv[1], "ispBeam") == 0) {
        integration_type = "Lobatto";
      } else {
        integration_type = "Legendre";
      }
    }

    if ((beamIntegr = GetBeamIntegration(integration_type, section_tags.size())) == nullptr) {
      opserr << OpenSees::PromptValueError << "invalid integration type or size\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    deleteBeamIntegr = true;
  }

  //
  //
  options.use_mass = use_mass;
  {
    Element *theElement = ndm == 2 
                        ? CreateFrame<2, CrdTransf, FrameSection>(*builder, argv[1], tag, multi_nodes, transfTag, 
                                                              section_tags, *beamIntegr, mass, max_iter, tol, options)
                        : CreateFrame<3, CrdTransf, FrameSection>(*builder, argv[1], tag, multi_nodes, transfTag, 
                                                                        section_tags, *beamIntegr, mass, max_iter, tol, options);

                                                                        
    if (theElement == nullptr) {
      status = TCL_ERROR;
      goto clean_up;
    }
    if (domain->addElement(theElement) == false) {
      opserr << OpenSees::PromptValueError 
             << "could not add element to the domain"
             << OpenSees::SignalMessageEnd;
      delete theElement;
      status = TCL_ERROR;
      goto clean_up;
    }
  }


clean_up:
  //
  // Clean up
  //
  if (deleteBeamIntegr && beamIntegr != nullptr)
    delete beamIntegr;

  if (removeHingeIntegr) {
    builder->removeObject<BeamIntegrationRule>(itg_tag);
    delete theRule;
  }

  return status;
}


//
//  BeamWithHinges
//
//     element beamWithHinges tag? ndI? ndJ? secTagI? lenI? secTagJ? lenJ? 
//        E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> 
//        <-iter maxIters tolerance>
//
int
TclBasicBuilder_addBeamWithHinges(ClientData clientData, Tcl_Interp *interp,
                                  int argc, TCL_Char ** const argv)
{
  BasicModelBuilder *builder = (BasicModelBuilder*)clientData;

  int NDM = builder->getNDM();
  int NDF = builder->getNDF();

  // Plane frame element
  if (NDM == 2 && NDF == 3) {
    if (argc < 13) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> "
                "<-iter maxIters tolerance>"
             << "\n";
      return TCL_ERROR;
    }

    double massDens = 0.0;
    int max_iters = 10;
    double tol = 1.0e-10;
    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, I;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &I) != TCL_OK) {
      opserr << "WARNING invalid I\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[12], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      return TCL_ERROR;
    }

    bool isShear = false;
    int shearTag = 0;

    if (argc > 13) {
      for (int i = 13; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-constHinge") == 0 && ++i < argc) {
          if (Tcl_GetInt(interp, argv[i], &shearTag) != TCL_OK) {
            opserr << "WARNING invalid constHinge tag\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
          isShear = true;
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &max_iters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    FrameSection *sectionI = builder->getTypedObject<FrameSection>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    FrameSection *sectionJ = builder->getTypedObject<FrameSection>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;

    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection2d elastic(8, E, A, I);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;

      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == nullptr) {
      opserr << "Unknown element type: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    if (isShear) {
      FrameSection *sectionL = builder->getTypedObject<FrameSection>(shearTag);
      if (sectionL == nullptr)
        return TCL_ERROR;

      sections[numSections++] = sectionL;
    }

    theElement = new ForceBeamColumn2d(tag, ndI, ndJ, numSections, 
                                       sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       max_iters, tol);

    delete theBeamIntegr;

    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add element to domain.\n";
      return TCL_ERROR;
    }
  }

  else if (NDM == 3 && NDF == 6) {
    if (argc < 16) {
      opserr << "WARNING insufficient arguments\n";
      opserr << "Want: element beamWithHinges tag? ndI? ndJ? secTagI? lenI? "
                "secTagJ? lenJ? ";
      opserr << "E? A? Iz? Iy? G? J? transfTag? <-shear shearLength?> <-mass "
                "massDens?> <-iter maxIters tolerance>"
             << "\n";
      return TCL_ERROR;
    }

    int tag, ndI, ndJ, secTagI, secTagJ, transfTag;
    double lenI, lenJ, E, A, Iz, Iy, G, J;
    double massDens = 0.0;
    int max_iters = 10;
    double tol = 1.0e-10;
    double shearLength = 1.0;

    if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
      opserr << "WARNING invalid beamWithHinges tag" << "\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[3], &ndI) != TCL_OK) {
      opserr << "WARNING invalid ndI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &ndJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[5], &secTagI) != TCL_OK) {
      opserr << "WARNING invalid secTagI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[6], &lenI) != TCL_OK) {
      opserr << "WARNING invalid lenI\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[7], &secTagJ) != TCL_OK) {
      opserr << "WARNING invalid ndJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[8], &lenJ) != TCL_OK) {
      opserr << "WARNING invalid lenJ\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[9], &E) != TCL_OK) {
      opserr << "WARNING invalid E\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[10], &A) != TCL_OK) {
      opserr << "WARNING invalid A\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[11], &Iz) != TCL_OK) {
      opserr << "WARNING invalid Iz\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[12], &Iy) != TCL_OK) {
      opserr << "WARNING invalid Iy\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[13], &G) != TCL_OK) {
      opserr << "WARNING invalid G\n";
      return TCL_ERROR;
    }

    if (Tcl_GetDouble(interp, argv[14], &J) != TCL_OK) {
      opserr << "WARNING invalid J\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[15], &transfTag) != TCL_OK) {
      opserr << "WARNING invalid transfTag\n";
      return TCL_ERROR;
    }


    if (argc > 16) {
      for (int i = 16; i < argc; ++i) {
        if (strcmp(argv[i], "-mass") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &massDens) != TCL_OK) {
            opserr << "WARNING invalid massDens\n";
            opserr << "BeamWithHinges: " << tag << "\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-shear") == 0 && ++i < argc) {
          if (Tcl_GetDouble(interp, argv[i], &shearLength) != TCL_OK) {
            opserr << "WARNING invalid shear\n";
            return TCL_ERROR;
          }
        }

        if (strcmp(argv[i], "-iter") == 0 && i + 2 < argc) {
          if (Tcl_GetInt(interp, argv[++i], &max_iters) != TCL_OK) {
            opserr << "WARNING invalid maxIters\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[++i], &tol) != TCL_OK) {
            opserr << "WARNING invalid tolerance\n";
            return TCL_ERROR;
          }
        }
      }
    }

    // Retrieve section I from the model builder
    SectionForceDeformation *sectionI = builder->getTypedObject<FrameSection>(secTagI);
    if (sectionI == nullptr)
      return TCL_ERROR;

    // Retrieve section J from the model builder
    SectionForceDeformation *sectionJ = builder->getTypedObject<FrameSection>(secTagJ);
    if (sectionJ == nullptr)
      return TCL_ERROR;


    CrdTransf *theTransf = builder->getTypedObject<CrdTransf>(transfTag);
    if (theTransf == nullptr)
      return TCL_ERROR;


    Element *theElement = nullptr;
    int numSections = 0;
    SectionForceDeformation *sections[10];
    BeamIntegration *theBeamIntegr = nullptr;

    ElasticSection3d elastic(0, E, A, Iz, Iy, G, J);

    if (strcmp(argv[1], "beamWithHinges1") == 0) {
      theBeamIntegr = new HingeMidpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges2") == 0) {
      theBeamIntegr = new HingeRadauTwoBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = sectionI;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = sectionJ;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges3") == 0 ||
               strcmp(argv[1], "beamWithHinges") == 0) {
      theBeamIntegr = new HingeRadauBeamIntegration(lenI, lenJ);

      numSections = 6;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = &elastic;
      sections[4] = &elastic;
      sections[5] = sectionJ;

    } else if (strcmp(argv[1], "beamWithHinges4") == 0) {
      theBeamIntegr = new HingeEndpointBeamIntegration(lenI, lenJ);

      numSections = 4;
      sections[0] = sectionI;
      sections[1] = &elastic;
      sections[2] = &elastic;
      sections[3] = sectionJ;
    }

    if (theBeamIntegr == nullptr) {
      opserr << "Unknown element type: " << argv[1] << "\n";
      return TCL_ERROR;
    }

    // TODO fix shear for beamWithHinges
    /*
    if (isShear) {
      SectionForceDeformation *sectionL = builder->getTypedObject<SectionForceDeformation>(shearTag);

      if (sectionL == 0) {
        opserr << "WARNING section L does not exist\n";
        opserr << "section: " << shearTag;
        opserr << "\nBeamWithHinges: " << tag << "\n";
        return TCL_ERROR;
      }
      sections[numSections++] = sectionL;
    }
    */

    theElement = new ForceBeamColumn3d(tag, ndI, ndJ, numSections, sections,
                                       *theBeamIntegr, *theTransf, massDens,
                                       max_iters, tol);

    delete theBeamIntegr;

    // Add to the domain
    if (builder->getDomain()->addElement(theElement) == false) {
      opserr << "WARNING could not add "
                "element to domain ";
      opserr << tag << "\n";
      return TCL_ERROR;
    }
  }

  else {
    opserr << "ERROR -- model dimension: " << NDM
           << " and nodal degrees of freedom: " << NDF
           << " are incompatible for BeamWithHinges element" << "\n";
    return TCL_ERROR;
  }

  return TCL_OK;
}
