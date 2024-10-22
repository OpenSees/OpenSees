//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation    
//
//===----------------------------------------------------------------------===//
//
// Written: cmp, mhs, rms, fmk
//
// Created: Feb 2023
//
//  BeamWithHinges
//
//     element beamWithHinges tag? ndI? ndJ? secTagI? lenI? secTagJ? lenJ? 
//        E? A? I? transfTag? <-shear shearLength?> <-mass massDens?> 
//        <-iter maxIters tolerance>
#if 1
  #include <array>

  #include <string>
  #include <array>
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
  
    // Geometry
  #include <FrameTransform.h>
  #include <LinearFrameTransf3d.h>
  
    // Elements
  #include "ElasticBeam2d.h"
  #include "ElasticBeam2d.h"
  #include "ElasticBeam3d.h"
  #include "ElasticBeam3d.h"
  #include "PrismFrame2d.h"
  #include "PrismFrame2d.h"
  #include "PrismFrame3d.h"
  #include "PrismFrame3d.h"
  
  #include <element/Frame/Basic/CubicFrame3d.h>
  #include <element/Frame/Basic/ForceFrame3d.h>
  #include <element/Frame/Basic/ForceDeltaFrame3d.h>
  #include <element/Frame/Basic/EulerFrame3d.h>
  #include <element/Frame/Basic/EulerDeltaFrame3d.h>
  #include <element/Frame/Basic/ExactFrame3d.h>
  
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

#else 

// Standard library
#include <string>
#include <array>
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

// Geometry
#include <FrameTransform.h>
#include <LinearFrameTransf3d.h>

// Elements
#include "ElasticBeam2d.h"
#include "ElasticBeam3d.h"
#include "PrismFrame2d.h"
#include "PrismFrame3d.h"

#include <element/Frame/Basic/CubicFrame3d.h>
#include <element/Frame/Basic/ForceFrame3d.h>
#include <element/Frame/Basic/ForceDeltaFrame3d.h>
#include <element/Frame/Basic/EulerFrame3d.h>
#include <element/Frame/Basic/EulerDeltaFrame3d.h>
#include <element/Frame/Basic/ExactFrame3d.h>

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
#endif

using namespace OpenSees;

struct Options {
  int mass_flag;
  int shear_flag;
  int geom_flag;
};


extern BeamIntegration*     GetBeamIntegration(TCL_Char* type);
extern BeamIntegrationRule* GetHingeStencil(int argc, TCL_Char ** const argv);


#if 0
Element*
CreateInelasticFrame(std::string, std::vector<int>& nodes,
                                  std::vector<FrameSection>&, 
                                  BeamIntegration&, 
                              //  FrameQuadrature&,
                                  FrameTransform&,
                                  Options&);
Element*
CreatePrismaticFrame(std::string);
#endif

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
//            ParseHingeScheme(pos[1])
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
    if (ndm == 3 && ndf == 6)
      ok = 1;

    if (ok == 0) {
      opserr << G3_ERROR_PROMPT << "ndm = " << ndm << " and ndf = " << ndf
             << " not compatible with Frame element" << "\n";
      return TCL_ERROR;
    }
  }

  if (argc < 6) {
    opserr << G3_ERROR_PROMPT << "insufficient arguments\n";
    opserr << "Want: element " << argv[1]
           << " tag? iNode? jNode? transfTag? ...\n";
    return TCL_ERROR;
  }

  //
  // Essential positional arguments
  //
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "invalid " << " tag " << tag << "\n";
    return TCL_ERROR;
  }

  int iNode, jNode;
  {
    if (Tcl_GetInt(interp, argv[3], &iNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid iNode\n";
      return TCL_ERROR;
    }

    if (Tcl_GetInt(interp, argv[4], &jNode) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid jNode\n";
      return TCL_ERROR;
    }
  }

  struct Options options;
  options.mass_flag  = 0;
  options.shear_flag = 1;
  options.geom_flag  = 0;

  int max_iter = 10;
  double tol  = 1.0e-12;
  double mass = 0.0;
  bool use_mass = false;
  int transfTag;
  std::vector<int> section_tags;
  BeamIntegration   *beamIntegr   = nullptr;
  BeamIntegrationRule  *theRule   = nullptr;
  int itg_tag;

  // If we get a BeamIntegration from a BeamIntegrationRule
  // then we dont own it and can't delete it
  bool deleteBeamIntegr = true;
  bool removeHingeIntegr = false;


  {
    int argi = 5;
    while (argi < argc) {
      // Shear
      if (strcmp(argv[argi], "-shear") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -shear $flag\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.shear_flag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid shear_flag, expected integer\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 2;
      }

      // Shear
      if (strcmp(argv[argi], "-order") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -order $flag\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &options.geom_flag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid geom_flag, expected integer\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 2;
      }

      // -iter $max_iter $tol 
      else if (strcmp(argv[argi], "-iter") == 0) {
        if (argc < argi + 3) {
          opserr << G3_ERROR_PROMPT << "not enough -iter args need -iter max_iter? tol?\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetInt(interp, argv[argi + 1], &max_iter) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid max_iter\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetDouble(interp, argv[argi + 2], &tol) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid tol\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi += 3;

      // mass
      } else if (strcmp(argv[argi], "-mass") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -mass $mass\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        if (Tcl_GetDouble(interp, argv[argi + 1], &mass) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid mass\n";
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
      } else if ((strcmp(argv[argi], "-cMass") == 0) ||
                 (strcmp(argv[argi], "cMass") == 0)) {
        options.mass_flag = 1;
        argi++;

      }

      // Quadrature
      else if (strcmp(argv[argi], "-integration") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -integration $integration\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        beamIntegr = GetBeamIntegration(argv[argi]);

        if (beamIntegr == nullptr) {
          opserr << G3_ERROR_PROMPT << "invalid integration type\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi++;
      }

      // Transform
      else if (strcmp(argv[argi], "-transform") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -transform $transform\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        if (Tcl_GetInt(interp, argv[argi], &transfTag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        argi++;
      }


      // Section
      else if (strcmp(argv[argi], "-section") == 0) {
        if (argc < argi + 2) {
          opserr << G3_ERROR_PROMPT << "not enough arguments, expected -section $section\n";
          status = TCL_ERROR;
          goto clean_up;
        }

        argi++;
        int sec_tag;
        if (Tcl_GetInt(interp, argv[argi], &sec_tag) != TCL_OK) {
          opserr << G3_ERROR_PROMPT << "invalid sec_tag\n";
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
    //  opserr << G3_ERROR_PROMPT << "problem merging list\n";
    //  return TCL_ERROR;
    //}
    //  int secc;
    //  TCL_Char ** secv;
    //  if (Tcl_SplitList(interp, argv[positions[2]], &secc, &secv) != TCL_OK) {
    //    opserr << G3_ERROR_PROMPT << "problem splitting list\n";
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
  if (positions.size() > 1 && strcmp(argv[positions[1]], "-sections") == 0) {

    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    // TODO: Make sure 2+nIP < positions.size()
  
    // Get section tags
    for (int i = 0; i < nIP; i++) {
      int secTag;
      if (Tcl_GetInt(interp, argv[positions[2+i]], &secTag) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "invalid secTag\n";
        status = TCL_ERROR;
        goto clean_up;
      }
      section_tags.push_back(secTag);
    }

    if (Tcl_GetInt(interp, argv[positions[2+nIP]], &transfTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

  }

  else if (positions.size() == 1) {
    if (section_tags.empty()) {
      status = TCL_ERROR;
      goto clean_up;
    }
  }

  else if (positions.size() == 2 || positions.size() > 3) {
    // Here we create a BeamIntegrationRule (theRule) which is a pair of
    // section tags and a BeamIntegration. In this case we do not
    // delete the BeamIntegration because it is owned by theRule.

    // Geometric transformation
    if (Tcl_GetInt(interp, argv[positions[0]], &transfTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    if (Tcl_GetInt(interp, argv[positions[1]], &itg_tag) == TCL_OK) {
      deleteBeamIntegr = false;
      removeHingeIntegr = false;
    }
    else {
      // If we fail to parse an integer tag, treat it like an inline definition
      builder->findFreeTag<BeamIntegrationRule>(itg_tag);
      std::string integrCommand{argv[positions[1]]};
      integrCommand.insert(integrCommand.find(" "), " "+std::to_string(itg_tag)+" ");
      integrCommand.insert(0, "beamIntegration ");
      if (Tcl_Eval(interp, integrCommand.c_str()) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to parse integration\n";
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

  //
  //
  else if (positions.size() == 3) {
    // Transform
    if (Tcl_GetInt(interp, argv[positions[2]], &transfTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid transfTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    int nIP;
    if (Tcl_GetInt(interp, argv[positions[0]], &nIP) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid nIP\n";
      status = TCL_ERROR;
      goto clean_up;
    }
    if (nIP <= 0) {
      opserr << G3_ERROR_PROMPT << "invalid nIP, must be > 0\n";
      status = TCL_ERROR;
      goto clean_up;
    }


    //
    int secTag;
    if (Tcl_GetInt(interp, argv[positions[1]], &secTag) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "invalid secTag\n";
      status = TCL_ERROR;
      goto clean_up;
    }

    for (int i=0; i < nIP; i++)
      section_tags.push_back(secTag);
  }

  
  // TODO
  if (section_tags.size() == 1 && theRule == nullptr) {
    if (strstr(argv[1], "isp") == 0) {
      section_tags.resize(5, section_tags[0]);
    } else {
      section_tags.resize(3, section_tags[0]);
    }
  }

  // Finalize the quadrature
  if (beamIntegr == nullptr) {
    if (strstr(argv[1], "ispBeam") == 0) {
      beamIntegr = new LobattoBeamIntegration();
    } else {
      beamIntegr = new LegendreBeamIntegration();
    }
    deleteBeamIntegr = true;
  }

  //
  //
  {
    Element          *theElement   = nullptr;
    std::vector<FrameSection*> sections;
    SectionForceDeformation** secptrs = nullptr;
    int nIP = 0;
    CrdTransf        *theTransf2d  = nullptr;
    FrameTransform3d *theTransf3d  = nullptr;

    // Finalize sections
    assert(section_tags.size() != 0);
    for (int tag : section_tags) {
      FrameSection *section = builder->getTypedObject<FrameSection>(tag);
      if (section == nullptr) {
        status = TCL_ERROR;
        goto clean_up;
      }
      sections.push_back(section);
    }
    nIP = sections.size();
    secptrs = (SectionForceDeformation**)(sections.data());


    // Finalize the coordinate transform
    switch (ndm) {
      case 2:
        theTransf2d = builder->getTypedObject<FrameTransform2d>(transfTag);
        if (theTransf2d == nullptr) {
          opserr << G3_ERROR_PROMPT << "transformation not found with tag " << transfTag << "\n";
          status = TCL_ERROR;
          goto clean_up;
        }
        break;

      case 3:
        theTransf3d = builder->getTypedObject<FrameTransform3d>(transfTag);
        if (theTransf3d == nullptr) {
          opserr << G3_ERROR_PROMPT << "transformation not found with tag " << transfTag << "\n";
          status = TCL_ERROR;
          goto clean_up;
        }
    }

    if (ndm == 2) {
      if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
        theElement = new ElasticForceBeamColumn2d(
            tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass);
      else if (strcmp(argv[1], "timoshenkoBeamColumn") == 0)
        theElement =
            new TimoshenkoBeamColumn2d(tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass);
      else if (strcmp(argv[1], "dispBeamColumn") == 0)
        theElement =
            new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass, options.mass_flag);
      else if (strcmp(argv[1], "dispBeamColumnNL") == 0)
        theElement =
            new DispBeamColumnNL2d(tag, iNode, jNode, nIP, secptrs,
                                   *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "dispBeamColumnThermal") == 0)
        theElement = new DispBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "dispBeamColumnWithSensitivity") == 0)
        theElement = new DispBeamColumn2d(tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass);


      // Force formulations
      else if (strcmp(argv[1], "forceBeamColumnCBDI") == 0)
        theElement = new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                               *beamIntegr, *theTransf2d, mass, false);

      else if (strcmp(argv[1], "forceBeamColumnCSBDI") == 0)
        theElement =  new ForceBeamColumnCBDI2d(tag, iNode, jNode, nIP, secptrs,
                                                *beamIntegr, *theTransf2d, mass, true);

      else if (strcmp(argv[1], "forceBeamColumnWarping") == 0)
        theElement =
            new ForceBeamColumnWarping2d(tag, iNode, jNode, nIP,
                                         secptrs, *beamIntegr, *theTransf2d);

      else if (strcmp(argv[1], "forceBeamColumnThermal") == 0)
        theElement = new ForceBeamColumn2dThermal(tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d, mass);

      else if (strcmp(argv[1], "elasticForceBeamColumnWarping") == 0)
        theElement = new ElasticForceBeamColumnWarping2d(
            tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf2d);
      else
        theElement =
            new ForceBeamColumn2d(tag, iNode, jNode, nIP, secptrs,
                                  *beamIntegr, *theTransf2d, mass, max_iter, tol);


    } else {
      // ndm == 3

      if (strcmp(argv[1], "CubicFrame") == 0) {
        std::array<int, 2> nodes {iNode, jNode};

        theElement = new EulerFrame3d(tag, nodes, nIP, sections.data(),
                                      *beamIntegr, *theTransf3d, mass, options.mass_flag);
      } 

      else if (strstr(argv[1], "Force") != 0) {
        std::array<int, 2> nodes {iNode, jNode};

        if (strcmp(argv[1], "ForceDeltaFrame") == 0 || options.geom_flag) {

          theElement = new ForceDeltaFrame3d(tag, nodes, sections,
                                        *beamIntegr, *theTransf3d, 
                                        mass, options.mass_flag, use_mass,
                                        max_iter, tol,
                                        options.shear_flag
                                        );
        } else {
          theElement = new ForceFrame3d(tag, nodes, sections,
                                        *beamIntegr, *theTransf3d,
                                        mass, options.mass_flag, use_mass,
                                        max_iter, tol
                                        );
        }
      }

      else if (strcmp(argv[1], "ExactFrame") == 0) {
        std::array<int, 2> nodes {iNode, jNode};
        theElement = new ExactFrame3d<2,1>(tag, nodes, sections.data(), *theTransf3d);
      }


      else if (strcmp(argv[1], "DisplFrame") == 0) {
        std::array<int, 2> nodes {iNode, jNode};
        theElement =
            new EulerDeltaFrame3d(tag, nodes, sections,
                                  *beamIntegr, *theTransf3d, 
                                  mass, options.mass_flag, use_mass);

      }

      else if (strcmp(argv[1], "elasticForceBeamColumn") == 0)
        theElement = new ElasticForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs, 
                                                  *beamIntegr, *theTransf3d, mass);

      else if (strcmp(argv[1], "dispBeamColumn") == 0)
        theElement = new DispBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                          *beamIntegr, *theTransf3d, 
                                          mass, options.mass_flag);

      else if (strcmp(argv[1], "dispBeamColumnWithSensitivity") == 0)
        theElement = new DispBeamColumn3d(
            tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf3d, mass);

      else if (strcmp(argv[1], "dispBeamColumnThermal") == 0)
        theElement = new DispBeamColumn3dThermal(
            tag, iNode, jNode, nIP, secptrs, *beamIntegr, *theTransf3d, mass);

      else if (strcmp(argv[1], "forceBeamColumnCBDI") == 0)
        theElement = new ForceBeamColumnCBDI3d(tag, iNode, jNode, nIP, secptrs,
                                               *beamIntegr, *theTransf3d, 
                                               mass, false, max_iter, tol);
      else
        theElement = new ForceBeamColumn3d(tag, iNode, jNode, nIP, secptrs,
                                           *beamIntegr, *theTransf3d, mass, max_iter, tol);
    }



    if (domain->addElement(theElement) == false) {
      opserr << G3_ERROR_PROMPT << "could not add element to the domain\n";
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


#if 1
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


    FrameTransform2d *theTransf = builder->getTypedObject<FrameTransform2d>(transfTag);
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
            opserr << "WARNING invalid shearLength\n";
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


    FrameTransform3d *theTransf = builder->getTypedObject<FrameTransform3d>(transfTag);
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
#endif



