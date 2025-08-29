//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
#include <tcl.h>
#include <Parsing.h>
#include <ArgumentTracker.h>
#include <Logging.h>
#include <BasicModelBuilder.h>
#include <Domain.h>
#include <vector>
#include <cstddef>

#include <ElementalLoad.h>
#include <FrameLoad.h>
#include <Beam2dPointLoad.h>
#include <Beam2dUniformLoad.h>
#include <Beam2dPartialUniformLoad.h>
#include <Beam2dTempLoad.h>
#include <Beam2dThermalAction.h>  //L.Jiang [SIF]
#include <Beam3dThermalAction.h>  //L.Jiang [SIF]
#include <ShellThermalAction.h>   //L.Jiang [SIF]
#include <ThermalActionWrapper.h> //L.Jiang [SIF]

#include <Beam3dPointLoad.h>
#include <Beam3dUniformLoad.h>
#include <Beam3dPartialUniformLoad.h>
#include <BrickSelfWeight.h>
#include <SurfaceLoader.h>
#include <SelfWeight.h>
#include <LoadPattern.h>

int
TclCommand_addFrameLoad(ClientData clientData, Tcl_Interp *interp, int argc,
                        TCL_Char **const argv)
{
  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *domain       = builder->getDomain();

  std::vector<int>   tags;
  std::vector<Vector3D> n(1);
  std::vector<Vector3D> m(1);
  std::vector<Vector3D> r(1);
  int shape, basis=FrameLoad::Reference;
  enum class Position : int {
    Force, End
  };
  ArgumentTracker<Position> tracker;


  LoadPattern *pattern = builder->getEnclosingPattern();

  // eleLoad FrameForce $shape -n $n -offset $r -pattern $pattern -basis $basis -ele $ele

  if (argc < 3) {
    opserr << "WARNING eleLoad FrameLoad $shape -force $n -couple $m -offset r -pattern pattern -basis $basis\n";
    return TCL_ERROR;
  }

  if (strcmp(argv[2], "Dirac") == 0)
    shape = FrameLoad::Dirac;
  else if (strcmp(argv[2], "Heaviside") == 0)
    shape = FrameLoad::Heaviside;
  else if (strcmp(argv[2], "Lagrange") == 0) {
    shape = FrameLoad::Lagrange;
    opserr << "Lagrange shape not yet implemented\n";
    return TCL_ERROR;
  }
  else {
    opserr << "WARNING unknown shape for FrameLoad " << argv[2] << "\n";
    return TCL_ERROR;
  }

  // Keywords
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i], "-pattern") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -pattern paramter missing required argument\n";
        return TCL_ERROR;
      }
      int ptag;
      if (Tcl_GetInt(interp, argv[i+1], &ptag) != TCL_OK) {
        opserr << "WARNING pattern parameter expected integer\n";
        return TCL_ERROR;
      }
      pattern = builder->getDomain()->getLoadPattern(ptag);
      if (pattern == nullptr) {
        opserr << "WARNING pattern " << argv[i+1] << " not found\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(argv[i], "-basis") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -basis paramter missing required argument\n";
        return TCL_ERROR;
      }
      if (strcmp(argv[i+1], "global") == 0)
        basis = FrameLoad::Embedding;
      else if ((strcmp(argv[i+1], "reference") == 0) || 
               (strcmp(argv[i+1], "local") == 0))
        basis = FrameLoad::Reference;
      else if (strcmp(argv[i+1], "director") == 0)
        basis = FrameLoad::Director;
      else {
        opserr << "WARNING unknown basis for FrameLoad " << argv[i+1] << "\n";
        return TCL_ERROR;
      }
      i++;
    }
    else if (strcmp(argv[i], "-force") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -force paramter missing required argument\n";
        return TCL_ERROR;
      }
      int list_argc;
      TCL_Char **list_argv;
      if (Tcl_SplitList(interp, argv[i+1], &list_argc, &list_argv) != TCL_OK) {
        opserr << "WARNING force parameter expected list of floats\n";
        return TCL_ERROR;
      }
      if (list_argc != 3) {
        opserr << "WARNING force parameter expected list of 3 floats\n";
        Tcl_Free((char *) list_argv);
        return TCL_ERROR;
      }
      Vector3D force;
      for (int j = 0; j < 3; j++) {
        if (Tcl_GetDouble(interp, list_argv[j], &force[j]) != TCL_OK) {
          opserr << "WARNING force parameter expected list of 3 floats\n";
          Tcl_Free((char *) list_argv);
          return TCL_ERROR;
        }
      }
      // n.push_back(force);
      n[0] = force;
      Tcl_Free((char *) list_argv);

      tracker.consume(Position::Force);
    }
    else if (strcmp(argv[i], "-couple") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -couple paramter missing required argument\n";
        return TCL_ERROR;
      }
      int list_argc;
      TCL_Char **list_argv;
      if (Tcl_SplitList(interp, argv[i+1], &list_argc, &list_argv) != TCL_OK) {
        opserr << "WARNING couple parameter expected list of floats\n";
        return TCL_ERROR;
      }
      if (list_argc != 3) {
        opserr << "WARNING couple parameter expected list of 3 floats\n";
        Tcl_Free((char *) list_argv);
        return TCL_ERROR;
      }
      Vector3D couple;
      for (int j = 0; j < 3; j++) {
        if (Tcl_GetDouble(interp, list_argv[j], &couple[j]) != TCL_OK) {
          opserr << "WARNING couple parameter expected list of 3 floats\n";
          Tcl_Free((char *) list_argv);
          return TCL_ERROR;
        }
      }
      // m.push_back(couple);
      m[0] = couple;
      Tcl_Free((char *) list_argv);

      tracker.consume(Position::Force);
    }
    else if (strcmp(argv[i], "-offset") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -offset paramter missing required argument\n";
        return TCL_ERROR;
      }
      int list_argc;
      TCL_Char **list_argv;
      if (Tcl_SplitList(interp, argv[i+1], &list_argc, &list_argv) != TCL_OK) {
        opserr << "WARNING offset parameter expected list of floats\n";
        return TCL_ERROR;
      }
      if (list_argc != 3) {
        opserr << "WARNING offset parameter expected list of 3 floats\n";
        Tcl_Free((char *) list_argv);
        return TCL_ERROR;
      }
      Vector3D offset;
      for (int j = 0; j < 3; j++) {
        if (Tcl_GetDouble(interp, list_argv[j], &offset[j]) != TCL_OK) {
          opserr << "WARNING offset parameter expected list of 3 floats\n";
          Tcl_Free((char *) list_argv);
          return TCL_ERROR;
        }
      }
      // r.push_back(offset);
      r[0] = offset;
      Tcl_Free((char *) list_argv);
    }
    else if (strcmp(argv[i], "-elements") == 0) {
      if (i == argc-1) {
        opserr << "WARNING -elements paramter missing required argument\n";
        return TCL_ERROR;
      }
      int list_argc;
      TCL_Char **list_argv;
      if (Tcl_SplitList(interp, argv[i+1], &list_argc, &list_argv) != TCL_OK) {
        opserr << "WARNING elements parameter expected list of integers\n";
        return TCL_ERROR;
      }
      for (int j = 0; j < list_argc; j++) {
        int tag;
        if (Tcl_GetInt(interp, list_argv[j], &tag) != TCL_OK) {
          opserr << "WARNING elements parameter expected list of integers\n";
          Tcl_Free((char *) list_argv);
          return TCL_ERROR;
        }
        tags.push_back(tag);
      }
      Tcl_Free((char *) list_argv);
    }
  }

  // Make sure we got everything we need.
  if (tracker.current() != Position::End) {
    opserr << OpenSees::PromptParseError
           << "missing required arguments: ";
    while (tracker.current() != Position::End) {
      switch (tracker.current()) {
        case Position::Force :
          opserr << "force ";
          break;
        case Position::End:
          break;
      }
      if (tracker.current() == Position::End)
        break;
      tracker.consume(tracker.current());
    }
    opserr << "\n";
    return TCL_ERROR;
  }

  if (pattern == nullptr) {
    opserr << "WARNING no current load pattern\n";
    return TCL_ERROR;
  }

  FrameLoad *load = new FrameLoad(basis, shape, n, m, r, *pattern);

  for (int i : tags) {
    Element *elem = domain->getElement(i);
    if (elem == nullptr) {
      opserr << "WARNING eleLoad - no element with tag " << i << "\n";
      delete load;
      return TCL_ERROR;
    }
    if (load->addElement(*elem) != 0) {
      opserr << "WARNING eleLoad - could not add load to element\n";
      delete load;
      return TCL_ERROR;
    }
  }

  if (domain->addElementalLoad(load, pattern->getTag()) == false) {
    opserr
        << "WARNING eleLoad - could not add load to domain\n ";
    delete load;
    return TCL_ERROR;
  }

  return TCL_OK;
}


int
TclCommand_addElementalLoad(ClientData clientData, Tcl_Interp *interp, int argc_main,
                            TCL_Char **const argv_main)
{
  if (argc_main < 2) {
    opserr << "WARNING eleLoad - expecting eleLoad type\n";
    return TCL_ERROR;
  }

  if (strcmp(argv_main[1], "Frame") == 0) {
    return TclCommand_addFrameLoad(clientData, interp, argc_main, argv_main);
  }


  BasicModelBuilder *builder = static_cast<BasicModelBuilder*>(clientData);
  Domain *domain       = builder->getDomain();
  static int eleLoadTag      = 0; // TODO: this is ugly


  int ndm = builder->getNDM();
  bool explicitPatternPassed = false;

  std::vector<int> element_tags;
  int loadPatternTag = 0;

  // Arguments argv_main will be copied into argv. 
  // We initialize with two placeholders for "-type" and $typeName.
  // Everything else will be appended.
  std::vector<const char *> argv{nullptr, nullptr};

  int typeIndex = -1;
  // Create an ID containing the ele tags of all elements
  // for which the load applies.
  int count    = 1;
  while (count < argc_main) {

    // Element tags
    if (strcmp(argv_main[count], "-ele") == 0) {
      count++;
      int eleStart = count;
      int eleEnd   = 0;
      int eleID;
      while (count < argc_main && eleEnd == 0) {
        if (Tcl_GetInt(interp, argv_main[count], &eleID) != TCL_OK)
          eleEnd = count;
        else
          count++;
      }
      if (eleStart != eleEnd) {
        for (int i = eleStart; i < eleEnd; ++i) {
          Tcl_GetInt(interp, argv_main[i], &eleID);
          element_tags.push_back(eleID);
        }
      }
    }

    else if (strcmp(argv_main[count], "-range") == 0) {
      count++;
      int eleStart, eleEnd;
      if (Tcl_GetInt(interp, argv_main[count], &eleStart) != TCL_OK) {
        opserr << "WARNING eleLoad -range invalid eleStart " << argv_main[count]
               << "\n";
        return TCL_ERROR;
      }
      count++;
      if (Tcl_GetInt(interp, argv_main[count], &eleEnd) != TCL_OK) {
        opserr << "WARNING eleLoad -range invalid eleEnd " << argv_main[count] << "\n";
        return TCL_ERROR;
      }
      count++;
      for (int i = eleStart; i <= eleEnd; ++i)
        element_tags.push_back(i);
    }

    else if (strcmp(argv_main[count], "-pattern") == 0) {
      if (count == argc_main - 1) {
        opserr << "WARNING eleLoad -pattern paramter missing required argument\n";
        return TCL_ERROR;
      }
      int ptag;
      if (Tcl_GetInt(interp, argv_main[++count], &ptag) != TCL_OK) {
        opserr << "WARNING eleLoad -pattern parameter expected integer\n";
        return TCL_ERROR;
      }
      explicitPatternPassed = true;
      loadPatternTag = ptag;
      count++;
    }

    else if (strcmp(argv_main[count], "-type") == 0) {
      argv[0] = argv_main[count++];
      argv[1] = argv_main[count++];
      if (count >= argc_main) {
        opserr << "WARNING eleLoad -type paramter missing required argument\n";
        return TCL_ERROR;
      }
      typeIndex = 0;

    } else {
      argv.push_back(argv_main[count++]);
    }
  }
  const int argc = static_cast<int>(argv.size());

  // If  -pattern  wasnt given explicitly, see if there is one
  // activated in the builder
  if (explicitPatternPassed == false) {
    LoadPattern *theTclLoadPattern = builder->getEnclosingPattern();
    if (theTclLoadPattern == nullptr) {
      opserr << "WARNING no current load pattern\n";
      return TCL_ERROR;

    } else {
      loadPatternTag = theTclLoadPattern->getTag();
    }
  }


  if (typeIndex == -1) {
    opserr << "WARNING missing required -type option"
           << "\n";
    return TCL_ERROR;
  }

  //
  // Create the load
  //
  count = typeIndex+1;

  if ((strcmp(argv[count], "-beamUniform") == 0) ||
      (strcmp(argv[count], "BeamUniform") == 0) ||
      (strcmp(argv[count], "beamUniform") == 0)) {
    //
    // see https://portwooddigital.com/2021/05/05/trapezoidal-beam-loads
    //
    count++;
    if (ndm == 2) {
      // wy wx a/L  b/L wyb wxb
      double wta;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wta) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      double waa = 0.0;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &waa) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
        return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
        return TCL_ERROR;
      }
      double wab = waa;
      double wtb = wta;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wtb) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wt for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wab) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wa for beamUniform \n";
        return TCL_ERROR;
      }

      for (int tag : element_tags) {
        ElementalLoad *theLoad = nullptr;
        if (aL > 0.0 || bL < 1.0 || wta != wtb || waa != wab)
          theLoad = new Beam2dPartialUniformLoad(eleLoadTag, wta, wtb, waa, wab,
                                                 aL, bL, tag);
        else
          theLoad = new Beam2dUniformLoad(eleLoadTag, wta, waa, tag);

        // add the load to the domain
        if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr
              << "WARNING eleLoad - could not add following load to domain:\n ";
          opserr << theLoad;
          delete theLoad;
          return TCL_ERROR;
        }
        eleLoadTag++;
      }

      return TCL_OK;

    } else if (ndm == 3) {
      // wy wz wx a/L b/L wyb wzb wxb
      double wy, wz;
      double wx  = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wy) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wy for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &wz) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wz for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wx) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wx for beamUniform \n";
        return TCL_ERROR;
      }
      double aL = 0.0;
      double bL = 1.0;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &aL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid aOverL for beamUniform \n";
        return TCL_ERROR;
      }
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &bL) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid bOverL for beamUniform \n";
        return TCL_ERROR;
      }
      
      //
      // Parse values at end "b"; when not supplied, set to value
      // at end "a"
      //

      double wyb = wy;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wyb) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wy for beamUniform \n";
        return TCL_ERROR;
      }

      double wzb = wz;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wzb) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wz for beamUniform \n";
        return TCL_ERROR;
      }

      double wxb = wx;
      count++;
      if (count < argc && Tcl_GetDouble(interp, argv[count], &wxb) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid wx for beamUniform \n";
        return TCL_ERROR;
      }

      for (int tag : element_tags) {
        ElementalLoad *theLoad = nullptr;
        if (aL > 0.0 || bL < 1.0)
          theLoad = new Beam3dPartialUniformLoad(eleLoadTag, wy, wz, wx, aL, bL, wyb, wzb, wxb, tag);
        else
          theLoad = new Beam3dUniformLoad(eleLoadTag, wy, wz, wx, tag);

        // add the load to the domain
        if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add following load to domain:\n ";
          delete theLoad;
          return TCL_ERROR;
        }
        eleLoadTag++;
      }
      return TCL_OK;
    }

    else {
      opserr << "WARNING eleLoad beamUniform currently only valid only for "
                "ndm=2 or 3\n";
      return TCL_ERROR;
    }

  } else if (strcmp(argv[count], "-beamPoint") == 0 ||
             strcmp(argv[count], "beamPoint") == 0) {
    count++;
    if (ndm == 2) {
      double P, x;
      double N = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &P) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid P for beamPoint\n";
        return TCL_ERROR;
      }
      if (count + 1 >= argc ||
          Tcl_GetDouble(interp, argv[count + 1], &x) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n";
        return TCL_ERROR;
      }
      if (count + 2 < argc &&
          Tcl_GetDouble(interp, argv[count + 2], &N) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid N for beamPoint\n";
        return TCL_ERROR;
      }

      if (x < 0.0 || x > 1.0) {
        opserr << "WARNING eleLoad - invalid xDivL of " << x;
        opserr << " for beamPoint (valid range [0.0, 1.0]\n";
        return TCL_ERROR;
      }

      for (int tag : element_tags) {
        ElementalLoad *theLoad = nullptr;
        theLoad = new Beam2dPointLoad(eleLoadTag, P, x, tag, N);

        // add the load to the domain
        if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr << "WARNING eleLoad - could not add load to domain:\n ";
          delete theLoad;
          return TCL_ERROR;
        }
        eleLoadTag++;
      }

      return TCL_OK;

    } else if (ndm == 3) {
      double Py, Pz, x;
      double N = 0.0;
      if (count >= argc || Tcl_GetDouble(interp, argv[count], &Py) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid Py for beamPoint\n";
        return TCL_ERROR;
      }
      if (count + 1 >= argc ||
          Tcl_GetDouble(interp, argv[count + 1], &Pz) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid Pz  for beamPoint\n";
        return TCL_ERROR;
      }
      if (count + 2 >= argc ||
          Tcl_GetDouble(interp, argv[count + 2], &x) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid xDivL for beamPoint\n";
        return TCL_ERROR;
      }
      if (count + 3 < argc &&
          Tcl_GetDouble(interp, argv[count + 3], &N) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid N for beamPoint\n";
        return TCL_ERROR;
      }

      if (x < 0.0 || x > 1.0) {
        opserr << "WARNING eleLoad - invalid xDivL of " << x;
        opserr << " for beamPoint (valid range [0.0, 1.0]\n";
        return TCL_ERROR;
      }

      for (int tag : element_tags) {
        ElementalLoad *theLoad = 
          new Beam3dPointLoad(eleLoadTag, Py, Pz, x, tag, N);

        // add the load to the domain
        if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
          opserr
              << "WARNING eleLoad - could not add following load to domain:\n ";
          opserr << theLoad;
          delete theLoad;
          return TCL_ERROR;
        }
        eleLoadTag++;
      }
      return 0;

    } else {
      opserr << "WARNING eleLoad beamPoint type currently only valid only for "
                "ndm=2 or 3\n";
      return TCL_ERROR;
    }
  }

  // Added Joey Yang UC Davis
  else if (strcmp(argv[count], "-BrickW") == 0) {

    for (int tag : element_tags) {
        ElementalLoad *theLoad = new BrickSelfWeight(eleLoadTag, tag);

      // add the load to the domain
      if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
        opserr
            << "WARNING eleLoad - could not add following load to domain:\n ";
        opserr << theLoad;
        delete theLoad;
        return TCL_ERROR;
      }
      eleLoadTag++;
    }
    return TCL_OK;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count], "-surfaceLoad") == 0) ||
           (strcmp(argv[count], "-SurfaceLoad") == 0)) {
    count++;
    for (int tag : element_tags) {
      ElementalLoad *theLoad = 
                new SurfaceLoader(eleLoadTag, tag);

      // add the load to the domain
      if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
        opserr
            << "WARNING eleLoad - could not add following load to domain:\n ";
        opserr << theLoad;
        delete theLoad;
        return TCL_ERROR;
      }
      eleLoadTag++;
    }
    return TCL_OK;
  }
  // Added: C.McGann, U.Washington
  else if ((strcmp(argv[count], "-selfWeight") == 0) ||
           (strcmp(argv[count], "-SelfWeight") == 0)) {
    count++;

    double xf = 0.0, yf = 0.0, zf = 0.0;
    if (Tcl_GetDouble(interp, argv[count], &xf) != TCL_OK) {
      opserr << "WARNING eleLoad - invalid xFactor " << argv[count]
             << " for -selfWeight\n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[count + 1], &yf) != TCL_OK) {
      opserr << "WARNING eleLoad - invalid yFactor " << argv[count + 1]
             << " for -selfWeight\n";
      return TCL_ERROR;
    }
    if (count + 2 < argc) { // adding to stop seg faults
      if (Tcl_GetDouble(interp, argv[count + 2], &zf) != TCL_OK) {
        opserr << "WARNING eleLoad - invalid zFactor " << argv[count + 2]
               << " for -selfWeight\n";
        return TCL_ERROR;
      }
    }

    for (std::size_t i=0; i< element_tags.size(); ++i) {
      ElementalLoad *theLoad = 
                new SelfWeight(eleLoadTag, xf, yf, zf, element_tags[i]);

      // add the load to the domain
      if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
        opserr << "WARNING eleLoad - could not add load to domain:\n ";
        delete theLoad;
        return TCL_ERROR;
      }
      eleLoadTag++;
    }
    return TCL_OK;
  }
  //--- Adding identifier for ThermalAction : [END] [SIF]---//

  //--- Adding tcl command for shell thermal action, 2013..[Begin]------
  else if (strcmp(argv[count], "-shellThermal") == 0) {
    count++;
    //so far three kinds of temperature distribution
    //(1) 9 temperature points, i.e. 8 layers
    //(2) 5 temperature points, i.e. 4 layers
    //(3) 2 temperature points, i.e. 1 layers: linear or uniform

    double t1, locY1, t2,
        locY2; //t3, locY3, t4, locY4, t5, locY5, t6, locY6, t7, locY7, t8, locY8, t9, locY9;
    // 9 temperature points are given,i.e. 8 layers are defined; Also the 9 corresponding vertical coordinate is given.
    // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
    //Start to add source file
    if (strcmp(argv[count], "-source") == 0) {
      if (strcmp(argv[count + 1], "-node") != 0) {
        count++;
#ifdef SIMULATION_INFO
        const char *pwd = getInterpPWD(interp);
        simulationInfo.addInputFile(argv[count], pwd);
#endif
        TimeSeries *theSeries =
            new PathTimeSeriesThermal(eleLoadTag, argv[count]);

        count++;

        double RcvLoc1, RcvLoc2;
        if (argc - count == 2) {

          if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid single loc  " << argv[count]
                   << " for -beamThermal\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid single loc  "
                   << argv[count + 1] << " for -beamThermal\n";
            return TCL_ERROR;
          }

        } else {
          opserr << "WARNING eleLoad - invalid input for -shellThermal\n";
        }

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = 
                    new ShellThermalAction(eleLoadTag, RcvLoc1, RcvLoc2,
                                           theSeries, element_tags[i]);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
      }
      // if not using nodal thermal action input
      else {
        for (int tag: element_tags) {
          ElementalLoad *theLoad = new ShellThermalAction(eleLoadTag, tag);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==
              false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
        return TCL_OK;
      } // end of <if(strcmp(argv[count+1],"-node") = 0)>
    }
    // end of the interface for importing temperature data from external file

    else {
      if (argc - count == 18) {
        double indata[18];
        double BufferData;

        for (int i = 0; i < 18; ++i) {
          if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid data " << argv[count]
                   << " for -beamThermal 3D\n";
            return TCL_ERROR;
          }
          indata[i] = BufferData;
          count++;
        }

        //temp1,loc1,temp2,loc2...temp9,loc9

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad =  new ShellThermalAction(
              eleLoadTag, indata[0], indata[1], indata[2], indata[3], indata[4],
              indata[5], indata[6], indata[7], indata[8], indata[9], indata[10],
              indata[11], indata[12], indata[13], indata[14], indata[15],
              indata[16], indata[17], element_tags[i]);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==
              false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
        return 0;
      }
      // 5 temperatures are given, i.e. 4 layers are defined.
      else if (argc - count == 10) {
        double indata[10];
        double BufferData;

        for (int i = 0; i < 10; ++i) {
          if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid data " << argv[count]
                   << " for -beamThermal 3D\n";
            return TCL_ERROR;
          }
          indata[i] = BufferData;
          count++;
        }

        //temp1,loc1,temp2,loc2...temp5,loc5

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = 
                    new ShellThermalAction(eleLoadTag, indata[0], indata[1],
                                           indata[2], indata[3], indata[4],
                                           indata[5], indata[6], indata[7],
                                           indata[8], indata[9], element_tags[i]);


          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==
              false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
        return 0;
      }

      // Two temperatures are given,
      // if the two temperatures are equal,i.e. uniform Temperature change in element
      // if the two temperatures are different,i.e. linear Temperature change in element
      else if (argc - count == 4) {
        if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid T1 " << argv[count]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 2], &t2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid T2 " << argv[count]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count + 3], &locY2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid LocY2 " << argv[count + 1]
                 << OpenSees::SignalMessageEnd;
          return TCL_ERROR;
        }

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = 
                    new ShellThermalAction(eleLoadTag, t1, locY1, t2, locY2,
                                           element_tags[i]);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==  false) {
            opserr << "WARNING eleLoad - could not add load to domain\n ";
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
        return 0;
      }
      //finish the temperature arguments
      else {
        opserr
            << "WARNING eleLoad -shellThermalLoad invalid number of "
               "temperature arguments,/n looking for 0, 1, 2 or 4 arguments.\n";
      }
    }
    //end of if(received argument is not "source" or direct temperature input)//Liming,2014
  }


  else if (strcmp(argv[count], "-ThermalWrapper") == 0 ||
           strcmp(argv[count], "-thermalWrapper") == 0) {

    count++;
    Vector loc = 0;

    ID NodalThermal = 0;
    int numNodal    = 0;

    if (strcmp(argv[count], "-nodeLoc") == 0 ||
        strcmp(argv[count], "nodeLoc") == 0) {
      count++;
      numNodal = (argc - count) / 2;
      loc.resize(numNodal);
      NodalThermal.resize(numNodal);

      for (int i = 0; i < numNodal; ++i) {

        double Dblloc;
        int NodalTtag;

        if (Tcl_GetDouble(interp, argv[count + i * 2 + 1], &Dblloc) != TCL_OK) {
          opserr << "WARNING NodalLoad - invalid loc  " << argv[count]
                 << " for NodalThermalAction\n";
          return TCL_ERROR;
        }

        if (Tcl_GetInt(interp, argv[count + 2 * i], &NodalTtag) != TCL_OK) {
          opserr << "WARNING invalid nodeId: " << argv[1];
          return TCL_ERROR;
        }

        loc(i)          = Dblloc;
        NodalThermal(i) = NodalTtag;
      }
      //end of for loop over numNodal
    }
    //end of nodelLoc

    else if (strcmp(argv[count], "-node") == 0 ||
             strcmp(argv[count], "node") == 0) {
      count++;
      numNodal = argc - count;
      NodalThermal.resize(numNodal);

      for (int i = 0; i < numNodal; ++i) {

        int NodalTtag;

        if (Tcl_GetInt(interp, argv[count + i], &NodalTtag) != TCL_OK) {
          opserr << "WARNING invalid nodeId " << argv[1];
          return TCL_ERROR;
        }

        NodalThermal(i) = NodalTtag;
      }
      //end of for loop over numNodal
    }
    //end of node tag

    //Obtain Pointers to NodalThermalAction;
    Node *theNode                         = 0;
    NodalThermalAction **theNodalThermals = 0;
    theNodalThermals                      = new NodalThermalAction *[numNodal];

    for (int i = 0; i < numNodal; ++i) {
      theNode             = domain->getNode(NodalThermal(i));
      theNodalThermals[i] = theNode->getNodalThermalActionPtr();
      if (theNodalThermals[i] == 0) {
        opserr << "WARNING:: An empty nodalThermalAction detected for "
                  "ThermalActionWrapper"
               << "\n";
        return TCL_ERROR;
      }
    }

    for (std::size_t i=0; i< element_tags.size(); ++i) {
      ElementalLoad *theLoad = nullptr;
      if (numNodal == 2) {
        theLoad =
            new ThermalActionWrapper(eleLoadTag, element_tags[i],
                                     theNodalThermals[0], theNodalThermals[1]);
      } else if (numNodal == 3) {
        theLoad = new ThermalActionWrapper(
            eleLoadTag, element_tags[i], theNodalThermals[0], theNodalThermals[1],
            theNodalThermals[2]);
      } else if (numNodal == 4) {
        theLoad = new ThermalActionWrapper(
            eleLoadTag, element_tags[i], theNodalThermals[0], theNodalThermals[1],
            theNodalThermals[2], theNodalThermals[3]);
      } else if (numNodal == 5) {
        theLoad = new ThermalActionWrapper(
            eleLoadTag, element_tags[i], theNodalThermals[0], theNodalThermals[1],
            theNodalThermals[2], theNodalThermals[3], theNodalThermals[4]);
      } else if (numNodal == 6) {
        theLoad = new ThermalActionWrapper(
            eleLoadTag, element_tags[i], theNodalThermals[0], theNodalThermals[1],
            theNodalThermals[2], theNodalThermals[3], theNodalThermals[4],
            theNodalThermals[5]);
      }


      if (loc != 0 && theLoad != nullptr)
        ((ThermalActionWrapper *)theLoad)->setRatios(loc);

      // add the load to the domain
      if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
        opserr
            << "WARNING eleLoad - could not add following load to domain:\n ";
        opserr << theLoad;
        delete theLoad;
        return TCL_ERROR;
      }
      eleLoadTag++;
    } //end of for loop
    return 0;
  }

  //-----------------Adding tcl command for beam thermal action(2D&3D), 2013..[Begin]---------------
  else if (strcmp(argv[count], "-beamThermal") == 0) {
    count++;
    //For two dimensional model
    if (ndm == 2) {
      // The thermal action can be defined with external data file, which is identified with "-source"
      // The external file itself can either be elemental data or nodal data, the latter is identified with "-node"
      if (strcmp(argv[count], "-source") == 0) {

        if (strcmp(argv[count + 1], "-node") == 0) {
          for (std::size_t i=0; i< element_tags.size(); ++i) {
            ElementalLoad *theLoad =         
                      new Beam2dThermalAction(eleLoadTag, element_tags[i]);

            // add the load to the domain
            if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                false) {
              opserr << "WARNING eleLoad - could not add following load to "
                        "domain:\n ";
              opserr << theLoad;
              delete theLoad;
              return TCL_ERROR;
            }
            eleLoadTag++;
          } //end of for loop
          return 0;

        }
        //end of <if(strcmp(argv[count+1],"-node") != 0)>
        else {
          count++;
#ifdef SIMULATION_INFO
          const char *pwd = getInterpPWD(interp);
          simulationInfo.addInputFile(argv[count], pwd);
#endif
          TimeSeries *theSeries =
              new PathTimeSeriesThermal(eleLoadTag, argv[count]);

          count++;
          Vector locs(9);
          //---------------for receiving 2 arguments
          if (argc - count == 2) {
            double RcvLoc1, RcvLoc2;
            if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  " << argv[count]
                     << " for -beamThermal\n";
              return TCL_ERROR;
            }
            if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  "
                     << argv[count + 1] << " for -beamThermal\n";
              return TCL_ERROR;
            }

            locs(0) = RcvLoc1;
            locs(8) = RcvLoc2;
            for (int i = 1; i < 8; ++i) {
              locs(i) = locs(0) - i * (locs(0) - locs(8)) / 8;
            }
          }
          //----------------for receiving 9 arguments
          else if (argc - count == 9) {

            int ArgStart = count;
            int ArgEnd   = argc;
            double data;

            if (ArgStart != ArgEnd) {
              for (int i = ArgStart; i < ArgEnd; ++i) {
                Tcl_GetDouble(interp, argv[i], &data);
                locs(i - ArgStart) = data;
              }
            }

          }
          //end of receiving 9 arguments
          else {
            opserr << "WARNING eleLoad - invalid input for -beamThermal\n";
          }

          for (std::size_t i=0; i< element_tags.size(); ++i) {
            ElementalLoad *theLoad =         
                      new Beam2dThermalAction(eleLoadTag, locs, theSeries,
                                              element_tags[i]);

            // add the load to the domain
            if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                false) {
              opserr << "WARNING eleLoad - could not add following load to "
                        "domain:\n ";
              opserr << theLoad;
              delete theLoad;
              return TCL_ERROR;
            }
            eleLoadTag++;
          } //end of for loop
          return 0;
        } // end of <if(strcmp(argv[count+1],"-node") = 0)>
        //--------------------------end for beam2DThermalAction with time series ----------------------------------------
      } else {
        //(1) 9 temperature points, i.e. 8 layers
        //(2) 5 temperature points, i.e. 4 layers
        //(3) 2 temperature points, i.e. 1 layers: linear or uniform
        //double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5, t6, locY6, t7, locY7, t8, locY8, t9, locY9;
        double Temp[9];
        double Loc[9];
        // 9 temperature points are given,i.e. 8 layers are defined; Also the 9 corresponding vertical coordinate is given.
        // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
        if (argc - count == 18) {
          double indata[18];
          double BufferData;

          for (int i = 0; i < 18; ++i) {
            if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid data " << argv[count]
                     << " for -beamThermal 3D\n";
              return TCL_ERROR;
            }
            indata[i] = BufferData;
            count++;
          }

          for (int i = 0; i < 9; ++i) {
            Temp[i] = indata[2 * i];
            Loc[i]  = indata[2 * i + 1];
          }

        }

        // 5 temperatures are given, i.e. 4 layers are defined.
        else if (argc - count == 10) {

          double indata[10];
          double BufferData;

          for (int i = 0; i < 10; ++i) {
            if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid data " << argv[count]
                     << " for -beamThermal 3D\n";
              return TCL_ERROR;
            }
            indata[i] = BufferData;
            count++;
          }

          Temp[0] = indata[0];
          Temp[2] = indata[2];
          Temp[4] = indata[4];
          Temp[6] = indata[6];
          Temp[8] = indata[8];
          Loc[0]  = indata[1];
          Loc[2]  = indata[3];
          Loc[4]  = indata[5];
          Loc[6]  = indata[7];
          Loc[8]  = indata[9];

          for (int i = 1; i < 5; ++i) {
            Temp[2 * i - 1] = (Temp[2 * i - 2] + Temp[2 * i]) / 2;
            Loc[2 * i - 1]  = (Loc[2 * i - 2] + Loc[2 * i]) / 2;
          }
        }
        //End for 5 inputs
        // two temperature is given,
        //if the two temperatures are equal,i.e. uniform Temperature change in element
        //if the two temperatures are different,i.e. linear Temperature change in element
        else if (argc - count == 4) {
          double indata[4];
          double BufferData;

          for (int i = 0; i < 4; ++i) {
            if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid data " << argv[count]
                     << " for -beamThermal 3D\n";
              return TCL_ERROR;
            }
            indata[i] = BufferData;
            count++;
          }

          Temp[0] = indata[0];
          Temp[8] = indata[2];
          Loc[0]  = indata[1];
          Loc[8]  = indata[3];
          for (int i = 1; i < 8; ++i) {
            Temp[i] = Temp[0] - i * (Temp[0] - Temp[8]) / 8;
            Loc[i]  = Loc[0] - i * (Loc[0] - Loc[8]) / 8;
          }
        }
        //end for 2 inputs

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad =         
                    new Beam2dThermalAction(
              eleLoadTag, Temp[0], Loc[0], Temp[1], Loc[1], Temp[2], Loc[2],
              Temp[3], Loc[3], Temp[4], Loc[4], Temp[5], Loc[5], Temp[6],
              Loc[6], Temp[7], Loc[7], Temp[8], Loc[8], element_tags[i]);


          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }

        return 0;
        //-End of Adding BeamThermalAction defined with direct input.
      }
      //end for no sourced pattern
    }
    // End of for the if (ndm==2)
    else if (ndm == 3) {
      //so far three kinds of temperature distribution
      double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5, t6, t7,
          locZ1, t8, t9, locZ2, t10, t11, locZ3, t12, t13, locZ4, t14, t15,
          locZ5;

      // the temperature at each fiber is obtained by interpolating of temperatures at the nearby temperature points.
      if (strcmp(argv[count], "-source") == 0) {
        count++;
        if (strcmp(argv[count], "-node") == 0) {
          for (std::size_t i=0; i< element_tags.size(); ++i) {
            ElementalLoad *theLoad =         
                      new Beam3dThermalAction(eleLoadTag, element_tags[i]);

            // add the load to the domain
            if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                false) {
              opserr << "WARNING eleLoad - could not add following load to "
                        "domain:\n ";
              opserr << theLoad;
              delete theLoad;
              return TCL_ERROR;
            }
            eleLoadTag++;
          } //end of loop tf all elements defined
          return 0;

        } //end for defing thermal action with nodal input
        else {
#ifdef SIMULATIN_INFO
          const char *pwd = getInterpPWD(interp);
          simulationInfo.addInputFile(argv[count], pwd);
#endif
          count++;
          // bool using2Ddata = false;

          double RcvLoc1, RcvLoc2, RcvLoc3, RcvLoc4;
          TimeSeries *theSeries;

          if (argc - count == 4) {
            theSeries =
                new PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 15);
            // using2Ddata = false;

            if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  " << argv[count]
                     << " for -beamThermal\n";
              return TCL_ERROR;
            }
            if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  "
                     << argv[count + 1] << " for -beamThermal\n";
              return TCL_ERROR;
            }
            if (Tcl_GetDouble(interp, argv[count + 2], &RcvLoc3) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  "
                     << argv[count + 2] << " for -beamThermal\n";
              return TCL_ERROR;
            }
            if (Tcl_GetDouble(interp, argv[count + 3], &RcvLoc4) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid single loc  "
                     << argv[count + 3] << " for -beamThermal\n";
              return TCL_ERROR;
            }

            // Create the loads
            for (int tag : element_tags) {

              ElementalLoad *theLoad =         
                  new Beam3dThermalAction(eleLoadTag, RcvLoc1, RcvLoc2, RcvLoc3,
                                          RcvLoc4, theSeries, tag);

              // add the load to the domain
              if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                  false) {
                opserr << "WARNING eleLoad - could not add following load to "
                          "domain:\n ";
                opserr << theLoad;
                delete theLoad;
                return TCL_ERROR;
              }
              eleLoadTag++;
            }
            return 0;

          }
          //end of defining 15 data points with external file
          else if (argc - count == 2 || argc - count == 9) {
            // for receiving data which has the similar structure as 2D beam section
            Vector locs(9);
            // using2Ddata = true;
            TimeSeries *theSeries =
                new PathTimeSeriesThermal(eleLoadTag, argv[count - 1], 9);
            if (argc - count == 2) {

              double RcvLoc1, RcvLoc2;
              if (Tcl_GetDouble(interp, argv[count], &RcvLoc1) != TCL_OK) {
                opserr << "WARNING eleLoad - invalid single loc  "
                       << argv[count] << " for -beamThermal\n";
                return TCL_ERROR;
              }
              if (Tcl_GetDouble(interp, argv[count + 1], &RcvLoc2) != TCL_OK) {
                opserr << "WARNING eleLoad - invalid single loc  "
                       << argv[count + 1] << " for -beamThermal\n";
                return TCL_ERROR;
              }

              locs(0) = RcvLoc1;
              locs(8) = RcvLoc2;
              for (int i = 1; i < 8; ++i) {
                locs(i) = locs(0) - i * (locs(0) - locs(8)) / 8;
              }
            }
            //end of receiving 2 data points
            else {
              double indata[9];
              double BufferData;
              for (int i = 0; i < 9; ++i) {
                if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
                  opserr << "WARNING eleLoad - invalid data " << argv[count]
                         << " for -beamThermal 3D\n";
                  return TCL_ERROR;
                }
                indata[i] = BufferData;
                count++;
              }

              for (int i = 0; i < 9; ++i) {
                locs(i) = indata[i];
              }
            }
            //end of receiving 9 data points

            for (std::size_t i=0; i< element_tags.size(); ++i) {
              ElementalLoad *theLoad =         
                        new Beam3dThermalAction(eleLoadTag, locs, theSeries,
                                                element_tags[i]);


              // add the load to the domain
              if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                  false) {
                opserr << "WARNING eleLoad - could not add following load to "
                          "domain:\n ";
                opserr << theLoad;
                delete theLoad;
                return TCL_ERROR;
              }
              eleLoadTag++;
            } //end of for loop
            return 0;
          }

          else {
            opserr << "WARNING eleLoad - invalid input for -beamThermal\n";
          }

        } //end for source beam element temperature data  -source -filePath

      }
      //--- end for importing temp data from external files --source -----

      else {

        //double t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5, locY5,
        //t6, t7, locZ1, t8, t9, locZ2, t10,t11, locZ3, t12, t13, locZ4, t14,t15, locZ5;

        if (argc - count == 25) {
          double indata[25];
          double BufferData;

          for (int i = 0; i < 25; ++i) {
            if (Tcl_GetDouble(interp, argv[count], &BufferData) != TCL_OK) {
              opserr << "WARNING eleLoad - invalid data " << argv[count]
                     << " for -beamThermal 3D\n";
              return TCL_ERROR;
            }
            indata[i] = BufferData;
            count++;
          }

          for (std::size_t i=0; i< element_tags.size(); ++i) {
            ElementalLoad *theLoad =         
                      new Beam3dThermalAction(
                eleLoadTag, indata[0], indata[1], indata[2], indata[3],
                indata[4], indata[5], indata[6], indata[7], indata[8],
                indata[9], indata[10], indata[11], indata[12], indata[13],
                indata[14], indata[15], indata[16], indata[17], indata[18],
                indata[19], indata[20], indata[21], indata[22], indata[23],
                indata[24], element_tags[i]);


            // add the load to the domain
            if (domain->addElementalLoad(theLoad, loadPatternTag) ==
                false) {
              opserr << "WARNING eleLoad - could not add following load to "
                        "domain:\n ";
              opserr << theLoad;
              delete theLoad;
              return TCL_ERROR;
            }
            eleLoadTag++;
          }
          return 0;
        } // end of  if (argc-count == 25){

        else if (argc - count == 4) {

          if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid T1 " << argv[count]
                   << " for -beamThermal\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1]
                   << " for -beamThermal\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[count + 2], &t5) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid T1 " << argv[count]
                   << " for -beamThermal\n";
            return TCL_ERROR;
          }
          if (Tcl_GetDouble(interp, argv[count + 3], &locY5) != TCL_OK) {
            opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1]
                   << " for -beamThermal\n";
            return TCL_ERROR;
          }

          locY2 = locY1 + (locY5 - locY1) / 4;
          locY3 = locY1 + (locY5 - locY1) / 2;
          locY4 = locY1 + 3 * (locY5 - locY1) / 4;
          t2    = t1 + (t5 - t1) / 4;
          t3    = t1 + (t5 - t1) / 2;
          t4    = t1 + 3 * (t5 - t1) / 4;
          locZ1 = locZ2 = locZ3 = locZ4 = locZ5 = 0;
          t6 = t7 = t8 = t9 = t10 = 0;
          t11 = t12 = t13 = t14 = t15 = 0;

          for (int tag : element_tags) {
            ElementalLoad *theLoad = new Beam3dThermalAction(
                eleLoadTag, t1, locY1, t2, locY2, t3, locY3, t4, locY4, t5,
                locY5, t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3, t12, t13,
                locZ4, t14, t15, locZ5, tag);


            // add the load to the domain
            if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
              opserr << "WARNING eleLoad - could not add following load to "
                        "domain:\n ";
              opserr << theLoad;
              delete theLoad;
              return TCL_ERROR;
            }
            eleLoadTag++;
          }
          return TCL_OK;
        }
        //end of  if (argc-count == 4){
			  else if (argc - count == 8) {

				  if (Tcl_GetDouble(interp, argv[count], &t1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 1], &locY1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocY1 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 2], &t5) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T5 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 3], &locY5) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocY5 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }

				  if (Tcl_GetDouble(interp, argv[count + 4], &t6) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T1 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 5], &locZ1) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocZ1 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 6], &t10) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid T10 " << argv[count] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }
				  if (Tcl_GetDouble(interp, argv[count + 7], &locZ5) != TCL_OK) {
					  opserr << "WARNING eleLoad - invalid LocZ5 " << argv[count + 1] << " for -beamThermal\n";
					  return TCL_ERROR;
				  }

				  locY2 = locY1 + (locY5 - locY1) / 4;
				  locY3 = locY1 + (locY5 - locY1) / 2;
				  locY4 = locY1 + 3 * (locY5 - locY1) / 4;
				  t2 = t1 + (t5 - t1) / 4;
				  t3 = t1 + (t5 - t1) / 2;
				  t4 = t1 + 3 * (t5 - t1) / 4;

				  locZ2 = locZ1 + (locZ5 - locZ1) / 4;
				  locZ3 = locZ1 + (locZ5 - locZ1) / 2;
				  locZ4 = locZ1 + 3 * (locZ5 - locZ1) / 4;
				  t11 = t6; t15 = t10;
				  t7 = t6 + (t10 - t6) / 4; t12 = t11 + (t15 - t11) / 4;
				  t8 = t6 + (t10 - t6) / 2; t13 = t11 + (t15 - t11) / 2;
				  t9 = t6 + 3*(t10 - t6) / 4; t14 = t11 + 3*(t15 - t11) / 4;

          for (int tag : element_tags) {
            ElementalLoad *theLoad = new Beam3dThermalAction(eleLoadTag,
						  t1, locY1, t2, locY2, t3, locY3, t4, locY4,
						  t5, locY5, t6, t7, locZ1, t8, t9, locZ2, t10, t11, locZ3,
						  t12, t13, locZ4, t14, t15, locZ5, tag);

					  // add the load to the domain
					  if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
						  opserr << "WARNING eleLoad - could not add following load to domain:\n ";
						  opserr << theLoad;
						  delete theLoad;
						  return TCL_ERROR;
					  }
					  eleLoadTag++;
				  }
				  return TCL_OK;
			  }
        else {
          opserr << "WARNING eleLoad Beam3dThermalAction: invalid number of "
                    "temperature arguments,/n looking for arguments for "
                    "Temperatures and coordinates.\n";
        }
      } //end for not sourced pattern
    }   //else if ndm=3
  }     //else if '-beamThermal'

  //--Adding identifier for ThermalAction:[END] [SIF]--//

  // Added by Scott R. Hamilton   - Stanford
  else if (strcmp(argv[count], "-beamTemp") == 0) {
    count++;
    if (ndm == 2) {
      double temp1, temp2, temp3, temp4;

      // Four temps given, Temp change at top node 1, bottom node 1, top node 2, bottom node 2.
      if (argc - count == 4) {
        if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Ttop1 " << argv[count]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count + 1], &temp2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot1 " << argv[count + 1]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 2], &temp3) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Ttop2 " << argv[count + 1]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }
        if (Tcl_GetDouble(interp, argv[count + 3], &temp4) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot2 " << argv[count + 1]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, temp3, temp4,
                                       element_tags[i]);


          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==
              false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }

        return TCL_OK;

      }
      // Two temps given, temp change at top, temp at bottom of element
      else if (argc - count == 2) {
        if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Ttop " << argv[count]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }

        if (Tcl_GetDouble(interp, argv[count + 1], &temp2) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot " << argv[count + 1]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = new Beam2dTempLoad(eleLoadTag, temp1, temp2, element_tags[i]);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) ==
              false) {
            opserr << "WARNING eleLoad - could not add following load to "
                      "domain:\n ";
            opserr << theLoad;
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }
      }

      // One twmp change give, uniform temp change in element
      else if (argc - count == 1) {
        if (Tcl_GetDouble(interp, argv[count], &temp1) != TCL_OK) {
          opserr << "WARNING eleLoad - invalid Tbot " << argv[count + 1]
                 << " for -beamTemp\n";
          return TCL_ERROR;
        }

        for (std::size_t i=0; i< element_tags.size(); ++i) {
          ElementalLoad *theLoad = new Beam2dTempLoad(eleLoadTag, temp1, element_tags[i]);

          // add the load to the domain
          if (domain->addElementalLoad(theLoad, loadPatternTag) == false) {
            opserr << "WARNING eleLoad - could not add load to "
                      "domain:\n ";
            delete theLoad;
            return TCL_ERROR;
          }
          eleLoadTag++;
        }

        return TCL_OK;

      }

      else {
        opserr << "WARNING eleLoad -beamTempLoad invalid number of temperature "
                  "arguments,/n looking for 0, 1, 2 or 4 arguments.\n";
      }
    } else {
      opserr << "WARNING eleLoad -beamTempLoad type currently only valid only "
                "for ndm=2\n";
      return TCL_ERROR;
    }
  }

  return TCL_ERROR;
}
