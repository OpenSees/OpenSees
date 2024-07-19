//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Description: This file implements selection of an integrator object.
//
#include <G3_Logging.h>
#include "integrator.h"
#include <assert.h>
#include <tcl.h>
#include <runtimeAPI.h>
#include <Domain.h>
#include <Node.h>
#include "BasicAnalysisBuilder.h"

// integrators
#include <LoadControl.h>
#include <ArcLength1.h>
#include <DisplacementControl.h>

#include <Newmark.h>
#include <TRBDF2.h>
#include <TRBDF3.h>
#include <Houbolt.h>
#include <ParkLMS3.h>
#include <BackwardEuler.h>

class G3_Runtime;

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);


StaticIntegrator*
G3Parse_newHSIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
StaticIntegrator*
G3Parse_newLoadControl(ClientData, Tcl_Interp*, int argc, TCL_Char *argv[]);
StaticIntegrator*
G3Parse_newEQPathIntegrator(ClientData, Tcl_Interp*, int argc, TCL_Char ** const);
StaticIntegrator*
G3Parse_newArcLengthIntegrator(ClientData, Tcl_Interp*, int argc, TCL_Char ** const);
StaticIntegrator*
G3Parse_newMinUnbalDispNormIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
StaticIntegrator*
G3Parse_newDisplacementControlIntegrator(ClientData, Tcl_Interp*, int, TCL_Char** const);
StaticIntegrator*
G3Parse_newStaticIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);

TransientIntegrator*
G3Parse_newNewmark1Integrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);
TransientIntegrator*
G3Parse_newNewmarkIntegrator(ClientData, Tcl_Interp*, int, TCL_Char** const);
TransientIntegrator*
G3Parse_newTransientIntegrator(ClientData, Tcl_Interp*, int, TCL_Char ** const);

//
// command invoked to select and construct an integrator
//
int
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  // make sure at least one other argument to select integrator
  if (argc < 2) {
    opserr << "WARNING need to specify an Integrator type \n";
    return TCL_ERROR;
  }

  StaticIntegrator* static_integrator = 
    G3Parse_newStaticIntegrator(clientData, interp, argc, argv);

  TransientIntegrator* transient_integrator = 
    G3Parse_newTransientIntegrator(clientData, interp, argc, argv);

  if (static_integrator != nullptr) {
    opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
    static_integrator->Print(opsdbg);
    builder->set(*static_integrator);

  } else if (transient_integrator != nullptr) {
    opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
    transient_integrator->Print(opsdbg);
    builder->set(*transient_integrator);
  }
  else
    return TCL_ERROR;

  return TCL_OK;
}

StaticIntegrator*
G3Parse_newStaticIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{

  StaticIntegrator* the_static_integrator = nullptr;

  // check argv[1] for type of integrator and create the object

  if (strcmp(argv[1], "LoadControl") == 0) {
    the_static_integrator = G3Parse_newLoadControl(clientData, interp, argc, argv);
  } 

  else if (strcmp(argv[1], "EQPath") == 0) {
    the_static_integrator = G3Parse_newEQPathIntegrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "MinUnbalDispNorm") == 0) {
    the_static_integrator = G3Parse_newMinUnbalDispNormIntegrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "DisplacementControl") == 0) {
    the_static_integrator = G3Parse_newDisplacementControlIntegrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ArcLength") == 0) {
    the_static_integrator = G3Parse_newArcLengthIntegrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "ArcLength1") == 0) {
    double arcLength;
    double alpha;
    if (argc != 4) {
      opserr << "WARNING integrator ArcLength1 arcLength alpha \n";
      return nullptr;
    }
    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
      return nullptr;
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)
      return nullptr;
    the_static_integrator = new ArcLength1(arcLength, alpha);
  }


//else if (strcmp(argv[1], "HSConstraint") == 0) {
//  the_static_integrator = G3Parse_newHSIntegrator(clientData, interp, argc, argv);
//}

  return the_static_integrator;
}

TransientIntegrator*
G3Parse_newTransientIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  G3_Runtime* rt = G3_getRuntime(interp);
  TransientIntegrator * theTransientIntegrator = nullptr;

  if ((strcmp(argv[1], "TRBDF2") == 0) ||
      (strcmp(argv[1], "Bathe") == 0)) {
    theTransientIntegrator = new TRBDF2();
  }

  else if ((strcmp(argv[1], "TRBDF3") == 0) ||
           (strcmp(argv[1], "Bathe3") == 0)) {
    theTransientIntegrator = new TRBDF3();
  }

  else if (strcmp(argv[1], "Houbolt") == 0) {
    theTransientIntegrator = new Houbolt();
  }

  /*else if (strcmp(argv[1],"ParkLMS3") == 0) {
      theTransientIntegrator = new ParkLMS3();
  }*/

  else if (strcmp(argv[1], "BackwardEuler") == 0) {
    int optn = 0;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &optn) != TCL_OK) {
        opserr << "WARNING integrator BackwardEuler <option> - undefined "
                  "option specified\n";
        return nullptr;
      }
    }
    theTransientIntegrator = new BackwardEuler(optn);
  }

  else if (strcmp(argv[1], "Newmark") == 0) {
    theTransientIntegrator = (TransientIntegrator *)G3Parse_newNewmarkIntegrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "GimmeMCK") == 0 || strcmp(argv[1], "MCK") == 0 || strcmp(argv[1], "ZZTop") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GimmeMCK(rt, argc, argv);
  } 
  else if (strcmp(argv[1], "NewmarkExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkExplicit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrReduct(rt, argc, argv);
  }

  else if (strcmp(argv[1], "NewmarkHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSIncrLimit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "NewmarkHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_NewmarkHSFixedNumIter(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHT") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHT_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHT_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralized_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTExplicit_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTGeneralizedExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTGeneralizedExplicit_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSIncrLimit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrLimit_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSIncrReduct_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSIncrReduct_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter(rt, argc, argv);
  }

  else if (strcmp(argv[1], "HHTHSFixedNumIter_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_HHTHSFixedNumIter_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "GeneralizedAlpha") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_GeneralizedAlpha(rt, argc, argv);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "KRAlphaExplicit_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_KRAlphaExplicit_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "AlphaOS") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS(rt, argc, argv);
  }

  else if (strcmp(argv[1], "AlphaOS_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOS_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized(rt, argc, argv);
  }

  else if (strcmp(argv[1], "AlphaOSGeneralized_TP") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_AlphaOSGeneralized_TP(rt, argc, argv);
  }

  else if (strcmp(argv[1], "Collocation") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_Collocation(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CollocationHSIncrReduct") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrReduct(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CollocationHSIncrLimit") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSIncrLimit(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CollocationHSFixedNumIter") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CollocationHSFixedNumIter(rt, argc, argv);
  }

  else if (strcmp(argv[1], "Newmark1") == 0) {
    theTransientIntegrator = G3Parse_newNewmark1Integrator(clientData, interp, argc, argv);
  }

  else if (strcmp(argv[1], "WilsonTheta") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_WilsonTheta(rt, argc, argv);
  }

  else if (strcmp(argv[1], "ExplicitDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_ExplicitDifference(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CentralDifference") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifference(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CentralDifferenceAlternative") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceAlternative(rt, argc, argv);
  }

  else if (strcmp(argv[1], "CentralDifferenceNoDamping") == 0) {
    theTransientIntegrator = (TransientIntegrator *)OPS_CentralDifferenceNoDamping(rt, argc, argv);
  }

  return theTransientIntegrator;
}

#include <HSConstraint.h>
StaticIntegrator*
G3Parse_newHSIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
  double arcLength, psi_u, psi_f, u_ref;

  if (argc < 3) {
    opserr << "WARNING integrator HSConstraint <arcLength> <psi_u> <psi_f> "
              "<u_ref> \n";
    return nullptr;
  }
  if (argc >= 3 && Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
    return nullptr;
  if (argc >= 4 && Tcl_GetDouble(interp, argv[3], &psi_u) != TCL_OK)
    return nullptr;
  if (argc >= 5 && Tcl_GetDouble(interp, argv[4], &psi_f) != TCL_OK)
    return nullptr;
  if (argc == 6 && Tcl_GetDouble(interp, argv[5], &u_ref) != TCL_OK)
    return nullptr;

  switch (argc) {
  case 3:
    return new HSConstraint(arcLength);
  case 4:
    return new HSConstraint(arcLength, psi_u);
  case 5:
    return new HSConstraint(arcLength, psi_u, psi_f);
  case 6:
    return new HSConstraint(arcLength, psi_u, psi_f, u_ref);
  default:
    return nullptr;
  }
}

StaticIntegrator*
G3Parse_newLoadControl(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator LoadControl dlam <Jd "
                "dlamMin dlamMax>\n";
      return nullptr;
    }
    if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)
      return nullptr;
    if (argc > 5) {
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)
        return nullptr;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    return new LoadControl(dLambda, numIter, minIncr, maxIncr);
}

#include <EQPath.h>
StaticIntegrator*
G3Parse_newEQPathIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, const char *argv[])
{
    double arcLength;
    int type;
    if (argc != 4) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << "REFS : \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849         \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879         \n";
      return nullptr;
    }

    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.849         \n";
      opserr << " https://doi.org/10.12989/sem.2013.48.6.879         \n";
      return nullptr;
    }

    if (Tcl_GetInt(interp, argv[3], &type) != TCL_OK) {
      opserr << "WARNING integrator EQPath $arc_length $type \n";
      opserr << "$type = 1 Minimum Residual Displacement \n";
      opserr << "$type = 2 Normal Plain \n";
      opserr << "$type = 3 Update Normal Plain \n";
      opserr << "$type = 4 Cylindrical Arc-Length \n";
      return nullptr;
    }

    return new EQPath(arcLength, type);
}

#include <ArcLength.h>
StaticIntegrator *
G3Parse_newArcLengthIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  double arcLength;
  double alpha   = 1.0;
  double expon   = 0.0;
  bool   use_det = false;
  double numIter = 5;
  ReferencePattern reference = ReferencePattern::Full;
  
  // parser variables
  bool got_alpha = false;
  
  // Begin parse
  if (argc < 3) {
    opserr << G3_ERROR_PROMPT << "integrator ArcLength arcLength alpha \n";
    return nullptr;
  }

  if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
    return nullptr;

  for (int i=3; i<argc; ++i) {
    if (strcmp(argv[i], "-det") == 0) {
      use_det = true;
    }
    else if ((strcmp(argv[i], "-exp") == 0)) {
      if (++i >= argc || Tcl_GetDouble(interp, argv[i], &expon) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read expon\n";
        return nullptr;
      }
    }
    else if ((strcmp(argv[i], "-j") == 0)) {
      if (++i >= argc || Tcl_GetDouble(interp, argv[i], &numIter) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read iter\n";
        return nullptr;
      }
    }
    else if ((strcmp(argv[i], "-reference") == 0)) {
      if (++i < argc && strcmp(argv[i], "point")==0) {
        reference = ReferencePattern::Point;
      }
    } else if (!got_alpha) {
      if (Tcl_GetDouble(interp, argv[i], &alpha) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "failed to read alpha, got " << argv[i] << "\n";
        return nullptr;
      }
      got_alpha = true;
    }

  }

  return new ArcLength(arcLength, alpha, numIter, expon, use_det, reference);
}


#include <MinUnbalDispNorm.h>
StaticIntegrator*
G3Parse_newMinUnbalDispNormIntegrator(ClientData clientData, Tcl_Interp* interp, int argc, TCL_Char ** const argv)
{
    if (argc < 3) {
      opserr << "WARNING integrator MinUnbalDispNorm lambda11 <Jd minLambda1j "
                "maxLambda1j>\n";
      return nullptr;
    }

    double lambda11;
    if (Tcl_GetDouble(interp, argv[2], &lambda11) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "expected float for lambda11 but got " << argv[2] << "\n";
      return nullptr;
    }

    // 
    // Optional Arguments
    //
    enum {
      NumIter = 1<<1,
      MinLamb = 1<<2,
      MaxLamb = 1<<3
    };
    int recvd = 0;

    // set defaults
    double minlambda = lambda11;
    double maxlambda = lambda11;
    int numIter = 1;
    int signFirstStepMethod = MinUnbalDispNorm::SIGN_LAST_STEP;

    for (int i=3; i < argc; ++i) {
      if ((strcmp(argv[i], "-determinant") == 0) || 
          (strcmp(argv[i], "-det") == 0)) {
          signFirstStepMethod = MinUnbalDispNorm::CHANGE_DETERMINANT;

      } else if ((recvd&NumIter) == 0) {
        if (Tcl_GetInt(interp, argv[i], &numIter) != TCL_OK)
          return nullptr;
        recvd |= NumIter;

      } else if ((recvd&MinLamb) == 0) {
        if (Tcl_GetDouble(interp, argv[i], &minlambda) != TCL_OK)
          return nullptr;

      } else if ((recvd&MaxLamb) == 0) {
        if (Tcl_GetDouble(interp, argv[i], &maxlambda) != TCL_OK)
          return nullptr;
      }
    }

    return new MinUnbalDispNorm(lambda11, numIter, minlambda,
                                maxlambda, signFirstStepMethod);
}

StaticIntegrator*
G3Parse_newDisplacementControlIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
    BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;
    Domain *domain = builder->getDomain();

    int node, dof, numIter;
    double increment, minIncr, maxIncr;

    if (argc < 5) {
      opserr << "WARNING integrator DisplacementControl node dof dU \n";
      opserr << "<Jd minIncrement maxIncrement>\n";
      return nullptr;
    }
    int tangFlag = 0;

    if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK) {
      return nullptr;
    }
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)
      return nullptr;
    if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)
      return nullptr;

    if (argc == 6 || argc == 9) {
      if (argc == 6) {
        if (strcmp(argv[5], "-initial") == 0)
          tangFlag = 1;
      } else if (strcmp(argv[8], "-initial") == 0)
        tangFlag = 1;
    }

    if (argc > 6) {
      if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK) {
        opserr << "WARNING failed to read numIter\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK) {
        opserr << "WARNING failed to read minIncr\n";
        return nullptr;
      }
      if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK) {
        opserr << "WARNING failed to read maxIncr\n";
        return nullptr;
      }
    } else {
      minIncr = increment;
      maxIncr = increment;
      numIter = 1;
    }

#ifdef _PARALLEL_PROCESSING
    the_static_integrator = new DistributedDisplacementControl(
        node, dof - 1, increment, numIter, minIncr, maxIncr);
#else
    Node *theNode = domain->getNode(node);
    if (theNode == nullptr) {
      opserr << "WARNING Node " << node << " does not exist\n";
      return nullptr;
    }

    int numDOF = theNode->getNumberDOF();
    if (dof <= 0 || dof > numDOF) {
      opserr << "WARNING invalid dof " << dof << "\n";
      return nullptr;
    }

    return
        new DisplacementControl(node, dof-1, increment, domain, numIter,
                                minIncr, maxIncr, tangFlag);

#endif
}

#if 0
StaticIntegrator*
G3Parse_newStagedLoadControlIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
    double dLambda;
    double minIncr, maxIncr;
    int numIter;
    if (argc < 3) {
      opserr << "WARNING incorrect # args - integrator StagedLoadControl dlam "
                "<Jd dlamMin dlamMax>\n";
      return nullptr;
    }
    if (Tcl_GetDouble(interp, argv[2], &dLambda) != TCL_OK)
      return nullptr;
    if (argc > 5) {
      if (Tcl_GetInt(interp, argv[3], &numIter) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[4], &minIncr) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[5], &maxIncr) != TCL_OK)
        return nullptr;
    } else {
      minIncr = dLambda;
      maxIncr = dLambda;
      numIter = 1;
    }
    return  new StagedLoadControl(dLambda, numIter, minIncr, maxIncr);
}
#endif


#ifdef _PARALLEL_INTERPRETERS
  else if ((strcmp(argv[1], "ParallelDisplacementControl") == 0) ||
           (strcmp(argv[1], "ParallelDisplacementControl") == 0)) {
    int node;
    int dof;
    double increment, minIncr, maxIncr;
    int numIter;
    if (argc < 5) {
      opserr << "WARNING integrator DisplacementControl node dof dU \n";
      opserr << "<Jd minIncrement maxIncrement>\n";
      return nullptr;
    }
    if (Tcl_GetInt(interp, argv[2], &node) != TCL_OK)
      return nullptr;
    if (Tcl_GetInt(interp, argv[3], &dof) != TCL_OK)
      return nullptr;
    if (Tcl_GetDouble(interp, argv[4], &increment) != TCL_OK)
      return nullptr;
    if (argc > 7) {
      if (Tcl_GetInt(interp, argv[5], &numIter) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[6], &minIncr) != TCL_OK)
        return nullptr;
      if (Tcl_GetDouble(interp, argv[7], &maxIncr) != TCL_OK)
        return nullptr;
    } else {
      minIncr = increment;
      maxIncr = increment;
      numIter = 1;
    }

    DistributedDisplacementControl *theDDC = new DistributedDisplacementControl(
        node, dof - 1, increment, numIter, minIncr, maxIncr);

    theDDC->setProcessID(OPS_rank);
    theDDC->setChannels(numChannels, theChannels);
    the_static_integrator = theDDC;
  }
#endif
