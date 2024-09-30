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
#include <BackwardEuler.h>

extern "C" int OPS_ResetInputNoBuilder(ClientData clientData,
                                       Tcl_Interp *interp, int cArg, int mArg,
                                       TCL_Char ** const argv, Domain *domain);


Tcl_CmdProc TclCommand_newStaticIntegrator;
Tcl_CmdProc TclCommand_newTransientIntegrator;

//
// Command invoked to select and construct an integrator
//
int
specifyIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  OPS_ResetInputNoBuilder(clientData, interp, 2, argc, argv, nullptr);

  // Make sure at least one other argument to select integrator
  if (argc < 2) {
    opserr << "WARNING need to specify an Integrator type \n";
    return TCL_ERROR;
  }

  if (TclCommand_newStaticIntegrator(clientData, interp, argc, argv) == TCL_OK)
    return TCL_OK;

  else if (TclCommand_newTransientIntegrator(clientData, interp, argc, argv) == TCL_OK)
    return TCL_OK;

  else
    return TCL_ERROR;

}

int
TclCommand_newStaticIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  auto tcl_cmd = StaticIntegratorLibrary.find(std::string(argv[1]));
  if (tcl_cmd != StaticIntegratorLibrary.end())
    return (*tcl_cmd->second)(clientData, interp, argc, &argv[0]);

  //
  StaticIntegrator* theStaticIntegrator = nullptr;

  // Check argv[1] for type of integrator and create the object

  if (strcmp(argv[1], "ArcLength1") == 0) {
    double arcLength;
    double alpha;
    if (argc != 4) {
      opserr << "WARNING integrator ArcLength1 arcLength alpha \n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[2], &arcLength) != TCL_OK)
      return TCL_ERROR;
    if (Tcl_GetDouble(interp, argv[3], &alpha) != TCL_OK)
      return TCL_ERROR;
    theStaticIntegrator = new ArcLength1(arcLength, alpha);
  }

  //
  // Done parsing; set the integrator
  //
  if (theStaticIntegrator != nullptr) {
    opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
    theStaticIntegrator->Print(opsdbg);
    builder->set(*theStaticIntegrator);
    return TCL_OK;
  } else
    return TCL_ERROR;
}

int
TclCommand_newTransientIntegrator(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  BasicAnalysisBuilder *builder = (BasicAnalysisBuilder*)clientData;

  auto tcl_cmd = TransientIntegratorLibrary.find(std::string(argv[1]));
  if (tcl_cmd != TransientIntegratorLibrary.end())
    return (*tcl_cmd->second)(clientData, interp, argc, &argv[0]);


  TransientIntegrator * theTransientIntegrator = nullptr;


  if (strcmp(argv[1], "BackwardEuler") == 0) {
    int optn = 0;
    if (argc == 3) {
      if (Tcl_GetInt(interp, argv[2], &optn) != TCL_OK) {
        opserr << "WARNING integrator BackwardEuler <option> - undefined "
                  "option specified\n";
        return TCL_ERROR;
      }
    }
    theTransientIntegrator = new BackwardEuler(optn);
  }

  //
  //
  //
  if (theTransientIntegrator == nullptr) {
    return TCL_ERROR;
  }

  opsdbg << G3_DEBUG_PROMPT << "Set integrator to \n";
  theTransientIntegrator->Print(opsdbg);
  builder->set(*theTransientIntegrator);
  return TCL_OK;
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
    theStaticIntegrator = new DistributedDisplacementControl(
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
    theStaticIntegrator = theDDC;
  }
#endif
