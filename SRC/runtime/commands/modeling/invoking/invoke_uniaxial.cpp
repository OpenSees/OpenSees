//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
// Description: This file contains the implementaion of functions
// used to directly invoke methods of a UniaxialMaterial from a
// Tcl interpreter.
//
// Written: cmp
//
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <tcl.h>
#include <G3_Logging.h>
#include <UniaxialMaterial.h>
#include <BasicModelBuilder.h>

typedef const char TCL_Char;

static Tcl_CmdProc TclCommand_setStrainUniaxialMaterial;
static Tcl_CmdProc TclCommand_commitState;
static Tcl_CmdProc TclCommand_getStressUniaxialMaterial;
static Tcl_CmdProc TclCommand_getTangUniaxialMaterial;
static Tcl_CmdProc TclCommand_integrateUniaxialMaterial;

const struct {const char*name; const Tcl_CmdProc*func;} command_table[] = {
  {"strain",      TclCommand_setStrainUniaxialMaterial },
  {"commit",      TclCommand_commitState               },
  {"stress",      TclCommand_getStressUniaxialMaterial },

  {"tangent",     TclCommand_getTangUniaxialMaterial   },
  {"stiffness",   TclCommand_getTangUniaxialMaterial   },

  {"integrate",   TclCommand_integrateUniaxialMaterial },
  // {"uniaxialTest",     TclCommand_setUniaxialMaterial}
};

//
// THE FUNCTIONS INVOKED BY THE INTERPRETER
//
int
TclCommand_useUniaxialMaterial(ClientData clientData,
                                          Tcl_Interp *interp, int argc,
                                          TCL_Char ** const argv)
{

  // Get the tag of the material to invoke
  int tag;
  if (Tcl_GetInt(interp, argv[2], &tag) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "could not read tag";
    return TCL_ERROR;
  }

  // Directly calling `getRegistryObject` as opposed to `getUniaxialMaterial`
  // prevents a copy from being created. This way, the state of the material
  // is preserved between invocations.
  UniaxialMaterial *theMaterial =
    ((BasicModelBuilder*)clientData)->getTypedObject<UniaxialMaterial>(tag);
  // UniaxialMaterial *theMaterial = 
  // ((BasicModelBuilder*)clientData)->getTypedObject<UniaxialMaterial>(argv[2]);

  if (theMaterial == nullptr) {
    opserr << G3_ERROR_PROMPT << "no material found with tag '" << tag << "'\n";
    return TCL_ERROR;

  } else {
    // theMaterial = theOrigMaterial->getCopy();
  }


  //
  // Add commands
  //
  const int ncmd = sizeof(command_table)/sizeof(command_table[0]);
  for (int i=0; i<ncmd; ++i)
    Tcl_CreateCommand(interp,
                      command_table[i].name,
                      command_table[i].func,
                      (ClientData)theMaterial, 
                      nullptr);
  //
  // Evaluate material commands
  //
  Tcl_Eval(interp, argv[3]);

  //
  // Clean up
  //
  Tcl_DeleteCommand(interp, "uniaxialTest");
  Tcl_DeleteCommand(interp, "strainUniaxialTest");
  Tcl_DeleteCommand(interp, "strain");
  Tcl_DeleteCommand(interp, "stress");
  Tcl_DeleteCommand(interp, "commit");
  Tcl_DeleteCommand(interp, "tangent");
  Tcl_DeleteCommand(interp, "stiffness");
  Tcl_DeleteCommand(interp, "integrate");

  return TCL_OK;
}

static int
TclCommand_setStrainUniaxialMaterial(ClientData clientData,
                                     Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;

  if (argc < 2) {
    opserr << G3_ERROR_PROMPT 
           << "bad arguments - want: strainUniaxialTest strain? <temp?>\n";
    return TCL_ERROR;
  }

  // get the tag from command line
  double strain;
  if (Tcl_GetDouble(interp, argv[1], &strain) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "could not read strain: strainUniaxialTest strain? "
              "<temp?>\n";
    return TCL_ERROR;
  }

  bool use_temp = false;
  bool commit = false;
  double temperature = 0.0;
  for (int i=2; i < argc; ++i) {
    if (strcmp(argv[i], "-commit")==0){
      commit = true;
    } else if (Tcl_GetDouble(interp, argv[2], &temperature) != TCL_OK) {
      opserr << G3_ERROR_PROMPT << "could not read strain: strainUniaxialTest strain? "
                "<temp?>\n";
      return TCL_ERROR;
    }
  }

  if (use_temp)
    theMaterial->setTrialStrain(strain, temperature, 0.0);
  else
    theMaterial->setTrialStrain(strain);

  if (commit)
    theMaterial->commitState();

  return TCL_OK;
}

int TclCommand_commitState(ClientData clientData, Tcl_Interp *interp, int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;
  theMaterial->commitState();
  return TCL_OK;
}

static int
TclCommand_getStressUniaxialMaterial(ClientData clientData,
                                     Tcl_Interp *interp,
                                     int argc, TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(theMaterial->getStress()));
  return TCL_OK;
}

static int
TclCommand_getTangUniaxialMaterial(ClientData clientData,
                                   Tcl_Interp *interp, int argc,
                                   TCL_Char ** const argv)
{
  assert(clientData != nullptr);
  UniaxialMaterial* theMaterial = (UniaxialMaterial*)clientData;

  double tangent = 0.0;

  tangent = theMaterial->getTangent();
  Tcl_SetObjResult(interp, Tcl_NewDoubleObj(tangent));
  return TCL_OK;

}

struct generalized_alpha {
  double alpha_m,
         alpha_f,
         beta,
         gamma;
};

static inline int
fsdof_integrate(struct generalized_alpha* conf,
    UniaxialMaterial& material, double M, double C,
    double scale, int n, double *p, double dt,
    double *response)
{ 
    const int max_iter   = 10;
    const double gamma   = conf->gamma;
    const double beta    = conf->beta;
    const double alpha_m = conf->alpha_m;
    const double alpha_f = conf->alpha_f;

    const double c1 = 1.0;
    const double c2 = gamma/(beta*dt);
    const double c3 = 1.0/(beta*dt*dt);

    const double a1 =     (1.0 -     gamma/beta);
    const double a2 =  dt*(1.0 - 0.5*gamma/beta);
    const double a3 = -1.0/(beta*dt);
    const double a4 =  1.0 - 0.5/beta;

    double  ua,
            va,
            aa,
            *u = &response[0],
            *v = &response[1],
            *a = &response[2];

    int i = 0;
    const int past = -3,
              pres =  0;

    // NOTE: The first row of the response array
    // is expected to be initialized!
    //  u[pres] = 0.0; v[pres] = 0.0;
    // const double M = material.getRho()   ;
    // const double C = material.getDampTangent() ;
    a[pres] = (p[i] - C*v[pres] - material.getStress())/M;

    enum {Accel, Veloc, Displ} init = Accel, 
                               form = Veloc;

    for (i = 1; i < n; ++i) {
      u += 3; v += 3; a += 3;

      // Predictor step
      u[pres] = u[past];
      v[pres] = a1*v[past] + a2*a[past];
      a[pres] = a4*a[past] + a3*v[past];

      ua = (1.0 - alpha_f)*u[past] + alpha_f*u[pres];
      va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
      aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];

      double du = 0.0;
      switch (init) {
        case Displ:   du = 0.0; break;
        case Veloc:   du = 0.0; break;
        case Accel:   du = 0.0; break;
      }

      for (int j=0; j<max_iter; j++) {
        //
        u[pres] += du;
        v[pres] += c2*du;
        a[pres] += c3*du;

        ua = (1.0 - alpha_f)*u[past] + alpha_f*u[pres];
        va = (1.0 - alpha_f)*v[past] + alpha_f*v[pres];
        aa = (1.0 - alpha_m)*a[past] + alpha_m*a[pres];

        material.setTrialStrain(ua); //u[pres]);
        // STATE DETERMINATION
        // double C = material.getDampTangent();
        double ri = (scale*p[i] - C*va - M*aa - material.getStress());
        double K  = material.getTangent();
        double ki = alpha_f*c1*K + alpha_f*c2*C + alpha_m*c3*M;

        // SOLVE
        du = ri / ki;

        if ((j > 0) && (fabs(ri) < 1e-7)) {
          material.commitState();
          break;
        }
      }
    }
    return 1;
}

static int
TclCommand_integrateUniaxialMaterial(ClientData clientData,
    Tcl_Interp* interp,
    int argc, const char** const argv)
{
  int n;
  double dt;
  double M=1.0, C=0.0;
  const char** str_values;
  struct generalized_alpha conf = {1.0, 1.0, 0.25, 0.5};
  UniaxialMaterial* material = static_cast<UniaxialMaterial*>(clientData);
  
  if (Tcl_GetDouble(interp, argv[1], &dt) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "problem reading time step, got '" << argv[1] << "'\n";
    return TCL_ERROR;
  }

  if (Tcl_SplitList(interp, argv[2], &n, &str_values) != TCL_OK) {
    opserr << G3_ERROR_PROMPT << "problem splitting path list " << argv[2] << "\n";
    return TCL_ERROR;
  }

  int argi = 3;
  while (argi < argc) {
    if (strcmp(argv[argi], "-alphaf") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &conf.alpha_f) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
    else if (strcmp(argv[argi], "-alpham") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &conf.alpha_m) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
    else if (strcmp(argv[argi], "-beta") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &conf.beta) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
    else if (strcmp(argv[argi], "-gamma") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &conf.gamma) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
    else if (strcmp(argv[argi], "-mass") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &M) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
    else if (strcmp(argv[argi], "-damp") == 0) {
      if (Tcl_GetDouble(interp, argv[1+argi], &C) != TCL_OK) {
        opserr << G3_ERROR_PROMPT << "problem reading " << argv[argi] << ", got '" << argv[argi+1] << "'\n";
        return TCL_ERROR;
      }
      argi += 2;
    }
  }
  
  // copy load vector strings into a double array
  double *load = new double[n];
  for (int i=0; i<n; ++i)
    Tcl_GetDouble(interp, str_values[i], &load[i]);

  auto  resp = new double[n][3];
  resp[0][0] = 0.0;
  resp[0][1] = 0.0;
  resp[0][2] = 0.0;

  fsdof_integrate(&conf, *material, M, C, 1.0, n, load, dt, (double*)resp);

  Tcl_Obj *displ = Tcl_NewListObj(0, NULL);
  for (int i=0; i<n; ++i)
    Tcl_ListObjAppendElement(interp, displ, Tcl_NewDoubleObj(resp[i][0]));

  Tcl_SetObjResult(interp, displ);

  // TODO: check that this frees everything
  Tcl_Free((char*)str_values);
  delete[] resp;

  return TCL_OK;
}
