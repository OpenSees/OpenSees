//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//
#include <map>
#include <vector>
#include <assert.h>
#include <elementAPI.h>
#include <stdlib.h>
#include <packages.h>
#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <runtimeAPI.h>
#include <G3_Runtime.h>
#include <G3_Logging.h>
#include <BasicModelBuilder.h>

#include <UniaxialMaterial.h>
#include <NDMaterial.h>
#include <FrameSection.h>
#include <SectionForceDeformation.h>
#include <CrdTransf.h>
#include <FrictionModel.h>

#include <DirectIntegrationAnalysis.h>
#include <StaticAnalysis.h>

#include <TimeSeries.h>

static Tcl_Interp *theInterp       = nullptr;
static TCL_Char **currentArgv      = nullptr;
static int currentArg = 0;
static int maxArg     = 0;


extern const char *getInterpPWD(Tcl_Interp *interp);

struct cmp_str {
  bool
  operator()(const char *a, const char *b)
  {
    return strcmp(a, b) < 0;
  }
};


extern "C" int
OPS_Error(const char *errorMessage, int length)
{
  opserr << errorMessage;
  opserr << endln;

  return 0;
}

extern "C" int
OPS_GetNumRemainingInputArgs()
{
  return maxArg - currentArg;
}

extern "C" int
OPS_ResetCurrentInputArg(int cArg)
{
  if (cArg < 0)
    currentArg += cArg;
  else
    currentArg = cArg;

  return 0;
}


// extern "C"
#if 0
int
OPS_ResetInput(ClientData clientData, Tcl_Interp *interp, int cArg, int mArg,
               TCL_Char ** const argv, void*, void*)
{
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;
  return 0;
}
#endif


extern "C" int
OPS_ResetInputNoBuilder(ClientData clientData, Tcl_Interp *interp, int cArg,
                        int mArg, TCL_Char ** const argv, Domain *domain)
{
  currentArgv = argv;
  currentArg = cArg;
  maxArg = mArg;
  return 0;
}

extern "C" int
OPS_GetIntInput(int *numData, int *data)
{
  int size = *numData;

  for (int i = 0; i < size; ++i) {
    if ((currentArg >= maxArg) ||
        (Tcl_GetInt(theInterp, currentArgv[currentArg], &data[i]) != TCL_OK)) {
      return -1;
    } else
      currentArg++;
  }

  return 0;
}

extern "C" int
OPS_GetDoubleInput(int *numData, double *data)
{
  int size = *numData;
  for (int i = 0; i < size; ++i) {
    if ((currentArg >= maxArg) ||
        (Tcl_GetDouble(theInterp, currentArgv[currentArg], &data[i]) != TCL_OK)) {
      return -1;
    } else
      currentArg++;
  }
  return 0;
}

extern "C" const char *
OPS_GetString(void)
{
  const char *res = 0;
  if (currentArg >= maxArg) {
    return res;
  }
  res = currentArgv[currentArg];

  currentArg++;

  return res;
}

extern "C" const char *
OPS_GetStringFromAll(char *buffer, int len)
{
  return OPS_GetString();
}

extern "C" int
OPS_SetString(const char *str)
{
  Tcl_SetResult(theInterp, (char *)str, TCL_VOLATILE);
  return 0;
}

int
OPS_GetStringCopy(char **arrayData)
{
  if (currentArg >= maxArg) {
    opserr << "OPS_GetStringInput -- error reading " << currentArg << endln;
    return -1;
  }
  char *newData = new char[strlen(currentArgv[currentArg]) + 1];
  strcpy(newData, currentArgv[currentArg]);
  *arrayData = newData;
  currentArg++;

  return 0;
}


extern "C"
int OPS_GetDoubleListInput(int* size, Vector* data)
{
    TCL_Char** strings;

    if (Tcl_SplitList(theInterp, currentArgv[currentArg],
        size, &strings) != TCL_OK) {
        opserr << "ERROR problem splitting list " << currentArgv[currentArg] << " \n";
        return -1;
    }

    data->resize(*size);
    for (int i = 0; i < *size; i++) {
        double value;
        if (Tcl_GetDouble(theInterp, strings[i], &value) != TCL_OK) {
            opserr << "ERROR problem reading data value " << strings[i] << " \n";
            // free up the array of strings .. see tcl man pages as to why
            Tcl_Free((char*)strings);
            return -1;
        }
        (*data)(i) = value;
    }
    // free up the array of strings .. see tcl man pages as to why
    Tcl_Free((char*)strings);

    currentArg++;

    return 0;
}

extern "C"
int OPS_EvalDoubleStringExpression(const char* theExpression, double& current_val) {
    if (Tcl_ExprDouble(theInterp, theExpression, &current_val) != TCL_OK) {
        opserr << "OPS_EvalDoubleStringExpression::evaluateExpression -- expression \"" << theExpression;
        opserr << "\" caused error:" << endln << Tcl_GetStringResult(theInterp) << endln;
        return -1;
    }
    return 0;
}


extern "C" int OPS_SetDoubleListsOutput(
    std::vector<std::vector<double>>& data) {
    // a vector holds tcl objects of lists
    std::vector<Tcl_Obj*> tclData(data.size());

    // for each sublist
    for (int i = 0; i < (int)data.size(); ++i) {
        std::vector<Tcl_Obj*> sublist(data[i].size());
        for (int j = 0; j < (int)data[i].size(); ++j) {
            sublist[j] = Tcl_NewDoubleObj(data[i][j]);
        }
        tclData[i] = Tcl_NewListObj((int)sublist.size(), &sublist[0]);
    }

    // Tcl object for list of list
    Tcl_Obj* lists = Tcl_NewListObj((int)tclData.size(), &tclData[0]);

    // set result
    Tcl_SetObjResult(theInterp, lists);

    return 0;
}

extern "C" int OPS_SetDoubleDictOutput(
    std::map<const char*, double>& data) {
    // dict object
    auto* dict = Tcl_NewDictObj();

    // for each item
    for (auto& item : data) {
        Tcl_DictObjPut(
            theInterp, dict,
            Tcl_NewStringObj(item.first, strlen(item.first)),
            Tcl_NewDoubleObj(item.second));
    }

    // set result
    Tcl_SetObjResult(theInterp, dict);

    return 0;
}

extern "C" int OPS_SetDoubleDictListOutput(
    std::map<const char*, std::vector<double>>& data) {
    // dict object
    auto* dict = Tcl_NewDictObj();

    // for each item
    for (auto& item : data) {
        // sublist
        std::vector<Tcl_Obj*> sublist(item.second.size());
        for (int j = 0; j < (int)item.second.size(); ++j) {
            sublist[j] = Tcl_NewDoubleObj(item.second[j]);
        }
        auto* obj = Tcl_NewListObj((int)sublist.size(), &sublist[0]);

        Tcl_DictObjPut(
            theInterp, dict,
            Tcl_NewStringObj(item.first, strlen(item.first)), obj);
    }

    // set result
    Tcl_SetObjResult(theInterp, dict);

    return 0;
}

// 
// END INTERPRETER STUFF
//


extern bool builtModel;
extern FE_Datastore *theDatabase;
static BasicModelBuilder *theModelBuilder = nullptr;

G3_Runtime *
G3_getRuntime(Tcl_Interp *interp)
{
  G3_Runtime *rt = (G3_Runtime*)Tcl_GetAssocData(interp, "G3_Runtime", nullptr);
  if (!rt)
    opserr << G3_WARN_PROMPT << " No runtime\n";;
  return rt;
}

Tcl_Interp *
G3_getInterpreter(G3_Runtime* rt) {return rt->m_interp;}

BasicModelBuilder *
G3_getModelBuilder(G3_Runtime *rt) {return rt->m_builder;}

int
G3_setModelBuilder(G3_Runtime *rt, BasicModelBuilder* builder)
{
  theModelBuilder = builder;
  rt->m_builder = builder;
  return 1;
}

BasicModelBuilder *
G3_getSafeBuilder(G3_Runtime *rt)
{
  return (BasicModelBuilder*)G3_getModelBuilder(rt);
}

int
G3_setDomain(G3_Runtime *rt, Domain* domain){
  int exists = rt->m_domain ? 1 : 0;
  rt->m_domain = domain;
  return exists;
}


Domain *
G3_getDomain(G3_Runtime *rt)
{
  return rt->m_domain;
}

int G3_addTimeSeries(G3_Runtime *rt, TimeSeries *series)
{
  assert(series != nullptr);
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
  return builder->addTaggedObject<TimeSeries>(*series);
}


TimeSeries *G3_getTimeSeries(G3_Runtime *rt, int tag)
{
  TimeSeries *series;
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
  if (builder) {
     series = builder->getTypedObject<TimeSeries>(tag);
     // TODO
#if 1
     if (series)
       return series->getCopy();
     else
      return nullptr;
#endif
  } else {
    return nullptr;
  }

  return series;
}


CrdTransf *
G3_getCrdTransf(G3_Runtime *rt, G3_Tag tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  if (!builder) {
    return nullptr;
  }
  return builder->getTypedObject<CrdTransf>(tag);
}

SectionForceDeformation*
G3_getSectionForceDeformation(G3_Runtime* rt, int tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  assert(builder);
  SectionForceDeformation* theSection = builder->getTypedObject<FrameSection>(tag);

  if (theSection != nullptr)
    return theSection;
  else
    return builder->getTypedObject<SectionForceDeformation>(tag);

}

UniaxialMaterial *
G3_getUniaxialMaterialInstance(G3_Runtime *rt, int tag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  assert(builder != nullptr);
  return builder->getTypedObject<UniaxialMaterial>(tag);
}

int G3_addUniaxialMaterial(G3_Runtime *rt, UniaxialMaterial *mat) {
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  assert(builder != nullptr);
  assert(mat != nullptr);

  return builder->addTaggedObject<UniaxialMaterial>(*mat);
}

NDMaterial *
G3_GetNDMaterial(G3_Runtime* rt, int matTag)
{
  BasicModelBuilder* builder = G3_getSafeBuilder(rt);
  assert(builder != nullptr);
  return builder->getTypedObject<NDMaterial>(matTag);
}

int
OPS_GetNDF()
{
  return theModelBuilder->getNDF();
}

int
G3_getNDM(G3_Runtime *rt)
{
  BasicModelBuilder *builder = G3_getSafeBuilder(rt);
  assert(builder != nullptr);
  return builder->getNDM();
}

int
OPS_GetNDM(void)
{
  return theModelBuilder->getNDM();
}

bool *
OPS_builtModel(void) 
{
  return &builtModel;
}

AnalysisModel **
G3_getAnalysisModelPtr(G3_Runtime *rt){return rt->m_analysis_model_ptr;}

FE_Datastore *
OPS_GetFEDatastore() {return theDatabase;}

const char *
OPS_GetInterpPWD() {return getInterpPWD(theInterp);}

#if 0
EquiSolnAlgo **
OPS_GetAlgorithm(void) {return &theAlgorithm;}

EigenSOE **
OPS_GetEigenSOE(void) {return &theEigenSOE;}

LinearSOE **
OPS_GetSOE(void) {return &theSOE;}

StaticAnalysis **
OPS_GetStaticAnalysis(void) {return &theStaticAnalysis;}
#endif

#if 0 && !defined(OPS_USE_RUNTIME)

modelState theModelState;

UniaxialMaterial *
OPS_GetUniaxialMaterial(int matTag)
{
  return OPS_getUniaxialMaterial(matTag);
}

Domain *
OPS_GetDomain(void) {return theDomain;}

AnalysisModel **
OPS_GetAnalysisModel(void){return &theAnalysisModel;}

CrdTransf *
OPS_GetCrdTransf(int crdTag) {return OPS_getCrdTransf(crdTag);}

StaticIntegrator **
OPS_GetStaticIntegrator(void) {return &theStaticIntegrator;}

TransientIntegrator **
OPS_GetTransientIntegrator(void) {return &theTransientIntegrator;}

DirectIntegrationAnalysis **
OPS_GetTransientAnalysis(void) {return &theTransientAnalysis;}

ConvergenceTest **
OPS_GetTest(void) {return &theTest;}

ConstraintHandler **
OPS_GetHandler(void) {return &theHandler;}

DOF_Numberer **
OPS_GetNumberer(void) {return &theGlobalNumberer;}

extern "C" int
OPS_InvokeMaterialDirectly2(matObject *theMat, modelState *model,
                            double *strain, double *stress, double *tang,
                            int *isw)
{
  int error = 0;
  if (theMat != nullptr)
    theMat->matFunctPtr(theMat, model, strain, tang, stress, isw, &error);
  else
    error = -1;

  return error;
}
static void
OPS_InvokeMaterialObject(struct matObject *theMat, modelState *theModel,
                         double *strain, double *tang, double *stress, int *isw,
                         int *result)
{
  int matType = (int)theMat->theParam[0];

  if (matType == 1) {
    //  UniaxialMaterial *theMaterial = theUniaxialMaterials[matCount];
    UniaxialMaterial *theMaterial = (UniaxialMaterial *)theMat->matObjectPtr;
    if (theMaterial == 0) {
      *result = -1;
      return;
    }

    if (*isw == ISW_COMMIT) {
      *result = theMaterial->commitState();
      return;
    } else if (*isw == ISW_REVERT) {
      *result = theMaterial->revertToLastCommit();
      return;
    } else if (*isw == ISW_REVERT_TO_START) {
      *result = theMaterial->revertToStart();
      return;
    } else if (*isw == ISW_FORM_TANG_AND_RESID) {
      double matStress = 0.0;
      double matTangent = 0.0;
      int res = theMaterial->setTrial(strain[0], matStress, matTangent);
      stress[0] = matStress;
      tang[0] = matTangent;
      *result = res;
      return;
    }
  }

  return;
}

#endif
