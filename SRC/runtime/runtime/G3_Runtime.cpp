//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// written: cmp
//
#include <string>
#include <vector>
#include "G3_Runtime.h"

#include <StaticAnalysis.h>
#include <DirectIntegrationAnalysis.h>
#include <LinearSOE.h>
#include <EigenSOE.h>
#include <DOF_Numberer.h>
#include <ConstraintHandler.h>
#include <StaticIntegrator.h>
#include <TransientIntegrator.h>
#include <EquiSolnAlgo.h>
#include <VariableTimeStepDirectIntegrationAnalysis.h>

// DEFAULTS
#include <AnalysisModel.h>
#include <ProfileSPDLinSolver.h>
#include <ProfileSPDLinDirectSolver.h>
#include <NewtonRaphson.h>
#include <TransformationConstraintHandler.h>
#include <CTestNormUnbalance.h>
#include <ProfileSPDLinSOE.h>
#include <PlainHandler.h>
#include <Newmark.h>
#include <RCM.h>
#include <LoadControl.h>

class StaticIntegrator;

#define G3Config_keyExists(conf, key) ((conf).find((key)) != (conf).end())

template <typename T>
using G3_Parse = T* (*)(G3_Runtime*, int, const char **const);


// Wrap a function with signature
//     T* G3Parse_newT(G3_Runtime*, int argc, G3_Char** argv)
// so that it works with a std::vector<std::string>
template<typename T, T* (*fn)(G3_Runtime*, int, G3_Char **)>
T* G3Object_newParsed(G3_Runtime *rt, G3_Char* command, std::vector<std::string> args) {
    std::vector<G3_Char *> cstrs;
    cstrs.reserve(args.size()+1);
    cstrs.push_back(command);
    for (auto &s : args)
      cstrs.push_back(const_cast<char *>(s.c_str()));
    return (*fn)(rt, cstrs.size(), cstrs.data());
}

DOF_Numberer*        G3Parse_newNumberer(G3_Runtime*, int, G3_Char**const);
EquiSolnAlgo*        G3Parse_newEquiSolnAlgo(G3_Runtime*, int, G3_Char **const);
TransientIntegrator* G3Parse_newTransientIntegrator(ClientData, Tcl_Interp*, int, G3_Char**const);
StaticIntegrator*    G3Parse_newStaticIntegrator(G3_Runtime*, int, G3_Char**const);



void *
G3_Runtime::newStaticAnalysis(G3_Config conf)
{
  StaticIntegrator* sintegrator = nullptr;

  // INTEGRATOR
  // if (G3Config_keyExists(conf, "integrator"))
  //   sintegrator = 
  //     G3Object_newParsed<StaticIntegrator, G3Parse_newStaticIntegrator>(this, "integrator", conf["integrator"]);
  // else
    sintegrator = new LoadControl(1, 1, 1, 1);

  // CONVERGENCE TEST
  ConvergenceTest *test = nullptr;
//   if (G3Config_keyExists(conf, "test"))
//     test = 
//       G3Object_newParsed<ConvergenceTest, TclDispatch_newConvergenceTest>(this, "test", conf["test"]);
//   else
    test = new CTestNormUnbalance(1.0e-6,25,0);

  // ALGORITHM
  EquiSolnAlgo* the_algorithm = nullptr;
  if (the_algorithm == nullptr)
    the_algorithm = new NewtonRaphson(*test);
  else
    the_algorithm->setConvergenceTest(test);


  // NUMBERER
  DOF_Numberer* the_numberer = nullptr;
  if (G3Config_keyExists(conf, "numberer"))
    the_numberer = 
      G3Object_newParsed<DOF_Numberer, G3Parse_newNumberer>(this, "numberer", conf["numberer"]);
  else
    the_numberer = this->m_global_strategy.m_numberer;

  if (the_numberer == nullptr)  {
    RCM *rcm  = new RCM(false);
    if (rcm)
      the_numberer = new DOF_Numberer(*rcm);
  }
    
  // CONSTRAINT HANDLER
  ConstraintHandler *the_handler = new TransformationConstraintHandler();

  // LINEAR SYSTEM
  LinearSOE* the_soe = nullptr;
//  if (G3Config_keyExists(conf, "system"))
//    the_soe = 
//      G3Object_newParsed<LinearSOE, G3Parse_newLinearSOE>(this, "system", conf["system"]);
//  else
    the_soe = this->m_global_strategy.m_linear_soe;

  if (the_soe == nullptr) 
      the_soe = new ProfileSPDLinSOE(*new ProfileSPDLinDirectSolver());


  if (m_analysis_model == nullptr)
    m_analysis_model = new AnalysisModel();

  return  new StaticAnalysis(*m_domain,
                             *the_handler,
                             *the_numberer,
                             *m_analysis_model,
                             *the_algorithm,
                             *the_soe,
                             *sintegrator,
                             test);
}
#if 1
void *
G3_Runtime::newTransientAnalysis(G3_Config conf)
{
  // NUMBERER
  DOF_Numberer* the_numberer = nullptr;
  if (G3Config_keyExists(conf, "numberer"))
    the_numberer = 
      G3Object_newParsed<DOF_Numberer, G3Parse_newNumberer>(this, "numberer", conf["numberer"]);
  else
    the_numberer = this->m_global_strategy.m_numberer;

  if (the_numberer == nullptr)  {
    RCM *rcm  = new RCM(false);
    if (rcm)
      the_numberer = new DOF_Numberer(*rcm);
  }

  // CONSTRAINT HANDLER
  ConstraintHandler *the_handler = new TransformationConstraintHandler();

  // CONVERGENCE TEST
  ConvergenceTest *test = nullptr;
//  if (G3Config_keyExists(conf, "test"))
//    test = 
//      G3Object_newParsed<ConvergenceTest, TclDispatch_newConvergenceTest>(this, "test", conf["test"]);
//  else
    test = new CTestNormUnbalance(1.0e-6,25,0);

  // ALGORITHM
  EquiSolnAlgo* the_algorithm = nullptr;
  /* TODO
  if (G3Config_keyExists(conf, "algorithm"))
    the_algorithm = 
      G3Object_newParsed<EquiSolnAlgo, G3Parse_newEquiSolnAlgo>(this, "algorithm", conf["algorithm"]);
  else
    the_algorithm = this->m_global_strategy.m_algorithm;
  */
  if (the_algorithm == nullptr)
    the_algorithm = new NewtonRaphson(*test);
  else
    the_algorithm->setConvergenceTest(test);


  // LINEAR SYSTEM
  LinearSOE* the_soe = nullptr;
//  if (G3Config_keyExists(conf, "system"))
//    the_soe = 
//      G3Object_newParsed<LinearSOE, G3Parse_newLinearSOE>(this, "system", conf["system"]);
//  else
    the_soe = this->m_global_strategy.m_linear_soe;

  if (the_soe == nullptr) 
      the_soe = new ProfileSPDLinSOE(*new ProfileSPDLinDirectSolver());


  // ANALYSIS MODEL
  if (m_analysis_model == nullptr)
    m_analysis_model = new AnalysisModel();
 

  TransientIntegrator* tintegrator = nullptr;
//   if (G3Config_keyExists(conf, "integrator"))
//     tintegrator = 
//       G3Object_newParsed<TransientIntegrator, G3Parse_newTransientIntegrator>(this, "integrator", conf["integrator"]);
//   else
      tintegrator = new Newmark(0.5, 0.25);



  if (G3Config_keyExists(conf, "analysis")) {
    if (!conf["analysis"].empty() && (conf["analysis"][0] == "Variable"))
      return new VariableTimeStepDirectIntegrationAnalysis(
                                       *m_domain,
                                       *the_handler,
                                       *the_numberer,
                                       *m_analysis_model,
                                       *the_algorithm,
                                       *the_soe,
                                       *tintegrator,
                                       test);
  }
    return new DirectIntegrationAnalysis(
                                       *m_domain,
                                       *the_handler,
                                       *the_numberer,
                                       *m_analysis_model,
                                       *the_algorithm,
                                       *the_soe,
                                       *tintegrator,
                                       test);
}

#endif
