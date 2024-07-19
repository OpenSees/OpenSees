//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// written: cmp
//
#include <stdio.h>
#include <unordered_map>
#include <string>
#include <vector>

#include <tcl.h>
#include <runtimeAPI.h>

typedef std::unordered_map<std::string, std::vector<std::string>> G3_Config;

class Domain;
class BasicModelBuilder;

class AnalysisModel;
class ConstraintHandler;
class LinearSOE;
class EigenSOE;
class DOF_Numberer;
class ConvergenceTest;
class StaticIntegrator;
class TransientIntegrator;


// #define G3_Builder G3_Runtime

class G3_Runtime {
public:

  Tcl_Interp     *m_interp = nullptr;
// MODEL BUILDING
  BasicModelBuilder *m_builder = nullptr;

  Domain         *m_domain  = nullptr;

  bool            model_is_built=false;

// ANALYSIS
  AnalysisModel  *m_analysis_model     = nullptr;
  AnalysisModel **m_analysis_model_ptr = &m_analysis_model;

  LinearSOE      *m_sys_of_eqn = nullptr;

  struct G3_Strategy {
    EquiSolnAlgo      *m_algorithm  = nullptr;
    ConvergenceTest   *m_convergence_test = nullptr;
    LinearSOE         *m_linear_sys = nullptr;
    EigenSOE          *m_eigen_sys  = nullptr;
    ConstraintHandler *m_handler    = nullptr;
    DOF_Numberer      *m_numberer   = nullptr;
    LinearSOE         *m_linear_soe = nullptr;
  } m_global_strategy;


  void *newStaticAnalysis(G3_Config);
  void *newTransientAnalysis(G3_Config);

// IO
  FILE* streams[3] = {stdin,stdout,stderr};
};


class G3_ParallelRuntime : public G3_Runtime {
  bool is_partitioned=false;
  int num_subdomains = 0;
  bool flag_MPID_SOE = false;
};


