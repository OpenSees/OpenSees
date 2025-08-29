//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// written: cmp
//
#pragma once
#include <stdio.h>
#include <tcl.h>

class Domain;
class BasicModelBuilder;

class AnalysisModel;
class ConstraintHandler;
class LinearSOE;
class EigenSOE;
class DOF_Numberer;


class G3_Runtime {
public:

  Tcl_Interp     *m_interp = nullptr;

// MODEL BUILDING
  BasicModelBuilder *m_builder = nullptr;
  Domain            *m_domain  = nullptr;

// ANALYSIS
  AnalysisModel  *m_analysis_model     = nullptr;
  AnalysisModel **m_analysis_model_ptr = &m_analysis_model;

// IO
  FILE* streams[3] = {stdin,stdout,stderr};
};



