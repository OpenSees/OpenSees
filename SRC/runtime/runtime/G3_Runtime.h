//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, OpenSees/Xara Developers
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
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



