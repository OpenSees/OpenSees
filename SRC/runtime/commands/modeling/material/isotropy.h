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
#pragma once
#include <set>
#include <tcl.h>
#include <Parsing.h>

struct IsotropicConstants {
  double E;      // Young's modulus
  double G;      // Shear modulus
  double K;      // Bulk modulus
  double nu;     // Poisson's ratio
  double lambda; // Lamé's first parameter
};

struct IsotropicParse {
  IsotropicConstants &constants;
  int required;
  std::set<int>       positions;
};

//---------------------------------------------------------------------
// Conversion routine:
//   Given two independent isotropic elastic constants—specified by the pair
//   (flag1, in1) and (flag2, in2)—compute the requested output property,
//   as indicated by flag_out. If one of the outputs is already among the
//   inputs its value is returned; otherwise the function first converts
//   the pair to the canonical pair (E, ν) and then computes the desired value.
//
//   The function returns 0 on success and a nonzero value if the conversion
//   cannot be performed (for example, due to an unrecognized combination or
//   division by zero).
//
//   Example: To compute the shear modulus (G) from a pair of inputs,
//     flag_out should be set to the value corresponding to Isotropy::Parameter::ShearModulus
//

int
TclCommand_setIsotropicParameters(ClientData, Tcl_Interp *, Tcl_Size argc, TCL_Char **argv);

int isotropic_convert(int flag1, double in1,
                         int flag2, double in2,
                         int flag_out,
                         double & out);

int isotropic_constants(int flag1, double in1, int flag2, double in2, IsotropicConstants &iso);
