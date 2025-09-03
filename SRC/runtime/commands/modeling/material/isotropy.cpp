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

#include "isotropy.h"
#include <string.h>
#include <cmath>
#include <assert.h>
#include <Logging.h>
#include <Parsing.h>

namespace Isotropy {
  enum class Parameter : int {
      YoungModulus  = 1 << 0,  // E
      ShearModulus  = 1 << 1,  // G
      BulkModulus   = 1 << 2,  // K
      PoissonsRatio = 1 << 3,  // ν
      LameLambda    = 1 << 4   // λ, Lame's first parameter
  };
}

namespace {
  const int E_FLAG       = static_cast<int>(Isotropy::Parameter::YoungModulus);
  const int G_FLAG       = static_cast<int>(Isotropy::Parameter::ShearModulus);
  const int K_FLAG       = static_cast<int>(Isotropy::Parameter::BulkModulus);
  const int NU_FLAG      = static_cast<int>(Isotropy::Parameter::PoissonsRatio);
  const int LAMBDA_FLAG  = static_cast<int>(Isotropy::Parameter::LameLambda);

  const double TOL = 1e-12;

  // Convert any input pair (flag1,in1) and (flag2,in2) to the canonical pair (E, ν).
  // Returns 0 on success, nonzero on failure.
  int convertToEN(int flag1, double in1,
                  int flag2, double in2,
                  double &E, double &nu)
  {
    // Case A: Already (E, ν)
    if ((flag1 == E_FLAG && flag2 == NU_FLAG) ||
        (flag1 == NU_FLAG && flag2 == E_FLAG)) {
      E  = (flag1 == E_FLAG)  ? in1 : in2;
      nu = (flag1 == NU_FLAG) ? in1 : in2;
      return 0;
    }
    // Case B: (G, ν): E = 2G(1+ν)
    if ((flag1 == G_FLAG && flag2 == NU_FLAG) ||
        (flag1 == NU_FLAG && flag2 == G_FLAG)) {
      double G = (flag1 == G_FLAG)  ? in1 : in2;
      nu       = (flag1 == NU_FLAG) ? in1 : in2;
      E = 2.0 * G * (1.0 + nu);
      return 0;
    }
    // Case C: (K, ν): E = 3K(1-2ν)
    if ((flag1 == K_FLAG && flag2 == NU_FLAG) ||
        (flag1 == NU_FLAG && flag2 == K_FLAG)) {
      double K = (flag1 == K_FLAG)  ? in1 : in2;
      nu       = (flag1 == NU_FLAG) ? in1 : in2;
      E = 3.0 * K * (1.0 - 2.0 * nu);
      return 0;
    }
    // Case D: (λ, ν): E = λ (1+ν)(1-2ν)/ν   (if ν != 0)
    if ((flag1 == LAMBDA_FLAG && flag2 == NU_FLAG) ||
        (flag1 == NU_FLAG && flag2 == LAMBDA_FLAG)) {
      double lambda = (flag1 == LAMBDA_FLAG) ? in1 : in2;
      nu            = (flag1 == NU_FLAG) ? in1 : in2;
      if (std::fabs(nu) < TOL)
        return -1;
      E = lambda * (1.0 + nu) * (1.0 - 2.0 * nu) / nu;
      return 0;
    }
    // Case E: (G, K): Use standard formulas:
    //    ν = (3K - 2G) / (2(3K+G))
    //    E = 9K·G/(3K+G)
    if ((flag1 == G_FLAG && flag2 == K_FLAG) ||
        (flag1 == K_FLAG && flag2 == G_FLAG)) {
      double G = (flag1 == G_FLAG) ? in1 : in2;
      double K = (flag1 == K_FLAG) ? in1 : in2;
      if (std::fabs(3.0*K + G) < TOL)
        return -1;
      nu = (3.0 * K - 2.0 * G) / (2.0 * (3.0 * K + G));
      E = 9.0 * K * G / (3.0 * K + G);
      return 0;
    }
    // Case F: (G, λ): From known relations for isotropic materials:
    //    ν = λ / [2(λ+G)]
    //    E = G(3λ+2G)/(λ+G)
    if ((flag1 == G_FLAG && flag2 == LAMBDA_FLAG) ||
        (flag1 == LAMBDA_FLAG && flag2 == G_FLAG)) {
      double G = (flag1 == G_FLAG) ? in1 : in2;
      double lambda = (flag1 == LAMBDA_FLAG) ? in1 : in2;
      if (std::fabs(lambda + G) < TOL) return -1;
      nu = lambda / (2.0 * (lambda + G));
      E = G * (3.0 * lambda + 2.0 * G) / (lambda + G);
      return 0;
    }
    // Case G: (K, λ): Use the relation K = λ + 2G/3.
    //    Hence G = 3(K - λ)/2, then
    //    ν = λ / (3K - λ)
    //    E = 9K(K - λ)/(3K - λ)
    if ((flag1 == K_FLAG && flag2 == LAMBDA_FLAG) ||
        (flag1 == LAMBDA_FLAG && flag2 == K_FLAG)) {
      double K = (flag1 == K_FLAG) ? in1 : in2;
      double lambda = (flag1 == LAMBDA_FLAG) ? in1 : in2;
      if (std::fabs(3.0*K - lambda) < TOL) return -1;
      nu = lambda / (3.0 * K - lambda);
      E = 9.0 * K * (K - lambda) / (3.0 * K - lambda);
      return 0;
    }
    // Case H: (E, G): ν = E/(2G) - 1.
    if ((flag1 == E_FLAG && flag2 == G_FLAG) ||
        (flag1 == G_FLAG && flag2 == E_FLAG)) {
      E = (flag1 == E_FLAG) ? in1 : in2;
      double G = (flag1 == G_FLAG) ? in1 : in2;
      if (std::fabs(G) < TOL)
        return -1;
      nu = E / (2.0 * G) - 1.0;
      return 0;
    }
    // Case I: (E, K): ν = (3K - E)/(6K).
    if ((flag1 == E_FLAG && flag2 == K_FLAG) ||
        (flag1 == K_FLAG && flag2 == E_FLAG)) {
      E = (flag1 == E_FLAG) ? in1 : in2;
      double K = (flag1 == K_FLAG) ? in1 : in2;
      if (std::fabs(K) < TOL) return -1;
      nu = (3.0 * K - E) / (6.0 * K);
      return 0;
    }
    // Case J: (E, λ): This leads to a quadratic in ν.
    if ((flag1 == E_FLAG && flag2 == LAMBDA_FLAG) ||
        (flag1 == LAMBDA_FLAG && flag2 == E_FLAG)) {
      E = (flag1 == E_FLAG) ? in1 : in2;
      double lambda = (flag1 == LAMBDA_FLAG) ? in1 : in2;
      // Equation:  λ(1+ν)(1-2ν) = Eν   →  2λ ν² + (λ+E)ν - λ = 0.
      double a = 2.0 * lambda;
      double b = lambda + E;
      double c = -lambda;
      double disc = b*b - 4.0 * a * c;
      if (disc < 0.0)
        return -1;
      double sqrt_disc = std::sqrt(disc);
      double nu1 = (-b + sqrt_disc) / (2.0 * a);
      double nu2 = (-b - sqrt_disc) / (2.0 * a);
      // Choose the solution in the physical range (-1, 0.5).
      if (nu1 > -1.0 && nu1 < 0.5) {
        nu = nu1; 
        return 0; 
      }
      if (nu2 > -1.0 && nu2 < 0.5) {
        nu = nu2; 
        return 0; 
      }
      return -1;
    }

    // Any other combination is not supported.
    return -1;
  }
} // namespace


int
isotropic_constants(int flag1, double in1, int flag2, double in2, IsotropicConstants &iso)
{
  // Convert to canonical (E, ν)
  double E, nu;
  if (convertToEN(flag1, in1, flag2, in2, E, nu) != 0)
    return -1;

  iso.E  = E;
  iso.nu = nu;
  iso.G  = E / (2.0 * (1.0 + nu));
  iso.K  = E / (3.0 * (1.0 - 2.0 * nu));
  iso.lambda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
  return 0;
}


int
isotropic_convert(int flag1, double in1,
                  int flag2, double in2,
                  int flag_out,
                  double & out)
{

  // First converts the given input pair into (E, ν) and then computes the requested
  // property (if not already provided as one of the inputs).

  // If the output is already one of the inputs, return it.
  if (flag_out == flag1) {
    out = in1;
    return 0;
  }
  if (flag_out == flag2) {
    out = in2;
    return 0;
  }

  // Convert to canonical (E, ν)
  double E, nu;
  if (convertToEN(flag1, in1, flag2, in2, E, nu) != 0)
    return -1;

  // Now, based solely on (E, ν), compute the desired property.
  if (flag_out == E_FLAG) {
    out = E;
    return 0;
  }
  if (flag_out == NU_FLAG) {
    out = nu;
    return 0;
  }
  if (flag_out == G_FLAG) {
    // G = E / [2(1+ν)]
    out = E / (2.0 * (1.0 + nu));
    return 0;
  }
  if (flag_out == K_FLAG) {
    // K = E / [3(1-2ν)]
    if (std::fabs(1.0 - 2.0*nu) < TOL) return -1;
    out = E / (3.0 * (1.0 - 2.0*nu));
    return 0;
  }
  if (flag_out == LAMBDA_FLAG) {
    // λ = E·ν / [(1+ν)(1-2ν)]
    if (std::fabs((1.0+nu)*(1.0-2.0*nu)) < TOL) return -1;
    out = E * nu / ((1.0 + nu) * (1.0 - 2.0*nu));
    return 0;
  }
  return -1;
}


int
TclCommand_setIsotropicParameters(
  ClientData clientData,
  Tcl_Interp *interp,
  Tcl_Size argc,
  TCL_Char ** const argv)
{
  // We expect the syntax:
  //    material Isotropic $tag <options>
  // The tag is still the first non-option argument (argv[2]).
  // Among the options, exactly two independent elastic parameters must be provided.
  // Valid elastic options: -E, -G, -K, -nu, -lambda.

  bool gotParam1 = false,
        gotParam2 = false;

  int flag1 = 0, 
      flag2 = 0;
  double val1 = 0.0, 
         val2 = 0.0;

  IsotropicParse* data = static_cast<IsotropicParse*>(clientData);
  IsotropicConstants* iso  = &data->constants;
  std::set<int>& positions =  data->positions;

  assert(iso != nullptr);

  // Process the remaining arguments.
  for (int i = 0; i < argc; i++) {
    if ((strcmp(argv[i], "-E") == 0) || 
        (strcmp(argv[i], "-youngs-modulus") == 0)) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (Tcl_GetDouble(interp, argv[i], &iso->E) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (!gotParam1) { 
          gotParam1 = true; 
          val1 = iso->E; 
          flag1 = static_cast<int>(Isotropy::Parameter::YoungModulus);
      }
      else if (!gotParam2) { 
          gotParam2 = true; 
          val2 = iso->E; 
          flag2 = static_cast<int>(Isotropy::Parameter::YoungModulus); 
      }
      else {
          opserr << "Too many elastic parameter options provided.\n";
          return TCL_ERROR;
      }
      positions.insert(i-1);
      positions.insert(i);
    }

    else if (strcmp(argv[i], "-G") == 0 || 
             strcmp(argv[i], "-shear-modulus") == 0 ||
             strcmp(argv[i], "-mu") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }

      if (Tcl_GetDouble(interp, argv[i], &iso->G) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (!gotParam1) {
        gotParam1 = true;
        val1 = iso->G;
        flag1 = static_cast<int>(Isotropy::Parameter::ShearModulus);
      }
      else if (!gotParam2) {
        gotParam2 = true;
        val2 = iso->G;
        flag2 = static_cast<int>(Isotropy::Parameter::ShearModulus);
      }
      else {
          opserr << "Too many elastic parameter options provided.\n";
          return TCL_ERROR;
      }
      positions.insert(i-1);
      positions.insert(i);
    }
    else if (strcmp(argv[i], "-K") == 0 || strcmp(argv[i], "-bulk-modulus") == 0) {
      if (++i >= argc) {
          opserr << "Missing value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      double val;
      if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
          opserr << "Invalid value for option " << argv[i-1] << "\n";
          return TCL_ERROR;
      }
      if (!gotParam1) {
        gotParam1 = true;
        val1 = val;
        flag1 = static_cast<int>(Isotropy::Parameter::BulkModulus);
      }
      else if (!gotParam2) {
        gotParam2 = true;
        val2 = val;
        flag2 = static_cast<int>(Isotropy::Parameter::BulkModulus);
      }
      else {
          opserr << "Too many elastic parameter options provided.\n";
          return TCL_ERROR;
      }
      positions.insert(i-1);
      positions.insert(i);
    }
    else if (strcmp(argv[i], "-nu") == 0 || strcmp(argv[i], "-poissons-ratio") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      double val;
      if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (!gotParam1) { 
        gotParam1 = true; 
        val1 = val; 
        flag1 = static_cast<int>(Isotropy::Parameter::PoissonsRatio);
      }
      else if (!gotParam2) {
        gotParam2 = true; 
        val2 = val; 
        flag2 = static_cast<int>(Isotropy::Parameter::PoissonsRatio);
      }
      else {
        opserr << "Too many elastic parameter options provided.\n";
        return TCL_ERROR;
      }
      positions.insert(i-1);
      positions.insert(i);
    }
    else if (strcmp(argv[i], "-lambda") == 0 || strcmp(argv[i], "-lame-lambda") == 0) {
      if (++i >= argc) {
        opserr << "Missing value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      double val;
      if (Tcl_GetDouble(interp, argv[i], &val) != TCL_OK) {
        opserr << "Invalid value for option " << argv[i-1] << "\n";
        return TCL_ERROR;
      }
      if (!gotParam1) {
        gotParam1 = true;
        val1 = val;
        flag1 = static_cast<int>(Isotropy::Parameter::LameLambda);
      }
      else if (!gotParam2) {
        gotParam2 = true;
        val2 = val;
        flag2 = static_cast<int>(Isotropy::Parameter::LameLambda);
      }
      else {
        opserr << "Too many elastic parameter options provided.\n";
        return TCL_ERROR;
      }
      positions.insert(i-1);
      positions.insert(i);
    }
  }
  
  // Compute canonical Young's modulus and Poisson's ratio.
  int ret = isotropic_constants(flag1, val1, flag2, val2, *iso);
  if (data->required == 1 && gotParam1)
    return TCL_OK;

  if (!gotParam1 || !gotParam2) {
    // opserr << "Must specify exactly two independent elastic parameters.\n";
    return TCL_ERROR;
  }

  if (ret != 0)
    return TCL_ERROR;

  return TCL_OK;
}
