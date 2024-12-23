//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
//  =============   ====================================================
//  m               mass
//  zeta            damping ratio
//  k               stiffness
//  Fy              yielding strength
//  alpha           strain-hardening ratio
//  dtF             time step for input data
//  filename        input data file, one force per line
//  dt              time step for analysis
//  uresidual       residual displacement at the end of previous analysis
//                             (optional, default=0)
//  umaxprev        previous displacement (optional, default=0)
//  =============   ====================================================

// The command returns a list of five response quantities.

// =============   =====================================================
// umax            maximum displacement during analysis
// u               displacement at end of analysis
// up              permanent residual displacement at end of analysis
// amax            maximum acceleration during analysis
// tamax           time when maximum accleration occurred
// =============   =====================================================
//
// https://portwooddigital.com/2021/02/14/how-many-clicks-does-it-take/
//
//
#include <tcl.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>

#ifndef TclPackage
#  define opserr std::cerr
#endif

#ifdef _WIN32
const char *STDIN_FILE_NAME = "CON";
#else
const char *STDIN_FILE_NAME = "/dev/stdin";
#endif

const char *Usage =
  "usage: sdof <m> <zeta> <k> <fy> <alpha> <dtF> <dt>\n\n"
  "   m               mass\n"
  "   zeta            damping ratio\n"
  "   k               stiffness\n"
  "   Fy              yielding strength\n"
  "   alpha           strain-hardening ratio\n"
  "   dtF             time step for input data\n"
  "   filename        input data file, one force per line\n"
  "   dt              time step for analysis\n"
  "   uresidual       residual displacement at the end of previous analysis\n"
  "                              (optional, default=0)\n"
  "   umaxprev        previous displacement (optional, default=0)\n";


struct SDOF_Response {
    double max_displ, u, up, max_accel, time_max_accel;
};


int
sdof_response(
    double m,
    double zeta,
    double k,
    double Fy,
    double alpha,
    double dtF,
    double dt,
    double uresidual,
    double max_prev_displ,
    const char* filename,
    struct SDOF_Response *result
    )
{
    if (filename == nullptr)
      filename = STDIN_FILE_NAME;

    std::ifstream infile(filename);
    
    const double gamma = 0.5;
    const double beta  = 0.25;
    const double tol   = 1.0e-8;
    const int maxIter  = 10;
 
    double c = zeta*2*sqrt(k*m);
    double Hkin = alpha/(1.0-alpha)*k;

    double p0  = 0.0;
    double u0  = uresidual;
    double v0  = 0.0;
    double fs0 = 0.0;
    double a0  = (p0-c*v0-fs0)/m;

    double a1 = m/(beta*dt*dt) + (gamma/(beta*dt))*c;
    double a2 = m/(beta*dt) + (gamma/beta-1.0)*c;
    double a3 = (0.5/beta-1.0)*m + dt*(0.5*gamma/beta-1.0)*c;

    double au = 1.0/(beta*dt*dt);
    double av = 1.0/(beta*dt);
    double aa = 0.5/beta-1.0;

    double vu = gamma/(beta*dt);
    double vv = 1.0-gamma/beta;
    double va = dt*(1.0-0.5*gamma/beta);
    
    double kT0 = k;

    double max_displ = fabs(max_prev_displ);
    double max_accel = 0.0;
    double time_max_accel = 0.0;
    double up = uresidual;
    double up0 = up;

    int i = 0;
    double ft, u=0, du, v, a, fs, zs, ftrial, kT, kTeff, dg, phat, R, R0;

    while (infile >> ft) {
        i++;
    
        u  = u0;
      
        fs = fs0;
        kT = kT0;
        up = up0;
      
        phat = ft + a1*u0 + a2*v0 + a3*a0;
      
        R = phat - fs - a1*u;
        R0 = R;
        if (R0 == 0.0) {
            R0 = 1.0;
        }
    
        int iter = 0;

        while (iter < maxIter && fabs(R/R0) > tol) {
            iter++;

            kTeff = kT + a1;

            du = R/kTeff;

            u = u + du;

            fs = k*(u-up0);
            zs = fs-Hkin*up0;
            ftrial = fabs(zs)-Fy;
            if (ftrial > 0) {
                dg = ftrial/(k+Hkin);
                if (fs < 0) {
                    fs = fs + dg*k;
                    up = up0 - dg;
                } else {
                    fs = fs - dg*k;
                    up = up0 + dg;
                }
                kT = k*Hkin/(k+Hkin);
            } else {
                kT = k;
            }
      
            R = phat - fs - a1*u;
        }

        v = vu*(u-u0) + vv*v0 + va*a0;
        a = au*(u-u0) - av*v0 - aa*a0;

        u0 = u;
        v0 = v;
        a0 = a;
        fs0 = fs;
        kT0 = kT;
        up0 = up;

        if (fabs(u) > max_displ) {
            max_displ = fabs(u);
        }
        if (fabs(a) > max_accel) {
            max_accel = fabs(a);
            time_max_accel = i*dt;
        }
    }
  
    infile.close();

    *result = (SDOF_Response){max_displ, u, up, max_accel, time_max_accel};

    return TCL_OK;
}

int
plastic_sdof(ClientData clientData, Tcl_Interp *interp, int argc, char **argv)
{
  const char* positional_arguments[] = {
    "m",
    "zeta",
    "k",
    "Fy",
    "alpha",
    "dtF",
    "dt",
    "uresidual",
    "max_prev_displ",
    "filename",
  };

  double m, zeta, k, Fy, alpha, dtF, dt;
  char *filename = nullptr;

  if (argc < 9) {
    opserr << "Insufficient arguments to sdfResponse\n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[1], &m) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read mass \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[2], &zeta) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read zeta \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[3], &k) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read k \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[4], &Fy) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read Fy \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[5], &alpha) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read alpha \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[6], &dtF) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dtF \n";
    return TCL_ERROR;
  }
  if (Tcl_GetDouble(interp, argv[8], &dt) != TCL_OK) {
    opserr << "WARNING sdfResponse -- could not read dt \n";
    return TCL_ERROR;
  }

  double uresidual = 0.0;
  double max_prev_displ = 0.0;
  if (argc > 9) {
    if (Tcl_GetDouble(interp, argv[9], &uresidual) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read uresidual \n";
      return TCL_ERROR;
    }
    if (Tcl_GetDouble(interp, argv[10], &max_prev_displ) != TCL_OK) {
      opserr << "WARNING sdfResponse -- could not read max_prev_displ \n";
      return TCL_ERROR;
    }
  }

  struct SDOF_Response res;
  sdof_response(m, zeta, k, Fy, alpha, dtF, dt, uresidual, max_prev_displ, filename, &res);

  char buffer[80];
  sprintf(buffer, "%f %f %f %f %f", res.max_displ, res.u, res.up, res.max_accel, res.time_max_accel);

  Tcl_SetResult(interp, buffer, TCL_VOLATILE);

  return TCL_OK;
}

int
main(int argc, char** argv)
{
  for (int i=0; i < argc; i++) {
    if (strcmp(argv[i], "-h") == 0 ||
        strcmp(argv[i], "--help") == 0) {
      puts(Usage);
      exit(0);
    }
  }

  Tcl_Interp *interp = Tcl_CreateInterp();

  plastic_sdof(nullptr, interp, argc, argv);

  fprintf(stdout, "max_displ   u    up    max_accel     time_max_accel\n");
  fprintf(stdout, "%s\n", Tcl_GetStringResult(interp));
}

