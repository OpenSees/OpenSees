#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include "catenary_solver.h"
#include "minpack.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdarg.h>
#include <float.h>

/* ========== Global Variables ========== */
typedef struct {
    double L, deltaZ, sag, EA, w;
    double tol;
    int maxIter, verbose;
    int iterCount;
} ProblemData;

static ProblemData g_solver_data;
static int g_solver_busy = 0;

/* ========== Internal Constants ========== */
#define MIN_H       1.0e-6
#define DEFAULT_TOL         1.0e-6
#define DEFAULT_MAXITER     200
#define MAX_STRETCH_RATIO   1.1
#define MIN_STRETCH_RATIO  -1.1

/* ========== Utility Functions ========== */

/**
 * @brief Safe inverse hyperbolic sine computation
 * @param x Input value
 * @return asinh(x) computed safely for all real x
 */
static double safe_asinh(double x) {
    return log(x + sqrt(x * x + 1.0));
}

/**
 * @brief Safe division with fallback value
 * @param a Numerator
 * @param b Denominator
 * @param default_val Value to return if denominator is too small
 * @return a/b if |b| > threshold, otherwise default_val
 */
static double safe_divide(double a, double b, double default_val) {
    return (fabs(b) > DBL_MIN * 10.0) ? a / b : default_val;
}

/* ========== Core Catenary Formulas ========== */

/**
 * @brief Calculate x-coordinate at arc length s
 * @param s Arc length coordinate (0 to L0)
 * @param L0 Total cable length
 * @param H Horizontal tension
 * @param V Vertical tension at left end
 * @param w Weight per unit length
 * @param EA Axial stiffness
 * @return x-coordinate at arc length s
 */
static double x_func(double s, double L0, double H, double V, double w, double EA) {
    double total_weight = w * L0;
    double term1 = H * s / EA;
    double term2 = (H * L0 / total_weight) *
        (safe_asinh(V / H) - safe_asinh((V - w * s) / H));
    return term1 + term2;
}

/**
 * @brief Calculate z-coordinate at arc length s
 * @param s Arc length coordinate (0 to L0)
 * @param L0 Total cable length
 * @param H Horizontal tension
 * @param V Vertical tension at left end
 * @param w Weight per unit length
 * @param EA Axial stiffness
 * @return z-coordinate at arc length s
 */
static double z_func(double s, double L0, double H, double V, double w, double EA) {
    double total_weight = w * L0;
    double term1 = -(total_weight * s / EA) * (V / total_weight - s / (2.0 * L0));
    double sqrt1 = sqrt(1.0 + (V / H) * (V / H));
    double tmp = (V - w * s) / H;
    double sqrt2 = sqrt(1.0 + tmp * tmp);
    double term2 = -(H * L0 / total_weight) * (sqrt1 - sqrt2);
    return term1 + term2;
}

/**
 * @brief Calculate sag at point s relative to chord line
 * @param s Arc length coordinate
 * @param L0 Total cable length
 * @param H Horizontal tension
 * @param V Vertical tension at left end
 * @param w Weight per unit length
 * @param EA Axial stiffness
 * @param L Horizontal span
 * @param deltaZ Height difference (right minus left)
 * @return Sag (vertical distance below chord line) at point s
 */
static double sag_at_point(double s, double L0, double H, double V,
    double w, double EA, double L, double deltaZ) {
    double chord_slope = safe_divide(deltaZ, L, 0.0);
    double x_s = x_func(s, L0, H, V, w, EA);
    double z_s = z_func(s, L0, H, V, w, EA);
    return chord_slope * x_s - z_s;
}

/**
 * @brief Find arc length coordinate of maximum sag using golden section search
 * @param L0 Total cable length
 * @param H Horizontal tension
 * @param V Vertical tension at left end
 * @param w Weight per unit length
 * @param EA Axial stiffness
 * @param L Horizontal span
 * @param deltaZ Height difference (right minus left)
 * @return Arc length coordinate s where sag is maximized
 */
static double find_max_sag_point(double L0, double H, double V,
    double w, double EA, double L, double deltaZ) {

    // Phase 1: Uniform sampling to find approximate maximum location
    const int COARSE_SAMPLES = 21;  // Includes endpoints, odd number ensures midpoint
    double max_sag = -INFINITY;
    double max_s = 0.0;

    for (int i = 0; i < COARSE_SAMPLES; i++) {
        double s = (double)i / (COARSE_SAMPLES - 1) * L0;
        double sag = sag_at_point(s, L0, H, V, w, EA, L, deltaZ);

        if (sag > max_sag) {
            max_sag = sag;
            max_s = s;
        }
    }

    // Phase 2: Refined search near maximum (adaptive bisection)
    const double TOL = 1e-8;      // Position tolerance
    const int MAX_REFINE = 30;    // Maximum refinement iterations

    // Determine search interval (ensure within [0, L0])
    double left = max(0.0, max_s - L0 / (COARSE_SAMPLES - 1));
    double right = min(L0, max_s + L0 / (COARSE_SAMPLES - 1));

    // If maximum is at endpoint, return directly
    if (max_s <= TOL || max_s >= L0 - TOL) {
        return max_s;
    }

    // Bisection refinement
    for (int refine = 0; refine < MAX_REFINE; refine++) {
        double mid = (left + right) * 0.5;
        double left_mid = (left + mid) * 0.5;
        double right_mid = (mid + right) * 0.5;

        double sag_left_mid = sag_at_point(left_mid, L0, H, V, w, EA, L, deltaZ);
        double sag_right_mid = sag_at_point(right_mid, L0, H, V, w, EA, L, deltaZ);

        if (sag_left_mid > sag_right_mid) {
            right = mid;  // Maximum is on left side
        }
        else {
            left = mid;   // Maximum is on right side
        }

        // Check convergence
        if (right - left < TOL) {
            break;
        }
    }

    return (left + right) * 0.5;
}

/* ========== Residual Function for MINPACK ========== */

/**
 * @brief Residual function for MINPACK's HYBRD solver
 *
 * This function computes the residuals for the three equations:
 * 1. x(L0) = L  (horizontal span constraint)
 * 2. z(L0) = deltaZ  (height difference constraint)
 * 3. sag_max = sag_target  (maximum sag constraint)
 *
 * @param n Number of variables (always 3)
 * @param x Variables [L0, H, V]
 * @param fvec Residual vector [f1, f2, f3]
 * @param iflag Error flag (not used)
 * @return Always returns 0
 */
static int catenary_residuals(int* n, double* x, double* fvec, int* iflag) {
    ProblemData* data = &g_solver_data;

    double L0 = x[0];
    double H = x[1];
    double V = x[2];

    double L = data->L;
    double deltaZ = data->deltaZ;
    double sag = data->sag;
    double EA = data->EA;
    double w = data->w;

    /* Physical constraints penalty */
    double penalty = 0.0;
//    if (H < MIN_H) penalty += 1.0e10 * (MIN_H - H);

    double chord = sqrt(L * L + deltaZ * deltaZ);
    double stretch = (L0 - chord) / chord;
 //   if (stretch < MIN_STRETCH_RATIO) penalty += 1.0e10 * (MIN_STRETCH_RATIO - stretch);
 //   if (stretch > MAX_STRETCH_RATIO) penalty += 1.0e10 * (stretch - MAX_STRETCH_RATIO);

    /* Calculate endpoint coordinates */
    double xL0 = x_func(L0, L0, H, V, w, EA);
    double zL0 = z_func(L0, L0, H, V, w, EA);

    /* Find maximum sag point */
    double s_min = find_max_sag_point(L0, H, V, w, EA, L, deltaZ);
    double sag_calc = sag_at_point(s_min, L0, H, V, w, EA, L, deltaZ);

    /* Penalize too small sag */
//    if (sag_calc < sag * 0.5) {
 //       penalty += 1e6 * (sag * 0.5 - sag_calc);
 //   }

    /* Compute residuals */
    fvec[0] = safe_divide(xL0 - L, L, 0.0);
    fvec[1] = (zL0 - deltaZ) / (fabs(deltaZ) + 1.0);
    fvec[2] = safe_divide(sag_calc - sag, sag, 0.0);

    /* Add penalty term */
    if (penalty > 0.0) {
        for (int i = 0; i < 3; i++) fvec[i] += penalty;
    }

    data->iterCount++;

    return 0;
}

/* ========== Initial Guess Functions ========== */

/**
 * @brief Compute initial guess for inclined cable (deltaZ != 0)
 * @param L Horizontal span
 * @param deltaZ Height difference
 * @param sag Target sag
 * @param w Weight per unit length
 * @param H0 [out] Initial guess for horizontal tension
 * @param V0 [out] Initial guess for vertical tension
 * @param L0_guess [out] Initial guess for cable length
 * @param verbose Verbosity level
 * @return Error code (always CAT_SUCCESS for this function)
 */
static int compute_inclined_initial(double L, double deltaZ, double sag, double w,
    double* H0, double* V0, double* L0_guess, int verbose) {
    
    *L0_guess = L + 8 * sag * sag / 3 / L + deltaZ * deltaZ / 2 / L;
    *H0 = w * L * L / 8 / sag;
    if (sag / L > 0.15) {
        double k = L / (2 * sag);
        *H0 = w * L / (2 * asinh(k));
    }
    *V0 = w * *L0_guess / 2 - *H0 * (deltaZ / L);
    
    if (verbose >= CAT_PRINT_DEBUG) {
        printf("[DEBUG] cable initial guess:\n");
        printf("  V0 = %.3f\n", *H0);
        printf("  V0 = %.3f\n", *V0);
        printf("  L0 = %.3f\n", *L0_guess);
    }
    return CAT_SUCCESS;
}

/**
 * @brief Compute initial guess for cable (general case)
 * @param L Horizontal span
 * @param deltaZ Height difference
 * @param sag Target sag
 * @param w Weight per unit length
 * @param H0 [out] Initial guess for horizontal tension
 * @param V0 [out] Initial guess for vertical tension
 * @param L0_guess [out] Initial guess for cable length
 * @param verbose Verbosity level
 * @return Error code
 */
static int compute_rigid_initial(double L, double deltaZ, double sag, double w,
    double* H0, double* V0, double* L0_guess,
    int verbose) {

    *L0_guess = L + 8 * sag * sag / 3 / L + deltaZ * deltaZ / 2 / L;
    *H0 = w * L * L / 8 / sag;
    if (sag / L > 0.15) {
        double k = L / (2 * sag);
        *H0 = w * L / (2 * asinh(k));
    }
    *V0 = w * *L0_guess / 2 - *H0 * (deltaZ / L);

    if (verbose >= CAT_PRINT_DEBUG) {
        printf("[DEBUG] cable initial guess:\n");
        printf("  V0 = %.3f\n", *H0);
        printf("  V0 = %.3f\n", *V0);
        printf("  L0 = %.3f\n", *L0_guess);
    }
    return CAT_SUCCESS;
}

/* ========== Public Function Implementations ========== */

/**
 * @brief Get required workspace size for the solver
 * @return Required workspace size in bytes
 */
size_t catenaryWorkspaceSize(void) {
    int n = 3;                     // Number of variables: L0, H, V
    int lr = n * (n + 1) / 2;      // Size of R matrix in QR decomposition
    return (n * n + lr + n * 7) * sizeof(double);  // All internal arrays
}

/**
 * @brief Get default parameter value
 * @param index Parameter index (PARAM_TOL, PARAM_MAXITER, PARAM_VERBOSE)
 * @return Default value for the specified parameter
 */
double getDefaultParam(int index) {
    switch (index) {
    case PARAM_TOL:     return DEFAULT_TOL;
    case PARAM_MAXITER: return DEFAULT_MAXITER;
    case PARAM_VERBOSE: return CAT_PRINT_NONE;
    default: return 0.0;
    }
}

/**
 * @brief Get error code description
 * @param error Error code
 * @return String describing the error
 */
const char* getCatenaryErrorString(int error) {
    switch (error) {
    case CAT_SUCCESS:              return "Success";
    case CAT_INVALID_INPUT:        return "Invalid input";
    case CAT_NO_CONVERGENCE:       return "No convergence";
    case CAT_PHYSICAL_VIOLATION:   return "Physical violation";
    case CAT_NUMERICAL_ERROR:      return "Numerical error";
    case CAT_INSUFFICIENT_MEMORY:  return "Insufficient memory";
    case CAT_INTERNAL_ERROR:       return "Internal error";
    default:                       return "Unknown error";
    }
}

/**
 * @brief Elastic catenary solver (simplified interface)
 *
 * Wrapper function that calls the extended interface with default parameters.
 */
int elasticCatenarySolve(double L, double deltaZ, double sag,
    double EA, double w,
    double* L0, double* H, double* V,
    int verbose) {

    double params[PARAM_COUNT];
    double results[RESULT_COUNT];

    /* Fill parameter array */
    params[PARAM_L] = L;
    params[PARAM_DELTAZ] = deltaZ;
    params[PARAM_SAG] = sag;
    params[PARAM_EA] = EA;
    params[PARAM_W] = w;
    params[PARAM_TOL] = getDefaultParam(PARAM_TOL);
    params[PARAM_MAXITER] = getDefaultParam(PARAM_MAXITER);
    params[PARAM_VERBOSE] = (double)verbose;

    /* Call extended interface */
    int error = elasticCatenarySolveEx(params, results, NULL, 0, NULL, 0);

    /* Extract results if successful */
    if (error == CAT_SUCCESS) {
        if (L0) *L0 = results[RESULT_L0];
        if (H)  *H = results[RESULT_H];
        if (V)  *V = results[RESULT_V];
    }

    return error;
}

/**
 * @brief Elastic catenary solver (extended interface)
 *
 * Main solver function that implements the nonlinear solution using MINPACK.
 */
int elasticCatenarySolveEx(const double* params,
    double* results,
    void* work,
    size_t workSize,
    char* errorMsg,
    size_t errorSize) {

    /* Initialize results */
    for (int i = 0; i < RESULT_COUNT; i++) {
        results[i] = 0.0;
    }
    results[RESULT_ERROR] = CAT_INTERNAL_ERROR;

    /* Validate input pointers */
    if (!params || !results) {
        if (errorMsg && errorSize > 0)
            snprintf(errorMsg, errorSize, "NULL pointer in params or results");
        return CAT_INVALID_INPUT;
    }

    /* Check if solver is already busy */
    if (g_solver_busy) {
        if (errorMsg && errorSize > 0)
            snprintf(errorMsg, errorSize, "Solver is already busy");
        return CAT_INTERNAL_ERROR;
    }
    g_solver_busy = 1;

    /* Set global solver data */
    g_solver_data.L = params[PARAM_L];
    g_solver_data.deltaZ = params[PARAM_DELTAZ];
    g_solver_data.sag = params[PARAM_SAG];
    g_solver_data.EA = params[PARAM_EA];
    g_solver_data.w = params[PARAM_W];
    g_solver_data.tol = (PARAM_COUNT > PARAM_TOL) ?
        params[PARAM_TOL] : getDefaultParam(PARAM_TOL);
    g_solver_data.maxIter = (PARAM_COUNT > PARAM_MAXITER) ?
        (int)params[PARAM_MAXITER] : (int)getDefaultParam(PARAM_MAXITER);
    g_solver_data.verbose = (PARAM_COUNT > PARAM_VERBOSE) ?
        (int)params[PARAM_VERBOSE] : (int)getDefaultParam(PARAM_VERBOSE);
    g_solver_data.iterCount = 0;

    /* Validate input parameters */
    if (g_solver_data.L <= 0.0 || g_solver_data.sag <= 0.0 ||
        g_solver_data.EA <= 0.0 || g_solver_data.w <= 0.0) {
        g_solver_busy = 0;
        if (errorMsg && errorSize > 0)
            snprintf(errorMsg, errorSize, "Invalid parameter values (must be positive)");
        return CAT_INVALID_INPUT;
    }

    /* Compute initial guess */
    double H0, V0, L0_guess;
    compute_rigid_initial(g_solver_data.L, g_solver_data.deltaZ,
        g_solver_data.sag, g_solver_data.w,
        &H0, &V0, &L0_guess, g_solver_data.verbose);

    /* Prepare variables for HYBRD */
    int n = 3;                          // Number of variables: L0, H, V
    double x[3] = { L0_guess, H0, V0 }; // Initial guess
    double fvec[3];                     // Residual vector
    double xtol = g_solver_data.tol;    // Convergence tolerance
    int maxfev = g_solver_data.maxIter * (n + 1); // Max function evaluations
    int ml = n - 1;                     // Number of subdiagonals in Jacobian
    int mu = n - 1;                     // Number of superdiagonals in Jacobian
    double epsfcn = 0.0;                // Step length for forward differences
    double factor = 100.0;              // Initial step bound
    int nprint = 0;                     // Print control (0 = no printing)
    double diag[3] = { 1.0, 1.0, 1.0 }; // Scaling factors
    int mode = 1;                       // Scaling mode
    int info, nfev;                     // Output from HYBRD

    const int ldfjac = n;               // Leading dimension of fjac
    const int lr = n * (n + 1) / 2;     // Size of R matrix

    /* Workspace arrays */
    double* fjac = NULL;
    double* r = NULL;
    double* qtf = NULL;
    double* wa1 = NULL, * wa2 = NULL, * wa3 = NULL, * wa4 = NULL;

    /* Allocate or use provided workspace */
    size_t requiredSize = catenaryWorkspaceSize();

    if (work != NULL && workSize >= requiredSize) {
        /* Use provided workspace */
        char* mem = (char*)work;
        fjac = (double*)mem; mem += n * n * sizeof(double);
        r = (double*)mem;    mem += lr * sizeof(double);
        qtf = (double*)mem;  mem += n * sizeof(double);
        wa1 = (double*)mem;  mem += n * sizeof(double);
        wa2 = (double*)mem;  mem += n * sizeof(double);
        wa3 = (double*)mem;  mem += n * sizeof(double);
        wa4 = (double*)mem;
    }
    else {
        /* Allocate workspace internally */
        fjac = (double*)malloc(n * n * sizeof(double));
        r = (double*)malloc(lr * sizeof(double));
        qtf = (double*)malloc(n * sizeof(double));
        wa1 = (double*)malloc(n * sizeof(double));
        wa2 = (double*)malloc(n * sizeof(double));
        wa3 = (double*)malloc(n * sizeof(double));
        wa4 = (double*)malloc(n * sizeof(double));

        if (!fjac || !r || !qtf || !wa1 || !wa2 || !wa3 || !wa4) {
            /* Cleanup on allocation failure */
            free(fjac); free(r); free(qtf);
            free(wa1); free(wa2); free(wa3); free(wa4);
            g_solver_busy = 0;
            if (errorMsg && errorSize > 0)
                snprintf(errorMsg, errorSize, "Memory allocation failed");
            return CAT_INSUFFICIENT_MEMORY;
        }
    }

    /* Call MINPACK's HYBRD solver */
    HYBRD(catenary_residuals, n, x, fvec, xtol, maxfev, ml, mu, epsfcn,
        diag, mode, factor, nprint, &info, &nfev,
        fjac, ldfjac, r, lr, qtf, wa1, wa2, wa3, wa4);

    /* Store solution */
    results[RESULT_L0] = x[0];
    results[RESULT_H] = x[1];
    results[RESULT_V] = x[2];
    results[RESULT_ITER] = g_solver_data.iterCount;
    results[RESULT_NFEV] = nfev;

    /* Calculate maximum sag point and actual sag */
    double s_min = find_max_sag_point(x[0], x[1], x[2],
        g_solver_data.w, g_solver_data.EA,
        g_solver_data.L, g_solver_data.deltaZ);

    double actualSag = sag_at_point(s_min, x[0], x[1], x[2],
        g_solver_data.w, g_solver_data.EA,
        g_solver_data.L, g_solver_data.deltaZ);

    results[RESULT_SMIN] = s_min;
    results[RESULT_ACTUALSAG] = actualSag;

    /* Convert MINPACK info to our error codes */
    if (info == 1) {
        results[RESULT_ERROR] = CAT_SUCCESS;
    }
    else if (info == 2) {
        results[RESULT_ERROR] = CAT_NO_CONVERGENCE;
    }
    else {
        results[RESULT_ERROR] = CAT_NUMERICAL_ERROR;
    }

    /* Print summary if requested */
    if (g_solver_data.verbose >= CAT_PRINT_SUMMARY) {
        printf("\n========================================\n");
        printf("Elastic Catenary Solution\n");
        printf("========================================\n");
        printf(" Cable length (L0):  %.6f m\n", x[0]);
        printf(" Horizontal tension (H): %.2f N\n", x[1]);
        printf(" Vertical tension (V):   %.2f N\n", x[2]);
        printf(" Max sag point (s_min):  %.6f m\n", s_min);
        printf(" Calculated sag:         %.6f m\n", actualSag);
        printf(" Target sag:             %.6f m\n", g_solver_data.sag);
        printf(" Iterations:             %d\n", g_solver_data.iterCount);
        printf(" Function evaluations:   %d\n", nfev);
        printf(" Status:                 %s\n",
            getCatenaryErrorString(results[RESULT_ERROR]));
        printf("========================================\n");
    }

    /* Cleanup if we allocated workspace internally */
    if (work == NULL) {
        free(fjac); free(r); free(qtf);
        free(wa1); free(wa2); free(wa3); free(wa4);
    }

    g_solver_busy = 0;
    return results[RESULT_ERROR];
}