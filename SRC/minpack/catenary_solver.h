#ifndef CATENARY_SOLVER_H
#define CATENARY_SOLVER_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * @defgroup error_codes Error Codes
     * @brief Error codes returned by the catenary solver
     * @{
     */
#define CAT_SUCCESS               0   ///< Calculation successful
#define CAT_INVALID_INPUT         1   ///< Invalid input parameters (zero, negative, or unreasonable)
#define CAT_NO_CONVERGENCE        2   ///< Iteration did not converge
#define CAT_PHYSICAL_VIOLATION    3   ///< Solution violates physical constraints
#define CAT_NUMERICAL_ERROR       4   ///< Numerical calculation error
#define CAT_INSUFFICIENT_MEMORY   5   ///< Insufficient memory
#define CAT_INTERNAL_ERROR       99   ///< Internal error
     /** @} */

     /**
      * @defgroup verbose_levels Verbose Levels
      * @brief Control the verbosity of solver output
      * @{
      */
#define CAT_PRINT_NONE     0   ///< No output
#define CAT_PRINT_SUMMARY  1   ///< Print final summary only
#define CAT_PRINT_ITER     2   ///< Print iteration information (for debugging)
#define CAT_PRINT_DEBUG    3   ///< Print detailed debug information
      /** @} */

      /**
       * @defgroup param_indices Parameter Indices
       * @brief Indices for parameter array passed to elasticCatenarySolveEx
       *
       * Use these indices to access specific parameters in the params array.
       * @{
       */
#define PARAM_L        0   ///< Cable horizontal span (m)
#define PARAM_DELTAZ   1   ///< Endpoint height difference (m), right minus left
#define PARAM_SAG      2   ///< Maximum sag (m)
#define PARAM_EA       3   ///< Axial stiffness (N) = Young's modulus ¡Á cross-sectional area
#define PARAM_W        4   ///< Weight per unit length (N/m)
#define PARAM_TOL      5   ///< Convergence tolerance (default: 1.0e-8)
#define PARAM_MAXITER  6   ///< Maximum number of iterations (default: 200)
#define PARAM_VERBOSE  7   ///< Verbosity level (CAT_PRINT_*)
#define PARAM_COUNT    8   ///< Total number of parameters
       /** @} */

       /**
        * @defgroup result_indices Result Indices
        * @brief Indices for result array returned by elasticCatenarySolveEx
        * @{
        */
#define RESULT_L0         0   ///< Actual cable length (m), including elastic elongation
#define RESULT_H          1   ///< Horizontal tension (N)
#define RESULT_V          2   ///< Vertical tension (N), vertical component at left end
#define RESULT_SMIN       3   ///< Arc length coordinate of maximum sag point (m)
#define RESULT_ACTUALSAG  4   ///< Calculated actual sag (m)
#define RESULT_ITER       5   ///< Actual number of iterations performed
#define RESULT_NFEV       6   ///< Number of function evaluations
#define RESULT_ERROR      7   ///< Error code
#define RESULT_COUNT      8   ///< Total number of results
        /** @} */

        /**
         * @brief Elastic catenary solver (simplified interface)
         *
         * Core function for solving elastic catenary problems. Uses Newton-Raphson
         * iteration to solve the nonlinear equations.
         *
         * @param L Cable horizontal span (m), must be > 0
         * @param deltaZ Endpoint height difference (m), right minus left
         * @param sag Target maximum sag (m), must be > 0
         * @param EA Axial stiffness (N) = Young's modulus ¡Á cross-sectional area, must be > 0
         * @param w Weight per unit length (N/m), must be > 0
         * @param L0 [out] Actual cable length (m), including elastic elongation
         * @param H [out] Horizontal tension (N)
         * @param V [out] Vertical tension (N), vertical component at left end
         * @param verbose Verbosity level (CAT_PRINT_NONE, CAT_PRINT_SUMMARY, etc.)
         *
         * @return Error code (CAT_SUCCESS indicates success)
         *
         * @note For large sag values or extreme parameters, convergence may fail.
         * @note Sign of vertical tension V: positive means upward, negative means downward.
         *
         * @example Basic usage
         * @code
         * double L0, H, V;
         * int error = elasticCatenarySolve(100.0, 5.0, 10.0,
         *                                  6.123e6, 10.0,
         *                                  &L0, &H, &V,
         *                                  CAT_PRINT_SUMMARY);
         * if (error == CAT_SUCCESS) {
         *     printf("L0 = %.3f m, H = %.1f N, V = %.1f N\n", L0, H, V);
         * }
         * @endcode
         */
    int elasticCatenarySolve(
        double L, double deltaZ, double sag,
        double EA, double w,
        double* L0, double* H, double* V,
        int verbose
    );

    /**
     * @brief Elastic catenary solver (extended interface)
     *
     * Extended interface providing more control options and detailed results.
     *
     * @param params Parameter array, must contain PARAM_COUNT elements
     *               (use PARAM_* indices for access)
     * @param results Result array, will be filled with RESULT_COUNT results
     *               (use RESULT_* indices for access)
     * @param work Optional workspace buffer, can be NULL
     * @param workSize Size of workspace buffer in bytes
     * @param errorMsg Optional error message buffer, can be NULL
     * @param errorSize Size of error message buffer in bytes
     *
     * @return Error code (CAT_SUCCESS indicates success)
     *
     * @note If work is NULL, the function will allocate memory internally.
     * @note Error messages will be truncated to fit errorSize-1 characters.
     *
     * @example Extended interface usage
     * @code
     * double params[PARAM_COUNT] = {
     *     [PARAM_L] = 100.0,       // horizontal span
     *     [PARAM_DELTAZ] = 5.0,    // height difference
     *     [PARAM_SAG] = 10.0,      // sag
     *     [PARAM_EA] = 6.123e6,    // axial stiffness
     *     [PARAM_W] = 10.0,        // weight per unit length
     *     [PARAM_TOL] = 1e-10,     // tolerance
     *     [PARAM_MAXITER] = 500,   // max iterations
     *     [PARAM_VERBOSE] = CAT_PRINT_SUMMARY
     * };
     *
     * double results[RESULT_COUNT];
     * int error = elasticCatenarySolveEx(params, results,
     *                                   NULL, 0, NULL, 0);
     * if (error == CAT_SUCCESS) {
     *     double L0 = results[RESULT_L0];
     *     double H = results[RESULT_H];
     *     double sag_actual = results[RESULT_ACTUALSAG];
     * }
     * @endcode
     */
    int elasticCatenarySolveEx(
        const double* params,
        double* results,
        void* work,
        size_t workSize,
        char* errorMsg,
        size_t errorSize
    );

    /**
     * @brief Get required workspace size for the solver
     *
     * Used to pre-allocate workspace memory to avoid repeated allocations.
     *
     * @return Required workspace size in bytes
     *
     * @note The returned size includes space for all internal arrays.
     * @note If the work parameter in elasticCatenarySolveEx is NULL,
     *       the function will internally allocate memory of the same size.
     */
    size_t catenaryWorkspaceSize(void);

    /**
     * @brief Get default parameter value
     *
     * @param index Parameter index (PARAM_TOL, PARAM_MAXITER, PARAM_VERBOSE)
     * @return Default value of the parameter
     *
     * @note Only applicable for PARAM_TOL, PARAM_MAXITER, and PARAM_VERBOSE.
     * @note Returns 0.0 for other indices.
     */
    double getDefaultParam(int index);

    /**
     * @brief Get error code description
     *
     * @param error Error code (CAT_SUCCESS, CAT_INVALID_INPUT, etc.)
     * @return Error description string
     *
     * @note The returned string is static and should not be freed.
     * @note Returns "Unknown error" for unknown error codes.
     */
    const char* getCatenaryErrorString(int error);

#ifdef __cplusplus
}
#endif

#endif /* CATENARY_SOLVER_H */