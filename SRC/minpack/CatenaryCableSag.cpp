// CatenaryCableMaxsag.cpp
// Factory function for CatenaryCable element with Maxsag support

#define _CRT_SECURE_NO_WARNINGS
#define _USE_MATH_DEFINES

#include <OPS_Globals.h>
#include <Domain.h>
#include <Node.h>
#include <CatenaryCable.h>
#include "catenary_solver.h"

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <algorithm>
#include <string>

/**
 * @brief Process catenary cable with specified maximum sag
 *
 * Computes the unstressed length L0 for a cable given a target sag value.
 * Uses parabolic approximation for initial guess followed by elastic
 * catenary solver iteration.
 *
 * @param theDomain Pointer to OpenSees domain
 * @param eleTag Element tag
 * @param nodeI First node ID
 * @param nodeJ Second node ID
 * @param weight Weight per unit length (N/m)
 * @param E Young's modulus (Pa)
 * @param A Cross-sectional area (m²)
 * @param pL0 IN: Target sag (m), OUT: Computed unstressed length (m)
 * @param alpha Thermal expansion coefficient (1/°C)
 * @param tempChange Temperature change (°C)
 * @param rho Material density (kg/m³)
 * @param errorTol Solution tolerance
 *
 * @return 0 on success, negative on error
 */
int processCatenaryWithMaxsag(
    Domain* theDomain,
    int eleTag, int nodeI, int nodeJ,
    double weight, double E, double A,
    double* pL0,
    double alpha, double tempChange, double rho, double errorTol)
{
//    opserr << "\n=== Processing Maxsag ===\n";
//    opserr << "Element tag: " << eleTag << "\n";
//    opserr << "Nodes: " << nodeI << " -> " << nodeJ << "\n";
//    opserr << "Target sag: " << *pL0 << " m\n";

    // 1. Validate domain
    if (theDomain == 0) {
        opserr << "Error: Domain is null\n";
        return -1;
    }

    // 2. Get node coordinates
    Node* node1 = theDomain->getNode(nodeI);
    Node* node2 = theDomain->getNode(nodeJ);

    if (node1 == 0 || node2 == 0) {
        opserr << "Error: Cannot find nodes\n";
        return -2;
    }

    const Vector& pos1 = node1->getCrds();
    const Vector& pos2 = node2->getCrds();

    if (pos1.Size() < 3 || pos2.Size() < 3) {
        opserr << "Error: Need 3D coordinates\n";
        return -3;
    }

    // 3. Compute geometric parameters
    double dx = pos2(0) - pos1(0);
    double dy = pos2(1) - pos1(1);
    double dz = pos2(2) - pos1(2);

    double L = sqrt(dx * dx + dy * dy);      // Horizontal span
    double Lspan = sqrt(L * L + dz * dz);    // Chord length
    double deltaZ = dz;                  // Height difference

 //   opserr << "Node " << nodeI << ": (" << pos1(0) << ", " << pos1(1) << ", " << pos1(2) << ")\n";
 //   opserr << "Node " << nodeJ << ": (" << pos2(0) << ", " << pos2(1) << ", " << pos2(2) << ")\n";
 //   opserr << "Horizontal span L = " << L << " m\n";
 //   opserr << "Height difference Δz = " << deltaZ << " m\n";
 //   opserr << "Chord length = " << Lspan << " m\n";

    // 4. Validate basic mathematical constraints
    double targetSag = *pL0;

    // Mathematical validity checks only
    if (L <= 1e-12) {
        opserr << "Error: Span too small or zero\n";
        return -10;
    }

    if (targetSag <= 0.0) {
        opserr << "Error: Target sag must be positive\n";
        return -4;
    }

    if (E <= 0.0 || A <= 0.0) {
        opserr << "Error: E and A must be positive\n";
        return -6;
    }

    // 5. Material properties
    double w = -weight;  // Sign correction for downward weight
    double EA = E * A;

    // 6. Generate initial guess using parabolic approximation
 //   opserr << "\n--- Generating initial guess ---\n";

    double L0_guess, H0_guess, V0_guess;

    // Parabolic approximation formulas (validated through testing)
/*    L0_guess = L + 8.0 * targetSag * targetSag / (3.0 * L)
        + deltaZ * deltaZ / (2.0 * L);

    H0_guess = w * L * L / (8.0 * targetSag);

    V0_guess = w * L0_guess / 2.0 - H0_guess * (deltaZ / L);

    double sagRatio = targetSag / L;  // For information only

    opserr << "Initial guess:\n";
    opserr << "  L0 = " << L0_guess << " m\n";
    opserr << "  H = " << H0_guess << " N\n";
    opserr << "  V = " << V0_guess << " N\n";
    opserr << "  Sag/span ratio = " << sagRatio * 100 << "%\n";

    // 7. Call elastic catenary solver
    opserr << "\n--- Solving elastic catenary ---\n";
*/
    double computedL0, computedH, computedV;
    int solver_status = elasticCatenarySolve(
        L, deltaZ, targetSag, EA, w,
        &computedL0, &computedH, &computedV,
        CAT_PRINT_NONE
    );

    // 8. Process solver results
//    if (solver_status != CAT_SUCCESS) {
//        opserr << "\n=== Solver failed ===\n";
//        opserr << "Error code: " << solver_status << "\n";

        // Fallback: use parabolic approximation
//        opserr << "Using parabolic approximation as fallback\n";

//        double strain = H0_guess / EA;
//        double unstressed_length = L0_guess / (1.0 + strain);

//        *pL0 = unstressed_length;

//        opserr << "Fallback unstressed length: " << unstressed_length << " m\n";

//        return -5;  // Warning: fallback used
//    }

    // 9. Success - compute unstressed length
    double strain = computedH / EA;
    double unstressed_length = computedL0 / (1.0 + strain);

    // 10. Output complete results for user evaluation
//    opserr << "\n=== Solution converged ===\n";
//    opserr << "Cable length (stressed): " << computedL0 << " m\n";
//    opserr << "Horizontal tension H: " << computedH << " N\n";
//    opserr << "Vertical tension V: " << computedV << " N\n";
//    opserr << "Resultant tension: " << sqrt(computedH * computedH + computedV * computedV) << " N\n";
//    opserr << "Unstressed length L0: " << unstressed_length << " m\n";
//    opserr << "Sag/span ratio: " << sagRatio * 100 << "%\n";
 //   opserr << "Strain: " << strain * 100 << "%\n";
 //   opserr << "L0/Chord ratio: " << unstressed_length / Lspan << "\n";

    // 11. Update output
    *pL0 = unstressed_length;

    return 0;  // Success
}
