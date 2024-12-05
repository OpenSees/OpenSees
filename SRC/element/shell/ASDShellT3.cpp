/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// $Revision: 1.10 $
// $Date: 2024/03

// Original implementation: Massimo Petracca (ASDEA)
//
// A 3-node general shear-deformable (thick) shell element based on 
// the MITC3 (Mixed Interpolation of Tensorial Components) formulation 
// to avoid transverse-shear-locking in the thin shell limit.
// The membrane behavior is based on the ANDeS formulation using the
// drilling DOF to improve the membrane behavior.
// 
// It supports both linear and corotational kinematics.

#include <ASDShellT3.h>
#include <ASDShellT3CorotationalTransformation.h>

#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <Renderer.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void*
OPS_ASDShellT3(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDShellT3 - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 5) {
        opserr << "Want: element ASDShellT3 $tag $iNode $jNode $kNode $secTag "
            "<-corotational> <-reducedIntegration> <-drillingNL> <-damp $dampTag>"
            "<-local $x1 $x2 $x3>";
        return 0;
    }

    int iData[5];
    int numData = 5;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: element ASDShellT3 \n";
        return 0;
    }

    bool corotational = false;
    bool reduced_int = false;
    ASDShellT3::DrillingDOFMode drill_mode = ASDShellT3::DrillingDOF_Elastic;
    int dampingTag = 0;
    Damping* m_damping = 0;
    Vector local_x(3);

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if ((strcmp(type, "-corotational") == 0) || (strcmp(type, "-Corotational") == 0)) {
            corotational = true;
        }
        else if (strcmp(type, "-drillingNL") == 0) {
            drill_mode = ASDShellT3::DrillingDOF_NonLinear;
        }
        else if (strcmp(type, "-reducedIntegration") == 0) {
            reduced_int = true;
        }
        else if (strcmp(type, "-local") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "Error: element ASDShellT3: not enough arguments for -local options (3 components are required)\n";
                return 0;
            }
            numData = 1;
            for (int i = 0; i < 3; ++i) {
                double local_x_com;
                if (OPS_GetDoubleInput(&numData, &local_x_com) == 0) {
                    local_x(i) = local_x_com;
                }
                else {
                    opserr << "Error: element ASDShellT3: cannot get the component " << i + 1 << " for the local X axis\n";
                    return 0;
                }
            }
        }
        else if (strcmp(type, "-damp") == 0) {
            if (OPS_GetNumRemainingInputArgs() > 0) {
                numData = 1;
                if (OPS_GetIntInput(&numData, &dampingTag) < 0) return 0;
                m_damping = OPS_getDamping(dampingTag);
                if (m_damping == 0) {
                    opserr << "Error: element ASDShellT3: damping not found\n";
                    return 0;
                }
            }
        }
    }
    SectionForceDeformation* section = OPS_getSectionForceDeformation(iData[4]);
    if (section == 0) {
        opserr << "ERROR:  element ASDShellT3 " << iData[0] << " section " << iData[4] << " not found\n";
        return 0;
    }

    if (!reduced_int) {
        drill_mode = ASDShellT3::DrillingDOF_Elastic;
    }

    return new ASDShellT3(iData[0], iData[1], iData[2], iData[3], section, local_x, corotational, reduced_int, drill_mode, m_damping);
}

// anonymous namespace for utilities
namespace
{
    // some typedefs
    typedef ASDVector3<double> Vector3Type;

    // calculation options
    constexpr int OPT_NONE = 0;
    constexpr int OPT_UPDATE = (1 << 0);
    constexpr int OPT_LHS = (1 << 1);
    constexpr int OPT_RHS = (1 << 2);
    constexpr int OPT_LHS_IS_INITIAL = (1 << 3);

    // gauss quadrature data (mid-edge)
    const static std::vector<double> XI = { 0.5, 0.0, 0.5 };
    const static std::vector<double> ETA = { 0.5, 0.5, 0.0 };
    const static std::vector<double> WTS = { 1.0 / 6.0, 1.0 / 6.0, 1.0 / 6.0 };

    // gauss quadrature data (reduced)
    const static std::vector<double> XI0 = { 1.0 / 3.0 };
    const static std::vector<double> ETA0 = { 1.0 / 3.0 };
    const static std::vector<double> WTS0 = { 1.0 / 2.0 };

    // drilling curvature scale factor.
    // make sure it is not so high to over-stiffen the element
    constexpr double DH_SCALE = 1.0e-1;

    // shape functions
    inline void shapeFunctions(double xi, double eta, Vector& N)
    {
        N(0) = 1.0 - xi - eta;
        N(1) = xi;
        N(2) = eta;
    }

    // shape function derivatives in iso-parametric space
    inline void shapeFunctionsNaturalDerivatives(double xi, double eta, Matrix& dN)
    {
        dN(0, 0) = -1.0;
        dN(1, 0) = 1.0;
        dN(2, 0) = 0.0;
        
        dN(0, 1) = -1.0;
        dN(1, 1) = 0.0;
        dN(2, 1) = 1.0;
    }

    /** \brief JacobianOperator
     *
     * This class is a utility to compute at a given integration point,
     * the Jacobian, its inverse, its determinant
     * and the derivatives of the shape functions in the local
     * cartesian coordinate system.
     */
    struct JacobianOperator
    {
        // Jacobian matrix
        Matrix J = Matrix(2, 2);
        // Jacobian inverse
        Matrix invJ = Matrix(2, 2);
        // Determinant of the Jacobian matrix
        double detJ = 0.0;

        void calculate(const ASDShellT3LocalCoordinateSystem& CS, const Matrix& dN)
        {
            // jacobian
            J(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() + dN(2, 0) * CS.X3();
            J(1, 0) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() + dN(2, 0) * CS.Y3();
            J(0, 1) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() + dN(2, 1) * CS.X3();
            J(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() + dN(2, 1) * CS.Y3();

            // determinant
            detJ = J(0, 0) * J(1, 1) - J(1, 0) * J(0, 1);
            double mult = 1.0 / detJ;

            // inv(jacobian)
            invJ(0, 0) = J(1, 1) * mult;
            invJ(1, 1) = J(0, 0) * mult;
            invJ(0, 1) = -J(0, 1) * mult;
            invJ(1, 0) = -J(1, 0) * mult;
        }

    };

    /** \brief ASDShellT3Globals
     *
     * This singleton class stores some data for the shell calculations that
     * can be statically instantiated to avoid useless re-allocations
     *
     */
    class ASDShellT3Globals
    {
    private:
        ASDShellT3Globals() = default;

    public:
        JacobianOperator jac; // Jacobian

        Vector UG = Vector(18); // global displacements
        Vector UL = Vector(18); // local displacements

        Matrix B = Matrix(8, 18); // strain-displacement matrix
        Matrix B1 = Matrix(8, 18); // strain-displacement matrix with inverted bending terms
        Matrix B1TD = Matrix(18, 8); // holds the B1^T*D terms
        Vector Bd = Vector(18); // strain-displacement matrix for drilling
        Vector Bhx = Vector(18); // strain-displacement matrix for drilling (curv x)
        Vector Bhy = Vector(18); // strain-displacement matrix for drilling (curv x)
        Vector N = Vector(3); // shape functions
        Matrix dN = Matrix(3, 2); // shape functions derivatives in isoparametric space
        Matrix dNdX = Matrix(3, 2); // shape functions derivatives in cartesian space
        Vector E = Vector(8); // strain vector
        Vector S = Vector(8); // stress vector
        Vector Elocal = Vector(8); // strain vector in local element coordinate system
        Matrix Re = Matrix(8, 8); // rotation matrix for strains
        Matrix Rs = Matrix(8, 8); // rotation matrix for stresses
        Matrix RsT = Matrix(8, 8); // transpose of above
        Matrix Dsection = Matrix(8, 8);
        Matrix D = Matrix(8, 8); // section tangent
        Matrix DRsT = Matrix(8, 8); // holds the product D * Rs^T

        // matrices for ANDeS higher order stiffness
        Matrix Te = Matrix(3, 3);
        Matrix T0 = Matrix(3, 9);
        Matrix Q1 = Matrix(3, 3);
        Matrix Q2 = Matrix(3, 3);
        Matrix Q3 = Matrix(3, 3);
        Matrix Q = Matrix(3, 3);
        Matrix QT0 = Matrix(3, 9);
        Matrix Bopt = Matrix(3, 9);

        // matrices for MITC3
        Matrix BsO = Matrix(2, 2);
        Matrix BsC = Matrix(2, 2);
        Matrix BsT = Matrix(2, 9);
        Matrix BsOC = Matrix(2, 2);
        Matrix Bs = Matrix(2, 9);

        // element matrices
        Matrix LHS = Matrix(18, 18); // LHS matrix (tangent stiffness)
        Matrix LHS_initial = Matrix(18, 18); // LHS matrix (initial stiffness)
        Matrix LHS_mass = Matrix(18, 18); // LHS matrix (mass matrix)
        Vector RHS = Vector(18); // RHS vector (residual vector)
        Vector RHS_winertia = Vector(18); // RHS vector (residual vector with inertia terms)

    public:
        static ASDShellT3Globals& instance() {
            static ASDShellT3Globals _instance;
            return _instance;
        }
    };

    // computes the B matrix for drilling DOF according to Huges and Brezzi
    // used only in case of reduced integration
    void computeBdrilling(
        const Matrix& dNdX, const Vector& N,
        Vector& Bd, Vector& Bhx, Vector& Bhy)
    {
        Bd.Zero();
        Bhx.Zero();
        Bhy.Zero();
        for (int i = 0; i < 3; ++i) {
            int ii = i * 6;
            Bd(ii) = -0.5 * dNdX(i, 1);
            Bd(ii + 1) = 0.5 * dNdX(i, 0);
            Bd(ii + 5) = -N(i);
            Bhx(ii + 5) = dNdX(i, 0);
            Bhy(ii + 5) = dNdX(i, 1);
        }
    }

    // computes the complete B matrix (membrane, bending and shear)
    void computeBMatrix(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const Matrix& dNdX, const Vector& N, double xi, double eta,
        bool do_opt, Matrix& B
    )
    {
        // initialize
        B.Zero();

        // geometric data
        const auto& p1 = LCS.P1();
        const auto& p2 = LCS.P2();
        const auto& p3 = LCS.P3();
        double y12 = p1.y() - p2.y();
        double y23 = p2.y() - p3.y();
        double y31 = p3.y() - p1.y();
        double x23 = p2.x() - p3.x();
        double x31 = p3.x() - p1.x();
        double x12 = p1.x() - p2.x();
        double x32 = p3.x() - p2.x();
        double y32 = p3.y() - p2.y();
        double x13 = p1.x() - p3.x();
        double y13 = p1.y() - p3.y();
        double x21 = p2.x() - p1.x();
        double y21 = p2.y() - p1.y();

        // membrane part (ANDeS, basic)
        constexpr double alpha_membrane = 1.5;
        double A2 = LCS.Area() * 2.0;
        B(0, 0) = y23 / A2;
        B(2, 0) = -x23 / A2;
        B(1, 1) = -x23 / A2;
        B(2, 1) = y23 / A2;
        B(0, 5) = alpha_membrane * y23 * (-y31 + y12) / 6 / A2;
        B(1, 5) = alpha_membrane * (-x23) * (x31 - x12) / 6 / A2;
        B(2, 5) = alpha_membrane * (-x31 * y31 + x12 * y12) / 3 / A2;

        B(0, 6) = y31 / A2;
        B(2, 6) = -x31 / A2;
        B(1, 7) = -x31 / A2;
        B(2, 7) = y31 / A2;
        B(0, 11) = alpha_membrane * y31 * (-y12 + y23) / 6 / A2;
        B(1, 11) = alpha_membrane * (-x31) * (x12 - x23) / 6 / A2;
        B(2, 11) = alpha_membrane * (-x12 * y12 + x23 * y23) / 3 / A2;

        B(0, 12) = y12 / A2;
        B(2, 12) = -x12 / A2;
        B(1, 13) = -x12 / A2;
        B(2, 13) = y12 / A2;
        B(0, 17) = alpha_membrane * y12 * (-y23 + y31) / 6 / A2;
        B(1, 17) = alpha_membrane * (-x12) * (x23 - x31) / 6 / A2;
        B(2, 17) = alpha_membrane * (-x23 * y23 + x31 * y31) / 3 / A2;

        // membrane part (ANDeS, higher order)
        if (do_opt) {
            double b1 = 1.0;
            double b2 = 2.0;
            double b3 = 1.0;
            double b4 = 0.0;
            double b5 = 1.0;
            double b6 = -1.0;
            double b7 = -1.0;
            double b8 = -1.0;
            double b9 = -2.0;
            double b0 = 0.5;
            double opt_scale = std::sqrt(3.0 / 4.0 * b0);
            // geometric data
            double l21 = x21 * x21 + y21 * y21;
            double l32 = x32 * x32 + y32 * y32;
            double l13 = x13 * x13 + y13 * y13;
            double A4 = LCS.Area() * 4.0;
            double A23 = 2.0 * LCS.Area() / 3.0;
            // matrix T0
            auto& T0 = ASDShellT3Globals::instance().T0;
            T0(0, 0) = x32; T0(0, 1) = y32; T0(0, 2) = A4; T0(0, 3) = x13; T0(0, 4) = y13; T0(0, 5) = 0; T0(0, 6) = x21; T0(0, 7) = y21; T0(0, 8) = 0;
            T0(1, 0) = x32; T0(1, 1) = y32; T0(1, 2) = 0; T0(1, 3) = x13; T0(1, 4) = y13; T0(1, 5) = A4; T0(1, 6) = x21; T0(1, 7) = y21; T0(1, 8) = 0;
            T0(2, 0) = x32; T0(2, 1) = y32; T0(2, 2) = 0; T0(2, 3) = x13; T0(2, 4) = y13; T0(2, 5) = 0; T0(2, 6) = x21; T0(2, 7) = y21; T0(2, 8) = A4;
            T0 /= A4;
            // matrix Te
            auto& Te = ASDShellT3Globals::instance().Te;
            Te(0, 0) = y23 * y13 * l21; Te(0, 1) = y31 * y21 * l32; Te(0, 2) = y12 * y32 * l13;
            Te(1, 0) = x23 * x13 * l21; Te(1, 1) = x31 * x21 * l32; Te(1, 2) = x12 * x32 * l13;
            Te(2, 0) = (y23 * x31 + x32 * y13) * l21;
            Te(2, 1) = (y31 * x12 + x13 * y21) * l32;
            Te(2, 2) = (y12 * x23 + x21 * y32) * l13;
            Te /= A4 * LCS.Area();
            // Q1
            auto& Q1 = ASDShellT3Globals::instance().Q1;
            Q1(0, 0) = A23 * b1 / l21; Q1(0, 1) = A23 * b2 / l21; Q1(0, 2) = A23 * b3 / l21;
            Q1(1, 0) = A23 * b4 / l32; Q1(1, 1) = A23 * b5 / l32; Q1(1, 2) = A23 * b6 / l32;
            Q1(2, 0) = A23 * b7 / l13; Q1(2, 1) = A23 * b8 / l13; Q1(2, 2) = A23 * b9 / l13;
            // Q2
            auto& Q2 = ASDShellT3Globals::instance().Q2;
            Q2(0, 0) = A23 * b9 / l21; Q2(0, 1) = A23 * b7 / l21; Q2(0, 2) = A23 * b8 / l21;
            Q2(1, 0) = A23 * b3 / l32; Q2(1, 1) = A23 * b1 / l32; Q2(1, 2) = A23 * b2 / l32;
            Q2(2, 0) = A23 * b6 / l13; Q2(2, 1) = A23 * b4 / l13; Q2(2, 2) = A23 * b5 / l13;
            // Q3
            auto& Q3 = ASDShellT3Globals::instance().Q3;
            Q3(0, 0) = A23 * b5 / l21; Q3(0, 1) = A23 * b6 / l21; Q3(0, 2) = A23 * b4 / l21;
            Q3(1, 0) = A23 * b8 / l32; Q3(1, 1) = A23 * b9 / l32; Q3(1, 2) = A23 * b7 / l32;
            Q3(2, 0) = A23 * b2 / l13; Q3(2, 1) = A23 * b3 / l13; Q3(2, 2) = A23 * b1 / l13;
            // Q
            auto& Q = ASDShellT3Globals::instance().Q;
            Q.addMatrix(0.0, Q1, N(0));
            Q.addMatrix(1.0, Q2, N(1));
            Q.addMatrix(1.0, Q3, N(2));
            // Bopt = Te*Q*T0
            auto& QT0 = ASDShellT3Globals::instance().QT0;
            auto& Bopt = ASDShellT3Globals::instance().Bopt;
            QT0.addMatrixProduct(0.0, Q, T0, 1.0);
            Bopt.addMatrixProduct(0.0, Te, QT0, 1.0);
            // add it to the membrane part of B (note: the sqrt(3) because we use 3 gp for this)
            for (int i = 0; i < 3; ++i) {
                for (int node = 0; node < 3; ++node) {
                    int j = node * 6;
                    int k = node * 3;
                    B(i, j) += sqrt(3.0) * opt_scale * Bopt(i, k);
                    B(i, j + 1) += sqrt(3.0) * opt_scale * Bopt(i, k + 1);
                    B(i, j + 5) += sqrt(3.0) * opt_scale * Bopt(i, k + 2);
                }
            }
        }

        // bending part
        for (int i = 0; i < 3; ++i) {
            int ii = i * 6;
            B(3, ii + 4) = -dNdX(i, 0);
            B(4, ii + 3) = dNdX(i, 1);
            B(5, ii + 3) = dNdX(i, 0);
            B(5, ii + 4) = -dNdX(i, 1);
        }

        // shear part (MITC3 treatment of transverse shear locking)
        double phi1 = std::atan2(y21, x21);
        double phi2 = 3.141592653589793 * 0.5 - std::atan2(x31, y31);
        auto& BsO = ASDShellT3Globals::instance().BsO;
        BsO(0, 0) = std::sin(phi2);
        BsO(0, 1) = -std::sin(phi1);
        BsO(1, 0) = -std::cos(phi2);
        BsO(1, 1) = std::cos(phi1);
        auto& BsC = ASDShellT3Globals::instance().BsC;
        BsC(0, 0) = std::sqrt(x13 * x13 + y13 * y13) / A2;
        BsC(0, 1) = 0.0;
        BsC(1, 0) = 0.0;
        BsC(1, 1) = std::sqrt(x21 * x21 + y21 * y21) / A2;
        auto& BsT = ASDShellT3Globals::instance().BsT;
        BsT(0, 0) = -1.0;
        BsT(0, 1) = -(y21 + y32 * xi) / 2.0;
        BsT(0, 2) = (x21 + x32 * xi) / 2.0;
        BsT(0, 3) = 1.0;
        BsT(0, 4) = -(y21 + y13 * xi) / 2.0;
        BsT(0, 5) = (x21 + x13 * xi) / 2.0;
        BsT(0, 7) = -(y21 * xi) / 2.0;
        BsT(0, 8) = (x21 * xi) / 2.0;
        BsT(1, 0) = -1.0;
        BsT(1, 1) = (y13 + y32 * eta) / 2.0;
        BsT(1, 2) = -(x13 + x32 * eta) / 2.0;
        BsT(1, 4) = (y13 * eta) / 2.0;
        BsT(1, 5) = -(x13 * eta) / 2.0;
        BsT(1, 6) = 1.0;
        BsT(1, 7) = (y13 + y21 * eta) / 2.0;
        BsT(1, 8) = -(x13 + x21 * eta) / 2.0;
        // Bs (modified shear B matrix = O*C*T)
        auto& BsOC = ASDShellT3Globals::instance().BsOC;
        BsOC.addMatrixProduct(0.0, BsO, BsC, 1.0);
        auto& Bs = ASDShellT3Globals::instance().Bs;
        Bs.addMatrixProduct(0.0, BsOC, BsT, 1.0);
        // add it to B
        for (int i = 0; i < 2; ++i) {
            for (int node = 0; node < 3; ++node) {
                int j = node * 6;
                int k = node * 3;
                for (int dof = 0; dof < 3; ++dof) {
                    B(i + 6, j + dof + 2) = Bs(i, k + dof);
                }
            }
        }
    }

    // invert bending terms
    void invertBBendingTerms(
        const Matrix& B,
        Matrix& B1)
    {
        // due to the convention in the shell sections, we need to change the sign of the bending terms
        // for the B^T case.
        B1.addMatrix(0.0, B, 1.0);
        for (int i = 3; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                B1(i, j * 6 + 3) *= -1.0;
                B1(i, j * 6 + 4) *= -1.0;
            }
        }
    }

    // computes the transformation matrix for generalized strains
    inline void getRotationMatrixForGeneralizedStrains(double radians, Matrix& T)
    {
        double c = std::cos(radians);
        double s = std::sin(radians);

        T.Zero();

        T(0, 0) = c * c;			T(0, 1) = s * s;			T(0, 2) = -s * c;
        T(1, 0) = s * s;			T(1, 1) = c * c;			T(1, 2) = s * c;
        T(2, 0) = 2.0 * s * c;		T(2, 1) = -2.0 * s * c;		T(2, 2) = c * c - s * s;

        T(3, 3) = c * c;			T(3, 4) = s * s;			T(3, 5) = -s * c;
        T(4, 3) = s * s;			T(4, 4) = c * c;			T(4, 5) = s * c;
        T(5, 3) = 2.0 * s * c;		T(5, 4) = -2.0 * s * c;		T(5, 5) = c * c - s * s;

        T(6, 6) = c;		T(6, 7) = s;
        T(7, 6) = -s;		T(7, 7) = c;
    }

    // computes the transformation matrix for generalized stresses
    inline void getRotationMatrixForGeneralizedStresses(double radians, Matrix& T)
    {
        double c = std::cos(radians);
        double s = std::sin(radians);

        T.Zero();

        T(0, 0) = c * c;		T(0, 1) = s * s;		T(0, 2) = -2.0 * s * c;
        T(1, 0) = s * s;		T(1, 1) = c * c;		T(1, 2) = 2.0 * s * c;
        T(2, 0) = s * c;		T(2, 1) = -s * c;		T(2, 2) = c * c - s * s;

        T(3, 3) = c * c;		T(3, 4) = s * s;		T(3, 5) = -2.0 * s * c;
        T(4, 3) = s * s;		T(4, 4) = c * c;		T(4, 5) = 2.0 * s * c;
        T(5, 3) = s * c;		T(5, 4) = -s * c;		T(5, 5) = c * c - s * s;

        T(6, 6) = c;		T(6, 7) = s;
        T(7, 6) = -s;		T(7, 7) = c;
    }

}

ASDShellT3::ASDShellT3()
    : Element(0, ELE_TAG_ASDShellT3)
{
}

ASDShellT3::ASDShellT3(
    int tag,
    int node1,
    int node2,
    int node3,
    SectionForceDeformation* section,
    const Vector& local_x,
    bool corotational,
    bool reduced_integration,
    DrillingDOFMode drill_mode,
    Damping* theDamping)
    : Element(tag, ELE_TAG_ASDShellT3)
    , m_reduced_integration(reduced_integration)
    , m_transformation(corotational ? new ASDShellT3CorotationalTransformation() : new ASDShellT3Transformation())
    , m_drill_mode(drill_mode)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;

    // copy sections
    for (int i = 0; i < 3; i++) {
        m_sections[i] = section->getCopy();
        if (m_sections[i] == 0) {
            opserr << "ASDShellT3::constructor - failed to get a material of type: ShellSection\n";
            exit(-1);
        }
    }

    // allocate non-linear drilling data
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill = new NLDrillingData();
    }

    // allocate local axes if provided
    if (local_x.Size() == 3 && local_x.Norm() > 0.0) {
        m_local_x = new Vector(local_x);
        m_local_x->Normalize();
    }

    // damping
    if (theDamping) {
        for (int i = 0; i < 3; i++) {
            m_damping[i] = theDamping->getCopy();
            if (!m_damping[i]) {
                opserr << "ASDShellT3::constructor - failed to get copy of damping\n";
                exit(-1);
            }
        }
    }
}

ASDShellT3::~ASDShellT3()
{
    // clean up section
    for (int i = 0; i < 3; i++)
        delete m_sections[i];

    // clean up coordinate transformation
    if (m_transformation)
        delete m_transformation;

    // clean up load vectors
    if (m_load)
        delete m_load;

    // clean up non-linear drilling data
    if (m_nldrill)
        delete m_nldrill;

    // clean up local axes
    if (m_local_x)
        delete m_local_x;

    // clean up damping
    for (int i = 0; i < 3; i++)
        if (m_damping[i])
            delete m_damping[i];
}

void ASDShellT3::setDomain(Domain* theDomain)
{
    // if domain is null
    if (theDomain == nullptr) {
        for (int i = 0; i < 3; i++)
            nodePointers[i] = nullptr;
        // set domain on transformation
        m_transformation->setDomain(theDomain, m_node_ids, m_initialized);
        // call base class implementation
        DomainComponent::setDomain(theDomain);
        return;
    }

    // node pointers
    for (int i = 0; i < 3; i++)
        nodePointers[i] = theDomain->getNode(m_node_ids(i));

    // set domain on transformation
    m_transformation->setDomain(theDomain, m_node_ids, m_initialized);

    // only if not already initialized from recvSelf
    if (!m_initialized) {

        // compute drilling stiffness
        m_drill_stiffness = 0.0;
        for (int i = 0; i < 3; i++)
            m_drill_stiffness += m_sections[i]->getInitialTangent()(2, 2);
        m_drill_stiffness /= 3.0; // average in-plane shear modulus
        m_drill_stiffness *= 0.233; // scale it to get the in-plane roto-distortional modulus

        // compute section orientation angle
        ASDShellT3LocalCoordinateSystem reference_cs = m_transformation->createReferenceCoordinateSystem();
        Vector3Type e1_local = reference_cs.Vx();
        Vector3Type e1;
        if (m_local_x) {
            // user-defined (already normalized in c-tor)
            e1 = Vector3Type(*m_local_x);
            // make sure it's on the reference x-y plane
            Vector3Type e3 = reference_cs.Vz();
            Vector3Type e2 = e3.cross(e1);
            double e2_norm = e2.normalize();
            if (e2_norm == 0.0) {
                opserr << "ASDShellT3::setDomain Error: The provided local X axis cannot be aligned with the shell normal vector";
                exit(-1);
            }
            e1 = e2.cross(e3);
            e1.normalize();
        }
        else {
            // default one
            Vector3Type P1(m_transformation->getNodes()[0]->getCrds());
            Vector3Type P2(m_transformation->getNodes()[1]->getCrds());
            Vector3Type e1 = (P2 - P1) / 2.0;
            e1.normalize();
        }
        m_angle = std::acos(std::max(-1.0, std::min(1.0, e1.dot(e1_local))));
        if (m_angle != 0.0) {
            // if they are not counter-clock-wise, let's change the sign of the angle
            const Matrix& R = reference_cs.Orientation();
            if ((e1(0) * R(1, 0) + e1(1) * R(1, 1) + e1(2) * R(1, 2)) < 0.0)
                m_angle = -m_angle;
        }

        for (int i = 0; i < 3; i++) {
            if (m_damping[i] && m_damping[i]->setDomain(theDomain, 8)) {
                opserr << "ASDShellT3::setDomain -- Error initializing damping\n";
                exit(-1);
            }
        }

        // initialized
        m_initialized = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDShellT3::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDShellT3\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << m_node_ids(0) << "\t" << m_node_ids(1);
        s << "\t" << m_node_ids(2) << "\t0.00";
        s << endln;
        s << "PROP_3D\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << -1 << "\tSHELL\t1.0\0.0";
        s << endln;
    }

    if (flag < -1) {
        int counter = (flag + 1) * -1;
        int eleTag = this->getTag();
        int i, j;
        int imax = m_reduced_integration ? 1 : 3;
        for (i = 0; i < imax; i++) {
            const Vector& stress = m_sections[i]->getStressResultant();
            s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
            for (j = 0; j < 6; j++)
                s << "\t" << stress(j);
            s << endln;
        }
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "ASDShellT3 Non-Locking Three Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << m_node_ids(0) << endln;
        s << "Node 2 : " << m_node_ids(1) << endln;
        s << "Node 3 : " << m_node_ids(2) << endln;

        s << "Material Information : \n ";
        m_sections[0]->Print(s, flag);

        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDShellT3\", ";
        s << "\"nodes\": [" << m_node_ids(0) << ", " << m_node_ids(1) << ", ";
        s << m_node_ids(2) << "], ";
        s << "\"section\": \"" << m_sections[0]->getTag() << "\"}";
    }
}

int
ASDShellT3::setDamping(Domain* theDomain, Damping* damping)
{
    if (theDomain && damping)
    {
        for (int i = 0; i < 3; i++) {
            if (m_damping[i]) 
                delete m_damping[i];
            m_damping[i] = damping->getCopy();
            if (!m_damping[i]) {
                opserr << "ASDShellT3::setDamping - failed to get copy of damping\n";
                return -1;
            }
            if (m_damping[i] && m_damping[i]->setDomain(theDomain, 8)) {
                opserr << "ASDShellT3::setDamping -- Error initializing damping\n";
                return -2;
            }
        }
    }
    return 0;
}

int ASDShellT3::getNumExternalNodes() const
{
    return 3;
}

const ID& ASDShellT3::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDShellT3::getNodePtrs(void)
{
    return m_transformation->getNodes().data();
}

int ASDShellT3::getNumDOF()
{
    return 18;
}

int ASDShellT3::commitState()
{
    int success = 0;

    // transformation
    m_transformation->commit();

    // section
    // sections
    for (int i = 0; i < 3; i++)
        success += m_sections[i]->commitState();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm = m_sections[0]->getStressResultant();
        m_nldrill->strain_comm = m_sections[0]->getSectionDeformation();
        m_nldrill->damage_comm = m_nldrill->damage;
    }

    // damping
    for (int i = 0; i < 3; i++)
        if (m_damping[i]) 
            success += m_damping[i]->commitState();

    // done
    return success;
}

int ASDShellT3::revertToLastCommit()
{
    int success = 0;

    // transformation
    m_transformation->revertToLastCommit();

    // section
    // sections
    for (int i = 0; i < 3; i++)
        success += m_sections[i]->revertToLastCommit();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm = m_sections[0]->getStressResultant();
        m_nldrill->strain_comm = m_sections[0]->getSectionDeformation();
        m_nldrill->damage = m_nldrill->damage_comm;
    }

    // damping
    for (int i = 0; i < 3; i++)
        if (m_damping[i]) 
            success += m_damping[i]->revertToLastCommit();

    // done
    return success;
}

int  ASDShellT3::revertToStart()
{
    int success = 0;

    // transformation
    m_transformation->revertToStart();

    // section
    for (int i = 0; i < 3; i++)
        success += m_sections[i]->revertToStart();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm.Zero();
        m_nldrill->strain_comm.Zero();
        m_nldrill->damage = m_nldrill->damage_comm = 0.0;
    }

    // damping
    for (int i = 0; i < 3; i++)
        if (m_damping[i])
            success += m_damping[i]->revertToStart();

    return success;
}

int ASDShellT3::update()
{
    // calculate
    auto& LHS = ASDShellT3Globals::instance().LHS;
    auto& RHS = ASDShellT3Globals::instance().RHS;
    return calculateAll(LHS, RHS, (OPT_UPDATE));
}

const Matrix& ASDShellT3::getTangentStiff()
{
    // calculate
    auto& LHS = ASDShellT3Globals::instance().LHS;
    auto& RHS = ASDShellT3Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_LHS));
    return LHS;
}

const Matrix& ASDShellT3::getInitialStiff()
{
    // calculate
    auto& LHS = ASDShellT3Globals::instance().LHS_initial;
    auto& RHS = ASDShellT3Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_LHS | OPT_LHS_IS_INITIAL));
    return LHS;
}

const Matrix& ASDShellT3::getMass()
{
    // Output matrix
    auto& LHS = ASDShellT3Globals::instance().LHS_mass;
    LHS.Zero();

    // Compute the reference coordinate system
    ASDShellT3LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Section density (already integrated through the thickness)
    // note: sections are equal...
    double rho = m_sections[0]->getRho();

    // Lumped translational mass contribution at each node
    double Tmass = rho * reference_cs.Area() / 3.0;

    // For each node ...
    for (int j = 0; j < 3; j++)
    {
        int index = j * 6;

        // Translational mass contribution
        for (int q = 0; q < 3; q++)
            LHS(index + q, index + q) += Tmass;

        // Rotational mass neglected for the moment ...
    }

    // Done
    return LHS;
}

void  ASDShellT3::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int
ASDShellT3::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "ASDShellT3::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int
ASDShellT3::addInertiaLoadToUnbalance(const Vector& accel)
{
    // Allocate load vector if necessary
    if (m_load == nullptr)
        m_load = new Vector(18);
    auto& F = *m_load;

    // Get mass matrix
    const auto& M = getMass();

    // Add -M*R*acc to unbalance, taking advantage of lumped mass matrix
    for (int i = 0; i < 3; i++)
    {
        const auto& RV = m_transformation->getNodes()[i]->getRV(accel);
        int index = i * 6;
        for (int j = 0; j < 6; j++)
            F(index + j) -= M(index + j, index + j) * RV(j);
    }

    // Done
    return 0;
}

const Vector& ASDShellT3::getResistingForce()
{
    // calculate
    auto& LHS = ASDShellT3Globals::instance().LHS;
    auto& RHS = ASDShellT3Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_RHS));
    return RHS;
}

const Vector& ASDShellT3::getResistingForceIncInertia()
{
    // calculate
    auto& LHS = ASDShellT3Globals::instance().LHS;
    auto& RHS = ASDShellT3Globals::instance().RHS_winertia;
    calculateAll(LHS, RHS, (OPT_RHS));

    // Add damping terms
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        RHS.addVector(1.0, getRayleighDampingForces(), 1.0);

    // Compute mass
    const auto& M = getMass();

    // Add M*acc to unbalance, taking advantage of lumped mass matrix
    for (int i = 0; i < 3; i++)
    {
        const auto& A = m_transformation->getNodes()[i]->getTrialAccel();
        int index = i * 6;
        for (int j = 0; j < 6; j++)
            RHS(index + j) += M(index + j, index + j) * A(j);
    }

    // Done
    return RHS;
}

int  ASDShellT3::sendSelf(int commitTag, Channel& theChannel)
{
    int res = 0;

    // note: we don't check for dataTag == 0 for Element
    // objects as that is taken care of in a commit by the Domain
    // object - don't want to have to do the check if sending data
    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // has load flag
    bool has_load = m_load != nullptr;

    // INT data
    // 1 tag +
    // 3 node tags +
    // 1 transformation type
    // 1 initialization flag
    // 1 has_load_flag
    // 6 -> 3 pairs of (section class tag, mt db tag)
    // 2 -> damping tag + damping db tag
    // 1 -> reduced integration flag
    // 1 -> non-linear drilling flag
    // 1 -> local_x flag
    static ID idData(18);
    counter = 0;
    idData(counter++) = getTag();
    for (int i = 0; i < 3; ++i)
        idData(counter++) = m_node_ids(i);
    idData(counter++) = static_cast<int>(m_transformation->isLinear());
    idData(counter++) = static_cast<int>(m_initialized);
    idData(counter++) = static_cast<int>(has_load);
    for (int i = 0; i < 3; i++) {
        idData(counter++) = m_sections[i]->getClassTag();
        int matDbTag = m_sections[i]->getDbTag();
        // NOTE: we do have to ensure that the material has a database
        // tag if we are sending to a database channel.
        if (matDbTag == 0) {
            matDbTag = theChannel.getDbTag();
            if (matDbTag != 0)
                m_sections[i]->setDbTag(matDbTag);
        }
        idData(counter++) = matDbTag;
    }
    if (m_damping[0]) {
        idData(counter++) = m_damping[0]->getClassTag();
        int dbTag = m_damping[0]->getDbTag();
        if (dbTag == 0) {
            dbTag = theChannel.getDbTag();
            if (dbTag != 0)
                for (int i = 0; i < 3; i++)
                    m_damping[i]->setDbTag(dbTag);
        }
        idData(counter++) = dbTag;
    }
    else {
        idData(counter++) = 0;
        idData(counter++) = 0;
    }
    idData(counter++) = static_cast<int>(m_reduced_integration);
    idData(counter++) = static_cast<int>(static_cast<bool>(m_nldrill));
    idData(counter++) = static_cast<int>(static_cast<bool>(m_local_x));

    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 4 damping factors +
    // 3 drilling strain +
    // 1 drilling stiffness +
    // 1 angle +
    // (optional) 3 local_x +
    // (optional) 18 non-linear drilling data +
    // (optional) 18 load +
    // NT transformation data
    int Nlocalx = m_local_x ? 3 : 0;
    int NLoad = has_load ? 18 : 0;
    int Nnldrill = m_nldrill ? 18 : 0;
    int NT = m_transformation->internalDataSize();
    Vector vectData(4 + 3 + 1 + 1 + Nlocalx + Nnldrill + NLoad + NT);
    counter = 0;
    vectData(counter++) = alphaM;
    vectData(counter++) = betaK;
    vectData(counter++) = betaK0;
    vectData(counter++) = betaKc;
    vectData(counter++) = m_drill_strain;
    vectData(counter++) = m_drill_curvature_x;
    vectData(counter++) = m_drill_curvature_y;
    vectData(counter++) = m_drill_stiffness;
    vectData(counter++) = m_angle;
    if (m_local_x) {
        for (int i = 0; i < 3; ++i)
            vectData(counter++) = (*m_local_x)(i);
    }
    if (m_nldrill) {
        for (int i = 0; i < 8; ++i)
            vectData(counter++) = m_nldrill->strain_comm(i);
        for (int i = 0; i < 8; ++i)
            vectData(counter++) = m_nldrill->stress_comm(i);
        vectData(counter++) = m_nldrill->damage;
        vectData(counter++) = m_nldrill->damage_comm;
    }
    if (has_load) {
        for (int i = 0; i < 18; ++i)
            vectData(counter++) = (*m_load)(i);
    }
    m_transformation->saveInternalData(vectData, counter);

    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send all sections
    for (int i = 0; i < 3; i++) {
        res = m_sections[i]->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send its Material\n";
            return res;
        }
    }

    // Ask the Damping to send itself
    if (m_damping[0]) {
        for (int i = 0; i < 3; i++) {
            res += m_damping[i]->sendSelf(commitTag, theChannel);
            if (res < 0) {
                opserr << "ASDShellT3::sendSelf -- could not send Damping\n";
                return res;
            }
        }
    }

    // done
    return res;
}

int  ASDShellT3::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 3 node tags +
    // 1 transformation type
    // 1 initialization flag
    // 1 has_load_flag
    // 6 -> 3 pairs of (section class tag, mt db tag)
    // 2 -> damping tag + damping db tag
    // 1 -> reduced integration flag
    // 1 -> non-linear drilling flag
    // 1 -> local_x flag
    static ID idData(18);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDShellT3::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    counter = 0;
    setTag(idData(counter++));
    for (int i = 0; i < 3; ++i)
        m_node_ids(i) = idData(counter++);
    bool linear_transform = static_cast<bool>(idData(counter++));
    m_initialized = static_cast<bool>(idData(counter++));
    bool has_load = static_cast<bool>(idData(counter++));
    // create sections
    for (int i = 0; i < 3; i++) {
        int matClassTag = idData(counter++);
        int matDbTag = idData(counter++);
        if (m_sections[i])
            delete m_sections[i];
        m_sections[i] = theBroker.getNewSection(matClassTag);
        if (m_sections[i] == 0) {
            opserr << "ASDShellT3::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
            return -1;
        }
        m_sections[i]->setDbTag(matDbTag);
    }
    int dmpTag = idData(counter++);
    int dmpDbTag = idData(counter++);
    m_reduced_integration = static_cast<bool>(idData(counter++));
    bool use_nldrill = static_cast<bool>(idData(counter++));
    bool use_local_x = static_cast<bool>(idData(counter++));

    // create transformation
    if (m_transformation)
        delete m_transformation;
    m_transformation = linear_transform ? new ASDShellT3Transformation() : new ASDShellT3CorotationalTransformation();
    // create load
    if (has_load) {
        if (m_load == nullptr)
            m_load = new Vector(18);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // create non-linear drilling
    if (m_nldrill) {
        delete m_nldrill;
        m_nldrill = nullptr;
    }
    if (use_nldrill)
        m_nldrill = new NLDrillingData();

    // clear local x
    if (m_local_x) {
        delete m_local_x;
        m_local_x = nullptr;
    }
    if (use_local_x)
        m_local_x = new Vector(3);

    // DOUBLE data
    // 4 damping factors +
    // 3 drilling strain +
    // 1 drilling stiffness +
    // 1 angle +
    // (optional) 3 local_x +
    // (optional) 18 non-linear drilling data +
    // (optional) 18 load +
    // NT transformation data
    int Nlocalx = use_local_x ? 3 : 0;
    int NLoad = has_load ? 18 : 0;
    int Nnldrill = use_nldrill ? 18 : 0;
    int NT = m_transformation->internalDataSize();
    Vector vectData(4 + 3 + 1 + 1 + Nnldrill + NLoad + NT);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellT3::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    counter = 0;
    alphaM = vectData(counter++);
    betaK = vectData(counter++);
    betaK0 = vectData(counter++);
    betaKc = vectData(counter++);
    m_drill_strain = vectData(counter++);
    m_drill_curvature_x = vectData(counter++);
    m_drill_curvature_y = vectData(counter++);
    m_drill_stiffness = vectData(counter++);
    m_angle = vectData(counter++);
    if (use_local_x) {
        for (int i = 0; i < 3; ++i)
            (*m_local_x)(i) = vectData(counter++);
    }
    if (use_nldrill) {
        for (int i = 0; i < 8; ++i)
            m_nldrill->strain_comm(i) = vectData(counter++);
        for (int i = 0; i < 8; ++i)
            m_nldrill->stress_comm(i) = vectData(counter++);
        m_nldrill->damage = vectData(counter++);
        m_nldrill->damage_comm = vectData(counter++);
    }
    if (has_load) {
        for (int i = 0; i < 18; ++i)
            (*m_load)(i) = vectData(counter++);
    }
    m_transformation->restoreInternalData(vectData, counter);

    // all sections
    for (int i = 0; i < 3; i++) {
        res = m_sections[i]->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "ASDShellT3::recvSelf() - material " << i << "failed to recv itself\n";
            return res;
        }
    }

    if (dmpTag) {
        for (int i = 0; i < 3; i++) {
            // Check if the Damping is null; if so, get a new one
            if (m_damping[i] == 0) {
                m_damping[i] = theBroker.getNewDamping(dmpTag);
                if (m_damping[i] == 0) {
                    opserr << "ASDShellT3::recvSelf -- could not get a Damping\n";
                    exit(-1);
                }
            }

            // Check that the Damping is of the right type; if not, delete
            // the current one and get a new one of the right type
            if (m_damping[i]->getClassTag() != dmpTag) {
                delete m_damping[i];
                m_damping[i] = theBroker.getNewDamping(dmpTag);
                if (m_damping[i] == 0) {
                    opserr << "ASDShellT3::recvSelf -- could not get a Damping\n";
                    exit(-1);
                }
            }

            // Now, receive the Damping
            m_damping[i]->setDbTag(dmpDbTag);
            res += m_damping[i]->recvSelf(commitTag, theChannel, theBroker);
            if (res < 0) {
                opserr << "ASDShellT3::recvSelf -- could not receive Damping\n";
                return res;
            }
        }
    }
    else {
        for (int i = 0; i < 3; i++)
        {
            if (m_damping[i])
            {
                delete m_damping[i];
                m_damping[i] = 0;
            }
        }
    }

    // done
    return res;
}

Response*
ASDShellT3::setResponse(const char** argv, int argc, OPS_Stream& output)
{
    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDShellT3");
    output.attr("eleTag", this->getTag());
    int numNodes = this->getNumExternalNodes();
    const ID& nodes = this->getExternalNodes();
    static char nodeData[32];

    for (int i = 0; i < numNodes; i++) {
        sprintf(nodeData, "node%d", i + 1);
        output.attr(nodeData, nodes(i));
    }

    if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
        strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
        const Vector& force = this->getResistingForce();
        int size = force.Size();
        for (int i = 0; i < size; i++) {
            sprintf(nodeData, "P%d", i + 1);
            output.tag("ResponseType", nodeData);
        }
        theResponse = new ElementResponse(this, 1, this->getResistingForce());
    }

    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0 ||
	     strcmp(argv[0], "section") == 0) {
        if (argc < 2) {
            opserr << "ASDShellT3::setResponse() - need to specify more data\n";
            return 0;
        }
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= 3) {
            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("eta", XI[pointNum - 1]);
            output.attr("neta", ETA[pointNum - 1]);
            int section_pos = m_reduced_integration ? 0 : pointNum - 1;
            theResponse = m_sections[section_pos]->setResponse(&argv[2], argc - 2, output);
            output.endTag();
        }
    }
    else if (strcmp(argv[0], "stresses") == 0) {

        for (int i = 0; i < 3; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", XI[i]);
            output.attr("neta", ETA[i]);

            output.tag("SectionForceDeformation");
            int section_pos = m_reduced_integration ? 0 : i;
            output.attr("classType", m_sections[section_pos]->getClassTag());
            output.attr("tag", m_sections[section_pos]->getTag());

            output.tag("ResponseType", "p11");
            output.tag("ResponseType", "p22");
            output.tag("ResponseType", "p12");
            output.tag("ResponseType", "m11");
            output.tag("ResponseType", "m22");
            output.tag("ResponseType", "m12");
            output.tag("ResponseType", "q1");
            output.tag("ResponseType", "q2");

            output.endTag(); // GaussPoint
            output.endTag(); // NdMaterialOutput
        }

        theResponse = new ElementResponse(this, 2, Vector(24));
    }

    else if (strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < 3; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", XI[i]);
            output.attr("neta", ETA[i]);

            output.tag("SectionForceDeformation");
            int section_pos = m_reduced_integration ? 0 : i;
            output.attr("classType", m_sections[section_pos]->getClassTag());
            output.attr("tag", m_sections[section_pos]->getTag());

            output.tag("ResponseType", "eps11");
            output.tag("ResponseType", "eps22");
            output.tag("ResponseType", "gamma12");
            output.tag("ResponseType", "theta11");
            output.tag("ResponseType", "theta22");
            output.tag("ResponseType", "theta33");
            output.tag("ResponseType", "gamma13");
            output.tag("ResponseType", "gamma23");

            output.endTag(); // GaussPoint
            output.endTag(); // NdMaterialOutput
        }

        theResponse = new ElementResponse(this, 3, Vector(24));
    }

    else if (m_damping && strcmp(argv[0], "dampingStresses") == 0) {

        for (int i = 0; i < 3; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", XI[i]);
            output.attr("neta", ETA[i]);

            output.tag("SectionForceDeformation");
            int section_pos = m_reduced_integration ? 0 : i;
            output.attr("classType", m_sections[section_pos]->getClassTag());
            output.attr("tag", m_sections[section_pos]->getTag());

            output.tag("ResponseType", "p11");
            output.tag("ResponseType", "p22");
            output.tag("ResponseType", "p12");
            output.tag("ResponseType", "m11");
            output.tag("ResponseType", "m22");
            output.tag("ResponseType", "m12");
            output.tag("ResponseType", "q1");
            output.tag("ResponseType", "q2");

            output.endTag(); // GaussPoint
            output.endTag(); // NdMaterialOutput
        }

        theResponse = new ElementResponse(this, 4, Vector(24));
    }

    output.endTag();
    return theResponse;
}

int
ASDShellT3::getResponse(int responseID, Information& eleInfo)
{
    static Vector stresses(24);
    static Vector strains(24);

    switch (responseID) {
    case 1: // global forces
        return eleInfo.setVector(this->getResistingForce());
        break;
    case 2: // stresses
        for (int i = 0; i < 3; i++) {

            // Get material stress response
            const Vector& sigma = m_sections[m_reduced_integration ? 0 : i]->getStressResultant();
            int offset = i * 8;
            stresses(0 + offset) = sigma(0);
            stresses(1 + offset) = sigma(1);
            stresses(2 + offset) = sigma(2);
            stresses(3 + offset) = sigma(3);
            stresses(4 + offset) = sigma(4);
            stresses(5 + offset) = sigma(5);
            stresses(6 + offset) = sigma(6);
            stresses(7 + offset) = sigma(7);
        }
        return eleInfo.setVector(stresses);
        break;
    case 3: //strain
        for (int i = 0; i < 3; i++) {

            // Get section deformation
            const Vector& deformation = m_sections[m_reduced_integration ? 0 : i]->getSectionDeformation();
            int offset = i * 8;
            strains(0 + offset) = deformation(0);
            strains(1 + offset) = deformation(1);
            strains(2 + offset) = deformation(2);
            strains(3 + offset) = deformation(3);
            strains(4 + offset) = deformation(4);
            strains(5 + offset) = deformation(5);
            strains(6 + offset) = deformation(6);
            strains(7 + offset) = deformation(7);
        }
        return eleInfo.setVector(strains);
        break;
    case 4: // damping stresses
        for (int i = 0; i < 3; i++) {

            // Get material stress response
            const Vector& sigma = m_damping[m_reduced_integration ? 0 : i]->getDampingForce();
            int offset = i * 8;
            stresses(0 + offset) = sigma(0);
            stresses(1 + offset) = sigma(1);
            stresses(2 + offset) = sigma(2);
            stresses(3 + offset) = sigma(3);
            stresses(4 + offset) = sigma(4);
            stresses(5 + offset) = sigma(5);
            stresses(6 + offset) = sigma(6);
            stresses(7 + offset) = sigma(7);
        }
        return eleInfo.setVector(stresses);
        break;
    default:
        return -1;
    }
}

int ASDShellT3::setParameter(const char** argv, int argc, Parameter& param)
{
    int res = -1;
    int matRes = res;
    for (int i = 0; i < 3; i++)
    {
        matRes = m_sections[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
    return res;
}

int ASDShellT3::calculateAll(Matrix& LHS, Vector& RHS, int options)
{
    // Check options
    if (!m_transformation->isLinear()) {
        // corotational calculation of the tangent LHS requires the RHS!
        if (options & OPT_LHS)
            options |= OPT_RHS;
    }

    // Zero output
    int result = 0;
    if (options & OPT_RHS)
        RHS.Zero();
    if (options & OPT_LHS)
        LHS.Zero();

    // Global displacements
    auto& UG = ASDShellT3Globals::instance().UG;
    m_transformation->computeGlobalDisplacements(UG);

    // update transformation
    if (options & OPT_UPDATE)
        m_transformation->update(UG);

    // Compute the reference coordinate system
    ASDShellT3LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Compute the local coordinate system.
    ASDShellT3LocalCoordinateSystem local_cs =
        m_transformation->createLocalCoordinateSystem(UG);

    // Some matrices/vectors
    auto& N = ASDShellT3Globals::instance().N;
    auto& dN = ASDShellT3Globals::instance().dN;
    auto& dNdX = ASDShellT3Globals::instance().dNdX;
    auto& jac = ASDShellT3Globals::instance().jac;
    auto& B = ASDShellT3Globals::instance().B;
    auto& Bd = ASDShellT3Globals::instance().Bd;
    auto& Bhx = ASDShellT3Globals::instance().Bhx;
    auto& Bhy = ASDShellT3Globals::instance().Bhy;
    auto& B1 = ASDShellT3Globals::instance().B1;
    auto& B1TD = ASDShellT3Globals::instance().B1TD;
    auto& E = ASDShellT3Globals::instance().E;
    auto& S = ASDShellT3Globals::instance().S;
    auto& D = ASDShellT3Globals::instance().D;
    auto& Dsection = ASDShellT3Globals::instance().Dsection;

    // matrices for orienting strains in section coordinate system
    auto& Re = ASDShellT3Globals::instance().Re;
    auto& Rs = ASDShellT3Globals::instance().Rs;
    if (m_angle != 0.0) {
        if (options & OPT_UPDATE)
            getRotationMatrixForGeneralizedStrains(-m_angle, Re);
        if ((options & OPT_RHS) || (options & OPT_LHS))
            getRotationMatrixForGeneralizedStresses(m_angle, Rs);
    }

    // Local displacements
    auto& UL = ASDShellT3Globals::instance().UL;
    m_transformation->calculateLocalDisplacements(local_cs, UG, UL);

    // Drilling data for drilling damage (optional)
    static Vector drill_dstrain(8);
    static Vector drill_dstress(8);
    static Vector drill_dstress_el(8);

    // Stenberg shear stabilization coefficient
    // to avoid shear oscillations in the thin limit
    double thickness = evaluateSectionThickness();
    double stab_coef = 0.1;
    double mesh_size = std::max(std::max(
        (reference_cs.P2() - reference_cs.P1()).norm(),
        (reference_cs.P3() - reference_cs.P2()).norm()),
        (reference_cs.P1() - reference_cs.P3()).norm());
    double sh_stab = thickness * thickness /
        (thickness * thickness + stab_coef * mesh_size * mesh_size);
    sh_stab = std::sqrt(sh_stab); // use the sqrt for symmetry (as in shear correction factor)

    // drilling curvature scale
    // Note: use lch^2 to make it size-invariant.
    // Here we use the drilling curvature only to stabilize the 2 spurious zero-energy modes,
    // we are not actually using a cosserat continuum with an internal length-scale, so
    // we don't want size-effects here.
    double dh_scale_factor = DH_SCALE * std::pow(getCharacteristicLength(), 2.0);

    // gauss integration
    const auto& gx = m_reduced_integration ? XI0 : XI;
    const auto& gy = m_reduced_integration ? ETA0 : ETA;
    const auto& gw = m_reduced_integration ? WTS0 : WTS;
    for(int igauss = 0; igauss < gx.size(); ++igauss)
    {
        // Current integration point data
        double xi = gx[igauss];
        double eta = gy[igauss];
        double w = gw[igauss];
        shapeFunctions(xi, eta, N);
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;
        dNdX.addMatrixProduct(0.0, dN, jac.invJ, 1.0);

        // Strain-displacement matrix
        computeBMatrix(reference_cs, dNdX, N, xi, eta, !m_reduced_integration, B);
        if (m_reduced_integration) {
            computeBdrilling(dNdX, N, Bd, Bhx, Bhy);
        }

        // Update strain strain
        if (options & OPT_UPDATE)
        {
            // Section deformation
            if (m_angle != 0.0) {
                auto& Elocal = ASDShellT3Globals::instance().Elocal;
                Elocal.addMatrixVector(0.0, B, UL, 1.0);
                E.addMatrixVector(0.0, Re, Elocal, 1.0);
            }
            else {
                E.addMatrixVector(0.0, B, UL, 1.0);
            }

            // apply Stenberg stabilization
            E(6) *= sh_stab;
            E(7) *= sh_stab;

            // Update section
            result += m_sections[igauss]->setTrialSectionDeformation(E);

            // Drilling strain Ed = Bd*UL
            if (m_reduced_integration) {
                m_drill_strain = 0.0;
                m_drill_curvature_x = 0.0;
                m_drill_curvature_y = 0.0;
                for (int i = 0; i < 18; i++) {
                    m_drill_strain += Bd(i) * UL(i);
                    m_drill_curvature_x += Bhx(i) * UL(i);
                    m_drill_curvature_y += Bhy(i) * UL(i);
                }
            }
        }

        // Invert bending terms for correct statement of equilibrium
        if ((options & OPT_RHS) || (options & OPT_LHS))
            invertBBendingTerms(B, B1);

        // Integrate RHS
        if (options & OPT_RHS)
        {
            // Section force
            if (m_angle != 0.0) {
                auto& Ssection = m_sections[igauss]->getStressResultant();
                S.addMatrixVector(0.0, Rs, Ssection, 1.0);
                if (m_damping[igauss]) {
                    m_damping[igauss]->update(Ssection);
                    auto& Sdsection = m_damping[igauss]->getDampingForce();
                    S.addMatrixVector(1.0, Rs, Sdsection, 1.0);
                }
            }
            else {
                S = m_sections[igauss]->getStressResultant();
                if (m_damping[igauss]) {
                    m_damping[igauss]->update(S);
                    S += m_damping[igauss]->getDampingForce();
                }
            }

            // apply Stenberg stabilization
            S(6) *= sh_stab;
            S(7) *= sh_stab;

            // Add current integration point contribution (RHS)
            RHS.addMatrixTransposeVector(1.0, B1, S, dA);

            // Drilling
            if (m_reduced_integration) {
                // Compute drilling damages
                if (m_drill_mode == DrillingDOF_NonLinear) {
                    drill_dstrain = m_sections[igauss]->getSectionDeformation();
                    drill_dstrain.addVector(1.0, m_nldrill->strain_comm, -1.0);
                    if (drill_dstrain.Norm() > 1.0e-10) {
                        drill_dstress = m_sections[igauss]->getStressResultant();
                        drill_dstress.addVector(1.0, m_nldrill->stress_comm, -1.0);
                        const Matrix& C0 = m_sections[igauss]->getInitialTangent();
                        drill_dstress_el.addMatrixVector(0.0, C0, drill_dstrain, 1.0);
                        double drill_damage = m_nldrill->damage_comm;
                        for (int j = 0; j < 8; ++j) {
                            double dx = drill_dstrain(j);
                            if (dx > 0.0) {
                                double dy = drill_dstress(j);
                                double dy0 = drill_dstress_el(j);
                                drill_damage = std::max(drill_damage, 1.0 - dy / dy0);
                            }
                        }
                        drill_damage = std::max(0.0, std::min(0.999, drill_damage));
                        m_nldrill->damage = drill_damage;
                    }
                }

                // Add drilling stress = Bd'*Sd * dA (RHS)
                double Sd = m_drill_stiffness * m_drill_strain;
                double Shx = m_drill_stiffness * m_drill_curvature_x * dh_scale_factor;
                double Shy = m_drill_stiffness * m_drill_curvature_y * dh_scale_factor;
                if (m_drill_mode == DrillingDOF_NonLinear) {
                    Sd *= (1.0 - m_nldrill->damage_comm);
                    Shx *= (1.0 - m_nldrill->damage_comm);
                    Shy *= (1.0 - m_nldrill->damage_comm);
                }
                for (int i = 0; i < 18; i++) {
                    RHS(i) += Bd(i) * Sd * dA;
                    RHS(i) += Bhx(i) * Shx * dA;
                    RHS(i) += Bhy(i) * Shy * dA;
                }
            }
        }

        // Integrate LHS
        if (options & OPT_LHS)
        {
            // Section tangent
            if (m_angle != 0.0) {
                Dsection = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
                if (m_damping[igauss]) 
                    Dsection *= m_damping[igauss]->getStiffnessMultiplier();
                auto& RsT = ASDShellT3Globals::instance().RsT;
                RsT.addMatrixTranspose(0.0, Rs, 1.0);
                auto& DRsT = ASDShellT3Globals::instance().DRsT;
                DRsT.addMatrixProduct(0.0, Dsection, RsT, 1.0);
                D.addMatrixProduct(0.0, Rs, DRsT, 1.0);
            }
            else {
                D = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
                if (m_damping[igauss])
                    D *= m_damping[igauss]->getStiffnessMultiplier();
            }

            // apply Stenberg stabilization
            for (int i = 0; i < 8; ++i) {
                D(6, i) *= sh_stab;
                D(7, i) *= sh_stab;
                D(i, 6) *= sh_stab;
                D(i, 7) *= sh_stab;
            }

            // Add current integration point contribution (LHS)
            B1TD.addMatrixTransposeProduct(0.0, B1, D, dA);
            LHS.addMatrixProduct(1.0, B1TD, B, 1.0);

            // Add drilling stiffness = Bd'*Kd*Bd * dA (LHS)
            if (m_reduced_integration) {
                double drill_tang = m_drill_stiffness;
                if (m_drill_mode == DrillingDOF_NonLinear)
                    drill_tang *= (1.0 - m_nldrill->damage_comm);
                for (int i = 0; i < 18; i++) {
                    for (int j = 0; j < 18; j++) {
                        LHS(i, j) += Bd(i) * drill_tang * Bd(j) * dA;
                        LHS(i, j) += Bhx(i) * drill_tang * Bhx(j) * dA * dh_scale_factor;
                        LHS(i, j) += Bhy(i) * drill_tang * Bhy(j) * dA * dh_scale_factor;
                    }
                }
            }
        }
    }

    // Transform LHS to global coordinate system
    m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));

    // Subtract external loads if any
    if ((options & OPT_RHS) && m_load)
        RHS.addVector(1.0, *m_load, -1.0);

    // Done
    return result;
}

double ASDShellT3::evaluateSectionThickness()
{
    const Matrix& D = m_sections[0]->getInitialTangent();
    double k00 = D(0, 0); // k00 = h * E / ( 1.0 - nu*nu ) -> membrane modulus
    if (k00 == 0.0) return 0.0;
    double k33 = -D(3, 3); // k33 = -E * (h * h * h) / 12.0 / (1.0 - nu * nu) -> bending modulus
    // k00/k33 = h/I = h/(h^2/12) = 12/h^2
    double h = std::sqrt(12.0 * k33 / k00);
    return h;
}

int
ASDShellT3::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(3, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
    }

    // set the quantity to be displayed at the nodes;
    static Vector values(3);
    for (int i = 0; i < 3; i++)
        values(i) = 0.0;

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}
