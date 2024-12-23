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
// $Date: 2020/05/18 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
// A 4-node general shell element based on the AGQ formulation 
// for the in-plane behavior, 
// and the MITC4 formulation for the out-of-plane behavior.
// 
// It supports both linear and corotational kinematics. Warped geometries
// can be modelled since this element is not assumed flat.
//

#include <ASDShellQ4.h>
#include <ASDShellQ4CorotationalTransformation.h>

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

void *
OPS_ASDShellQ4(void)
{
    static bool first_done = false;
    if (!first_done) {
        opserr << "Using ASDShellQ4 - Developed by: Massimo Petracca, Guido Camata, ASDEA Software Technology\n";
        first_done = true;
    }

    int numArgs = OPS_GetNumRemainingInputArgs();
    if (numArgs < 6) {
        opserr << "Want: element ASDShellQ4 $tag $iNode $jNode $kNode $lNode $secTag "
            "<-corotational> <-noeas> <-drillingStab $drillingStab> <-drillingNL>"
            "<-local $x1 $x2 $x3>";
        return 0;
    }

    int iData[6];
    int numData = 6;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: element ASDShellQ4 \n";
        return 0;
    }

    bool corotational = false;
    bool use_eas = true;
    bool use_drill_stab = false;
    double drilling_stab = 0.01;
    ASDShellQ4::DrillingDOFMode drill_mode = ASDShellQ4::DrillingDOF_Elastic;
    int dampingTag = 0;
    Damping *m_damping = 0;
    Vector local_x(3);
  
    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* type = OPS_GetString();
        if ((strcmp(type, "-corotational") == 0) || (strcmp(type, "-Corotational") == 0)) {
            corotational = true;
        }
        else if (strcmp(type, "-noeas") == 0) {
            use_eas = false;
        }
        else if (strcmp(type, "-drillingStab") == 0) {
            if (drill_mode != ASDShellQ4::DrillingDOF_Elastic) {
                opserr << "Error: element ASDShellQ4: -drillingStab and -drillingNL options are mutually exclusive\n";
                return 0;
            }
            if (OPS_GetNumRemainingInputArgs() == 0) {
                opserr << "Error: element ASDShellQ4: drilling stabilization parameter not provided with -drillingStab option\n";
                return 0;
            }
            numData = 1;
            if (OPS_GetDoubleInput(&numData, &drilling_stab) == 0) {
                drilling_stab = std::max(0.0, std::min(1.0, drilling_stab));
            }
            else {
                opserr << "Error: element ASDShellQ4: cannot get drilling stabilization parameter with -drillingStab option\n";
                return 0;
            }
            drill_mode = ASDShellQ4::DrillingDOF_Elastic;
            use_drill_stab = true;
        }
        else if (strcmp(type, "-drillingNL") == 0) {
            if (use_drill_stab) {
                opserr << "Error: element ASDShellQ4: -drillingStab and -drillingNL options are mutually exclusive\n";
                return 0;
            }
            drill_mode = ASDShellQ4::DrillingDOF_NonLinear;
            drilling_stab = 1.0;
        }
        else if (strcmp(type, "-local") == 0) {
            if (OPS_GetNumRemainingInputArgs() < 3) {
                opserr << "Error: element ASDShellQ4: not enough arguments for -local options (3 components are required)\n";
                return 0;
            }
            numData = 1;
            for (int i = 0; i < 3; ++i) {
                double local_x_com;
                if (OPS_GetDoubleInput(&numData, &local_x_com) == 0) {
                    local_x(i) = local_x_com;
                }
                else {
                    opserr << "Error: element ASDShellQ4: cannot get the component " << i + 1 << " for the local X axis\n";
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
                    opserr << "Error: element ASDShellQ4: damping not found\n";
                    return 0;
                }
            }
        }
    }
    SectionForceDeformation* section = OPS_getSectionForceDeformation(iData[5]);

    if (section == 0) {
        opserr << "ERROR:  element ASDShellQ4 " << iData[0] << " section " << iData[5] << " not found\n";
        return 0;
    }
    
    return new ASDShellQ4(iData[0], iData[1], iData[2], iData[3], iData[4], section, local_x, corotational, use_eas, drill_mode, drilling_stab, m_damping);
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

    // gauss quadrature data
    constexpr double GLOC = 0.577350269189626;
    constexpr std::array<double, 4> XI = { -GLOC, GLOC, GLOC, -GLOC };
    constexpr std::array<double, 4> ETA = { -GLOC, -GLOC, GLOC, GLOC };
    constexpr std::array<double, 4> WTS = { 1.0, 1.0, 1.0, 1.0 };

    // shape functions
    inline void shapeFunctions(double xi, double eta, Vector& N)
    {
        N(0) = 0.25 * (1.0 - xi) * (1.0 - eta);
        N(1) = 0.25 * (1.0 + xi) * (1.0 - eta);
        N(2) = 0.25 * (1.0 + xi) * (1.0 + eta);
        N(3) = 0.25 * (1.0 - xi) * (1.0 + eta);
    }

    // shape function derivatives in iso-parametric space
    inline void shapeFunctionsNaturalDerivatives(double xi, double eta, Matrix& dN)
    {
        dN(0, 0) = -(1.0 - eta) * 0.25;
        dN(1, 0) = (1.0 - eta) * 0.25;
        dN(2, 0) = (1.0 + eta) * 0.25;
        dN(3, 0) = -(1.0 + eta) * 0.25;

        dN(0, 1) = -(1.0 - xi) * 0.25;
        dN(1, 1) = -(1.0 + xi) * 0.25;
        dN(2, 1) = (1.0 + xi) * 0.25;
        dN(3, 1) = (1.0 - xi) * 0.25;
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

        void calculate(const ASDShellQ4LocalCoordinateSystem& CS, const Matrix& dN)
        {
            // jacobian
            J(0, 0) = dN(0, 0) * CS.X1() + dN(1, 0) * CS.X2() + dN(2, 0) * CS.X3() + dN(3, 0) * CS.X4();
            J(1, 0) = dN(0, 0) * CS.Y1() + dN(1, 0) * CS.Y2() + dN(2, 0) * CS.Y3() + dN(3, 0) * CS.Y4();
            J(0, 1) = dN(0, 1) * CS.X1() + dN(1, 1) * CS.X2() + dN(2, 1) * CS.X3() + dN(3, 1) * CS.X4();
            J(1, 1) = dN(0, 1) * CS.Y1() + dN(1, 1) * CS.Y2() + dN(2, 1) * CS.Y3() + dN(3, 1) * CS.Y4();

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

    /** \brief MITC4Params
     *
     * This class performs some operations and stores some data to compute
     * the transverse shear contribution of the stiffness matrix using the
     * M.I.T.C. formulation.
     *
     */
    struct MITC4Params
    {
        double Ax = 0.0;
        double Ay = 0.0;
        double Bx = 0.0;
        double By = 0.0;
        double Cx = 0.0;
        double Cy = 0.0;
        Matrix transformation = Matrix(2, 2);
        Matrix shearStrains = Matrix(4, 24);

        void compute(const ASDShellQ4LocalCoordinateSystem& LCS)
        {
            double x21 = LCS.X2() - LCS.X1();
            double y21 = LCS.Y2() - LCS.Y1();
            double x34 = LCS.X3() - LCS.X4();
            double y34 = LCS.Y3() - LCS.Y4();
            double x41 = LCS.X4() - LCS.X1();
            double y41 = LCS.Y4() - LCS.Y1();
            double x32 = LCS.X3() - LCS.X2();
            double y32 = LCS.Y3() - LCS.Y2();

            Ax = -LCS.X1() + LCS.X2() + LCS.X3() - LCS.X4();
            Bx = LCS.X1() - LCS.X2() + LCS.X3() - LCS.X4();
            Cx = -LCS.X1() - LCS.X2() + LCS.X3() + LCS.X4();
            Ay = -LCS.Y1() + LCS.Y2() + LCS.Y3() - LCS.Y4();
            By = LCS.Y1() - LCS.Y2() + LCS.Y3() - LCS.Y4();
            Cy = -LCS.Y1() - LCS.Y2() + LCS.Y3() + LCS.Y4();

            double Alpha = std::atan(Ay / Ax);
            double Beta = 3.141592653589793 * 0.5 - std::atan(Cx / Cy);

            transformation(0, 0) = std::sin(Beta);
            transformation(0, 1) = -std::sin(Alpha);
            transformation(1, 0) = -std::cos(Beta);
            transformation(1, 1) = std::cos(Alpha);

            shearStrains.Zero();

            shearStrains(0, 2) = -0.5;
            shearStrains(0, 3) = -y41 * 0.25;
            shearStrains(0, 4) = x41 * 0.25;

            shearStrains(0, 20) = 0.5;
            shearStrains(0, 21) = -y41 * 0.25;
            shearStrains(0, 22) = x41 * 0.25;

            shearStrains(1, 2) = -0.5;
            shearStrains(1, 3) = -y21 * 0.25;
            shearStrains(1, 4) = x21 * 0.25;

            shearStrains(1, 8) = 0.5;
            shearStrains(1, 9) = -y21 * 0.25;
            shearStrains(1, 10) = x21 * 0.25;

            shearStrains(2, 8) = -0.5;
            shearStrains(2, 9) = -y32 * 0.25;
            shearStrains(2, 10) = x32 * 0.25;

            shearStrains(2, 14) = 0.5;
            shearStrains(2, 15) = -y32 * 0.25;
            shearStrains(2, 16) = x32 * 0.25;

            shearStrains(3, 14) = 0.5;
            shearStrains(3, 15) = -y34 * 0.25;
            shearStrains(3, 16) = x34 * 0.25;

            shearStrains(3, 20) = -0.5;
            shearStrains(3, 21) = -y34 * 0.25;
            shearStrains(3, 22) = x34 * 0.25;
        }

    };

    /** \brief AGQIParams
     *
     * This class performs some operations and stores some data to compute
     * the AGQI enhancement of the membrane part
     *
     */
    struct AGQIParams
    {
        std::array<double, 4> X = { {0.0, 0.0, 0.0, 0.0} };
        std::array<double, 4> Y = { {0.0, 0.0, 0.0, 0.0} };
        std::array<double, 4> b = { {0.0, 0.0, 0.0, 0.0} };
        std::array<double, 4> c = { {0.0, 0.0, 0.0, 0.0} };
        double A1 = 0.0;
        double A2 = 0.0;
        double A3 = 0.0;
        double A = 0.0;
        std::array<double, 4> g = { {0.0, 0.0, 0.0, 0.0} };

        void compute(const ASDShellQ4LocalCoordinateSystem& LCS)
        {
            // extract the x and y coordinates
            // and compute bi = yj-yk, and ci = xk-xj
            // Eq 3
            for (int i = 0; i < 4; i++) {
                X[i] = LCS.X(i);
                Y[i] = LCS.Y(i);
            }
            for (int i = 0; i < 4; i++) {
                int j = i + 1; if (j > 3) j = 0;
                int k = j + 1; if (k > 3) k = 0;
                b[i] = Y[j] - Y[k];
                c[i] = X[k] - X[j];
            }

            // areas of sub-triangles
            A1 = 0.5 * (X[1] * Y[3] + X[0] * Y[1] + X[3] * Y[0] - X[1] * Y[0] - X[3] * Y[1] - X[0] * Y[3]);
            A2 = 0.5 * (X[1] * Y[2] + X[0] * Y[1] + X[2] * Y[0] - X[1] * Y[0] - X[2] * Y[1] - X[0] * Y[2]);
            A3 = 0.5 * (X[2] * Y[3] + X[1] * Y[2] + X[3] * Y[1] - X[2] * Y[1] - X[3] * Y[2] - X[1] * Y[3]);
            A = A1 + A3;
            // characteristic parameters of a quadrilateral
            // Eq 4
            g[0] = A1 / A; g[1] = A2 / A; g[2] = 1.0 - g[0]; g[3] = 1.0 - g[1];
        }

    };

    /** \brief ASDShellQ4Globals
     *
     * This singleton class stores some data for the shell calculations that
     * can be statically instantiated to avoid useless re-allocations
     *
     */
    class ASDShellQ4Globals
    {
    private:
        ASDShellQ4Globals() = default;

    public:

        JacobianOperator jac; // Jacobian
        MITC4Params mitc; // MITC4 parameters
        AGQIParams agq; // GQ12 parameters

        Vector UG = Vector(24); // global displacements
        Vector UL = Vector(24); // local displacements

        Matrix B = Matrix(8, 24); // strain-displacement matrix
        Matrix B1 = Matrix(8, 24); // strain-displacement matrix with inverted bending terms
        Matrix B1TD = Matrix(24, 8); // holds the B1^T*D terms
        Vector Bd = Vector(24); // strain-displacement matrix for drilling
        Vector Bd0 = Vector(24); // strain-displacement matrix for drilling (reduced integration)
        Vector N = Vector(4); // shape functions
        Matrix dN = Matrix(4, 2); // shape functions derivatives in isoparametric space
        Matrix dNdX = Matrix(4, 2); // shape functions derivatives in cartesian space
        Vector E = Vector(8); // strain vector
        Vector S = Vector(8); // stress vector
        Vector Elocal = Vector(8); // strain vector in local element coordinate system
        Matrix Re = Matrix(8, 8); // rotation matrix for strains
        Matrix Rs = Matrix(8, 8); // rotation matrix for stresses
        Matrix RsT = Matrix(8, 8); // transpose of above
        Matrix Dsection = Matrix(8, 8);
        Matrix D = Matrix(8, 8); // section tangent
        Matrix DRsT = Matrix(8, 8); // holds the product D * Rs^T

        Matrix BQ = Matrix(8, 4); // B matrix for internal DOFs
        Matrix BQ_mean = Matrix(8, 4); // Average B matrix for internal DOFs
        Matrix BQTD = Matrix(4, 8); // holds the BQ^T*D terms
        Matrix DBQ = Matrix(8, 4); // holds the D*BQ^T terms

        Matrix LHS = Matrix(24, 24); // LHS matrix (tangent stiffness)
        Matrix LHS_initial = Matrix(24, 24); // LHS matrix (initial stiffness)
        Matrix LHS_mass = Matrix(24, 24); // LHS matrix (mass matrix)
        Vector RHS = Vector(24); // RHS vector (residual vector)
        Vector RHS_winertia = Vector(24); // RHS vector (residual vector with inertia terms)

    public:
        static ASDShellQ4Globals& instance() {
            static ASDShellQ4Globals _instance;
            return _instance;
        }
    };

    // computes only the drilling B matrix
    void computeBdrilling(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        double xi, double eta,
        const JacobianOperator& Jac,
        const AGQIParams& agq,
        const Vector& N,
        const Matrix& dN,
        Vector& Bd,
        bool use_eas
    )
    {
        // cartesian derivatives of standard shape function
        auto& dNdX = ASDShellQ4Globals::instance().dNdX;
        dNdX.addMatrixProduct(0.0, dN, Jac.invJ, 1.0);

        // initialize
        Bd.Zero();

        // AGQI proc ***********************************************************************************************

        if (use_eas) {

            // area coordinates of the gauss point (Eq 7)
            std::array<double, 4> L;
            L[0] = 0.25 * (1.0 - xi) * (agq.g[1] * (1.0 - eta) + agq.g[2] * (1.0 + eta));
            L[1] = 0.25 * (1.0 - eta) * (agq.g[3] * (1.0 - xi) + agq.g[2] * (1.0 + xi));
            L[2] = 0.25 * (1.0 + xi) * (agq.g[0] * (1.0 - eta) + agq.g[3] * (1.0 + eta));
            L[3] = 0.25 * (1.0 + eta) * (agq.g[0] * (1.0 - xi) + agq.g[1] * (1.0 + xi));

            // computed modified shape function gradients for the strain matrix for external dofs
            // using the area coordinate method as written in the reference paper of the AGQI element
            static constexpr std::array<double, 4> NXI = { -1.0, 1.0, 1.0, -1.0 };
            static constexpr std::array<double, 4> NETA = { -1.0, -1.0, 1.0, 1.0 };
            for (int i = 0; i < 4; i++)
            {
                int j = i + 1; if (j > 3) j = 0;
                int k = j + 1; if (k > 3) k = 0;
                double SX = 0.0;
                double SY = 0.0;
                for (int ii = 0; ii < 4; ii++)
                {
                    int jj = ii + 1; if (jj > 3) jj = 0;
                    int kk = jj + 1; if (kk > 3) kk = 0;
                    int mm = kk + 1; if (mm > 3) mm = 0;
                    SX = SX + agq.b[ii] * NXI[ii] * NETA[ii] * (3.0 * (L[jj] - L[mm]) + (agq.g[jj] - agq.g[kk]));
                    SY = SY + agq.c[ii] * NXI[ii] * NETA[ii] * (3.0 * (L[jj] - L[mm]) + (agq.g[jj] - agq.g[kk]));
                }
                dNdX(i, 0) = (agq.b[i] + agq.b[j]) / agq.A / 2.0 +
                    NXI[i] * NETA[i] * agq.g[k] * SX / 2.0 / agq.A / (1.0 + agq.g[0] * agq.g[2] + agq.g[1] * agq.g[3]);
                dNdX(i, 1) = (agq.c[i] + agq.c[j]) / agq.A / 2.0 +
                    NXI[i] * NETA[i] * agq.g[k] * SY / 2.0 / agq.A / (1.0 + agq.g[0] * agq.g[2] + agq.g[1] * agq.g[3]);
            }
        }

        // We use the drilling penalty as defined by hughes and brezzi,
        // where we link the independent rotation to the skew symmetric part of the in-plane displacement field.

        Bd(0) = -0.5 * dNdX(0, 1);
        Bd(1) = 0.5 * dNdX(0, 0);
        Bd(5) = -N(0);

        Bd(6) = -0.5 * dNdX(1, 1);
        Bd(7) = 0.5 * dNdX(1, 0);
        Bd(11) = -N(1);

        Bd(12) = -0.5 * dNdX(2, 1);
        Bd(13) = 0.5 * dNdX(2, 0);
        Bd(17) = -N(2);

        Bd(18) = -0.5 * dNdX(3, 1);
        Bd(19) = 0.5 * dNdX(3, 0);
        Bd(23) = -N(3);
    }

    // computes the complete B matrix
    void computeBMatrix(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        double xi, double eta,
        const JacobianOperator& Jac, 
        const AGQIParams& agq,
        const MITC4Params& mitc,
        const Vector& N,
        const Matrix& dN,
        const Matrix& BQ_mean,
        Matrix& B,
        Matrix& BQ,
        Vector& Bd,
        bool use_eas
    )
    {
        // cartesian derivatives of standard shape function
        auto& dNdX = ASDShellQ4Globals::instance().dNdX;
        dNdX.addMatrixProduct(0.0, dN, Jac.invJ, 1.0);

        // initialize
        B.Zero();
        Bd.Zero();

        // AGQI proc ***********************************************************************************************

        if (use_eas) {

            // area coordinates of the gauss point (Eq 7)
            std::array<double, 4> L;
            L[0] = 0.25 * (1.0 - xi) * (agq.g[1] * (1.0 - eta) + agq.g[2] * (1.0 + eta));
            L[1] = 0.25 * (1.0 - eta) * (agq.g[3] * (1.0 - xi) + agq.g[2] * (1.0 + xi));
            L[2] = 0.25 * (1.0 + xi) * (agq.g[0] * (1.0 - eta) + agq.g[3] * (1.0 + eta));
            L[3] = 0.25 * (1.0 + eta) * (agq.g[0] * (1.0 - xi) + agq.g[1] * (1.0 + xi));

            // computed modified shape function gradients for the strain matrix for external dofs
            // using the area coordinate method as written in the reference paper of the AGQI element
            static constexpr std::array<double, 4> NXI = { -1.0, 1.0, 1.0, -1.0 };
            static constexpr std::array<double, 4> NETA = { -1.0, -1.0, 1.0, 1.0 };
            for (int i = 0; i < 4; i++)
            {
                int j = i + 1; if (j > 3) j = 0;
                int k = j + 1; if (k > 3) k = 0;
                double SX = 0.0;
                double SY = 0.0;
                for (int ii = 0; ii < 4; ii++)
                {
                    int jj = ii + 1; if (jj > 3) jj = 0;
                    int kk = jj + 1; if (kk > 3) kk = 0;
                    int mm = kk + 1; if (mm > 3) mm = 0;
                    SX = SX + agq.b[ii] * NXI[ii] * NETA[ii] * (3.0 * (L[jj] - L[mm]) + (agq.g[jj] - agq.g[kk]));
                    SY = SY + agq.c[ii] * NXI[ii] * NETA[ii] * (3.0 * (L[jj] - L[mm]) + (agq.g[jj] - agq.g[kk]));
                }
                dNdX(i, 0) = (agq.b[i] + agq.b[j]) / agq.A / 2.0 +
                    NXI[i] * NETA[i] * agq.g[k] * SX / 2.0 / agq.A / (1.0 + agq.g[0] * agq.g[2] + agq.g[1] * agq.g[3]);
                dNdX(i, 1) = (agq.c[i] + agq.c[j]) / agq.A / 2.0 +
                    NXI[i] * NETA[i] * agq.g[k] * SY / 2.0 / agq.A / (1.0 + agq.g[0] * agq.g[2] + agq.g[1] * agq.g[3]);
            }

            // compute the BQ matrix for the internal DOFs
            BQ.Zero();
            for (int i = 0; i < 2; i++)
            {
                unsigned int j = i + 1; if (j > 3) j = 0;
                unsigned int k = j + 1; if (k > 3) k = 0;
                double NQX = (agq.b[i] * L[k] + agq.b[k] * L[i]) / agq.A / 2.0;
                double NQY = (agq.c[i] * L[k] + agq.c[k] * L[i]) / agq.A / 2.0;
                unsigned int index1 = i * 2;
                unsigned int index2 = index1 + 1;
                BQ(0, index1) += NQX;
                BQ(1, index2) += NQY;
                BQ(2, index1) += NQY;
                BQ(2, index2) += NQX;
            }
            BQ.addMatrix(1.0, BQ_mean, -1.0);
        }

        // membrane ************************************************************************************************
        // this is the memebrane part of the compatible standard in-plane displacement field

        B(0, 0) = dNdX(0, 0);   B(0, 6) = dNdX(1, 0);   B(0, 12) = dNdX(2, 0);   B(0, 18) = dNdX(3, 0);
        B(1, 1) = dNdX(0, 1);   B(1, 7) = dNdX(1, 1);   B(1, 13) = dNdX(2, 1);   B(1, 19) = dNdX(3, 1);
        B(2, 0) = dNdX(0, 1);   B(2, 6) = dNdX(1, 1);   B(2, 12) = dNdX(2, 1);   B(2, 18) = dNdX(3, 1);
        B(2, 1) = dNdX(0, 0);   B(2, 7) = dNdX(1, 0);   B(2, 13) = dNdX(2, 0);   B(2, 19) = dNdX(3, 0);

        // We use the drilling penalty as defined by hughes and brezzi,
        // where we link the independent rotation to the skew symmetric part of the in-plane displacement field.

        Bd(0) = -0.5 * dNdX(0, 1);
        Bd(1) = 0.5 * dNdX(0, 0);
        Bd(5) = -N(0);
        
        Bd(6) = -0.5 * dNdX(1, 1);
        Bd(7) = 0.5 * dNdX(1, 0);
        Bd(11) = -N(1);
        
        Bd(12) = -0.5 * dNdX(2, 1);
        Bd(13) = 0.5 * dNdX(2, 0);
        Bd(17) = -N(2);
        
        Bd(18) = -0.5 * dNdX(3, 1);
        Bd(19) = 0.5 * dNdX(3, 0);
        Bd(23) = -N(3);

        // bending *************************************************************************************************

        B(3, 4) = -dNdX(0, 0);   B(3, 10) = -dNdX(1, 0);   B(3, 16) = -dNdX(2, 0);   B(3, 22) = -dNdX(3, 0);
        B(4, 3) =  dNdX(0, 1);   B(4, 9)  =  dNdX(1, 1);   B(4, 15) =  dNdX(2, 1);   B(4, 21) =  dNdX(3, 1);
        B(5, 3) =  dNdX(0, 0);   B(5, 9)  =  dNdX(1, 0);   B(5, 15) =  dNdX(2, 0);   B(5, 21) =  dNdX(3, 0);
        B(5, 4) = -dNdX(0, 1);   B(5, 10) = -dNdX(1, 1);   B(5, 16) = -dNdX(2, 1);   B(5, 22) = -dNdX(3, 1);

        // shear ***************************************************************************************************

        // MITC modified shape functions
        static Matrix MITCShapeFunctions(2, 4);
        MITCShapeFunctions.Zero();
        MITCShapeFunctions(1, 0) = 1.0 - xi;
        MITCShapeFunctions(0, 1) = 1.0 - eta;
        MITCShapeFunctions(1, 2) = 1.0 + xi;
        MITCShapeFunctions(0, 3) = 1.0 + eta;

        // strain displacement matrix in natural coordinate system.
        // interpolate the shear strains given in MITC4Params
        // using the modified shape function
        static Matrix BN(2, 24);
        BN.addMatrixProduct(0.0, MITCShapeFunctions, mitc.shearStrains, 1.0);

        // modify the shear strain intensity in the tying points
        // to match the values that would be obtained using standard
        // interpolations
        double Temp1, Temp2, Temp3;
        Temp1 = mitc.Cx + xi * mitc.Bx;
        Temp3 = mitc.Cy + xi * mitc.By;
        Temp1 = Temp1 * Temp1 + Temp3 * Temp3;
        Temp1 = std::sqrt(Temp1) / (8.0 * Jac.detJ);
        Temp2 = mitc.Ax + eta * mitc.Bx;
        Temp3 = mitc.Ay + eta * mitc.By;
        Temp2 = Temp2 * Temp2 + Temp3 * Temp3;
        Temp2 = std::sqrt(Temp2) / (8.0 * Jac.detJ);
        for (int i = 0; i < 24; i++) {
            BN(0, i) *= Temp1;
            BN(1, i) *= Temp2;
        }
        
        // transform the strain-displacement matrix from natural
        // to local coordinate system taking into account the element distortion
        static Matrix TBN(2, 24);
        TBN.addMatrixProduct(0.0, mitc.transformation, BN, 1.0);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 24; j++)
                B(i + 6, j) = TBN(i, j);
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
            for (int j = 0; j < 4; j++) {
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

ASDShellQ4::ASDShellQ4() 
    : Element(0, ELE_TAG_ASDShellQ4)
    , m_transformation(nullptr)
{
    for (int i = 0; i < 4; i++) m_damping[i] = nullptr;
}

ASDShellQ4::ASDShellQ4(
    int tag,
    int node1,
    int node2,
    int node3,
    int node4,
    SectionForceDeformation* section,
    const Vector& local_x,
    bool corotational,
    bool use_eas,
    DrillingDOFMode drill_mode,
    double drilling_stab,
    Damping *damping)
    : Element(tag, ELE_TAG_ASDShellQ4)
    , m_transformation(corotational ? new ASDShellQ4CorotationalTransformation() : new ASDShellQ4Transformation())
    , m_drill_mode(drill_mode)
    , m_drill_stab(drilling_stab)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;
    m_node_ids(3) = node4;

    // copy sections
    for (int i = 0; i < 4; i++) {
        m_sections[i] = section->getCopy();
        if (m_sections[i] == 0) {
            opserr << "ASDShellQ4::constructor - failed to get a material of type: ShellSection\n";
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

    // allocate eas data
    if (use_eas) {
        m_eas = new EASData();
    }

    // damping
    if (damping) {
        for (int i = 0; i < 4; i++) {
            m_damping[i] = damping->getCopy();
            if (!m_damping[i]) {
                opserr << "ASDShellQ4::constructor - failed to get copy of damping\n";
            }
        }
    }
    else {
        for (int i = 0; i < 4; i++) m_damping[i] = nullptr;
    }
}

ASDShellQ4::~ASDShellQ4( )
{
    // clean up section
    for (int i = 0; i < 4; i++)
        delete m_sections[i];

    // clean up coordinate transformation
    if(m_transformation)
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

    // clean up eas data
    if (m_eas)
        delete m_eas;

    // clean up damping
    for (int i = 0; i < 4; i++) {
        if (m_damping[i]) {
            delete m_damping[i];
            m_damping[i] = 0;
        }
    }
}

void  ASDShellQ4::setDomain(Domain* theDomain)
{
    // if domain is null
    if (theDomain == nullptr) {
        for (int i = 0; i < 4; i++)
            nodePointers[i] = nullptr;
        // set domain on transformation
        m_transformation->setDomain(theDomain, m_node_ids, m_initialized);
        // call base class implementation
        DomainComponent::setDomain(theDomain);
        return;
    }

    // node pointers
    for (int i = 0; i < 4; i++)
        nodePointers[i] = theDomain->getNode(m_node_ids(i));

    // set domain on transformation
    m_transformation->setDomain(theDomain, m_node_ids, m_initialized);

    // only if not already initialized from recvSelf
    if (!m_initialized) {

        // compute drilling stiffness
        m_drill_stiffness = 0.0;
        for (int i = 0; i < 4; i++)
            m_drill_stiffness += m_sections[i]->getInitialTangent()(2, 2);
        m_drill_stiffness /= 4.0; // average in-plane shear modulus
        m_drill_stiffness *= 0.233; // scale it to get the in-plane roto-distortional modulus

        // compute section orientation angle
        ASDShellQ4LocalCoordinateSystem reference_cs = m_transformation->createReferenceCoordinateSystem();
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
                opserr << "ASDShellQ4::setDomain Error: The provided local X axis cannot be aligned with the shell normal vector";
                exit(-1);
            }
            e1 = e2.cross(e3);
            e1.normalize();
        }
        else {
            // default one
            Vector3Type P1(m_transformation->getNodes()[0]->getCrds());
            Vector3Type P2(m_transformation->getNodes()[1]->getCrds());
            Vector3Type P3(m_transformation->getNodes()[2]->getCrds());
            Vector3Type P4(m_transformation->getNodes()[3]->getCrds());
            Vector3Type e1 = (P2 + P3) / 2.0 - (P1 + P4) / 2.0;
            e1.normalize();
        }
        m_angle = std::acos(std::max(-1.0, std::min(1.0, e1.dot(e1_local))));
        if (m_angle != 0.0) {
            // if they are not counter-clock-wise, let's change the sign of the angle
            const Matrix& R = reference_cs.Orientation();
            if ((e1(0) * R(1, 0) + e1(1) * R(1, 1) + e1(2) * R(1, 2)) < 0.0)
                m_angle = -m_angle;
        }

        // AGQI internal DOFs
        if(m_eas)
            AGQIinitialize();

        for (int i = 0; i < 4; i++) {
          if (m_damping[i] && m_damping[i]->setDomain(theDomain, 8)) {
            opserr << "ASDShellQ4::setDomain -- Error initializing damping\n";
            exit(-1);
          }
        }

        // initialized
        m_initialized = true;
    }

    // call base class implementation
    DomainComponent::setDomain(theDomain);
}

void ASDShellQ4::Print(OPS_Stream& s, int flag)
{
    if (flag == -1) {
        int eleTag = this->getTag();
        s << "EL_ASDShellQ4\t" << eleTag << "\t";
        s << eleTag << "\t" << 1;
        s << "\t" << m_node_ids(0) << "\t" << m_node_ids(1);
        s << "\t" << m_node_ids(2) << "\t" << m_node_ids(3) << "\t0.00";
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
        for (i = 0; i < 4; i++) {
            const Vector& stress = m_sections[i]->getStressResultant();
            s << "STRESS\t" << eleTag << "\t" << counter << "\t" << i << "\tTOP";
            for (j = 0; j < 6; j++)
                s << "\t" << stress(j);
            s << endln;
        }
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "MITC4 Non-Locking Four Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << m_node_ids(0) << endln;
        s << "Node 2 : " << m_node_ids(1) << endln;
        s << "Node 3 : " << m_node_ids(2) << endln;
        s << "Node 4 : " << m_node_ids(3) << endln;

        s << "Material Information : \n ";
        m_sections[0]->Print(s, flag);

        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDShellQ4\", ";
        s << "\"nodes\": [" << m_node_ids(0) << ", " << m_node_ids(1) << ", ";
        s << m_node_ids(2) << ", " << m_node_ids(3) << "], ";
        s << "\"section\": \"" << m_sections[0]->getTag() << "\"}";
    }
}

int
ASDShellQ4::setDamping(Domain *theDomain, Damping *damping)
{
  if (theDomain && damping)
  {
    for (int i = 0; i < 4; i++) {
      if (m_damping[i]) delete m_damping[i];
      m_damping[i] = damping->getCopy();
      if (!m_damping[i]) {
        opserr << "ASDShellQ4::setDamping - failed to get copy of damping\n";
        return -1;
      }
      if (m_damping[i] && m_damping[i]->setDomain(theDomain, 8)) {
        opserr << "ASDShellQ4::setDamping -- Error initializing damping\n";
        return -2;
      }
    }
  }
  return 0;
}

int  ASDShellQ4::getNumExternalNodes() const
{
    return 4;
}

const ID& ASDShellQ4::getExternalNodes()
{
    return m_node_ids;
}

Node**
ASDShellQ4::getNodePtrs(void)
{
    return m_transformation->getNodes().data();
}

int  ASDShellQ4::getNumDOF()
{
    return 24;
}

int  ASDShellQ4::commitState()
{
    int success = 0;

    // transformation
    m_transformation->commit();

    // sections
    for (int i = 0; i < 4; i++)
        success += m_sections[i]->commitState();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        for (int i = 0; i < 4; i++) {
            m_nldrill->stress_comm[i] = m_sections[i]->getStressResultant();
            m_nldrill->strain_comm[i] = m_sections[i]->getSectionDeformation();
            m_nldrill->damage_comm[i] = m_nldrill->damage[i];
        }
    }

    // AGQI internal DOFs
    if(m_eas) {
        m_eas->U_converged = m_eas->U;
        m_eas->Q_converged = m_eas->Q;
    }

    // damping
    for (int i = 0; i < 4; i++ )
      if (m_damping[i]) success += m_damping[i]->commitState();

    // done
    return success;
}

int  ASDShellQ4::revertToLastCommit()
{
    int success = 0;

    // transformation
    m_transformation->revertToLastCommit();

    // sections
    for (int i = 0; i < 4; i++)
        success += m_sections[i]->revertToLastCommit();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        for (int i = 0; i < 4; i++) {
            m_nldrill->stress_comm[i] = m_sections[i]->getStressResultant();
            m_nldrill->strain_comm[i] = m_sections[i]->getSectionDeformation();
            m_nldrill->damage[i] = m_nldrill->damage_comm[i];
        }
    }

    // AGQI internal DOFs
    if (m_eas) {
        m_eas->U = m_eas->U_converged;
        m_eas->Q = m_eas->Q_converged;
    }

    // damping
    for (int i = 0; i < 4; i++ )
      if (m_damping[i]) success += m_damping[i]->revertToLastCommit();

    // done
    return success;
}

int  ASDShellQ4::revertToStart()
{
    int success = 0;

    // transformation
    m_transformation->revertToStart();

    // sections
    for (int i = 0; i < 4; i++)
        success += m_sections[i]->revertToStart();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        for (int i = 0; i < 4; i++) {
            m_nldrill->stress_comm[i].Zero();
            m_nldrill->strain_comm[i].Zero();
            m_nldrill->damage[i] = m_nldrill->damage_comm[i] = 0.0;
        }
    }

    // AGQI internal DOFs
    if(m_eas)
        AGQIinitialize();

    // damping
    for (int i = 0; i < 4; i++ )
      if (m_damping[i]) success += m_damping[i]->revertToStart();

    return success;
}

int ASDShellQ4::update()
{
    // calculate
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    return calculateAll(LHS, RHS, (OPT_UPDATE));
}

const Matrix& ASDShellQ4::getTangentStiff()
{
    // calculate
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_LHS));
    return LHS;
}

const Matrix& ASDShellQ4::getInitialStiff()
{
    // calculate
    auto& LHS = ASDShellQ4Globals::instance().LHS_initial;
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_LHS | OPT_LHS_IS_INITIAL));
    return LHS;
}

const Matrix& ASDShellQ4::getMass()
{
    // Output matrix
    auto& LHS = ASDShellQ4Globals::instance().LHS_mass;
    LHS.Zero();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& dN = ASDShellQ4Globals::instance().dN;
    auto& N = ASDShellQ4Globals::instance().N;

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctions(xi, eta, N);
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Section density (already integrated through the thickness)
        double rho = m_sections[i]->getRho();

        // Add current integration point contribution
        for (int j = 0; j < 4; j++)
        {
            int index = j * 6;

            // Translational mass contribution
            double Tmass = N(j) * rho * dA;
            for (int q = 0; q < 3; q++)
                LHS(index + q, index + q) += Tmass;

            // Rotational mass neglected for the moment ...
        }
    }

    // Done
    return LHS;
}

void  ASDShellQ4::zeroLoad()
{
    if (m_load)
        m_load->Zero();
}

int
ASDShellQ4::addLoad(ElementalLoad* theLoad, double loadFactor)
{
    opserr << "ASDShellQ4::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
    return -1;
}

int
ASDShellQ4::addInertiaLoadToUnbalance(const Vector& accel)
{
    // Allocate load vector if necessary
    if (m_load == nullptr)
        m_load = new Vector(24);
    auto& F = *m_load;

    // Get mass matrix
    const auto& M = getMass();

    // Add -M*R*acc to unbalance, taking advantage of lumped mass matrix
    for (int i = 0; i < 4; i++)
    {
        const auto& RV = m_transformation->getNodes()[i]->getRV(accel);
        int index = i * 6;
        for (int j = 0; j < 6; j++)
            F(index + j) -= M(index + j, index + j) * RV(j);
    }

    // Done
    return 0;
}

const Vector& ASDShellQ4::getResistingForce()
{
    // calculate
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    calculateAll(LHS, RHS, (OPT_RHS));
    return RHS;
}

const Vector& ASDShellQ4::getResistingForceIncInertia()
{
    // calculate
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    auto& RHS = ASDShellQ4Globals::instance().RHS_winertia;
    calculateAll(LHS, RHS, (OPT_RHS));

    // Add damping terms
    if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
        RHS.addVector(1.0, getRayleighDampingForces(), 1.0);

    // Compute mass
    const auto& M = getMass();

    // Add M*acc to unbalance, taking advantage of lumped mass matrix
    for (int i = 0; i < 4; i++)
    {
        const auto& A = m_transformation->getNodes()[i]->getTrialAccel();
        int index = i * 6;
        for (int j = 0; j < 6; j++)
            RHS(index + j) += M(index + j, index + j) * A(j);
    }

    // Done
    return RHS;
}

int  ASDShellQ4::sendSelf(int commitTag, Channel& theChannel)
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
    // 4 node tags +
    // 1 transformation type
    // 1 initialization flag
    // 1 has_load_flag
    // 8 -> 4 pairs of (section class tag, mt db tag)
    // 2 -> damping tag + damping db tag
    // 1 -> EAS flag
    // 1 -> non-linear drilling flag
    // 1 -> local_x flag
    static ID idData(21);
    counter = 0;
    idData(counter++) = getTag();
    for(int i = 0; i < 4; ++i)
        idData(counter++) = m_node_ids(i);
    idData(counter++) = static_cast<int>(m_transformation->isLinear());
    idData(counter++) = static_cast<int>(m_initialized);
    idData(counter++) = static_cast<int>(has_load);
    for (int i = 0; i < 4; i++) {
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
                for (int i = 0; i < 4; i++)
                    m_damping[i]->setDbTag(dbTag);
        }
        idData(counter++) = dbTag;
    }
    else {
        idData(counter++) = 0;
        idData(counter++) = 0;
    }
    idData(counter++) = static_cast<int>(static_cast<bool>(m_eas));
    idData(counter++) = static_cast<int>(static_cast<bool>(m_nldrill));
    idData(counter++) = static_cast<int>(static_cast<bool>(m_local_x));

    res = theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send ID\n";
        return res;
    }

    // DOUBLE data
    // 4 damping factors +
    // 4 drilling strain +
    // 1 drilling stiffness +
    // 1 drilling stabilization parameter +
    // 1 angle +
    // (optional) 3 local_x +
    // (optional) 268 AGQI internal data + 
    // (optional) 72 non-linear drilling data +
    // (optional) 24 load +
    // NT transformation data
    int Nlocalx = m_local_x ? 3 : 0;
    int NLoad = has_load ? 24 : 0;
    int Neas = m_eas ? 268 : 0;
    int Nnldrill = m_nldrill ? 72 : 0;
    int NT = m_transformation->internalDataSize();
    Vector vectData(4 + 4 + 1 + 1 + 1 + Nlocalx + Neas + Nnldrill + NLoad + NT);
    counter = 0;
    vectData(counter++) = alphaM;
    vectData(counter++) = betaK;
    vectData(counter++) = betaK0;
    vectData(counter++) = betaKc;
    for (int i = 0; i < 4; ++i)
        vectData(counter++) = m_drill_strain[i];
    vectData(counter++) = m_drill_stiffness;
    vectData(counter++) = m_drill_stab;
    vectData(counter++) = m_angle;
    if (m_local_x) {
        for (int i = 0; i < 3; ++i)
            vectData(counter++) = (*m_local_x)(i);
    }
    if (m_eas) {
        for (int i = 0; i < 4; ++i)
            vectData(counter++) = m_eas->Q(i);
        for (int i = 0; i < 4; ++i)
            vectData(counter++) = m_eas->Q_converged(i);
        for (int i = 0; i < 24; ++i)
            vectData(counter++) = m_eas->U(i);
        for (int i = 0; i < 24; ++i)
            vectData(counter++) = m_eas->U_converged(i);
        for (int i = 0; i < 4; ++i)
            vectData(counter++) = m_eas->Q_residual(i);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                vectData(counter++) = m_eas->KQQ_inv(i, j);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 24; ++j)
                vectData(counter++) = m_eas->KQU(i, j);
        for (int i = 0; i < 24; ++i)
            for (int j = 0; j < 4; ++j)
                vectData(counter++) = m_eas->KUQ(i, j);
    }
    if (m_nldrill) {
        for (const auto& item : m_nldrill->strain_comm)
            for (int i = 0; i < 8; ++i)
                vectData(counter++) = item(i);
        for (const auto& item : m_nldrill->stress_comm)
            for (int i = 0; i < 8; ++i)
                vectData(counter++) = item(i);
        for (auto item : m_nldrill->damage)
            vectData(counter++) = item;
        for (auto item : m_nldrill->damage_comm)
            vectData(counter++) = item;
    }
    if (has_load) {
        for (int i = 0; i < 24; ++i)
            vectData(counter++) = (*m_load)(i);
    }
    m_transformation->saveInternalData(vectData, counter);

    res = theChannel.sendVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send Vector\n";
        return res;
    }

    // send all sections
    for (int i = 0; i < 4; i++) {
        res = m_sections[i]->sendSelf(commitTag, theChannel);
        if (res < 0) {
            opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send its Material\n";
            return res;
        }
    }

    // Ask the Damping to send itself
    if (m_damping[0]) {
      for (int i = 0 ;  i < 4; i++) {
        res += m_damping[i]->sendSelf(commitTag, theChannel);
        if (res < 0) {
          opserr << "ASDShellQ4::sendSelf -- could not send Damping\n";
          return res;
        }
      }
    }

    // done
    return res;
}

int  ASDShellQ4::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    int dataTag = this->getDbTag();

    // a counter
    int counter;

    // INT data
    // 1 tag +
    // 4 node tags +
    // 1 transformation type
    // 1 initialization flag
    // 1 has_load_flag
    // 8 -> 4 pairs of (section class tag, mt db tag)
    // 2 -> damping tag + damping db tag
    // 1 -> EAS flag
    // 1 -> non-linear drilling flag
    // 1 -> local_x flag
    static ID idData(21);
    res = theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
        opserr << "WARNING ASDShellQ4::recvSelf() - " << this->getTag() << " failed to receive ID\n";
        return res;
    }

    counter = 0;
    setTag(idData(counter++));
    for (int i = 0; i < 4; ++i)
        m_node_ids(i) = idData(counter++);
    bool linear_transform = static_cast<bool>(idData(counter++));
    m_initialized = static_cast<bool>(idData(counter++));
    bool has_load = static_cast<bool>(idData(counter++));
    // create sections
    for (int i = 0; i < 4; i++) {
        int matClassTag = idData(counter++);
        int matDbTag = idData(counter++);
        if (m_sections[i])
            delete m_sections[i];
        m_sections[i] = theBroker.getNewSection(matClassTag);
        if (m_sections[i] == 0) {
            opserr << "ASDShellQ4::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
            return -1;
        }
        m_sections[i]->setDbTag(matDbTag);
    }
    int dmpTag = idData(counter++);
    int dmpDbTag = idData(counter++);
    bool use_eas = static_cast<bool>(idData(counter++));
    bool use_nldrill = static_cast<bool>(idData(counter++));
    bool use_local_x = static_cast<bool>(idData(counter++));
    
    // create transformation
    if (m_transformation)
        delete m_transformation;
    m_transformation = linear_transform ? new ASDShellQ4Transformation() : new ASDShellQ4CorotationalTransformation();
    // create load
    if (has_load) {
        if (m_load == nullptr)
            m_load = new Vector(24);
    }
    else {
        if (m_load) {
            delete m_load;
            m_load = nullptr;
        }
    }

    // create eas
    if (m_eas) {
        delete m_eas;
        m_eas = nullptr;
    }
    if (use_eas)
        m_eas = new EASData();

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
    // 4 drilling strain +
    // 1 drilling stiffness +
    // 1 drilling stabilization parameter +
    // 1 angle +
    // (optional) 3 local_x +
    // (optional) 268 AGQI internal data + 
    // (optional) 72 non-linear drilling data +
    // (optional) 24 load +
    // NT transformation data
    int Nlocalx = use_local_x ? 3 : 0;
    int NLoad = has_load ? 24 : 0;
    int Neas = use_eas ? 268 : 0;
    int Nnldrill = use_nldrill ? 72 : 0;
    int NT = m_transformation->internalDataSize();
    Vector vectData(4 + 4 + 1 + 1 + 1 + Neas + Nnldrill + NLoad + NT);
    res = theChannel.recvVector(dataTag, commitTag, vectData);
    if (res < 0) {
        opserr << "WARNING ASDShellQ4::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
        return res;
    }

    counter = 0;
    alphaM = vectData(counter++);
    betaK = vectData(counter++);
    betaK0 = vectData(counter++);
    betaKc = vectData(counter++);
    for (int i = 0; i < 4; ++i)
        m_drill_strain[i] = vectData(counter++);
    m_drill_stiffness = vectData(counter++);
    m_drill_stab = vectData(counter++);
    m_angle = vectData(counter++);
    if (use_local_x) {
        for (int i = 0; i < 3; ++i)
            (*m_local_x)(i) = vectData(counter++);
    }
    if (use_eas) {
        for (int i = 0; i < 4; ++i)
            m_eas->Q(i) = vectData(counter++);
        for (int i = 0; i < 4; ++i)
            m_eas->Q_converged(i) = vectData(counter++);
        for (int i = 0; i < 24; ++i)
            m_eas->U(i) = vectData(counter++);
        for (int i = 0; i < 24; ++i)
            m_eas->U_converged(i) = vectData(counter++);
        for (int i = 0; i < 4; ++i)
            m_eas->Q_residual(i) = vectData(counter++);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                m_eas->KQQ_inv(i, j) = vectData(counter++);
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 24; ++j)
                m_eas->KQU(i, j) = vectData(counter++);
        for (int i = 0; i < 24; ++i)
            for (int j = 0; j < 4; ++j)
                m_eas->KUQ(i, j) = vectData(counter++);
    }
    if (use_nldrill) {
        for (auto& item : m_nldrill->strain_comm)
            for (int i = 0; i < 8; ++i)
                item(i) = vectData(counter++);
        for (auto& item : m_nldrill->stress_comm)
            for (int i = 0; i < 8; ++i)
                item(i) = vectData(counter++);
        for (auto& item : m_nldrill->damage)
            item = vectData(counter++);
        for (auto& item : m_nldrill->damage_comm)
            item = vectData(counter++);
    }
    if (has_load) {
        for (int i = 0; i < 24; ++i)
            (*m_load)(i) = vectData(counter++);
    }
    m_transformation->restoreInternalData(vectData, counter);

    // all sections
    for (int i = 0; i < 4; i++) {
        res = m_sections[i]->recvSelf(commitTag, theChannel, theBroker);
        if (res < 0) {
            opserr << "ASDShellQ4::recvSelf() - material " << i << "failed to recv itself\n";
            return res;
        }
    }

  if (dmpTag) {
    for (int i = 0 ;  i < 4; i++) {
      // Check if the Damping is null; if so, get a new one
      if (m_damping[i] == 0) {
        m_damping[i] = theBroker.getNewDamping(dmpTag);
        if (m_damping[i] == 0) {
          opserr << "ASDShellQ4::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Check that the Damping is of the right type; if not, delete
      // the current one and get a new one of the right type
      if (m_damping[i]->getClassTag() != dmpTag) {
        delete m_damping[i];
        m_damping[i] = theBroker.getNewDamping(dmpTag);
        if (m_damping[i] == 0) {
          opserr << "ASDShellQ4::recvSelf -- could not get a Damping\n";
          exit(-1);
        }
      }
  
      // Now, receive the Damping
      m_damping[i]->setDbTag(dmpDbTag);
      res += m_damping[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "ASDShellQ4::recvSelf -- could not receive Damping\n";
        return res;
      }
    }
  }
  else {
    for (int i = 0; i < 4; i++)
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
ASDShellQ4::setResponse(const char **argv, int argc, OPS_Stream &output)
{
    Response* theResponse = 0;

    output.tag("ElementOutput");
    output.attr("eleType", "ASDShellQ4");
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
            opserr << "ASDShellQ4::setResponse() - need to specify more data\n";
            return 0;
        }
        int pointNum = atoi(argv[1]);
        if (pointNum > 0 && pointNum <= 4) {

            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("eta", XI[pointNum - 1]);
            output.attr("neta", ETA[pointNum - 1]);

            theResponse = m_sections[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

            output.endTag();
        }

    }
    else if (strcmp(argv[0], "stresses") == 0) {

        for (int i = 0; i < 4; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", XI[i]);
            output.attr("neta", ETA[i]);

            output.tag("SectionForceDeformation");
            output.attr("classType", m_sections[i]->getClassTag());
            output.attr("tag", m_sections[i]->getTag());

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

        theResponse = new ElementResponse(this, 2, Vector(32));
    }

    else if (strcmp(argv[0], "strains") == 0) {

        for (int i = 0; i < 4; i++) {
            output.tag("GaussPoint");
            output.attr("number", i + 1);
            output.attr("eta", XI[i]);
            output.attr("neta", ETA[i]);

            output.tag("SectionForceDeformation");
            output.attr("classType", m_sections[i]->getClassTag());
            output.attr("tag", m_sections[i]->getTag());

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

        theResponse = new ElementResponse(this, 3, Vector(32));
    }

    else if (m_damping[0] && strcmp(argv[0],"dampingStresses") ==0) {
  
      for (int i=0; i<4; i++) {
        output.tag("GaussPoint");
        output.attr("number",i+1);
        output.attr("eta",XI[i]);
        output.attr("neta",ETA[i]);
        
        output.tag("SectionForceDeformation");
        output.attr("classType", m_damping[i]->getClassTag());
        output.attr("tag", m_damping[i]->getTag());
        
        output.tag("ResponseType","p11");
        output.tag("ResponseType","p22");
        output.tag("ResponseType","p12");
        output.tag("ResponseType","m11");
        output.tag("ResponseType","m22");
        output.tag("ResponseType","m12");
        output.tag("ResponseType","q1");
        output.tag("ResponseType","q2");
        
        output.endTag(); // GaussPoint
        output.endTag(); // NdMaterialOutput
      }
      
      theResponse =  new ElementResponse(this, 4, Vector(32));
    }

    output.endTag();
    return theResponse;
}

int
ASDShellQ4::getResponse(int responseID, Information &eleInfo)
{
    static Vector stresses(32);
    static Vector strains(32);

    switch (responseID) {
    case 1: // global forces
        return eleInfo.setVector(this->getResistingForce());
        break;

    case 2: // stresses
        for (int i = 0; i < 4; i++) {

            // Get material stress response
            const Vector& sigma = m_sections[i]->getStressResultant();
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
        for (int i = 0; i < 4; i++) {

            // Get section deformation
            const Vector& deformation = m_sections[i]->getSectionDeformation();
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
      for (int i = 0; i < 4; i++) {
  
        // Get material stress response
        const Vector &sigma = m_damping[i]->getDampingForce();
        stresses(0) = sigma(0);
        stresses(1) = sigma(1);
        stresses(2) = sigma(2);
        stresses(3) = sigma(3);
        stresses(4) = sigma(4);
        stresses(5) = sigma(5);
        stresses(6) = sigma(6);
        stresses(7) = sigma(7);
      }
      return eleInfo.setVector(stresses);
      break;
    default:
        return -1;
    }
}

int ASDShellQ4::setParameter(const char** argv, int argc, Parameter& param)
{
    int res = -1;
    int matRes = res;
    for (int i = 0; i < 4; i++)
    {
        matRes = m_sections[i]->setParameter(argv, argc, param);
        if (matRes != -1)
            res = matRes;
    }
    return res;
}

double ASDShellQ4::getCharacteristicLength(void)
{
    // by default we return the min distance between 2 nodes, which should coincide
    // with lch = sqrt(A) in non-distorted elements.
    // however, when the EAS is used, the element localizes only in a row of gauss points
    // so the lch should be half the original one
    double lch = Element::getCharacteristicLength();
    if (m_eas)
        lch /= 2.0;
    return lch;
}

int ASDShellQ4::calculateAll(Matrix& LHS, Vector& RHS, int options)
{
    // Check options
    if (!m_transformation->isLinear()) {
        // corotational calculation of the tangent LHS requires the RHS!
        if (options & OPT_LHS)
            options |= OPT_RHS;
    }

    // Zero output
    int result = 0;
    if(options & OPT_RHS)
        RHS.Zero();
    if(options & OPT_LHS)
        LHS.Zero();

    // Global displacements
    auto& UG = ASDShellQ4Globals::instance().UG;
    m_transformation->computeGlobalDisplacements(UG);

    // update transformation
    if(options & OPT_UPDATE)
        m_transformation->update(UG);

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Compute the local coordinate system.
    ASDShellQ4LocalCoordinateSystem local_cs =
        m_transformation->createLocalCoordinateSystem(UG);

    // Prepare all the parameters needed for the MITC4
    // and AGQI formulations.
    // This is to be done here outside the Gauss Loop.
    auto& mitc = ASDShellQ4Globals::instance().mitc;
    auto& agq = ASDShellQ4Globals::instance().agq;
    mitc.compute(reference_cs);
    if(m_eas)
        agq.compute(reference_cs);

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices/vectors
    auto& B = ASDShellQ4Globals::instance().B;
    auto& Bd = ASDShellQ4Globals::instance().Bd;
    auto& Bd0 = ASDShellQ4Globals::instance().Bd0;
    auto& BQ = ASDShellQ4Globals::instance().BQ;
    auto& BQ_mean = ASDShellQ4Globals::instance().BQ_mean;
    auto& BQTD = ASDShellQ4Globals::instance().BQTD;
    auto& DBQ = ASDShellQ4Globals::instance().DBQ;
    auto& B1 = ASDShellQ4Globals::instance().B1;
    auto& B1TD = ASDShellQ4Globals::instance().B1TD;
    auto& N = ASDShellQ4Globals::instance().N;
    auto& dN = ASDShellQ4Globals::instance().dN;
    auto& E = ASDShellQ4Globals::instance().E;
    auto& S = ASDShellQ4Globals::instance().S;
    auto& D = ASDShellQ4Globals::instance().D;
    auto& Dsection = ASDShellQ4Globals::instance().Dsection;

    // matrices for orienting strains in section coordinate system
    auto& Re = ASDShellQ4Globals::instance().Re;
    auto& Rs = ASDShellQ4Globals::instance().Rs;
    if (m_angle != 0.0) {
        if (options & OPT_UPDATE)
            getRotationMatrixForGeneralizedStrains(-m_angle, Re);
        if ((options & OPT_RHS) || (options & OPT_LHS))
            getRotationMatrixForGeneralizedStresses(m_angle, Rs);
    }

    // Local displacements
    auto& UL = ASDShellQ4Globals::instance().UL;
    m_transformation->calculateLocalDisplacements(local_cs, UG, UL);

    // Update AGQI internal DOFs
    if (options & OPT_UPDATE)
        if(m_eas)
            AGQIupdate(UL);

    // AGQI begin gauss loop (computes the BQ_mean correction matrix)
    if(m_eas)
        AGQIbeginGaussLoop(reference_cs);

    // Drilling strain-displacement matrix at center for reduced integration
    shapeFunctions(0.0, 0.0, N);
    shapeFunctionsNaturalDerivatives(0.0, 0.0, dN);
    computeBdrilling(reference_cs, 0.0, 0.0, jac, agq, N, dN, Bd0, m_eas);
    static Vector drill_dstrain(8);
    static Vector drill_dstress(8);
    static Vector drill_dstress_el(8);

    // Gauss loop
    for (int igauss = 0; igauss < 4; igauss++)
    {
        // Current integration point data
        double xi = XI[igauss];
        double eta = ETA[igauss];
        double w = WTS[igauss];
        shapeFunctions(xi, eta, N);
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Strain-displacement matrix
        computeBMatrix(reference_cs, xi, eta, jac, agq, mitc, N, dN, BQ_mean, B, BQ, Bd, m_eas);

        // The drilling according to Hughes-Brezzi plays and important role in warped shells and
        // the stiffness should be proportional to the initial (elastic) in-plane shear modulus.
        // If fully integrated however, and in case of strain-softening materials,
        // it can over-stiffen the membrane response. If under-integrated the element
        // has spurious zero-energy modes. Thus here we use the compute the drilling B matrix
        // as a weighted sum of the reduced integrated one (Bd0) and a scaled  version of the fully
        // integrated one, just to suppress spurious zero energy modes.
        // Kd0 = Bd0'*Bd0*kdrill*V
        // Kd = int(Bd'*Bd*kdrill*dV)
        // K = stab*Kd + (1-stab)*Kd0
        double drill_alpha = std::sqrt(m_drill_stab);
        double drill_beta = std::sqrt(1.0 - m_drill_stab);
        Bd.addVector(drill_alpha, Bd0, drill_beta);

        // Update strain strain
        if (options & OPT_UPDATE)
        {
            // Section deformation
            if (m_angle != 0.0) {
                auto& Elocal = ASDShellQ4Globals::instance().Elocal;
                Elocal.addMatrixVector(0.0, B, UL, 1.0);
                if (m_eas)
                    Elocal.addMatrixVector(1.0, BQ, m_eas->Q, 1.0);
                E.addMatrixVector(0.0, Re, Elocal, 1.0);
            }
            else {
                E.addMatrixVector(0.0, B, UL, 1.0);
                if (m_eas)
                    E.addMatrixVector(1.0, BQ, m_eas->Q, 1.0);
            }

            // Update section
            result += m_sections[igauss]->setTrialSectionDeformation(E);

            // Drilling strain Ed = Bd*UL
            double Ed = 0.0;
            for (int i = 0; i < 24; i++)
                Ed += Bd(i) * UL(i);
            m_drill_strain[igauss] = Ed;
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

            // Add current integration point contribution (RHS)
            RHS.addMatrixTransposeVector(1.0, B1, S, dA);

            // Compute drilling damages
            if (m_drill_mode == DrillingDOF_NonLinear) {
                drill_dstrain = m_sections[igauss]->getSectionDeformation();
                drill_dstrain.addVector(1.0, m_nldrill->strain_comm[igauss], -1.0);
                if (drill_dstrain.Norm() > 1.0e-10) {
                    drill_dstress = m_sections[igauss]->getStressResultant();
                    drill_dstress.addVector(1.0, m_nldrill->stress_comm[igauss], -1.0);
                    const Matrix& C0 = m_sections[igauss]->getInitialTangent();
                    drill_dstress_el.addMatrixVector(0.0, C0, drill_dstrain, 1.0);
                    double drill_damage = m_nldrill->damage_comm[igauss];
                    for (int j = 0; j < 8; ++j) {
                        double dx = drill_dstrain(j);
                        if (dx > 0.0) {
                            double dy = drill_dstress(j);
                            double dy0 = drill_dstress_el(j);
                            drill_damage = std::max(drill_damage, 1.0 - dy / dy0);
                        }
                    }
                    drill_damage = std::max(0.0, std::min(0.999, drill_damage));
                    m_nldrill->damage[igauss] = drill_damage;
                }
            }

            // Add drilling stress = Bd'*Sd * dA (RHS)
            double Sd = m_drill_stiffness * m_drill_strain[igauss];
            if (m_drill_mode == DrillingDOF_NonLinear)
                Sd *= (1.0 - m_nldrill->damage_comm[igauss]);
            for (int i = 0; i < 24; i++)
                RHS(i) += Bd(i) * Sd * dA;

            // AGQI: update residual for internal DOFs
            if(m_eas)
                m_eas->Q_residual.addMatrixTransposeVector(1.0, BQ, S, -dA);
        }

        // AGQI: due to the static condensation, the following items
        // must be computed if we need either the RHS or the LHS, or both of them
        if (((options & OPT_RHS) && m_eas) || (options & OPT_LHS))
        {
            // Section tangent
            if (m_angle != 0.0) {
                Dsection = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
                if(m_damping[igauss]) Dsection *= m_damping[igauss]->getStiffnessMultiplier();
                auto& RsT = ASDShellQ4Globals::instance().RsT;
                RsT.addMatrixTranspose(0.0, Rs, 1.0);
                auto& DRsT = ASDShellQ4Globals::instance().DRsT;
                DRsT.addMatrixProduct(0.0, Dsection, RsT, 1.0);
                D.addMatrixProduct(0.0, Rs, DRsT, 1.0);
            }
            else {
                D = (options & OPT_LHS_IS_INITIAL) ?
                    m_sections[igauss]->getInitialTangent() :
                    m_sections[igauss]->getSectionTangent();
                if(m_damping[igauss]) D *= m_damping[igauss]->getStiffnessMultiplier();
            }

            // Matrices for AGQI static condensation of internal DOFs
            if (m_eas) {
                BQTD.addMatrixTransposeProduct(0.0, BQ, D, dA);
                DBQ.addMatrixProduct(0.0, D, BQ, dA);
                m_eas->KQQ_inv.addMatrixProduct(1.0, BQTD, BQ, 1.0);
                m_eas->KQU.addMatrixProduct(1.0, BQTD, B, 1.0);
                m_eas->KUQ.addMatrixTransposeProduct(1.0, B1, DBQ, 1.0);
            }
        }

        // Integrate LHS
        if (options & OPT_LHS) 
        {
            // Add current integration point contribution (LHS)
            B1TD.addMatrixTransposeProduct(0.0, B1, D, dA);
            LHS.addMatrixProduct(1.0, B1TD, B, 1.0);

            // Add drilling stiffness = Bd'*Kd*Bd * dA (LHS)
            double drill_tang = m_drill_stiffness;
            if (m_drill_mode == DrillingDOF_NonLinear) 
                drill_tang *= (1.0 - m_nldrill->damage_comm[igauss]);
            for (int i = 0; i < 24; i++)
                for (int j = 0; j < 24; j++)
                    LHS(i, j) += Bd(i) * drill_tang * Bd(j) * dA;
        }

    } // End of Gauss Loop

    // AGQI: static condensation
    if (((options & OPT_RHS) || (options & OPT_LHS)) && m_eas)
    {
        static Matrix KQQ = Matrix(4, 4);
        static Matrix KUQ_KQQ_inv = Matrix(24, 4);
        KQQ = m_eas->KQQ_inv;
        int inv_res = KQQ.Invert(m_eas->KQQ_inv);
        KUQ_KQQ_inv.addMatrixProduct(0.0, m_eas->KUQ, m_eas->KQQ_inv, 1.0);
        if (options & OPT_RHS)
            RHS.addMatrixVector(1.0, KUQ_KQQ_inv, m_eas->Q_residual, 1.0);
        if (options & OPT_LHS)
            LHS.addMatrixProduct(1.0, KUQ_KQQ_inv, m_eas->KQU, -1.0);
    }

    // Transform LHS to global coordinate system
    m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));

    // Subtract external loads if any
    if ((options & OPT_RHS) && m_load)
        RHS.addVector(1.0, *m_load, -1.0);

    // Done
    return result;
}

void ASDShellQ4::AGQIinitialize()
{
    // Global displacements
    auto& UG = ASDShellQ4Globals::instance().UG;
    m_transformation->computeGlobalDisplacements(UG);

    ASDShellQ4LocalCoordinateSystem local_cs = m_transformation->createLocalCoordinateSystem(UG);

    // Local displacements
    auto& UL = ASDShellQ4Globals::instance().UL;
    m_transformation->calculateLocalDisplacements(local_cs, UG, UL);

    // Initialize internal DOFs members
    m_eas->Q.Zero();
    m_eas->Q_converged.Zero();
    m_eas->U = UL;
    m_eas->U_converged = UL;
}

void ASDShellQ4::AGQIupdate(const Vector& UL)
{
    // Compute incremental displacements
    static Vector dUL(24);
    dUL = UL;
    dUL.addVector(1.0, m_eas->U, -1.0);

    // Save current trial displacements
    m_eas->U = UL;

    // Update internal DOFs
    static Vector temp(4);
    temp.addMatrixVector(0.0, m_eas->KQU, dUL, 1.0);
    temp.addVector(1.0, m_eas->Q_residual, -1.0);
    m_eas->Q.addMatrixVector(1.0, m_eas->KQQ_inv, temp, -1.0);
}

void ASDShellQ4::AGQIbeginGaussLoop(const ASDShellQ4LocalCoordinateSystem& reference_cs)
{
    // set to zero vectors and matrices for the static condensation
    // of internal DOFs before proceeding with gauss integration
    m_eas->KQU.Zero();
    m_eas->KUQ.Zero();
    m_eas->KQQ_inv.Zero();
    m_eas->Q_residual.Zero();

    // Some matrices/vectors
    auto& N = ASDShellQ4Globals::instance().N;
    auto& dN = ASDShellQ4Globals::instance().dN;

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // this one is already computed in the caller method
    auto& agq = ASDShellQ4Globals::instance().agq;

    // The AGQI uses incompatible modes and only the weak patch test is passed.
    // Here we computed the mean Bq matrix for the enhanced strains. It will be subtracted
    // from the gauss-wise Bq matrices to make this element pass the strict patch test:
    // BQ = BQ - 1/A*int{BQ*dV}
    auto& BQ_mean = ASDShellQ4Globals::instance().BQ_mean;
    BQ_mean.Zero();
    double Atot = 0.0;

    // Gauss loop
    std::array<double, 4> L;
    for (int igauss = 0; igauss < 4; igauss++)
    {
        // Current integration point data
        double xi = XI[igauss];
        double eta = ETA[igauss];
        double w = WTS[igauss];
        shapeFunctions(xi, eta, N);
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;
        Atot += dA;

        // area coordinates of the gauss point (Eq 7)
        L[0] = 0.25 * (1.0 - xi) * (agq.g[1] * (1.0 - eta) + agq.g[2] * (1.0 + eta));
        L[1] = 0.25 * (1.0 - eta) * (agq.g[3] * (1.0 - xi) + agq.g[2] * (1.0 + xi));
        L[2] = 0.25 * (1.0 + xi) * (agq.g[0] * (1.0 - eta) + agq.g[3] * (1.0 + eta));
        L[3] = 0.25 * (1.0 + eta) * (agq.g[0] * (1.0 - xi) + agq.g[1] * (1.0 + xi));

        // strain matrix for internal dofs
        for (int i = 0; i < 2; i++)
        { 
            int j = i + 1; if (j > 3) j = 0;
            int k = j + 1; if (k > 3) k = 0;
            double NQX = (agq.b[i] * L[k] + agq.b[k] * L[i]) / agq.A / 2.0;
            double NQY = (agq.c[i] * L[k] + agq.c[k] * L[i]) / agq.A / 2.0;
            int index1 = i * 2;
            int index2 = index1 + 1;
            BQ_mean(0, index1) += NQX * dA;
            BQ_mean(1, index2) += NQY * dA;
            BQ_mean(2, index1) += NQY * dA;
            BQ_mean(2, index2) += NQX * dA;
        }
    }

    // Average
    BQ_mean /= Atot;
}

int
ASDShellQ4::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** modes, int numMode)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
    nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
    nodePointers[2]->getDisplayCrds(v3, fact, displayMode);
    nodePointers[3]->getDisplayCrds(v4, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(4, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
        coords(3, i) = v4(i);
    }

    // set the quantity to be displayed at the nodes;
    static Vector values(4);
    for (int i = 0; i < 4; i++)
        values(i) = 0.0;

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}
