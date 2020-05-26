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
        opserr << "Want: element ASDShellQ4 $tag $iNode $jNoe $kNode $lNode $secTag <-corotational>";
        return 0;
    }

    int iData[6];
    int numData = 6;
    if (OPS_GetInt(&numData, iData) != 0) {
        opserr << "WARNING invalid integer tag: element ASDShellQ4 \n";
        return 0;
    }
    bool corotational = false;

    if (numArgs == 7) {
        const char* type = OPS_GetString();
        if ((strcmp(type, "-corotational") == 0) || (strcmp(type, "-Corotational") == 0))
            corotational = true;
    }

    SectionForceDeformation* section = OPS_getSectionForceDeformation(iData[5]);

    if (section == 0) {
        opserr << "ERROR:  element ASDShellQ4 " << iData[0] << "section " << iData[5] << " not found\n";
        return 0;
    }

    return new ASDShellQ4(iData[0], iData[1], iData[2], iData[3], iData[4], section, corotational);
}

// anonymous namespace for utilities
namespace
{

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

    // shape functions in iso-parametric space
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
        double Ax;
        double Ay;
        double Bx;
        double By;
        double Cx;
        double Cy;
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

    /** \brief GQ12Params
     *
     * This class performs some operations and stores some data to compute
     * the higher order shape functions for the drilling dof
     *
     */
    struct GQ12Params
    {
        double a1 = 0.0;
        double a2 = 0.0;
        double a3 = 0.0;
        double b1 = 0.0;
        double b2 = 0.0;
        double b3 = 0.0;

        void compute(const ASDShellQ4LocalCoordinateSystem& LCS)
        {
            static constexpr std::array<double, 4> s = { -1.0, 1.0, 1.0, -1.0 };
            static constexpr std::array<double, 4> t = { -1.0, -1.0, 1.0, 1.0 };

            for (size_t i = 0; i < 4; i++) {
                a1 += s[i] * LCS.X(i) / 4.0;
                a2 += t[i] * LCS.X(i) / 4.0;
                a3 += s[i] * t[i] * LCS.X(i) / 4.0;
                b1 += s[i] * LCS.Y(i) / 4.0;
                b2 += t[i] * LCS.Y(i) / 4.0;
                b3 += s[i] * t[i] * LCS.Y(i) / 4.0;
            }
        }
    };

    /** \brief ASShellQ4Globals
     *
     * This singleton class stores some data for the shell calculations that
     * can be statically instanciated to avoid useless re-allocations
     *
     */
    class ASDShellQ4Globals
    {
    private:
        ASDShellQ4Globals() {}

    public:

        JacobianOperator jac; // Jacobian
        MITC4Params mitc; // MITC4 parameters
        GQ12Params gq12; // GQ12 parameters

        Vector UG = Vector(24); // global displacements
        Vector UL = Vector(24); // local displacements

        Matrix B = Matrix(8, 24); // strain-displacement matrix
        Matrix dN = Matrix(4, 2); // shape functions derivatives in isoparametric space
        Vector E = Vector(8); // strain vector

        Matrix LHS = Matrix(24, 24); // LHS matrix
        Vector RHS = Vector(24); // RHS vector

    public:
        static ASDShellQ4Globals& instance() {
            static ASDShellQ4Globals _instance;
            return _instance;
        }
    };

    static 

    // computes the complete B matrix
    void computeBMatrix(
        double xi, double eta,
        const JacobianOperator& Jac, 
        const GQ12Params& gq12,
        const MITC4Params& mitc,
        const Matrix& dN,
        Matrix& B)
    {
        // cartesian derivatives of standard shape function
        static Matrix dNdX(4, 2);
        dNdX.addMatrixProduct(0.0, dN, Jac.invJ, 1.0);

        // cartesian derivatives of the modified shape functions for the 
        // drilling DOF according to the GQ12 formulation
        static constexpr std::array<double, 4> s = { -1.0, 1.0, 1.0, -1.0 };
        static constexpr std::array<double, 4> t = { -1.0, -1.0, 1.0, 1.0 };
        static Matrix dNU(4, 2);
        static Matrix dNV(4, 2);
        static Matrix dNUdX(4, 2);
        static Matrix dNVdX(4, 2);
        for (size_t i = 0; i < 4; i++) {
            dNU(i, 0) = 1.0 / 8.0 * (-2.0 * s[i] * xi * (gq12.b1 + gq12.b3 * t[i]) * (1.0 + t[i] * eta) +
                s[i] * t[i] * (1.0 - eta * eta) * (gq12.b2 + gq12.b3 * s[i]));
            dNU(i, 1) = 1.0 / 8.0 * (s[i] * t[i] * (1.0 - xi * xi) * (gq12.b1 + gq12.b3 * t[i]) -
                2.0 * t[i] * eta * (gq12.b2 + gq12.b3 * s[i]) * (1.0 + s[i] * xi));
            dNV(i, 0) = -1.0 / 8.0 * (-2.0 * s[i] * xi * (gq12.a1 + gq12.a3 * t[i]) * (1.0 + t[i] * eta) +
                s[i] * t[i] * (1.0 - eta * eta) * (gq12.a2 + gq12.a3 * s[i]));
            dNV(i, 1) = -1.0 / 8.0 * (s[i] * t[i] * (1.0 - xi * xi) * (gq12.a1 + gq12.a3 * t[i]) -
                2.0 * t[i] * eta * (gq12.a2 + gq12.a3 * s[i]) * (1.0 + s[i] * xi));
        }
        dNUdX.addMatrixProduct(0.0, dNU, Jac.invJ, 1.0);
        dNVdX.addMatrixProduct(0.0, dNV, Jac.invJ, 1.0);

        // membrane ************************************************************************************************

        B(0, 0) = dNdX(0, 0);   B(0, 6) = dNdX(1, 0);   B(0, 12) = dNdX(2, 0);   B(0, 18) = dNdX(3, 0);
        B(1, 1) = dNdX(0, 1);   B(1, 7) = dNdX(1, 1);   B(1, 13) = dNdX(2, 1);   B(1, 19) = dNdX(3, 1);
        B(2, 0) = dNdX(0, 1);   B(2, 6) = dNdX(1, 1);   B(2, 12) = dNdX(2, 1);   B(2, 18) = dNdX(3, 1);
        B(2, 1) = dNdX(0, 0);   B(2, 7) = dNdX(1, 0);   B(2, 13) = dNdX(2, 0);   B(2, 19) = dNdX(3, 0);

        // drilling ************************************************************************************************

        B(0, 2) = dNUdX(0, 0);
        B(1, 2) = dNVdX(0, 1);
        B(2, 2) = dNUdX(0, 1) + dNVdX(0, 0);

        B(0, 8) = dNUdX(1, 0);
        B(1, 8) = dNVdX(1, 1);
        B(2, 8) = dNUdX(1, 1) + dNVdX(1, 0);

        B(0, 14) = dNUdX(2, 0);
        B(1, 14) = dNVdX(2, 1);
        B(2, 14) = dNUdX(2, 1) + dNVdX(2, 0);

        B(0, 20) = dNUdX(3, 0);
        B(1, 20) = dNVdX(3, 1);
        B(2, 20) = dNUdX(3, 1) + dNVdX(3, 0);

        // bending *************************************************************************************************

        B(3, 4) =  dNdX(0, 0);   B(3, 10) =  dNdX(1, 0);   B(3, 16) =  dNdX(2, 0);   B(3, 22) =  dNdX(3, 0);
        B(4, 3) = -dNdX(0, 1);   B(4, 9)  = -dNdX(1, 1);   B(4, 15) = -dNdX(2, 1);   B(4, 21) = -dNdX(3, 1);
        B(5, 3) = -dNdX(0, 0);   B(5, 9)  = -dNdX(1, 0);   B(5, 15) = -dNdX(2, 0);   B(5, 21) = -dNdX(3, 0);
        B(5, 4) =  dNdX(0, 1);   B(5, 10) =  dNdX(1, 1);   B(5, 16) =  dNdX(2, 1);   B(5, 22) =  dNdX(3, 1);

        // shear ***************************************************************************************************

        // MITC modified shape functions
        static Matrix MITCShapeFunctions(2, 4);
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
        // to local coordinate system taking into account the element distorsion
        static Matrix TBN(2, 24);
        TBN.addMatrixProduct(0.0, mitc.transformation, BN, 1.0);
        for (int i = 0; i < 2; i++)
            for (int j = 0; j < 24; j++)
                B(i + 6, j) = TBN(i, j);
    }

}

ASDShellQ4::ASDShellQ4() 
    : Element(0, ELE_TAG_ASDShellQ4)
    , m_transformation(new ASDShellQ4Transformation())
{
}

ASDShellQ4::ASDShellQ4(
    int tag,
    int node1,
    int node2,
    int node3,
    int node4,
    SectionForceDeformation* section,
    bool corotational)
    : Element(tag, ELE_TAG_ASDShellQ4)
    , m_transformation(corotational ? new ADShellQ4CorotationalTransformation() : new ASDShellQ4Transformation())
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
}

ASDShellQ4::~ASDShellQ4( )
{
    // clean up section
    for (int i = 0; i < 4; i++)
        delete m_sections[i];

    // clean up coordinate transformation
    delete m_transformation;

    // clean up load vectors
    if (m_load)
        delete m_load;
}

void  ASDShellQ4::setDomain(Domain* theDomain)
{
    // set domain on transformation
    m_transformation->setDomain(theDomain, m_node_ids);

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

    return success;
}

int ASDShellQ4::update()
{
    // Output code
    int result = 0;

    // Compute the local coordinate system.
    ASDShellQ4LocalCoordinateSystem local_cs = 
        m_transformation->createLocalCoordinateSystem();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs = 
        m_transformation->createReferenceCoordinateSystem();

    // Prepare all the parameters needed for the MITC4
    // and GQ12 formulations.
    // This is to be done here outside the Gauss Loop.
    auto& mitc = ASDShellQ4Globals::instance().mitc;
    auto& gq12 = ASDShellQ4Globals::instance().gq12;
    mitc.compute(reference_cs);
    gq12.compute(reference_cs);

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& B = ASDShellQ4Globals::instance().B;
    auto& dN = ASDShellQ4Globals::instance().dN;
    auto& E = ASDShellQ4Globals::instance().E;

    // Displacements
    auto& UG = ASDShellQ4Globals::instance().UG;
    auto& UL = ASDShellQ4Globals::instance().UL;
    m_transformation->computeGlobalDisplacements(UG);
    m_transformation->calculateLocalDisplacements(local_cs, UG, UL);

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);

        // Strain-displacement matrix
        computeBMatrix(xi, eta, jac, gq12, mitc, dN, B);

        // Section deformation
        E.addMatrixVector(0.0, B, UL, 1.0);
        
        // Update section
        result += m_sections[i]->setTrialSectionDeformation(E);
    }

    // Done
    return result;
}

const Matrix& ASDShellQ4::getTangentStiff()
{
    // Output matrix
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    LHS.Zero();

    // Compute the local coordinate system.
    ASDShellQ4LocalCoordinateSystem local_cs =
        m_transformation->createLocalCoordinateSystem();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Prepare all the parameters needed for the MITC4
    // and GQ12 formulations.
    // This is to be done here outside the Gauss Loop.
    auto& mitc = ASDShellQ4Globals::instance().mitc;
    auto& gq12 = ASDShellQ4Globals::instance().gq12;
    mitc.compute(reference_cs);
    gq12.compute(reference_cs);

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& B = ASDShellQ4Globals::instance().B;
    auto& dN = ASDShellQ4Globals::instance().dN;

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Strain-displacement matrix
        computeBMatrix(xi, eta, jac, gq12, mitc, dN, B);

        // Section tangent
        const auto& D = m_sections[i]->getSectionTangent();

        // Add current integration point contribution
        LHS.addMatrixTripleProduct(1.0, B, D, dA);
    }

    // Tranform LHS to global coordinate system
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    m_transformation->transformToGlobal(local_cs, LHS, RHS, false, true);

    // Done
    return LHS;
}

const Matrix& ASDShellQ4::getInitialStiff()
{
    // Output matrix
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    LHS.Zero();

    // Compute the local coordinate system.
    ASDShellQ4LocalCoordinateSystem local_cs =
        m_transformation->createLocalCoordinateSystem();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Prepare all the parameters needed for the MITC4
    // and GQ12 formulations.
    // This is to be done here outside the Gauss Loop.
    auto& mitc = ASDShellQ4Globals::instance().mitc;
    auto& gq12 = ASDShellQ4Globals::instance().gq12;
    mitc.compute(reference_cs);
    gq12.compute(reference_cs);

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& B = ASDShellQ4Globals::instance().B;
    auto& dN = ASDShellQ4Globals::instance().dN;

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Strain-displacement matrix
        computeBMatrix(xi, eta, jac, gq12, mitc, dN, B);

        // Section tangent
        const auto& D = m_sections[i]->getInitialTangent();

        // Add current integration point contribution
        LHS.addMatrixTripleProduct(1.0, B, D, dA);
    }

    // Tranform LHS to global coordinate system
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    m_transformation->transformToGlobal(local_cs, LHS, RHS, false, true);

    // Done
    return LHS;
}

const Matrix& ASDShellQ4::getMass()
{
    // Output matrix
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    LHS.Zero();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& dN = ASDShellQ4Globals::instance().dN;

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Section rho
        double rho = m_sections[i]->getRho();

        // Add current integration point contribution
        for (int j = 0; j < 4; j++)
        {
            int index = j * 6;

            // Translational mass contribution
            for (int q = 0; q < 3; q++)
                LHS(index + q, index + q) += rho * dA;

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
    // Output vector
    auto& RHS = ASDShellQ4Globals::instance().RHS;
    RHS.Zero();

    // Compute the local coordinate system.
    ASDShellQ4LocalCoordinateSystem local_cs =
        m_transformation->createLocalCoordinateSystem();

    // Compute the reference coordinate system
    ASDShellQ4LocalCoordinateSystem reference_cs =
        m_transformation->createReferenceCoordinateSystem();

    // Prepare all the parameters needed for the MITC4
    // and GQ12 formulations.
    // This is to be done here outside the Gauss Loop.
    auto& mitc = ASDShellQ4Globals::instance().mitc;
    auto& gq12 = ASDShellQ4Globals::instance().gq12;
    mitc.compute(reference_cs);
    gq12.compute(reference_cs);

    // Jacobian
    auto& jac = ASDShellQ4Globals::instance().jac;

    // Some matrices
    auto& B = ASDShellQ4Globals::instance().B;
    auto& dN = ASDShellQ4Globals::instance().dN;

    // Gauss loop
    for (int i = 0; i < 4; i++)
    {
        // Current integration point data
        double xi = XI[i];
        double eta = ETA[i];
        double w = WTS[i];
        shapeFunctionsNaturalDerivatives(xi, eta, dN);
        jac.calculate(reference_cs, dN);
        double dA = w * jac.detJ;

        // Strain-displacement matrix
        computeBMatrix(xi, eta, jac, gq12, mitc, dN, B);

        // Section force
        const auto& S = m_sections[i]->getStressResultant();

        // Add current integration point contribution
        RHS.addMatrixTransposeVector(1.0, B, S, dA);
    }

    // Tranform RHS to global coordinate system
    auto& LHS = ASDShellQ4Globals::instance().LHS;
    m_transformation->transformToGlobal(local_cs, LHS, RHS, true, false);

    // Subtract external loads if any
    if (m_load) {
        RHS.addVector(1.0, *m_load, -1.0);
    }

    // Done
    return RHS;
}

const Vector& ASDShellQ4::getResistingForceIncInertia()
{
    // Output vector
    auto& RHS = ASDShellQ4Globals::instance().RHS;

    // Compute residual (sets RHS above...)
    getResistingForce();

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
}

int  ASDShellQ4::sendSelf(int commitTag, Channel& theChannel)
{
    return -1;
    //int res = 0;

    //// note: we don't check for dataTag == 0 for Element
    //// objects as that is taken care of in a commit by the Domain
    //// object - don't want to have to do the check if sending data
    //int dataTag = this->getDbTag();

    //// send the ids of sections
    //int matDbTag;

    //static ID idData(14);

    //for (int i = 0; i < 4; i++) {
    //    idData(i) = m_sections[i]->getClassTag();
    //    matDbTag = m_sections[i]->getDbTag();
    //    // NOTE: we do have to ensure that the material has a database
    //    // tag if we are sending to a database channel.
    //    if (matDbTag == 0) {
    //        matDbTag = theChannel.getDbTag();
    //        if (matDbTag != 0)
    //            m_sections[i]->setDbTag(matDbTag);
    //    }
    //    idData(i + 4) = matDbTag;
    //}

    //idData(8) = getTag();
    //idData(9) = m_node_ids(0);
    //idData(10) = m_node_ids(1);
    //idData(11) = m_node_ids(2);
    //idData(12) = m_node_ids(3);
    //idData(13) = m_transformation->isLinear() ? 0 : 1;

    //res += theChannel.sendID(dataTag, commitTag, idData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    //    return res;
    //}

    //static Vector vectData(5 + 6 * 4);
    //vectData(0) = Ktt;
    //vectData(1) = alphaM;
    //vectData(2) = betaK;
    //vectData(3) = betaK0;
    //vectData(4) = betaKc;

    //int pos = 0;
    //for (int node = 0; node < 4; ++node)
    //{
    //    for (int dof = 0; dof < 6; ++dof)
    //    {
    //        vectData(5 + pos) = init_disp[node][dof];
    //        pos++;
    //    }
    //}

    //res += theChannel.sendVector(dataTag, commitTag, vectData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    //    return res;
    //}

    //// Finally, quad asks its material objects to send themselves
    //for (i = 0; i < 4; i++) {
    //    res += m_sections[i]->sendSelf(commitTag, theChannel);
    //    if (res < 0) {
    //        opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    //        return res;
    //    }
    //}

    //return res;
}

int  ASDShellQ4::recvSelf(int commitTag,
    Channel& theChannel,
    FEM_ObjectBroker& theBroker)
{
    return -1;
    //int res = 0;

    //int dataTag = this->getDbTag();

    //static ID idData(14);
    //// Quad now receives the tags of its four external nodes
    //res += theChannel.recvID(dataTag, commitTag, idData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellQ4::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    //    return res;
    //}

    //this->setTag(idData(8));
    //connectedExternalNodes(0) = idData(9);
    //connectedExternalNodes(1) = idData(10);
    //connectedExternalNodes(2) = idData(11);
    //connectedExternalNodes(3) = idData(12);
    //if (idData(13) == 0)
    //    doUpdateBasis = true;
    //else
    //    doUpdateBasis = false;

    //static Vector vectData(5 + 6 * 4);
    //res += theChannel.recvVector(dataTag, commitTag, vectData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellQ4::sendSelf() - " << this->getTag() << " failed to send ID\n";
    //    return res;
    //}

    //Ktt = vectData(0);
    //alphaM = vectData(1);
    //betaK = vectData(2);
    //betaK0 = vectData(3);
    //betaKc = vectData(4);


    //int pos = 0;
    //for (int node = 0; node < 4; ++node)
    //{
    //    for (int dof = 0; dof < 6; ++dof)
    //    {
    //        init_disp[node][dof] = vectData(5 + pos);
    //        pos++;
    //    }
    //}


    //int i;

    //if (materialPointers[0] == 0) {
    //    for (i = 0; i < 4; i++) {
    //        int matClassTag = idData(i);
    //        int matDbTag = idData(i + 4);
    //        // Allocate new material with the sent class tag
    //        materialPointers[i] = theBroker.getNewSection(matClassTag);
    //        if (materialPointers[i] == 0) {
    //            opserr << "ASDShellQ4::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
    //            return -1;
    //        }
    //        // Now receive materials into the newly allocated space
    //        materialPointers[i]->setDbTag(matDbTag);
    //        res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
    //        if (res < 0) {
    //            opserr << "ASDShellQ4::recvSelf() - material " << i << "failed to recv itself\n";
    //            return res;
    //        }
    //    }
    //}
    //// Number of materials is the same, receive materials into current space
    //else {
    //    for (i = 0; i < 4; i++) {
    //        int matClassTag = idData(i);
    //        int matDbTag = idData(i + 4);
    //        // Check that material is of the right type; if not,
    //        // delete it and create a new one of the right type
    //        if (materialPointers[i]->getClassTag() != matClassTag) {
    //            delete materialPointers[i];
    //            materialPointers[i] = theBroker.getNewSection(matClassTag);
    //            if (materialPointers[i] == 0) {
    //                opserr << "ASDShellQ4::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;
    //                exit(-1);
    //            }
    //        }
    //        // Receive the material
    //        materialPointers[i]->setDbTag(matDbTag);
    //        res += materialPointers[i]->recvSelf(commitTag, theChannel, theBroker);
    //        if (res < 0) {
    //            opserr << "ASDShellQ4::recvSelf() - material " << i << "failed to recv itself\n";
    //            return res;
    //        }
    //    }
    //}

    //return res;
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

    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0) {
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
            output.tag("ResponseType", "p1212");
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
    case 3: //strain
        for (int i = 0; i < 4; i++) {

            // Get section deformation
            const Vector& deformation = m_sections[i]->getSectionDeformation();
            strains(0) = deformation(0);
            strains(1) = deformation(1);
            strains(2) = deformation(2);
            strains(3) = deformation(3);
            strains(4) = deformation(4);
            strains(5) = deformation(5);
            strains(6) = deformation(6);
            strains(7) = deformation(7);
        }
        return eleInfo.setVector(strains);
        break;
    default:
        return -1;
    }
}









