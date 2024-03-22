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
// the Cell-based Smoothed Discrete Shear Gap (CS-DSG) formulation 
// to avoid transverse-shear-locking in the thin shell limit.
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
            "<-corotational> <-drillingNL> <-damp $dampTag>"
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

    return new ASDShellT3(iData[0], iData[1], iData[2], iData[3], section, local_x, corotational, drill_mode, m_damping);
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
    constexpr double XI = 1.0 / 3.0;
    constexpr double ETA = 1.0 / 3.0;
    constexpr double WTS = 0.5;

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

    // computes the complete B matrix
    void computeBMatrix(
        const Vector3Type& p1,
        const Vector3Type& p2,
        const Vector3Type& p3,
        double A,
        Matrix& B,
        Vector& Bd
    )
    {
        // some constants
        double a = p2.x() - p1.x();
        double b = p2.y() - p1.y();
        double c = p3.y() - p1.y();
        double d = p3.x() - p1.x();
        double w = 1.0 / (2.0 * A);

        // initialize
        B.Zero();
        Bd.Zero();

        // membrane part (compatible part)
        // N1
        B(0, 0) = (b - c) * w;
        B(1, 1) = (d - a) * w;
        B(2, 0) = (d - a) * w;
        B(2, 1) = (b - c) * w;
        // N2
        B(0, 6) =  c * w;
        B(1, 7) = -d * w;
        B(2, 6) = -d * w;
        B(2, 7) =  c * w;
        // N3
        B(0, 12) = -b * w;
        B(1, 13) =  a * w;
        B(2, 12) =  a * w;
        B(2, 13) = -b * w;

        // bending part
        // N1
        B(3, 4) = (c - b) * w;
        B(4, 3) = (d - a) * w;
        B(5, 3) = (b - c) * w;
        B(5, 4) = (a - d) * w;
        // N2
        B(3, 10) = -c * w;
        B(4, 9) = -d * w;
        B(5, 9) = c * w;
        B(5, 10) =  d * w;
        // N3
        B(3, 16) = b * w;
        B(4, 15) = a * w;
        B(5, 15) = -b * w;
        B(5, 16) = -a * w;

        // shear part (DSG treatment of transverse shear locking)
        // N1
        B(6, 2) = (b - c) * w;
        B(6, 4) = 0.5;
        B(7, 2) = (d - a) * w;
        B(7, 3) = -0.5;
        // N2
        B(6, 8) = c * w;
        B(6, 9) = -(b * c) / 2.0 * w;
        B(6, 10) = (a * c) / 2.0 * w;
        B(7, 8) = -d * w;
        B(7, 9) = (b * d) / 2.0 * w;
        B(7, 10) = -(a * d) / 2.0 * w;
        // N3
        B(6, 14) = -b * w;
        B(6, 15) = (b * c) / 2.0 * w;
        B(6, 16) = -(b * d) / 2.0 * w;
        B(7, 14) = a * w;
        B(7, 15) = -(a * c) / 2.0 * w;
        B(7, 16) = (a * d) / 2.0 * w;

        // We use the drilling penalty as defined by hughes and brezzi,
        // where we link the independent rotation to the skew symmetric part of the in-plane displacement field.

        //Bd(0) = -0.5 * dNdX(0, 1);
        //Bd(1) = 0.5 * dNdX(0, 0);
        //Bd(5) = -N(0);

        //Bd(6) = -0.5 * dNdX(1, 1);
        //Bd(7) = 0.5 * dNdX(1, 0);
        //Bd(11) = -N(1);

        //Bd(12) = -0.5 * dNdX(2, 1);
        //Bd(13) = 0.5 * dNdX(2, 0);
        //Bd(17) = -N(2);

        //Bd(18) = -0.5 * dNdX(3, 1);
        //Bd(19) = 0.5 * dNdX(3, 0);
        //Bd(23) = -N(3);
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
    DrillingDOFMode drill_mode,
    Damping* theDamping)
    : Element(tag, ELE_TAG_ASDShellT3)
    , m_transformation(corotational ? new ASDShellT3CorotationalTransformation() : new ASDShellT3Transformation())
    , m_drill_mode(drill_mode)
{
    // save node ids
    m_node_ids(0) = node1;
    m_node_ids(1) = node2;
    m_node_ids(2) = node3;

    // copy section
    m_section = section->getCopy();
    if (m_section == nullptr) {
        opserr << "ASDShellT3::constructor - failed to get a material of type: ShellSection\n";
        exit(-1);
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
        m_damping = theDamping->getCopy();
        if (m_damping == nullptr) {
            opserr << "ASDShellT3::constructor - failed to get copy of damping\n";
        }
    }
}

ASDShellT3::~ASDShellT3()
{
    // clean up section
    if(m_section)
        delete m_section;

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
    if (m_damping)
        delete m_damping;
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
        m_drill_stiffness = m_section->getInitialTangent()(2, 2);
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

        if (m_damping && m_damping->setDomain(theDomain, 8)) {
            opserr << "ASDShellT3::setDomain -- Error initializing damping\n";
            exit(-1);
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
        const Vector& stress = m_section->getStressResultant();
        s << "STRESS\t" << eleTag << "\t" << counter << "\t" << 0 << "\tTOP";
        for (int j = 0; j < 6; j++)
            s << "\t" << stress(j);
        s << endln;
    }

    if (flag == OPS_PRINT_CURRENTSTATE) {
        s << endln;
        s << "ASDShellT3 Non-Locking Three Node Shell \n";
        s << "Element Number: " << this->getTag() << endln;
        s << "Node 1 : " << m_node_ids(0) << endln;
        s << "Node 2 : " << m_node_ids(1) << endln;
        s << "Node 3 : " << m_node_ids(2) << endln;

        s << "Material Information : \n ";
        m_section->Print(s, flag);

        s << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": " << this->getTag() << ", ";
        s << "\"type\": \"ASDShellT3\", ";
        s << "\"nodes\": [" << m_node_ids(0) << ", " << m_node_ids(1) << ", ";
        s << m_node_ids(2) << "], ";
        s << "\"section\": \"" << m_section->getTag() << "\"}";
    }
}

int
ASDShellT3::setDamping(Domain* theDomain, Damping* damping)
{
    if (theDomain && damping)
    {
        if (m_damping)
            delete m_damping;
        m_damping = damping->getCopy();
        if (!m_damping) {
            opserr << "ASDShellT3::setDamping - failed to get copy of damping\n";
            return -1;
        }
        if (m_damping && m_damping->setDomain(theDomain, 8)) {
            opserr << "ASDShellT3::setDamping -- Error initializing damping\n";
            return -2;
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
    success += m_section->commitState();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm = m_section->getStressResultant();
        m_nldrill->strain_comm = m_section->getSectionDeformation();
        m_nldrill->damage_comm = m_nldrill->damage;
    }

    // damping
    if (m_damping) 
        success += m_damping->commitState();

    // done
    return success;
}

int ASDShellT3::revertToLastCommit()
{
    int success = 0;

    // transformation
    m_transformation->revertToLastCommit();

    // section
    success += m_section->revertToLastCommit();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm = m_section->getStressResultant();
        m_nldrill->strain_comm = m_section->getSectionDeformation();
        m_nldrill->damage = m_nldrill->damage_comm;
    }

    // damping
    if (m_damping) 
        success += m_damping->revertToLastCommit();

    // done
    return success;
}

int  ASDShellT3::revertToStart()
{
    int success = 0;

    // transformation
    m_transformation->revertToStart();

    // section
    success += m_section->revertToStart();
    if (m_drill_mode == DrillingDOF_NonLinear) {
        m_nldrill->stress_comm.Zero();
        m_nldrill->strain_comm.Zero();
        m_nldrill->damage = m_nldrill->damage_comm = 0.0;
    }


    // damping
    if (m_damping) 
        success += m_damping->revertToStart();

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
    double rho = m_section->getRho();

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

    //// note: we don't check for dataTag == 0 for Element
    //// objects as that is taken care of in a commit by the Domain
    //// object - don't want to have to do the check if sending data
    //int dataTag = this->getDbTag();

    //// a counter
    //int counter;

    //// has load flag
    //bool has_load = m_load != nullptr;

    //// INT data
    //// 1 tag +
    //// 4 node tags +
    //// 1 transformation type
    //// 1 initialization flag
    //// 1 has_load_flag
    //// 8 -> 4 pairs of (section class tag, mt db tag)
    //// 2 -> damping tag + damping db tag
    //// 1 -> EAS flag
    //// 1 -> non-linear drilling flag
    //// 1 -> local_x flag
    //static ID idData(21);
    //counter = 0;
    //idData(counter++) = getTag();
    //for (int i = 0; i < 4; ++i)
    //    idData(counter++) = m_node_ids(i);
    //idData(counter++) = static_cast<int>(m_transformation->isLinear());
    //idData(counter++) = static_cast<int>(m_initialized);
    //idData(counter++) = static_cast<int>(has_load);
    //for (int i = 0; i < 4; i++) {
    //    idData(counter++) = m_sections[i]->getClassTag();
    //    int matDbTag = m_sections[i]->getDbTag();
    //    // NOTE: we do have to ensure that the material has a database
    //    // tag if we are sending to a database channel.
    //    if (matDbTag == 0) {
    //        matDbTag = theChannel.getDbTag();
    //        if (matDbTag != 0)
    //            m_sections[i]->setDbTag(matDbTag);
    //    }
    //    idData(counter++) = matDbTag;
    //}
    //if (m_damping[0]) {
    //    idData(counter++) = m_damping[0]->getClassTag();
    //    int dbTag = m_damping[0]->getDbTag();
    //    if (dbTag == 0) {
    //        dbTag = theChannel.getDbTag();
    //        if (dbTag != 0)
    //            for (int i = 0; i < 4; i++)
    //                m_damping[i]->setDbTag(dbTag);
    //    }
    //    idData(counter++) = dbTag;
    //}
    //else {
    //    idData(counter++) = 0;
    //    idData(counter++) = 0;
    //}
    //idData(counter++) = static_cast<int>(static_cast<bool>(m_eas));
    //idData(counter++) = static_cast<int>(static_cast<bool>(m_nldrill));
    //idData(counter++) = static_cast<int>(static_cast<bool>(m_local_x));

    //res = theChannel.sendID(dataTag, commitTag, idData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send ID\n";
    //    return res;
    //}

    //// DOUBLE data
    //// 4 damping factors +
    //// 4 drilling strain +
    //// 1 drilling stiffness +
    //// 1 drilling stabilization parameter +
    //// 1 angle +
    //// (optional) 3 local_x +
    //// (optional) 268 AGQI internal data + 
    //// (optional) 72 non-linear drilling data +
    //// (optional) 24 load +
    //// NT transformation data
    //int Nlocalx = m_local_x ? 3 : 0;
    //int NLoad = has_load ? 24 : 0;
    //int Neas = m_eas ? 268 : 0;
    //int Nnldrill = m_nldrill ? 72 : 0;
    //int NT = m_transformation->internalDataSize();
    //Vector vectData(4 + 4 + 1 + 1 + 1 + Nlocalx + Neas + Nnldrill + NLoad + NT);
    //counter = 0;
    //vectData(counter++) = alphaM;
    //vectData(counter++) = betaK;
    //vectData(counter++) = betaK0;
    //vectData(counter++) = betaKc;
    //for (int i = 0; i < 4; ++i)
    //    vectData(counter++) = m_drill_strain[i];
    //vectData(counter++) = m_drill_stiffness;
    //vectData(counter++) = m_drill_stab;
    //vectData(counter++) = m_angle;
    //if (m_local_x) {
    //    for (int i = 0; i < 3; ++i)
    //        vectData(counter++) = (*m_local_x)(i);
    //}
    //if (m_eas) {
    //    for (int i = 0; i < 4; ++i)
    //        vectData(counter++) = m_eas->Q(i);
    //    for (int i = 0; i < 4; ++i)
    //        vectData(counter++) = m_eas->Q_converged(i);
    //    for (int i = 0; i < 24; ++i)
    //        vectData(counter++) = m_eas->U(i);
    //    for (int i = 0; i < 24; ++i)
    //        vectData(counter++) = m_eas->U_converged(i);
    //    for (int i = 0; i < 4; ++i)
    //        vectData(counter++) = m_eas->Q_residual(i);
    //    for (int i = 0; i < 4; ++i)
    //        for (int j = 0; j < 4; ++j)
    //            vectData(counter++) = m_eas->KQQ_inv(i, j);
    //    for (int i = 0; i < 4; ++i)
    //        for (int j = 0; j < 24; ++j)
    //            vectData(counter++) = m_eas->KQU(i, j);
    //    for (int i = 0; i < 24; ++i)
    //        for (int j = 0; j < 4; ++j)
    //            vectData(counter++) = m_eas->KUQ(i, j);
    //}
    //if (m_nldrill) {
    //    for (const auto& item : m_nldrill->strain_comm)
    //        for (int i = 0; i < 8; ++i)
    //            vectData(counter++) = item(i);
    //    for (const auto& item : m_nldrill->stress_comm)
    //        for (int i = 0; i < 8; ++i)
    //            vectData(counter++) = item(i);
    //    for (auto item : m_nldrill->damage)
    //        vectData(counter++) = item;
    //    for (auto item : m_nldrill->damage_comm)
    //        vectData(counter++) = item;
    //}
    //if (has_load) {
    //    for (int i = 0; i < 24; ++i)
    //        vectData(counter++) = (*m_load)(i);
    //}
    //m_transformation->saveInternalData(vectData, counter);

    //res = theChannel.sendVector(dataTag, commitTag, vectData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    //    return res;
    //}

    //// send all sections
    //for (int i = 0; i < 4; i++) {
    //    res = m_sections[i]->sendSelf(commitTag, theChannel);
    //    if (res < 0) {
    //        opserr << "WARNING ASDShellT3::sendSelf() - " << this->getTag() << " failed to send its Material\n";
    //        return res;
    //    }
    //}

    //// Ask the Damping to send itself
    //if (m_damping[0]) {
    //    for (int i = 0; i < 4; i++) {
    //        res += m_damping[i]->sendSelf(commitTag, theChannel);
    //        if (res < 0) {
    //            opserr << "ASDShellT3::sendSelf -- could not send Damping\n";
    //            return res;
    //        }
    //    }
    //}

    // done
    return res;
}

int  ASDShellT3::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
    int res = 0;

    //int dataTag = this->getDbTag();

    //// a counter
    //int counter;

    //// INT data
    //// 1 tag +
    //// 4 node tags +
    //// 1 transformation type
    //// 1 initialization flag
    //// 1 has_load_flag
    //// 8 -> 4 pairs of (section class tag, mt db tag)
    //// 2 -> damping tag + damping db tag
    //// 1 -> EAS flag
    //// 1 -> non-linear drilling flag
    //// 1 -> local_x flag
    //static ID idData(21);
    //res = theChannel.recvID(dataTag, commitTag, idData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellT3::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    //    return res;
    //}

    //counter = 0;
    //setTag(idData(counter++));
    //for (int i = 0; i < 4; ++i)
    //    m_node_ids(i) = idData(counter++);
    //bool linear_transform = static_cast<bool>(idData(counter++));
    //m_initialized = static_cast<bool>(idData(counter++));
    //bool has_load = static_cast<bool>(idData(counter++));
    //// create sections
    //for (int i = 0; i < 4; i++) {
    //    int matClassTag = idData(counter++);
    //    int matDbTag = idData(counter++);
    //    if (m_sections[i])
    //        delete m_sections[i];
    //    m_sections[i] = theBroker.getNewSection(matClassTag);
    //    if (m_sections[i] == 0) {
    //        opserr << "ASDShellT3::recvSelf() - Broker could not create NDMaterial of class type" << matClassTag << endln;;
    //        return -1;
    //    }
    //    m_sections[i]->setDbTag(matDbTag);
    //}
    //int dmpTag = idData(counter++);
    //int dmpDbTag = idData(counter++);
    //bool use_eas = static_cast<bool>(idData(counter++));
    //bool use_nldrill = static_cast<bool>(idData(counter++));
    //bool use_local_x = static_cast<bool>(idData(counter++));

    //// create transformation
    //if (m_transformation)
    //    delete m_transformation;
    //m_transformation = linear_transform ? new ASDShellT3Transformation() : new ASDShellT3CorotationalTransformation();
    //// create load
    //if (has_load) {
    //    if (m_load == nullptr)
    //        m_load = new Vector(24);
    //}
    //else {
    //    if (m_load) {
    //        delete m_load;
    //        m_load = nullptr;
    //    }
    //}

    //// create eas
    //if (m_eas) {
    //    delete m_eas;
    //    m_eas = nullptr;
    //}
    //if (use_eas)
    //    m_eas = new EASData();

    //// create non-linear drilling
    //if (m_nldrill) {
    //    delete m_nldrill;
    //    m_nldrill = nullptr;
    //}
    //if (use_nldrill)
    //    m_nldrill = new NLDrillingData();

    //// clear local x
    //if (m_local_x) {
    //    delete m_local_x;
    //    m_local_x = nullptr;
    //}
    //if (use_local_x)
    //    m_local_x = new Vector(3);

    //// DOUBLE data
    //// 4 damping factors +
    //// 4 drilling strain +
    //// 1 drilling stiffness +
    //// 1 drilling stabilization parameter +
    //// 1 angle +
    //// (optional) 3 local_x +
    //// (optional) 268 AGQI internal data + 
    //// (optional) 72 non-linear drilling data +
    //// (optional) 24 load +
    //// NT transformation data
    //int Nlocalx = use_local_x ? 3 : 0;
    //int NLoad = has_load ? 24 : 0;
    //int Neas = use_eas ? 268 : 0;
    //int Nnldrill = use_nldrill ? 72 : 0;
    //int NT = m_transformation->internalDataSize();
    //Vector vectData(4 + 4 + 1 + 1 + 1 + Neas + Nnldrill + NLoad + NT);
    //res = theChannel.recvVector(dataTag, commitTag, vectData);
    //if (res < 0) {
    //    opserr << "WARNING ASDShellT3::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
    //    return res;
    //}

    //counter = 0;
    //alphaM = vectData(counter++);
    //betaK = vectData(counter++);
    //betaK0 = vectData(counter++);
    //betaKc = vectData(counter++);
    //for (int i = 0; i < 4; ++i)
    //    m_drill_strain[i] = vectData(counter++);
    //m_drill_stiffness = vectData(counter++);
    //m_drill_stab = vectData(counter++);
    //m_angle = vectData(counter++);
    //if (use_local_x) {
    //    for (int i = 0; i < 3; ++i)
    //        (*m_local_x)(i) = vectData(counter++);
    //}
    //if (use_eas) {
    //    for (int i = 0; i < 4; ++i)
    //        m_eas->Q(i) = vectData(counter++);
    //    for (int i = 0; i < 4; ++i)
    //        m_eas->Q_converged(i) = vectData(counter++);
    //    for (int i = 0; i < 24; ++i)
    //        m_eas->U(i) = vectData(counter++);
    //    for (int i = 0; i < 24; ++i)
    //        m_eas->U_converged(i) = vectData(counter++);
    //    for (int i = 0; i < 4; ++i)
    //        m_eas->Q_residual(i) = vectData(counter++);
    //    for (int i = 0; i < 4; ++i)
    //        for (int j = 0; j < 4; ++j)
    //            m_eas->KQQ_inv(i, j) = vectData(counter++);
    //    for (int i = 0; i < 4; ++i)
    //        for (int j = 0; j < 24; ++j)
    //            m_eas->KQU(i, j) = vectData(counter++);
    //    for (int i = 0; i < 24; ++i)
    //        for (int j = 0; j < 4; ++j)
    //            m_eas->KUQ(i, j) = vectData(counter++);
    //}
    //if (use_nldrill) {
    //    for (auto& item : m_nldrill->strain_comm)
    //        for (int i = 0; i < 8; ++i)
    //            item(i) = vectData(counter++);
    //    for (auto& item : m_nldrill->stress_comm)
    //        for (int i = 0; i < 8; ++i)
    //            item(i) = vectData(counter++);
    //    for (auto& item : m_nldrill->damage)
    //        item = vectData(counter++);
    //    for (auto& item : m_nldrill->damage_comm)
    //        item = vectData(counter++);
    //}
    //if (has_load) {
    //    for (int i = 0; i < 24; ++i)
    //        (*m_load)(i) = vectData(counter++);
    //}
    //m_transformation->restoreInternalData(vectData, counter);

    //// all sections
    //for (int i = 0; i < 4; i++) {
    //    res = m_sections[i]->recvSelf(commitTag, theChannel, theBroker);
    //    if (res < 0) {
    //        opserr << "ASDShellT3::recvSelf() - material " << i << "failed to recv itself\n";
    //        return res;
    //    }
    //}

    //if (dmpTag) {
    //    for (int i = 0; i < 4; i++) {
    //        // Check if the Damping is null; if so, get a new one
    //        if (m_damping[i] == 0) {
    //            m_damping[i] = theBroker.getNewDamping(dmpTag);
    //            if (m_damping[i] == 0) {
    //                opserr << "ASDShellT3::recvSelf -- could not get a Damping\n";
    //                exit(-1);
    //            }
    //        }

    //        // Check that the Damping is of the right type; if not, delete
    //        // the current one and get a new one of the right type
    //        if (m_damping[i]->getClassTag() != dmpTag) {
    //            delete m_damping[i];
    //            m_damping[i] = theBroker.getNewDamping(dmpTag);
    //            if (m_damping[i] == 0) {
    //                opserr << "ASDShellT3::recvSelf -- could not get a Damping\n";
    //                exit(-1);
    //            }
    //        }

    //        // Now, receive the Damping
    //        m_damping[i]->setDbTag(dmpDbTag);
    //        res += m_damping[i]->recvSelf(commitTag, theChannel, theBroker);
    //        if (res < 0) {
    //            opserr << "ASDShellT3::recvSelf -- could not receive Damping\n";
    //            return res;
    //        }
    //    }
    //}
    //else {
    //    for (int i = 0; i < 4; i++)
    //    {
    //        if (m_damping[i])
    //        {
    //            delete m_damping[i];
    //            m_damping[i] = 0;
    //        }
    //    }
    //}

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

    else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0) {
        if (argc < 2) {
            opserr << "ASDShellT3::setResponse() - need to specify more data\n";
            return 0;
        }
        int pointNum = atoi(argv[1]);
        if (pointNum == 1) {

            output.tag("GaussPoint");
            output.attr("number", pointNum);
            output.attr("eta", XI);
            output.attr("neta", ETA);

            theResponse = m_section->setResponse(&argv[2], argc - 2, output);

            output.endTag();
        }

    }
    else if (strcmp(argv[0], "stresses") == 0) {

        output.tag("GaussPoint");
        output.attr("number", 1);
        output.attr("eta", XI);
        output.attr("neta", ETA);

        output.tag("SectionForceDeformation");
        output.attr("classType", m_section->getClassTag());
        output.attr("tag", m_section->getTag());

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

        theResponse = new ElementResponse(this, 2, Vector(8));
    }

    else if (strcmp(argv[0], "strains") == 0) {

        output.tag("GaussPoint");
        output.attr("number", 1);
        output.attr("eta", XI);
        output.attr("neta", ETA);

        output.tag("SectionForceDeformation");
        output.attr("classType", m_section->getClassTag());
        output.attr("tag", m_section->getTag());

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

        theResponse = new ElementResponse(this, 3, Vector(8));
    }

    else if (m_damping && strcmp(argv[0], "dampingStresses") == 0) {

        output.tag("GaussPoint");
        output.attr("number", 1);
        output.attr("eta", XI);
        output.attr("neta", ETA);

        output.tag("SectionForceDeformation");
        output.attr("classType", m_damping->getClassTag());
        output.attr("tag", m_damping->getTag());

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

        theResponse = new ElementResponse(this, 4, Vector(8));
    }

    output.endTag();
    return theResponse;
}

int
ASDShellT3::getResponse(int responseID, Information& eleInfo)
{
    static Vector stresses(8);
    static Vector strains(8);

    switch (responseID) {
    case 1: // global forces
        return eleInfo.setVector(this->getResistingForce());
        break;
    case 2: // stresses
        return eleInfo.setVector(m_section->getStressResultant());
        break;
    case 3: //strain
        return eleInfo.setVector(m_section->getSectionDeformation());
        break;
    case 4: // damping stresses
        return eleInfo.setVector(m_damping->getDampingForce());
        break;
    default:
        return -1;
    }
}

int ASDShellT3::setParameter(const char** argv, int argc, Parameter& param)
{
    return m_section->setParameter(argv, argc, param);;
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

    // Jacobian
    auto& jac = ASDShellT3Globals::instance().jac;

    // Some matrices/vectors
    auto& B = ASDShellT3Globals::instance().B;
    auto& Bd = ASDShellT3Globals::instance().Bd;
    auto& B1 = ASDShellT3Globals::instance().B1;
    auto& B1TD = ASDShellT3Globals::instance().B1TD;
    auto& N = ASDShellT3Globals::instance().N;
    auto& dN = ASDShellT3Globals::instance().dN;
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

    // Drilling strain-displacement matrix at center for reduced integration
    static Vector drill_dstrain(8);
    static Vector drill_dstress(8);
    static Vector drill_dstress_el(8);

    // Current integration point data
    shapeFunctions(XI, ETA, N);
    shapeFunctionsNaturalDerivatives(XI, ETA, dN);
    jac.calculate(reference_cs, dN);
    double dA = WTS * jac.detJ;

    // Strain-displacement matrix
    computeBMatrix(
        reference_cs.P1(),
        reference_cs.P2(),
        reference_cs.P3(),
        reference_cs.Area(),
        B, Bd
    );

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

        // Update section
        result += m_section->setTrialSectionDeformation(E);

        // Drilling strain Ed = Bd*UL
        double Ed = 0.0;
        for (int i = 0; i < 18; i++)
            Ed += Bd(i) * UL(i);
        m_drill_strain = Ed;
    }

    // Invert bending terms for correct statement of equilibrium
    if ((options & OPT_RHS) || (options & OPT_LHS))
        invertBBendingTerms(B, B1);

    // Integrate RHS
    if (options & OPT_RHS)
    {
        // Section force
        if (m_angle != 0.0) {
            auto& Ssection = m_section->getStressResultant();
            S.addMatrixVector(0.0, Rs, Ssection, 1.0);
            if (m_damping) {
                m_damping->update(Ssection);
                auto& Sdsection = m_damping->getDampingForce();
                S.addMatrixVector(1.0, Rs, Sdsection, 1.0);
            }
        }
        else {
            S = m_section->getStressResultant();
            if (m_damping) {
                m_damping->update(S);
                S += m_damping->getDampingForce();
            }
        }

        // Add current integration point contribution (RHS)
        RHS.addMatrixTransposeVector(1.0, B1, S, dA);

        // Compute drilling damages
        if (m_drill_mode == DrillingDOF_NonLinear) {
            drill_dstrain = m_section->getSectionDeformation();
            drill_dstrain.addVector(1.0, m_nldrill->strain_comm, -1.0);
            if (drill_dstrain.Norm() > 1.0e-10) {
                drill_dstress = m_section->getStressResultant();
                drill_dstress.addVector(1.0, m_nldrill->stress_comm, -1.0);
                const Matrix& C0 = m_section->getInitialTangent();
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
        if (m_drill_mode == DrillingDOF_NonLinear)
            Sd *= (1.0 - m_nldrill->damage_comm);
        for (int i = 0; i < 18; i++)
            RHS(i) += Bd(i) * Sd * dA;
    }

    // AGQI: due to the static condensation, the following items
    // must be computed if we need either the RHS or the LHS, or both of them
    if (options & OPT_LHS)
    {
        // Section tangent
        if (m_angle != 0.0) {
            Dsection = (options & OPT_LHS_IS_INITIAL) ?
                m_section->getInitialTangent() :
                m_section->getSectionTangent();
            if (m_damping) Dsection *= m_damping->getStiffnessMultiplier();
            auto& RsT = ASDShellT3Globals::instance().RsT;
            RsT.addMatrixTranspose(0.0, Rs, 1.0);
            auto& DRsT = ASDShellT3Globals::instance().DRsT;
            DRsT.addMatrixProduct(0.0, Dsection, RsT, 1.0);
            D.addMatrixProduct(0.0, Rs, DRsT, 1.0);
        }
        else {
            D = (options & OPT_LHS_IS_INITIAL) ?
                m_section->getInitialTangent() :
                m_section->getSectionTangent();
            if (m_damping) D *= m_damping->getStiffnessMultiplier();
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
            drill_tang *= (1.0 - m_nldrill->damage_comm);
        for (int i = 0; i < 18; i++)
            for (int j = 0; j < 18; j++)
                LHS(i, j) += Bd(i) * drill_tang * Bd(j) * dA;
    }

    // Transform LHS to global coordinate system
    m_transformation->transformToGlobal(local_cs, UG, UL, LHS, RHS, (options & OPT_LHS));

    // Subtract external loads if any
    if ((options & OPT_RHS) && m_load)
        RHS.addVector(1.0, *m_load, -1.0);

    // Done
    return result;
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
