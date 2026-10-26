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
// Implementation of a corotational coordinate transformation 4-node shells
//

#ifndef ASDShellQ4CorotationalTransformation_h
#define ASDShellQ4CorotationalTransformation_h

#include <ASDEICR.h>
#include <ASDShellQ4Transformation.h>

// this is experimental: it fits the corotational frame following the polar
// decomposition rather than the 1-2 side alignment as per Felippa's work
#define USE_POLAR_DECOMP_ALLIGN

/** \brief ASDShellQ4CorotationalTransformation
*
* This class represents a corotational (nonlinear) coordinate transformation
* that can be used by any element whose geometry is a QUAD 4 in 3D space,
* with 6 D.O.F.s per node.
* Its main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
*    removing rigid body displacements and rotations.
* 3) Transform the outgoing matrices and vectors in global coordinate system
*    with rigid body displacements and rotations.
*
* - Makes use of Quaternions (Euler Parameters) to parametrize finite rotations
*   in an efficient and robust way.
*
* References:
*
* - C.A.Felippa,B.Haugen, "Unified formulation of small-strain corotational
*   finite elements: I. Theory",
*   CU-CAS-05-02, January 2005
* - C.A.Felippa, AFEM.Ch.38, "Quadrilateral Shell Elements",
*   Chapter 5 of B.Haugen's Thesis.
*   link: http://www.colorado.edu/engineering/CAS/courses.d/AFEM.d/
*/
class ASDShellQ4CorotationalTransformation : public ASDShellQ4Transformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 4> NodeContainerType;

public:

    ASDShellQ4CorotationalTransformation()
        : ASDShellQ4Transformation()
    {
    }

    virtual ~ASDShellQ4CorotationalTransformation()
    {
    }

public:

    virtual ASDShellQ4Transformation* create()const
    {
        return new ASDShellQ4CorotationalTransformation();
    }

    virtual bool isLinear() const
    {
        return false;
    }

    virtual void revertToStart()
    {
        // create the reference (undeformed configuration) coordinate system
        ASDShellQ4LocalCoordinateSystem LCS = createReferenceCoordinateSystem();

        // save reference orientation and center
        m_Q0 = QuaternionType::FromRotationMatrix(LCS.Orientation());
        m_C0 = LCS.Center();

        // save initial rotations, no need to take current rotation
        // since we will remove the initial ones (themselves)...
        for (int i = 0; i < 4; i++)
        {
            m_RV[i] = Vector3Type(0.0, 0.0, 0.0);
            m_QN[i] = QuaternionType::FromRotationVector(m_RV[i]);

            m_RV_converged[i] = m_RV[i];
            m_QN_converged[i] = m_QN[i];
        }
    }

    virtual void setDomain(Domain* domain, const ID& node_ids, bool initialized)
    {
        // call base class setDomain to
        // get nodes and save initial displacements and rotations
        ASDShellQ4Transformation::setDomain(domain, node_ids, initialized);

        // quick return
        if (domain == nullptr || initialized)
            return;

        // init state variables
        revertToStart();
    }

    virtual void revertToLastCommit()
    {
        for (int i = 0; i < 4; i++)
        {
            m_RV[i] = m_RV_converged[i];
            m_QN[i] = m_QN_converged[i];
        }
    }

    virtual void commit()
    {
        for (int i = 0; i < 4; i++)
        {
            m_RV_converged[i] = m_RV[i];
            m_QN_converged[i] = m_QN[i];
        }
    }

    virtual void update(const VectorType& globalDisplacements)
    {
        for (int i = 0; i < 4; i++)
        {
            // compute current rotation vector removing initial rotations if any
            Vector3Type currentRotVec;
            int index = i * 6;
            currentRotVec(0) = globalDisplacements(index + 3) - m_U0(index + 3);
            currentRotVec(1) = globalDisplacements(index + 4) - m_U0(index + 4);
            currentRotVec(2) = globalDisplacements(index + 5) - m_U0(index + 5);

            // compute incremental rotation vector
            Vector3Type incrementalRotation = currentRotVec - m_RV[i];

            // save current rotation vector
            m_RV[i] = currentRotVec;

            // compute incremental quaternion from incremental rotation vector
            QuaternionType incrementalQuaternion = QuaternionType::FromRotationVector(incrementalRotation);

            // update nodal quaternion
            m_QN[i] = incrementalQuaternion * m_QN[i];
        }
    }

    virtual ASDShellQ4LocalCoordinateSystem createLocalCoordinateSystem(const VectorType& globalDisplacements)const
    {
        // reference coordinate system
        ASDShellQ4LocalCoordinateSystem a = createReferenceCoordinateSystem();

        // compute nodal positions at current configuration removing initial displacements if any
        std::array<Vector3Type, 4> def = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds()),
            Vector3Type(m_nodes[3]->getCrds())
        };
        for (int i = 0; i < 4; i++) {
            int index = i * 6;
            Vector3Type& iP = def[i];
            iP(0) += globalDisplacements(index) - m_U0(index);
            iP(1) += globalDisplacements(index + 1) - m_U0(index + 1);
            iP(2) += globalDisplacements(index + 2) - m_U0(index + 2);
        }

        // current coordinate system
        ASDShellQ4LocalCoordinateSystem b(def[0], def[1], def[2], def[3]);

#ifndef USE_POLAR_DECOMP_ALLIGN
        return b;
#endif // !USE_POLAR_DECOMP_ALLIGN

        double aX1 = a.X1(); double aY1 = a.Y1();
        double bX1 = b.X1(); double bY1 = b.Y1();
        double aX2 = a.X2(); double aY2 = a.Y2();
        double bX2 = b.X2(); double bY2 = b.Y2();
        double aX3 = a.X3(); double aY3 = a.Y3();
        double bX3 = b.X3(); double bY3 = b.Y3();
        double aX4 = a.X4(); double aY4 = a.Y4();
        double bX4 = b.X4(); double bY4 = b.Y4();

        // now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
        // which is the same for both coordinate systems.
        // now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

        double C1 = 1.0 / (aX1 * aY2 - aX2 * aY1 - aX1 * aY4 + aX2 * aY3 - aX3 * aY2 + aX4 * aY1 + aX3 * aY4 - aX4 * aY3);
        double C2 = bY1 / 4.0 + bY2 / 4.0 - bY3 / 4.0 - bY4 / 4.0;
        double C3 = bY1 / 4.0 - bY2 / 4.0 - bY3 / 4.0 + bY4 / 4.0;
        double C4 = bX1 / 4.0 + bX2 / 4.0 - bX3 / 4.0 - bX4 / 4.0;
        double C5 = bX1 / 4.0 - bX2 / 4.0 - bX3 / 4.0 + bX4 / 4.0;
        double C6 = aX1 + aX2 - aX3 - aX4;
        double C7 = aX1 - aX2 - aX3 + aX4;
        double C8 = aY1 + aY2 - aY3 - aY4;
        double C9 = aY1 - aY2 - aY3 + aY4;
        double f11 = 2.0 * C1 * C5 * C8 - 2.0 * C1 * C4 * C9;
        double f12 = 2.0 * C1 * C4 * C7 - 2.0 * C1 * C5 * C6;
        double f21 = 2.0 * C1 * C3 * C8 - 2.0 * C1 * C2 * C9;
        double f22 = 2.0 * C1 * C2 * C7 - 2.0 * C1 * C3 * C6;

        // now we can extrapolate the rotation angle that makes this deformation gradient symmetric.
        // F = R*U -> find R such that R'*F = U
        double alpha = std::atan2(f21 - f12, f11 + f22);

        // this final coordinate system is the one in which 
        // the deformation gradient is equal to the stretch tensor
        return ASDShellQ4LocalCoordinateSystem(def[0], def[1], def[2], def[3], alpha);
    }

    virtual void calculateLocalDisplacements(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        // orientation and center of current local coordinate system
        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());
        const Vector3Type& C = LCS.Center();

        for (int i = 0; i < 4; i++)
        {
            int index = i * 6;

            // centered undeformed position
            Vector3Type X0 = Vector3Type(m_nodes[i]->getCrds());
            X0 -= m_C0;

            // centered deformed position
            Vector3Type X = X0 + Vector3Type(globalDisplacements, index);
            X -= C;

            // get deformational displacements
            Q.rotateVector(X);
            m_Q0.rotateVector(X0);
            Vector3Type deformationalDisplacements = X - X0;

            localDisplacements[index] = deformationalDisplacements[0];
            localDisplacements[index + 1] = deformationalDisplacements[1];
            localDisplacements[index + 2] = deformationalDisplacements[2];

            // get deformational rotations
            QuaternionType Qd = Q * m_QN[i] * m_Q0.conjugate();
            Qd.toRotationVector(
                localDisplacements[index + 3],
                localDisplacements[index + 4],
                localDisplacements[index + 5]);
        }
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        // Get the total rotation matrix (local - to - global)
        // Note: do NOT include the warpage correction matrix!
        // Explanation:
        // The Warpage correction matrix computed by the LocalCoordinateSystem is a Linear Projector.
        // It should be used in a LinearCoordinateTransformation.
        // Here instead we already calculate a nonlinear Projector (P = Pu - S * G)!

        static MatrixType T(24, 24);
        LCS.ComputeTotalRotationMatrix(T);

        // Form all matrices:
        // S: Spin-Fitter matrix
        // G: Spin-Lever matrix
        // P: Projector (Translational & Rotational)
        static MatrixType P(24, 24);
        static MatrixType S(24, 3);
        static MatrixType G(3, 24);
        EICR::Compute_Pt(4, P);
        EICR::Compute_S(LCS.Nodes(), S);
        RotationGradient(LCS, globalDisplacements, G);
        P.addMatrixProduct(1.0, S, G, -1.0); // P -= S*G

        // Compute the projected local forces ( pe = P' * RHS ).
        // Note: here the RHS is already given as a residual vector -> - internalForces -> (pe = - Ke * U)
        // so projectedLocalForces = - P' * Ke * U

        static VectorType projectedLocalForces(24);
        projectedLocalForces.addMatrixTransposeVector(0.0, P, RHS, 1.0);

        // Compute the Right-Hand-Side vector in global coordinate system (- T' * P' * Km * U).
        // At this point the computation of the Right-Hand-Side is complete.

        RHS.addMatrixTransposeVector(0.0, T, projectedLocalForces, 1.0);

        // Begin the computation of the Left-Hand-Side Matrix :

        if (!LHSrequired) 
            return; // avoid useless calculations!

        // H: Axial Vector Jacobian
        static MatrixType H(24, 24);
        EICR::Compute_H(localDisplacements, H);

        // Step 1: ( K.M : Material Stiffness Matrix )
        // Apply the projector to the Material Stiffness Matrix (Ke = P' * Km * H * P)
        // At this point 'LHS' contains the 'projected' Material Stiffness matrix
        // in local corotational coordinate system

        static MatrixType temp(24, 24);
        temp.addMatrixProduct(0.0, LHS, H, 1.0);
        LHS.addMatrixProduct(0.0, temp, P, 1.0);
        temp.addMatrixTransposeProduct(0.0, P, LHS, 1.0);
        LHS = temp;

        // Step 2: ( K.GP: Equilibrium Projection Geometric Stiffness Matrix )
        // First assemble the 'Fnm' matrix with the Spins of the nodal forces.
        // Actually at this point the 'Fnm' Matrix is the 'Fn' Matrix,
        // because it only contains the spins of the 'translational' forces.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P))

        static MatrixType Fnm(24, 3);
        Fnm.Zero();
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 0);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 6);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 12);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 18);

        static MatrixType FnmT(3, 24);
        FnmT.addMatrixTranspose(0.0, Fnm, 1.0);

        temp.addMatrixTransposeProduct(0.0, G, FnmT, 1.0);
        // the element accumulates the internal force with a positive sign (+B'*S*dA),
        // so the geometric terms are subtracted
        LHS.addMatrixProduct(1.0, temp, P, -1.0);

        // Step 3: ( K.GR: Rotational Geometric Stiffness Matrix )
        // Add the Spins of the nodal moments to 'Fnm'.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G))

        EICR::Spin_AtRow(projectedLocalForces, Fnm, 3);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 9);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 15);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 21);

        LHS.addMatrixProduct(1.0, Fnm, G, -1.0);

        // Step 4: (Global Stiffness Matrix)
        // Transform the LHS to the Global coordinate system.
        // T' * [(P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G)] * T
        temp.addMatrixProduct(0.0, LHS, T, 1.0);
        LHS.addMatrixTransposeProduct(0.0, T, temp, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType globalDisplacements(24);
        static VectorType localDisplacements(24);
        computeGlobalDisplacements(globalDisplacements);
        calculateLocalDisplacements(LCS, globalDisplacements, localDisplacements);
        transformToGlobal(LCS, globalDisplacements, localDisplacements, LHS, RHS, LHSrequired);
    }

    virtual int internalDataSize() const
    {
        // 24 -> initial displacement +
        // 9*4 -> 9 quaternions +
        // 9*3 -> 9 3d vectors
        return 87;
    }

    virtual void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellQ4CorotationalTransformation - failed to save internal data: vector too small\n";
            exit(-1);
        }

        // 24 -> initial displacement +
        for (int i = 0; i < 24; i++)
            v(pos++) = m_U0(i);

        // 9*4 -> 9 quaternions +
        auto lamq = [&v, &pos](const QuaternionType& x) {
            v(pos++) = x.w();
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamq(m_Q0);
        for (int i = 0; i < 4; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 4; i++)
            lamq(m_QN_converged[i]);

        // 9*3 -> 9 3d vectors +
        auto lamv = [&v, &pos](const Vector3Type& x) {
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamv(m_C0);
        for (int i = 0; i < 4; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 4; i++)
            lamv(m_RV_converged[i]);
    }

    virtual void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellQ4CorotationalTransformation - failed to restore internal data: vector too small\n";
            exit(-1);
        }
        
        // 24 -> initial displacement +
        for (int i = 0; i < 24; i++)
            m_U0(i) = v(pos++);

        // 9*4 -> 9 quaternions +
        auto lamq = [&v, &pos](QuaternionType& x) {
            x = QuaternionType(v(pos), v(pos+1), v(pos+2), v(pos+3));
            pos += 4;
        };
        lamq(m_Q0);
        for (int i = 0; i < 4; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 4; i++)
            lamq(m_QN_converged[i]);

        // 9*3 -> 9 3d vectors +
        auto lamv = [&v, &pos](Vector3Type& x) {
            x = Vector3Type(v(pos++), v(pos++), v(pos++));
        };
        lamv(m_C0);
        for (int i = 0; i < 4; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 4; i++)
            lamv(m_RV_converged[i]);
    }

private:

    /**
    * Computes the Spin Fitter Matrix, i.e. the rotation gradient G in dw = G * dv,
    * where dw is the spin of the corotational frame in local components.
    * This is the only matrix not included in the EICR, because it depends on how
    * the corotational frame follows the element.
    * With the polar decomposition fitting, the frame is built in three steps,
    * and G follows them:
    * - e3 from the cross product of the deformed diagonals (rows 0 and 1);
    * - e1 from side 1-2 projected onto the plane normal to e3 (row 2, part 1);
    * - the in-plane polar angle alpha = atan2(f21 - f12, f11 + f22) (row 2, part 2).
    * Since alpha is a rotation about the e3 shared by both frames, row 2 is the
    * sum of the drilling spin of the side-aligned frame and dalpha.
    * The frame depends on the nodal translations only, so the rotational
    * columns of G are zero.
    * @return the Spin Fitter Matrix
    */
    inline void RotationGradient(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements, 
        MatrixType& G)
    {
        G.Zero();

#ifdef USE_POLAR_DECOMP_ALLIGN

        // deformed nodal positions, as in createLocalCoordinateSystem
        std::array<Vector3Type, 4> p = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds()),
            Vector3Type(m_nodes[3]->getCrds())
        };
        for (int i = 0; i < 4; i++) {
            int index = i * 6;
            for (int k = 0; k < 3; k++)
                p[i](k) += globalDisplacements(index + k) - m_U0(index + k);
        }

        // side-aligned current frame (b in createLocalCoordinateSystem).
        // alpha is defined from its local coordinates.
        ASDShellQ4LocalCoordinateSystem b(p[0], p[1], p[2], p[3]);
        const MatrixType& Ob = b.Orientation();
        const auto& xb = b.Nodes();

        // final frame (rotated by alpha): rows are e1, e2, e3
        const MatrixType& Of = LCS.Orientation();
        Vector3Type e1(Of(0, 0), Of(0, 1), Of(0, 2));
        Vector3Type e2(Of(1, 0), Of(1, 1), Of(1, 2));
        Vector3Type e3(Of(2, 0), Of(2, 1), Of(2, 2));
        // in-plane y axis of the side-aligned frame (e3 is the same)
        Vector3Type e2t(Ob(1, 0), Ob(1, 1), Ob(1, 2));

        // geometry
        Vector3Type d13 = p[2] - p[0];
        Vector3Type d24 = p[3] - p[1];
        double nn = 2.0 * b.Area(); // ||d13 x d24||
        Vector3Type v = p[1] - p[0];
        double v_e3 = v.dot(e3);
        Vector3Type Pv = v - v_e3 * e3;
        double lPv = Pv.norm();
        if (nn < 1.0e-300 || lPv < 1.0e-300)
            return; // degenerate element
        double inn = 1.0 / nn;

        // constants of the polar angle from the reference local coordinates
        ASDShellQ4LocalCoordinateSystem a = createReferenceCoordinateSystem();
        double aX1 = a.X1(); double aY1 = a.Y1();
        double aX2 = a.X2(); double aY2 = a.Y2();
        double aX3 = a.X3(); double aY3 = a.Y3();
        double aX4 = a.X4(); double aY4 = a.Y4();
        double C1 = 1.0 / (aX1 * aY2 - aX2 * aY1 - aX1 * aY4 + aX2 * aY3
            - aX3 * aY2 + aX4 * aY1 + aX3 * aY4 - aX4 * aY3);
        double C6 = aX1 + aX2 - aX3 - aX4;
        double C7 = aX1 - aX2 - aX3 + aX4;
        double C8 = aY1 + aY2 - aY3 - aY4;
        double C9 = aY1 - aY2 - aY3 + aY4;

        // current s = f21 - f12 and c = f11 + f22. Both are linear in the
        // local coordinates of the side-aligned frame, so their variations
        // are linear in the nodal increments.
        double C2 = 0.25 * (xb[0][1] + xb[1][1] - xb[2][1] - xb[3][1]);
        double C3 = 0.25 * (xb[0][1] - xb[1][1] - xb[2][1] + xb[3][1]);
        double C4 = 0.25 * (xb[0][0] + xb[1][0] - xb[2][0] - xb[3][0]);
        double C5 = 0.25 * (xb[0][0] - xb[1][0] - xb[2][0] + xb[3][0]);
        double ss = 2.0 * C1 * (C3 * C8 - C2 * C9 - C4 * C7 + C5 * C6);
        double cc = 2.0 * C1 * (C5 * C8 - C4 * C9 + C2 * C7 - C3 * C6);
        double den = ss * ss + cc * cc;
        if (den < 1.0e-300)
            return;
        double iden = 1.0 / den;

        // one column for each node i and local direction j
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {

                // local axis j in global components
                Vector3Type du(Of(j, 0), Of(j, 1), Of(j, 2));

                // variation of n = d13 x d24 due to du at node i
                Vector3Type dn;
                switch (i) {
                case 0:  dn = d24.cross(du); break;
                case 1:  dn = du.cross(d13); break;
                case 2:  dn = du.cross(d24); break;
                default: dn = d13.cross(du); break;
                }

                // tilt of e3: dw0 = -de3.e2, dw1 = de3.e1
                double e2tdn = e2t.dot(dn);
                double dw0 = -e2.dot(dn) * inn;
                double dw1 = e1.dot(dn) * inn;

                // drilling spin of the side-aligned frame
                double dvterm = 0.0;
                if (i == 0)      dvterm = -e2t.dot(du);
                else if (i == 1) dvterm = e2t.dot(du);
                double dw2t = (dvterm - v_e3 * e2tdn * inn) / lPv;

                // variation of the local coordinates in the side-aligned frame.
                // The tilt of the frame would add -dw1t*z and +dw0t*z, but these
                // terms cancel: since n is orthogonal to both diagonals, the local
                // z of the nodes is (h, -h, h, -h), which is orthogonal to the
                // nodal weights (+,+,-,-) and (+,-,-,+) used in C2..C5.
                double dbX[4], dbY[4];
                for (int k = 0; k < 4; k++) {
                    double f = (k == i) ? 0.75 : -0.25;
                    Vector3Type dk = f * du;
                    double ax = Ob(0, 0) * dk(0) + Ob(0, 1) * dk(1) + Ob(0, 2) * dk(2);
                    double ay = Ob(1, 0) * dk(0) + Ob(1, 1) * dk(1) + Ob(1, 2) * dk(2);
                    dbX[k] = ax + dw2t * xb[k][1];
                    dbY[k] = ay - dw2t * xb[k][0];
                }
                double dC2 = 0.25 * (dbY[0] + dbY[1] - dbY[2] - dbY[3]);
                double dC3 = 0.25 * (dbY[0] - dbY[1] - dbY[2] + dbY[3]);
                double dC4 = 0.25 * (dbX[0] + dbX[1] - dbX[2] - dbX[3]);
                double dC5 = 0.25 * (dbX[0] - dbX[1] - dbX[2] + dbX[3]);
                double ds = 2.0 * C1 * (dC3 * C8 - dC2 * C9 - dC4 * C7 + dC5 * C6);
                double dc = 2.0 * C1 * (dC5 * C8 - dC4 * C9 + dC2 * C7 - dC3 * C6);

                // variation of the polar angle
                double dalpha = (cc * ds - ss * dc) * iden;

                int col = i * 6 + j;
                G(0, col) = dw0;
                G(1, col) = dw1;
                G(2, col) = dw2t + dalpha;
            }
        }

#else // !USE_POLAR_DECOMP_ALLIGN

        const auto& P1 = LCS.P1();
        const auto& P2 = LCS.P2();
        const auto& P3 = LCS.P3();
        const auto& P4 = LCS.P4();

        double Ap = 2.0 * LCS.Area();
        double m = 1.0 / Ap;

        Vector3Type D12(P2 - P1);
        Vector3Type D24(P4 - P2);
        Vector3Type D13(P3 - P1);

        double x42 = D24(0);
        double x24 = -x42;
        double y42 = D24(1);
        double y24 = -y42;
        double x31 = D13(0);
        double x13 = -x31;
        double y31 = D13(1);
        double y13 = -y31;

        // Note, assuming the input vectors are in local CR, 
        // l12 is the length of the side 1-2 projected onto the xy plane.
        double l12 = std::sqrt(D12(0) * D12(0) + D12(1) * D12(1));

        // G1

        G(0, 2) = x42 * m;
        G(1, 2) = y42 * m;
        G(2, 1) = -1.0 / l12;

        // G2

        G(0, 8) = x13 * m;
        G(1, 8) = y13 * m;
        G(2, 7) = 1.0 / l12;

        // G3

        G(0, 14) = x24 * m;
        G(1, 14) = y24 * m;

        // G4

        G(0, 20) = x31 * m;
        G(1, 20) = y31 * m;

#endif // USE_POLAR_DECOMP_ALLIGN

    }

private:

    QuaternionType m_Q0;
    Vector3Type m_C0;
    std::array<QuaternionType, 4> m_QN;
    std::array<Vector3Type, 4> m_RV;
    std::array<QuaternionType, 4> m_QN_converged;
    std::array<Vector3Type, 4> m_RV_converged;
};

#endif // !ASDShellQ4CorotationalTransformation_h
