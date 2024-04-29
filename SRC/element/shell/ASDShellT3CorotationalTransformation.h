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
// Implementation of a corotational coordinate transformation 3-node shells
//

#ifndef ASDShellT3CorotationalTransformation_h
#define ASDShellT3CorotationalTransformation_h

#include <ASDEICR.h>
#include <ASDShellT3Transformation.h>

// this is experimental: it fits the corotational frame following the polar
// decomposition rather than the 1-2 side alignment as per Felippa's work
#define USE_POLAR_DECOMP_ALLIGN

/** \brief ASDShellT3CorotationalTransformation
*
* This class represents a corotational (nonlinear) coordinate transformation
* that can be used by any element whose geometry is a Triangle 3 in 3D space,
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
class ASDShellT3CorotationalTransformation : public ASDShellT3Transformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 3> NodeContainerType;

public:

    ASDShellT3CorotationalTransformation()
        : ASDShellT3Transformation()
    {
    }

    virtual ~ASDShellT3CorotationalTransformation()
    {
    }

public:

    virtual ASDShellT3Transformation* create()const
    {
        return new ASDShellT3CorotationalTransformation();
    }

    virtual bool isLinear() const
    {
        return false;
    }

    virtual void revertToStart()
    {
        // create the reference (undeformed configuration) coordinate system
        ASDShellT3LocalCoordinateSystem LCS = createReferenceCoordinateSystem();

        // save reference orientation and center
        m_Q0 = QuaternionType::FromRotationMatrix(LCS.Orientation());
        m_C0 = LCS.Center();

        // save initial rotations, no need to take current rotation
        // since we will remove the initial ones (themselves)...
        for (int i = 0; i < 3; i++)
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
        ASDShellT3Transformation::setDomain(domain, node_ids, initialized);

        // quick return
        if (domain == nullptr || initialized)
            return;

        // init state variables
        revertToStart();
    }

    virtual void revertToLastCommit()
    {
        for (int i = 0; i < 3; i++)
        {
            m_RV[i] = m_RV_converged[i];
            m_QN[i] = m_QN_converged[i];
        }
    }

    virtual void commit()
    {
        for (int i = 0; i < 3; i++)
        {
            m_RV_converged[i] = m_RV[i];
            m_QN_converged[i] = m_QN[i];
        }
    }

    virtual void update(const VectorType& globalDisplacements)
    {
        for (int i = 0; i < 3; i++)
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

    virtual ASDShellT3LocalCoordinateSystem createLocalCoordinateSystem(const VectorType& globalDisplacements)const
    {
        // reference coordinate system
        ASDShellT3LocalCoordinateSystem a = createReferenceCoordinateSystem();

        // compute nodal positions at current configuration removing initial displacements if any
        std::array<Vector3Type, 3> def = {
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds())
        };
        for (int i = 0; i < 3; i++) {
            int index = i * 6;
            Vector3Type& iP = def[i];
            iP(0) += globalDisplacements(index) - m_U0(index);
            iP(1) += globalDisplacements(index + 1) - m_U0(index + 1);
            iP(2) += globalDisplacements(index + 2) - m_U0(index + 2);
        }

        // current coordinate system
        ASDShellT3LocalCoordinateSystem b(def[0], def[1], def[2]);

#ifndef USE_POLAR_DECOMP_ALLIGN
        return b;
#endif // !USE_POLAR_DECOMP_ALLIGN

        double aX1 = a.X1(); double aY1 = a.Y1();
        double bX1 = b.X1(); double bY1 = b.Y1();
        double aX2 = a.X2(); double aY2 = a.Y2();
        double bX2 = b.X2(); double bY2 = b.Y2();
        double aX3 = a.X3(); double aY3 = a.Y3();
        double bX3 = b.X3(); double bY3 = b.Y3();

        // now we are in the local coordinate systems (reference and current), i.e. we are looking in the local Z direction
        // which is the same for both coordinate systems.
        // now we can compute the 2D deformation gradient between the 2 configurations, at the element center.

        double C1 = 1.0 / (aX1*aY2 - aX1*aY3 - aX2*aY1 + aX2*aY3 + aX3*aY1 - aX3*aY2);
        double f11 = C1*(-(aY1 - aY2)*(bX1 - bX3) + (aY1 - aY3)*(bX1 - bX2));
        double f12 = C1*((aX1 - aX2)*(bX1 - bX3) - (aX1 - aX3)*(bX1 - bX2));
        double f21 = C1*(-(aY1 - aY2)*(bY1 - bY3) + (aY1 - aY3)*(bY1 - bY2));
        double f22 = C1*((aX1 - aX2)*(bY1 - bY3) - (aX1 - aX3)*(bY1 - bY2));

        // now we can extrapolate the rotation angle that makes this deformation gradient symmetric.
        // F = R*U -> find R such that R'*F = U
        double alpha = std::atan2(f21 - f12, f11 + f22);

        // this final coordinate system is the one in which 
        // the deformation gradient is equal to the stretch tensor
        return ASDShellT3LocalCoordinateSystem(def[0], def[1], def[2], alpha);
    }

    virtual void calculateLocalDisplacements(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        // orientation and center of current local coordinate system
        QuaternionType Q = QuaternionType::FromRotationMatrix(LCS.Orientation());
        const Vector3Type& C = LCS.Center();

        for (int i = 0; i < 3; i++)
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
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        // Get the total rotation matrix (local - to - global)

        static MatrixType T(18, 18);
        LCS.ComputeTotalRotationMatrix(T);

        // Form all matrices:
        // S: Spin-Fitter matrix
        // G: Spin-Lever matrix
        // P: Projector (Translational & Rotational)
        static MatrixType P(18, 18);
        static MatrixType S(18, 3);
        static MatrixType G(3, 18);
        EICR::Compute_Pt(3, P);
        EICR::Compute_S(LCS.Nodes(), S);
        RotationGradient(LCS, globalDisplacements, G);
        P.addMatrixProduct(1.0, S, G, -1.0); // P -= S*G

        // Compute the projected local forces ( pe = P' * RHS ).
        // Note: here the RHS is already given as a residual vector -> - internalForces -> (pe = - Ke * U)
        // so projectedLocalForces = - P' * Ke * U

        static VectorType projectedLocalForces(18);
        projectedLocalForces.addMatrixTransposeVector(0.0, P, RHS, 1.0);

        // Compute the Right-Hand-Side vector in global coordinate system (- T' * P' * Km * U).
        // At this point the computation of the Right-Hand-Side is complete.

        RHS.addMatrixTransposeVector(0.0, T, projectedLocalForces, 1.0);

        // Begin the computation of the Left-Hand-Side Matrix :

        if (!LHSrequired) 
            return; // avoid useless calculations!

        // H: Axial Vector Jacobian
        static MatrixType H(18, 18);
        EICR::Compute_H(localDisplacements, H);

        // Step 1: ( K.M : Material Stiffness Matrix )
        // Apply the projector to the Material Stiffness Matrix (Ke = P' * Km * H * P)
        // At this point 'LHS' contains the 'projected' Material Stiffness matrix
        // in local corotational coordinate system

        static MatrixType temp(18, 18);
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

        static MatrixType Fnm(18, 3);
        Fnm.Zero();
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 0);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 6);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 12);

        static MatrixType FnmT(3, 18);
        FnmT.addMatrixTranspose(0.0, Fnm, 1.0);

        temp.addMatrixTransposeProduct(0.0, G, FnmT, 1.0);
        LHS.addMatrixProduct(1.0, temp, P, 1.0); // note: '+' not '-' because the RHS already has the negative sign

        // Step 3: ( K.GR: Rotational Geometric Stiffness Matrix )
        // Add the Spins of the nodal moments to 'Fnm'.
        // At this point 'LHS' contains also this term of the Geometric stiffness
        // (Ke = (P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G))

        EICR::Spin_AtRow(projectedLocalForces, Fnm, 3);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 9);
        EICR::Spin_AtRow(projectedLocalForces, Fnm, 15);

        LHS.addMatrixProduct(1.0, Fnm, G, 1.0); // note: '+' not '-' because the RHS already has the negative sign

        // Step 4: (Global Stiffness Matrix)
        // Transform the LHS to the Global coordinate system.
        // T' * [(P' * Km * H * P) - (G' * Fn' * P) - (Fnm * G)] * T
        temp.addMatrixProduct(0.0, LHS, T, 1.0);
        LHS.addMatrixTransposeProduct(0.0, T, temp, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellT3LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType globalDisplacements(18);
        static VectorType localDisplacements(18);
        computeGlobalDisplacements(globalDisplacements);
        calculateLocalDisplacements(LCS, globalDisplacements, localDisplacements);
        transformToGlobal(LCS, globalDisplacements, localDisplacements, LHS, RHS, LHSrequired);
    }

    virtual int internalDataSize() const
    {
        // 18 -> initial displacement +
        // 7*4 -> 7 quaternions +
        // 7*3 -> 7 3d vectors
        return 67;
    }

    virtual void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3CorotationalTransformation - failed to save internal data: vector too small\n";
            exit(-1);
        }

        // 18 -> initial displacement +
        for (int i = 0; i < 18; i++)
            v(pos++) = m_U0(i);

        // 7*4 -> 7 quaternions +
        auto lamq = [&v, &pos](const QuaternionType& x) {
            v(pos++) = x.w();
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamq(m_Q0);
        for (int i = 0; i < 3; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 3; i++)
            lamq(m_QN_converged[i]);

        // 7*3 -> 7 3d vectors +
        auto lamv = [&v, &pos](const Vector3Type& x) {
            v(pos++) = x.x();
            v(pos++) = x.y();
            v(pos++) = x.z();
        };
        lamv(m_C0);
        for (int i = 0; i < 3; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 3; i++)
            lamv(m_RV_converged[i]);
    }

    virtual void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3CorotationalTransformation - failed to restore internal data: vector too small\n";
            exit(-1);
        }
        
        // 18 -> initial displacement +
        for (int i = 0; i < 18; i++)
            m_U0(i) = v(pos++);

        // 7*4 -> 7 quaternions +
        auto lamq = [&v, &pos](QuaternionType& x) {
            x = QuaternionType(v(pos), v(pos+1), v(pos+2), v(pos+3));
            pos += 4;
        };
        lamq(m_Q0);
        for (int i = 0; i < 3; i++)
            lamq(m_QN[i]);
        for (int i = 0; i < 3; i++)
            lamq(m_QN_converged[i]);

        // 7*3 -> 7 3d vectors +
        auto lamv = [&v, &pos](Vector3Type& x) {
            x = Vector3Type(v(pos++), v(pos++), v(pos++));
        };
        lamv(m_C0);
        for (int i = 0; i < 3; i++)
            lamv(m_RV[i]);
        for (int i = 0; i < 3; i++)
            lamv(m_RV_converged[i]);
    }

private:

    /**
    * Computes the Spin Fitter Matrix.
    * This is the only matrix not included in the EICR, because it depends on how
    * the corotational frame follows the element.
    * This implementation works for the fitting based on polar decomposition.
    * For the moment the gradient is calculated numerically (perturbation).
    * TODO: Find a closed form of the derivative.
    * @return the Spin Fitter Matrix
    */
    inline void RotationGradient(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements, 
        MatrixType& G)
    {
        G.Zero();

#ifdef USE_POLAR_DECOMP_ALLIGN

        /**
        Note: when attaching the local system using the polar decomposition, it turns out
        that the G matrix should be exactly 0!
        */

        //// perturbation relative to initial area
        //double pert = std::sqrt(createReferenceCoordinateSystem().Area()) * 1.0e-8;

        //// un-perturbed coordinate system
        //auto CS0 = createLocalCoordinateSystem(globalDisplacements);
        //Vector3Type e30 = CS0.Vz();
        //Vector3Type e10 = CS0.Vx();
        //QuaternionType Q0 = QuaternionType::FromRotationMatrix(CS0.Orientation());

        //// perturbed displacement vector
        //static Vector UGP(24);
        //UGP = globalDisplacements;

        //// for each node...
        //for (int i = 0; i < 3; i++)
        //{
        //    int index = i * 6;

        //    // for each component in [x, y, z] ...
        //    for (int j = 0; j < 3; j++)
        //    {
        //        // apply perturbation
        //        UGP(index + j) += pert;

        //        // current coordinate system (perturbed at node i component j)
        //        auto CSi = createLocalCoordinateSystem(UGP);

        //        // save the (numerical) rotation gradient

        //        Vector3Type e3 = CSi.Vz();// - e30;
        //        Vector3Type e1 = CSi.Vx();// - e10;
        //        Q0.rotateVector(e3);
        //        Q0.rotateVector(e1);

        //        G(0, index + j) = e3(1) / pert; // - d(e3.y)/dxj , where xj is x,y,z for j=0,1,2
        //        G(1, index + j) = -e3(0) / pert; // + d(e3.x)/dxj , where xj is x,y,z for j=0,1,2
        //        G(2, index + j) = -e1(1) / pert; // + d(e1.y)/dxj , where xj is x,y,z for j=0,1,2

        //        UGP(index + j) = globalDisplacements(index + j); // restore the current coordinate
        //    }
        //}

#else // !USE_POLAR_DECOMP_ALLIGN

        const auto& P1 = LCS.P1();
        const auto& P2 = LCS.P2();
        const auto& P3 = LCS.P3();

        double Ap = 2.0 * LCS.Area();
        double m = 1.0 / Ap;

        double x13 = P1(0) - P3(0);
        double y13 = P1(1) - P3(1);
        double x32 = P3(0) - P2(0);
        double y32 = P3(1) - P2(1);

        double L3 = (P2 - P1).norm();
        double h3 = Ap/L3;

        // G1

        G(0, 2) = x32 * m;
        G(1, 2) = y32 * m;
        G(2, 1) = -h3;

        // G2

        G(0, 8) = x13 * m;
        G(1, 8) = y13 * m;
        G(2, 7) = h3;

        // G3

        G(0, 14) = x32 * m;
        G(1, 14) = y32 * m;

#endif // USE_POLAR_DECOMP_ALLIGN

    }

private:

    QuaternionType m_Q0;
    Vector3Type m_C0;
    std::array<QuaternionType, 3> m_QN;
    std::array<Vector3Type, 3> m_RV;
    std::array<QuaternionType, 3> m_QN_converged;
    std::array<Vector3Type, 3> m_RV_converged;
};

#endif // !ASDShellT3CorotationalTransformation_h
