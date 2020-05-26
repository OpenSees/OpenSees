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
// Implementation of a linear coordinate transformation 4-node shells
//

#ifndef ASDShellQ4Transformation_h
#define ASDShellQ4Transformation_h

#include <ASDShellQ4LocalCoordinateSystem.h>
#include <Node.h>
#include <Domain.h>

/** \brief ASDShellQ4Transformation
*
* This class represents a basic (linear) coordinate transformation that can be used
* by any element whose geometry is a QUAD 4 in 3D space, with 6 D.O.F.s per node.
* Its main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
* 3) Transform the outgoing matrices and vectors in global coordinate system
*/
class ASDShellQ4Transformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 4> NodeContainerType;

public:

    ASDShellQ4Transformation()
    {
    }

    virtual ~ASDShellQ4Transformation()
    {
    }

private:

    ASDShellQ4Transformation(const ASDShellQ4Transformation& other) = delete;

    ASDShellQ4Transformation& operator = (const ASDShellQ4Transformation& other) = delete;

public:

    virtual ASDShellQ4Transformation* create()const
    {
        return new ASDShellQ4Transformation();
    }

    virtual bool isLinear() const
    {
        return true;
    }

    virtual void revertToStart()
    {
    }

    virtual void setDomain(Domain* domain, const ID& node_ids)
    {
        // get nodes and save initial displacements and rotations
        for (size_t i = 0; i < 4; i++) {
            m_nodes[i] = domain->getNode(node_ids(i));
            if (m_nodes[i] == nullptr) {
                opserr << "ASDShellQ4Transformation::setDomain - no node " << node_ids(i)
                    << " exists in the model\n";
                exit(-1);
            }
            const Vector& iU = m_nodes[i]->getTrialDisp();
            if (iU.Size() != 6) {
                opserr << "ASDShellQ4Transformation::setDomain - node " << node_ids(i)
                    << " has " << iU.Size() << " DOFs, while 6 are expected\n";
                exit(-1);
            }
            size_t index = i * 6;
            for (size_t j = 0; j < 6; j++)
                m_U0(index + j) = iU(j);
        }
    }

    virtual void revertToLastCommit()
    {
    }

    virtual void commit()
    {
    }

    virtual void update()
    {
    }

    virtual ASDShellQ4LocalCoordinateSystem createReferenceCoordinateSystem()const
    {
        return ASDShellQ4LocalCoordinateSystem(
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds()),
            Vector3Type(m_nodes[3]->getCrds())
        );
    }

    virtual ASDShellQ4LocalCoordinateSystem createLocalCoordinateSystem()const
    {
        return createReferenceCoordinateSystem();
    }

    virtual void computeGlobalDisplacements(VectorType& globalDisplacements) const
    {
        for (int i = 0; i < 4; i++) {
            int index = i * 6;
            const VectorType& iU = m_nodes[i]->getTrialDisp();
            for (int j = 0; j < 6; j++) {
                globalDisplacements(index + j) = iU(j) - m_U0(index + j);
            }
        }
    }

    virtual const MatrixType& computeTransformationMatrix(const ASDShellQ4LocalCoordinateSystem& LCS) const
    {
        static MatrixType R(24, 24);
        static MatrixType T(24, 24);
        static MatrixType W(24, 24);
        if (LCS.IsWarped()) {
            LCS.ComputeTotalRotationMatrix(R);
            LCS.ComputeTotalWarpageMatrix(W);
            T.addMatrixProduct(0.0, W, R, 1.0);
        }
        else {
            LCS.ComputeTotalRotationMatrix(T);
        }
        return T;
    }

    virtual void calculateLocalDisplacements(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        const MatrixType& R = computeTransformationMatrix(LCS);
        localDisplacements.addMatrixVector(0.0, R, globalDisplacements, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool RHSrequired,
        bool LHSrequired)
    {
        static MatrixType temp(24, 24);
        const MatrixType& R = computeTransformationMatrix(LCS);
        if (LHSrequired) {
            temp.addMatrixTransposeProduct(0.0, R, LHS, 1.0);
            LHS.addMatrixProduct(0.0, temp, R, 1.0);
        }
        if (RHSrequired) {
            RHS.addMatrixTransposeVector(0.0, R, RHS, 1.0);
        }
    }

    virtual void transformToGlobal(
        const ASDShellQ4LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool RHSrequired,
        bool LHSrequired)
    {
        static VectorType dummy;
        transformToGlobal(LCS, dummy, dummy, LHS, RHS, RHSrequired, LHSrequired);
    }

public:

    inline const NodeContainerType& getNodes()const { return m_nodes; }
    inline NodeContainerType& getNodes() { return m_nodes; }

protected:

    NodeContainerType m_nodes = { nullptr, nullptr, nullptr, nullptr };
    Vector m_U0 = Vector(24);
};

#endif // !ASDShellQ4Transformation_h

