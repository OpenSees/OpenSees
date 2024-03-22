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
// Implementation of a linear coordinate transformation 3-node shells
//

#ifndef ASDShellT3Transformation_h
#define ASDShellT3Transformation_h

#include <ASDShellT3LocalCoordinateSystem.h>
#include <Node.h>
#include <Domain.h>

/** \brief ASDShellT3Transformation
*
* This class represents a basic (linear) coordinate transformation that can be used
* by any element whose geometry is a Triangle 3 in 3D space, with 6 D.O.F.s per node.
* Its main aim is to:
* 1) Create the local coordinate system
* 2) Transform the incoming global displacements in local coordinate system
* 3) Transform the outgoing matrices and vectors in global coordinate system
*/
class ASDShellT3Transformation
{

public:

    typedef ASDVector3<double> Vector3Type;

    typedef ASDQuaternion<double> QuaternionType;

    typedef Vector VectorType;

    typedef Matrix MatrixType;

    typedef std::array<Node*, 3> NodeContainerType;

public:

    ASDShellT3Transformation()
    {
    }

    virtual ~ASDShellT3Transformation()
    {
    }

private:

    ASDShellT3Transformation(const ASDShellT3Transformation& other) = delete;

    ASDShellT3Transformation& operator = (const ASDShellT3Transformation& other) = delete;

public:

    virtual ASDShellT3Transformation* create()const
    {
        return new ASDShellT3Transformation();
    }

    virtual bool isLinear() const
    {
        return true;
    }

    virtual void revertToStart()
    {
    }

    virtual void setDomain(Domain* domain, const ID& node_ids, bool initialized)
    {
        // if domain is null
        if (domain == nullptr) {
            for (size_t i = 0; i < 3; i++) {
                m_nodes[i] = nullptr;
            }
            return;
        }

        // get nodes and save initial displacements and rotations
        for (size_t i = 0; i < 3; i++) {
            m_nodes[i] = domain->getNode(node_ids(i));
            if (m_nodes[i] == nullptr) {
                opserr << "ASDShellT3Transformation::setDomain - no node " << node_ids(i)
                    << " exists in the model\n";
                exit(-1);
            }
            if (!initialized) {
                const Vector& iU = m_nodes[i]->getTrialDisp();
                if (iU.Size() != 6) {
                    opserr << "ASDShellT3Transformation::setDomain - node " << node_ids(i)
                        << " has " << iU.Size() << " DOFs, while 6 are expected\n";
                    exit(-1);
                }
                size_t index = i * 6;
                for (size_t j = 0; j < 6; j++)
                    m_U0(index + j) = iU(j);
            }
        }
    }

    virtual void revertToLastCommit()
    {
    }

    virtual void commit()
    {
    }

    virtual void update(const VectorType& globalDisplacements)
    {
    }

    virtual ASDShellT3LocalCoordinateSystem createReferenceCoordinateSystem()const
    {
        // the reference coordinate system in the underformed configuration
        // using the default alignment to the first column of the jacobian at center
        return ASDShellT3LocalCoordinateSystem(
            Vector3Type(m_nodes[0]->getCrds()),
            Vector3Type(m_nodes[1]->getCrds()),
            Vector3Type(m_nodes[2]->getCrds())
        );
    }

    virtual ASDShellT3LocalCoordinateSystem createLocalCoordinateSystem(const VectorType& globalDisplacements)const
    {
        // same as reference
        return createReferenceCoordinateSystem();
    }

    virtual void computeGlobalDisplacements(VectorType& globalDisplacements) const
    {
        for (int i = 0; i < 3; i++) {
            int index = i * 6;
            const VectorType& iU = m_nodes[i]->getTrialDisp();
            for (int j = 0; j < 6; j++) {
                globalDisplacements(index + j) = iU(j) - m_U0(index + j);
            }
        }
    }

    virtual const MatrixType& computeTransformationMatrix(const ASDShellT3LocalCoordinateSystem& LCS) const
    {
        static MatrixType T(18, 18);
        LCS.ComputeTotalRotationMatrix(T);
        return T;
    }

    virtual void calculateLocalDisplacements(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        VectorType& localDisplacements)
    {
        const MatrixType& R = computeTransformationMatrix(LCS);
        localDisplacements.addMatrixVector(0.0, R, globalDisplacements, 1.0);
    }

    virtual void transformToGlobal(
        const ASDShellT3LocalCoordinateSystem& LCS,
        const VectorType& globalDisplacements,
        const VectorType& localDisplacements,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static MatrixType RT_LHS(18, 18);
        static VectorType RHScopy(18);
        const MatrixType& R = computeTransformationMatrix(LCS);
        RHScopy = RHS;
        RHS.addMatrixTransposeVector(0.0, R, RHScopy, 1.0);
        if (LHSrequired) {
            RT_LHS.addMatrixTransposeProduct(0.0, R, LHS, 1.0);
            LHS.addMatrixProduct(0.0, RT_LHS, R, 1.0);
        }
    }

    virtual void transformToGlobal(
        const ASDShellT3LocalCoordinateSystem& LCS,
        MatrixType& LHS,
        VectorType& RHS,
        bool LHSrequired)
    {
        static VectorType dummy;
        transformToGlobal(LCS, dummy, dummy, LHS, RHS, LHSrequired);
    }

    virtual int internalDataSize() const
    {
        // just the size of the initial displacements
        return 18;
    }

    virtual void saveInternalData(VectorType& v, int pos) const
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3Transformation - failed to save internal data: vector too small\n";
            exit(-1);
        }
        for (int i = 0; i < 18; i++)
            v(pos++) = m_U0(i);
    }

    virtual void restoreInternalData(const VectorType& v, int pos)
    {
        if ((v.Size() - pos) < internalDataSize()) {
            opserr << "ASDShellT3Transformation - failed to restore internal data: vector too small\n";
            exit(-1);
        }
        for (int i = 0; i < 18; i++)
            m_U0(i) = v(pos++);
    }

public:

    inline const NodeContainerType& getNodes()const { return m_nodes; }
    inline NodeContainerType& getNodes() { return m_nodes; }

protected:

    NodeContainerType m_nodes = { {nullptr, nullptr, nullptr} };
    Vector m_U0 = Vector(18);
};

#endif // !ASDShellT3Transformation_h

