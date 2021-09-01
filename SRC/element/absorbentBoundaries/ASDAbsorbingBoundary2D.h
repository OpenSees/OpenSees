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
// $Date: 2021/08/30 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)

#ifndef ASDAbsorbingBoundary2D_h
#define ASDAbsorbingBoundary2D_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class ASDAbsorbingBoundary2D : public Element
{
public:
    enum StageType {
        Stage_StaticConstraint = 0,
        Stage_Absorbing = 1
    };
    enum BoundaryType {
        Boundary_Horizontal = 0,
        Boundary_Vertical = 1
    };

public:

    // life cycle
    ASDAbsorbingBoundary2D();
    ASDAbsorbingBoundary2D(
        int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        double G,
        double rho,
        double thickness);
    virtual ~ASDAbsorbingBoundary2D();

    // class type
    const char* getClassType(void) const;

    // domain
    void setDomain(Domain* theDomain);

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int revertToLastCommit();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff(void);
    const Matrix& getInitialStiff(void);
    const Matrix& getDamp(void);
    const Matrix& getMass(void);

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

private:

    // compute a consistent penalty value
    double penaltyFactor();
    // fills the penalty stiffness matrix in stage = 0
    void addKPenalty(Matrix& K);
    // fills the penalty resisting forces
    void addRPenalty(Vector& R);

private:

    // nodal ids
    ID m_node_ids = ID(4);
    // node pointers
    std::vector<Node*> m_nodes;
    // shear modulus
    double m_G = 0.0;
    // mass density
    double m_rho = 0.0;
    // thickness
    double m_thickness = 1.0;
    // stage
    StageType m_stage = Stage_StaticConstraint;
    // boundary type
    BoundaryType m_boundary = Boundary_Horizontal;
    // direction sign
    double m_n = 1.0;
    // length in X
    double m_lx = 0.0;
    // length in Y
    double m_ly = 0.0;
    // total number of dofs
    int m_num_dofs = 0;
    // a vector containing the local id mapping for assembling
    // into the element matrix and vectors
    ID m_mapping = ID(8);
    // initial displacement
    Vector m_U0 = Vector(8);
    // initial velocity
    Vector m_V0 = Vector(8);

};

#endif // ASDAbsorbingBoundary2D_h
