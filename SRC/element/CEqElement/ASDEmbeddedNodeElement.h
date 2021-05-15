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
// $Date: 2021/04/28 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
//
// Notes:
//
// 

#ifndef ASDEmbeddedNodeElement_h
#define ASDEmbeddedNodeElement_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class ASDEmbeddedNodeElement : public Element
{

public:

    // life cycle
    ASDEmbeddedNodeElement();
    ASDEmbeddedNodeElement(int tag, int cNode, int rNode1, int rNode2, int rNode3, bool rot_flag, double K);
    ASDEmbeddedNodeElement(int tag, int cNode, int rNode1, int rNode2, int rNode3, int rNode4, bool rot_flag, double K);
    virtual ~ASDEmbeddedNodeElement();

    // domain
    const char* getClassType(void) const;
    void setDomain(Domain* theDomain);

    // print
    void Print(OPS_Stream& s, int flag);

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID& getExternalNodes();
    Node** getNodePtrs();
    int getNumDOF();

    // methods dealing with committed state and update
    int update();
    int revertToLastCommit();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
    const Vector& getGlobalDisplacements() const;
    const Matrix& TRI_2D_U();
    const Matrix& TRI_2D_UR();
    const Matrix& TRI_3D_U();
    const Matrix& TRI_3D_UR();
    const Matrix& TET_3D_U();
    const Matrix& TET_3D_UR();

private:

    // the nodal ids, the first one is the constrained node,
    // the other 3 (triangle in 2D or shell triangle in 3D) or 4 (tetrahedron in 3D)
    // are the retained nodes
    ID m_node_ids;
    // the ndoes
    std::vector<Node*> m_nodes;
    // store the number of dimensions (2 or 3 are allowed)
    int m_ndm = 0;
    // total number of dofs
    int m_num_dofs = 0;
    // user input to constrained, if necessary, the rotation of the constrained node
    // if the constrained node has rotational DOFs
    bool m_rot_c_flag = false;
    // true if the constrained node has rotational DOFs and the user flag is true
    bool m_rot_c = false;
    // a vector containing the local id mapping for assembling
    // into the element matrix and vectors
    ID m_mapping;
    // stiffness penalty value to impose the constraint
    double m_K = 1.0e18;
    // initial displacements
    Vector m_U0;
    bool m_U0_computed = false;

};

#endif // ASDEmbeddedNodeElement_h
