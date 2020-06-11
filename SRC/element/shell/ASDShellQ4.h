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
// A 4-node general shell element based on the AGQI formulation 
// for the in-plane behavior, 
// and the MITC4 formulation for the out-of-plane behavior.
// 
// It supports both linear and corotational kinematics. Warped geometries
// can be modelled since this element is not assumed flat.
//
//  References:
// - Dvorkin, Bathe, "A continuum mechanics based four node shell
//   element for general nonlinear analysis",
//   Eng.Comput., vol. 1, 77 - 88, 1984
// - Bathe, Dvorkin, "Short communication A four-node plate bending element
//   based on Mindlin / Reissner plate theory and a Mixed Interpolation",
//   International Journal for Numerical Methods in Eng.,
//   vol. 21, 367 - 383, 1985
//

#ifndef ASDShellQ4_h
#define ASDShellQ4_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>

class SectionForceDeformation;
class ASDShellQ4Transformation;
class ASDShellQ4LocalCoordinateSystem;

class ASDShellQ4 : public Element
{

public:

    // life cycle
    ASDShellQ4();
    ASDShellQ4(
        int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        SectionForceDeformation* section,
        bool corotational = false);
    virtual ~ASDShellQ4();

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
    int commitState();
    int revertToLastCommit();
    int revertToStart();
    int update();

    // methods to return the current linearized stiffness,
    // damping and mass matrices
    const Matrix& getTangentStiff();
    const Matrix& getInitialStiff();
    const Matrix& getMass();

    // methods for applying loads
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

    int setParameter(const char** argv, int argc, Parameter& param);

private:

    // internal method to compute everything using switches...
    int calculateAll(Matrix& LHS, Vector& RHS, int options);

    void AGQIinitialize();
    void AGQIupdate(const Vector& UL);
    void AGQIbeginGaussLoop(const ASDShellQ4LocalCoordinateSystem& reference_cs);

private:

    // cross sections
    SectionForceDeformation* m_sections[4] = { nullptr, nullptr, nullptr, nullptr };

    // nodal ids
    ID m_node_ids = ID(4);

    // coordinate transformation
    ASDShellQ4Transformation* m_transformation = nullptr;

    // vectors for applying load (allocated only if necessary)
    Vector* m_load = nullptr;

    // drilling strain for the indipendent rotation field (Hughes-Brezzi)
    double m_drill_strain[4] = { 0.0, 0.0, 0.0, 0.0 };
    double m_drill_stiffness = 0.0;

    // section orientation with respect to the local coordinate system
    double m_angle = 0.0;

    // members for non-linear treatement of AGQI internal DOFs:
    // it has 24 displacement DOFs
    // and 4 internal DOFs for membrane enhancement
    Vector m_Q = Vector(4);
    Vector m_Q_converged = Vector(4);
    Vector m_U = Vector(24);
    Vector m_U_converged = Vector(24);
    Vector m_Q_residual = Vector(4);
    Matrix m_KQQ_inv = Matrix(4, 4);
    Matrix m_KQU = Matrix(4, 24); // L = G'*C*B
    Matrix m_KUQ = Matrix(24, 4); // L^T = B'*C'*G
};

#endif // ASDShellQ4_h
