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
// A 4-node general shell element based on 
// the AGQ6-I formulation for the membrane part, 
// and the MITC4 formulation for the shear-deformable bending part.
// 
// It supports both linear and corotational kinematics. Warped geometries
// can be modelled since this element is not assumed flat.
//
//  References:
// 
// 1 Chen, Xiao-Ming, et al. "Membrane elements insensitive to distortion 
//   using the quadrilateral area coordinate method." 
//   Computers & Structures 82.1 (2004): 35-54.
//   http://www.paper.edu.cn/scholar/showpdf/MUT2ANwINTT0Ax5h
//
// 2 Dvorkin, Eduardo N., and Klaus-Jurgen Bathe. 
//   "A continuum mechanics based four-node shell element for 
//   general non-linear analysis." Engineering computations (1984).
//   https://www.researchgate.net/profile/Eduardo_Dvorkin/publication/235313212_A_Continuum_mechanics_based_four-node_shell_element_for_general_nonlinear_analysis/links/00b7d52611d8813ffe000000.pdf
//
// 3 Bathe, Klaus-Jurgen, and Eduardo N. Dvorkin. 
//   "A four-node plate bending element based on Mindlin/Reissner plate theory 
//   and a mixed interpolation." 
//   International Journal for Numerical Methods in Engineering 21.2 (1985): 367-383.
//   http://www.simytec.com/docs/Short_communicaion_%20four_node_plate.pdf
//
// 4 Hughes, Thomas JR, and F. Brezzi. 
//   "On drilling degrees of freedom." 
//   Computer methods in applied mechanics and engineering 72.1 (1989): 105-121.
//   https://www.sciencedirect.com/science/article/pii/0045782589901242
//
// Notes:
//
// - The AGQ6-I formulation greatly enhances the membrane behavior of the
//   standard 4-node iso-parametric element. However it does not pass the constant-strain
//   patch test as most of the enhanced elements based on incompatible-modes.
//   In this implementation the enhancing BQ matrix is corrected to make the element
//   pass the constant strain patch test. At each gauss point the BQ = BQ - 1/V*BQ_mean
//
// - The drilling DOF is treated using the Hughes-Brezzi formulation. The drilling
//   stiffness is taken equal to the initial (elastic) in-plane section shear modulus.
//   Taking the real shear modulus (and not a scaled version of it, as suggested in many
//   articles) is mandatory to obtain a correct response for warped shells, i.e. when
//   the drilling rotation affects the bending rotation. However using the full shear modulus
//   deteriorates the element membrane behavior, making the element stiffer. For example,
//   when using softening material models, this drilling stiffness makes the element lock
//   failing in representing crack propagation correctly.
//   Here we solved this problem using a reduced integration for the drilling part. However
//   a pure reduced integration for the drilling DOFs, leads to 3 spurious zero-energy modes.
//   Therefore we also include a fully integrated contribution of the drilling DOFs
//   scaling the drilling stiffness by 1.0e-4, so just to suppress these spurious zero energy modes.
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
