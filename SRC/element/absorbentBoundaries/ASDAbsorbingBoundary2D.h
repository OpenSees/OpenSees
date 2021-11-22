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

class TimeSeries;

class ASDAbsorbingBoundary2D : public Element
{
public:
    enum StageType {
        Stage_StaticConstraint = 0,
        Stage_Absorbing = 1
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
        double v,
        double rho,
        double thickness,
        int btype,
        TimeSeries* actionx,
        TimeSeries* actiony);
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

    // computation of reactions
    int addResistingForceToNodalReaction(int flag);

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

    // response
    Response* setResponse(const char** argv, int argc, OPS_Stream& output);
    int getResponse(int responseID, Information& eleInfo);

private:

    // methods to get displacement/velocity
    void addDisplacement(Vector& U);
    const Vector& getDisplacement();
    const Vector& getVelocity();
    const Vector& getAcceleration();
    // get element sizes
    void getElementSizes(double& lx, double& ly, double& nx);
    // update stage
    void updateStage();
    // compute a consistent penalty value
    void penaltyFactor(double& sp, double& mp);
    // fills the penalty stiffness matrix in stage = 0
    void addKPenaltyStage0(Matrix& K);
    // fills the penalty resisting forces in stage = 0
    void addRPenaltyStage0(Vector& R);
    // fills the penalty stiffness matrix in stage = 1
    void addKPenaltyStage1(Matrix& K);
    // fills the penalty resisting forces in stage = 1
    void addRPenaltyStage1(Vector& R);
    // restore reactions
    void addRReactions(Vector& R);
    // fills the mass matrix of the free-field
    void addMff(Matrix& M, double scale = 1.0);
    // fills the force vector due to mass*acc
    void addRMff(Vector& R);
    // fills the stiffness matrix of the free-field
    void addKff(Matrix& K, double scale = 1.0);
    // fills the resisting forces of the free-field
    void addRff(Vector& R);
    // fills the stiffness matrix of the free-field forces tranfered to the soil domain
    void addKffToSoil(Matrix& K);
    // fills the forces transfered from the free-field to the soil domain
    void addRffToSoil(Vector& R);
    // compute damping parameters
    void getDampParam(double& alpha, double& beta);
    // fills the damping matrix of the free-field
    void addCff(Matrix& C);
    // fills the damping forces of the free-field
    void addRCff(Vector& R);
    // computes the dashpot coefficients
    void getLKcoeff(double& ap, double& as);
    // fills the damping matrix of the L-K dashpots
    void addClk(Matrix& C);
    // fills the damping forces of the L-K dashpots
    void addRlk(Vector& R);
    // fills the input action forces
    void addBaseActions(Vector& R);

private:

    // nodal ids
    ID m_node_ids = ID(4);
    // node pointers
    std::vector<Node*> m_nodes = std::vector<Node*>(4, nullptr);
    // shear modulus
    double m_G = 0.0;
    // poisson's ratio
    double m_v = 0.0;
    // mass density
    double m_rho = 0.0;
    // thickness
    double m_thickness = 1.0;
    // stage
    StageType m_stage = Stage_StaticConstraint;
    // boundary type
    int m_boundary = 0;
    // total number of dofs
    int m_num_dofs = 0;
    // a vector containing the local DOF mapping for assembling
    // into the element matrix and vectors
    ID m_dof_map = ID(8);
    // a vector containing the local node id mapping
    std::vector<std::size_t> m_node_map = std::vector<std::size_t>(4, 0);
    // initial displacement
    Vector m_U0;
    // reaction forces saved at the end of stage 0
    Vector m_R0;
    // flag for computation of reactions
    bool m_is_computing_reactions = false;
    // initialized flag
    bool m_initialized = false;
    // time series for base actions
    TimeSeries* m_tsx = nullptr;
    TimeSeries* m_tsy = nullptr;

};

#endif // ASDAbsorbingBoundary2D_h
