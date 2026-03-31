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
// $Date: 2023/12/18 22:51:21 $

// Implementation: Massimo Petracca (ASDEA)
// Based on feap fortran element elmt03 by Ushnish Basu

#ifndef FSIFluidElement2D_h
#define FSIFluidElement2D_h

#include <Element.h>
#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <vector>

class TimeSeries;

class FSIFluidElement2D : public Element
{

public:

    // life cycle
    FSIFluidElement2D();
    FSIFluidElement2D(
        int tag,
        int node1,
        int node2,
        int node3,
        int node4,
        double c,
        double thickness);
    virtual ~FSIFluidElement2D();

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
    void zeroLoad();
    int addLoad(ElementalLoad* theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector& accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector& getResistingForce();
    const Vector& getResistingForceIncInertia();

    // public methods for element output
    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

private:
    const Vector& getU() const;
    const Vector& getV() const;
    const Vector& getA() const;
    const Matrix& q4gauss() const;
    const Vector& q4N(double xi, double eta) const;
    const Matrix& q4dN(double xi, double eta) const;
    const Matrix& q4Jac(const Matrix& dN, double& detJ) const;

private:

    // nodal ids
    ID m_node_ids = ID(4);
    // node pointers
    std::vector<Node*> m_nodes = std::vector<Node*>(4, nullptr);
    // acousting wave velocity
    double m_c = 0.0;
    // thickness
    double m_thickness = 1.0;
    // vectors for applying load (allocated only if necessary)
    Vector* m_load = nullptr;

};

#endif // FSIFluidElement2D_h
