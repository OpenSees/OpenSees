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

// $Revision $
// $Date$
// $URL$

// Written: Minjie Zhu (zhum@oregonstate.edu)
//
// Description: MINI element which mass matrix can be lumped
//

#ifndef MINI_H
#define MINI_H

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <vector>

class MINI : public Element
{
public:

    MINI();
    MINI(int tag, int nd1, int nd2, int nd3,
	 double r, double m, double b1, double b2, 
	 double thk=1.0, double ka=2.15e9);
    MINI(int tag, int nd1, int nd2, int nd3, int nd4,
	 double r, double m, double b1, double b2, double b3,
	 double thk=1.0, double ka=2.15e9);
    ~MINI();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes() const;
    const ID &getExternalNodes();
    Node **getNodePtrs();
    int getNumDOF();

    // public methods to set the state of the element
    int revertToLastCommit();
    int update();
    int commitState();

    // public methods to obtain stiffness, mass, damping and residual information
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getDamp();
    const Matrix &getMass();

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector &accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector &getResistingForce();
    const Vector &getResistingForceIncInertia();

    // MovableObject
    const char *getClassType() const;
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel,
		 FEM_ObjectBroker &theBroker);

    // DomainComponent
    void setDomain(Domain *theDomain);

    // TaggedObject
    void Print(OPS_Stream &s, int flag =0);
    int displaySelf(Renderer &, int mode, float fact,
		    const char **displayModes=0, int numModes=0);

private:

    double det(const Matrix& m);

private:

    // nodes
    ID ntags;
    std::vector<Node*> nodes;

    // data: rho,mu,thk,kappa,bx,by,bz
    Vector data;

    // dof for each node
    ID dofs;

    // matrix
    static Matrix mat;
    static Vector vec;

    // Jacobian cofactor matrix
    // [a1,b1,c1,d1]
    // [a2,b2,c2,d2]
    // [a3,b3,c3,d3]
    // [a4,b4,c4,d4]
    Matrix Jco; 

};


#endif
