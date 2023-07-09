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
// Description: P2P1 element, quadratic velocity, linear pressure
//

#ifndef TaylorHood2D_H
#define TaylorHood2D_H

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <Pressure_Constraint.h>
#include <vector>

class TaylorHood2D : public Element
{
public:

    TaylorHood2D();
    TaylorHood2D(int tag, int nd1, int nd2, int nd3,
		 double r, double m, double b1, double b2, 
		 double thk=1.0, double ka=2.15e9);
    ~TaylorHood2D();

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

    // get current coordinates
    int getX(Matrix& x);
    
    // Jacobian matrix
    void Jacobian(Matrix& J, const Matrix& X,
		  const Matrix& dNdL);

    // shape functions
    void vshape(Vector& N, double l1, double l2, double l3);
    void pshape(Vector& N, double l1, double l2, double l3);
    void derishape(Matrix& dN, double l1, double l2, double l3);

    // matrices
    double Mab(int a, int b, double detJ, const Vector& N);
    double K11(int a, int b, double detJ,
	       const Matrix& dNdx);
    double K12(int a, int b, double detJ,
	       const Matrix& dNdx);
    double K21(int a, int b, double detJ,
	       const Matrix& dNdx);
    double K22(int a, int b, double detJ,
	       const Matrix& dNdx);
    double Gab(int a, int b, double detJ, const Matrix& dNdx,
	       int i, const Vector& Np);
    double Fa(int a, double detJ, const Vector& N, double bx);
    double Mpab(int a, int b, double detJ, const Vector& Np);

private:

    // nodes
    ID ntags;
    std::vector<Node*> nodes;
    std::vector<Pressure_Constraint*> pcs;

    // data
    double rho, mu, b1, b2, thk, kappa;

    // number of dofs
    int ndf;
    ID vxdof, vydof, pdof;

    // matrix
    static Matrix mat;
    static Vector vec;

};


#endif
