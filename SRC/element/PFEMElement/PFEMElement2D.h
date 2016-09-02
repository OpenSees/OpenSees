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
                                                                        
// $Revision: 1.00 $
// $Date: 2012/01/11 13:48:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement2D.h,v $
                                                                        
// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2012
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement2D.

#ifndef PFEMElement2D_h
#define PFEMElement2D_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>

class Pressure_Constraint;

class PFEMElement2D : public Element
{
public:
    PFEMElement2D();
    PFEMElement2D(int tag, int nd1, int nd2, int nd3,
                  double r, double m, double b1, double b2, double thk=1.0,
		  double ka=-1, bool lumped=true);
    
    ~PFEMElement2D();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes(void) const {return ntags.Size();}
    const ID &getExternalNodes(void) {return ntags;}
    Node **getNodePtrs(void) {return nodes;}
    int getNumDOF(void) {return ndf;}

    // public methods to set the state of the element    
    int revertToLastCommit(void) {return 0;}
    int revertToStart(void) {return Element::revertToStart();}
    int update(void);
    int commitState(void) {return Element::commitState();}

    // public methods to obtain stiffness, mass, damping and residual information    
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getDamp();
    const Matrix &getMass(void);    

    // methods for applying loads
    int addInertiaLoadToUnbalance(const Vector &accel);

    // methods for obtaining resisting force (force includes elemental loads)
    const Vector &getResistingForce(void);
    const Vector &getResistingForceIncInertia(void);   

    // MovableObject
    const char *getClassType(void) const;
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    // DomainComponent
    void setDomain(Domain *theDomain); 

    // TaggedObject
    void Print(OPS_Stream &s, int flag =0);
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

private:

    void inverse(Matrix& mat) const;
    
private:

    ID ntags; // Tags of nodes
    Node* nodes[6]; // pointers of nodes
    Pressure_Constraint* thePCs[3];
    double rho;  // density
    double mu;   // viscocity
    double b1, b2; // body force
    double thickness;
    double kappa;
    bool lumped;
    int ndf;
    int vxdof[3], vydof[3], pdof[3];
    double M, Mp;
    Matrix Km, S;
    Vector Gx, Gy, F, Fp;

    static Matrix K;
    static Vector P;
};

#endif


