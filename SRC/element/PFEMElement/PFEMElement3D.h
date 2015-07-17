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
// $Date: 2013/01/29 13:41:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/PFEMElement/PFEMElement3D.h,v $
                                                                        
// Written: Minjie Zhu (zhum@engr.orst.edu)
// Created: Jan 2013
// Revised: --------
//
// Description: This file contains the class definition for PFEMElement3D.

#ifndef PFEMElement3D_h
#define PFEMElement3D_h

#include <Matrix.h>
#include <Vector.h>
#include <Element.h>

class Pressure_Constraint;

class PFEMElement3D : public Element
{
public:
    PFEMElement3D(int tag, int nd1, int nd2, int nd3, int nd4,
                  double r, double m, double b1, double b2, double b3);
    
    ~PFEMElement3D();

    // methods dealing with nodes and number of external dof
    int getNumExternalNodes(void) const;
    const ID &getExternalNodes(void);
    Node **getNodePtrs(void);
    int getNumDOF(void);

    // public methods to set the state of the element    
    int revertToLastCommit(void);
    //int revertToStart(void);   
    int update(void);
    int commitState(void);    

    // public methods to obtain stiffness, mass, damping and residual information  
    const Matrix &getTangentStiff(void);
    const Matrix &getInitialStiff(void);    
    const Matrix &getDamp();
    const Matrix &getDampWithK();
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
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

    // TaggedObject
    void Print(OPS_Stream &s, int flag =0);

protected:
    
private:

    ID ntags; // Tags of nodes
    Node* nodes[8]; // pointers of nodes
    Pressure_Constraint* thePCs[4];
    double rho;  // density
    double mu;   // viscocity
    double bx;    // body force
    double by;    // body force
    double bz;    // body force
    double dNdx[4], dNdy[4], dNdz[4];
    double J;
    ID numDOFs;

    static Matrix K;
    static Vector P;

private:
    double det33(const Matrix& A);

};

#endif


