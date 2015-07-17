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

// by Jinchi Lu and Zhaohui Yang (May 2004)
//
// 20-8 Noded brick element: TwentyEightNodeBrickUP
//

#ifndef _TwentyEightNodeBrickUP_H
#define _TwentyEightNodeBrickUP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class TwentyEightNodeBrickUP : public Element {

public :

    //null constructor
    TwentyEightNodeBrickUP( ) ;

    //full constructor
    TwentyEightNodeBrickUP( int tag,
		int node1,
		int node2,
		int node3,
		int node4,
		int node5,
		int node6,
		int node7,
		int node8,
		int node9,
		int node10,
		int node11,
		int node12,
		int node13,
		int node14,
		int node15,
		int node16,
		int node17,
		int node18,
		int node19,
		int node20,
		NDMaterial &theMaterial,double bulk, double rhof,
		double perm1, double perm2, double perm3,
		double b1 = 0.0, double b2 = 0.0, double b3 = 0.0) ;

    //destructor
    virtual ~TwentyEightNodeBrickUP( ) ;
    const char *getClassType(void) const {return "TwentyEightNodeBrickUP";};
    //set domain
    void setDomain( Domain *theDomain ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;

    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs(void);

    //return number of dofs
    int getNumDOF( ) ;

    //commit state
    int commitState( ) ;

    //revert to last commit
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;

	int update(void);

    //return stiffness matrix
    const Matrix &getTangentStiff( ) ;
    const Matrix &getInitialStiff( ) ;
    const Matrix &getDamp(void);
    const Matrix &getMass( ) ;

    void zeroLoad( ) ;
    int addLoad(ElementalLoad *theLoad, double loadFactor);
    int addInertiaLoadToUnbalance(const Vector &accel);

    //get residual
    const Vector &getResistingForce( ) ;

    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker
		&theBroker);

    Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);




    //plotting
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    // RWB; PyLiq1 & TzLiq1 need to see the excess pore pressure and initial stresses.
    friend class PyLiq1;
    friend class TzLiq1;

private :
    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damp ;

    //quadrature data
    static const int nintu;
    static const int nintp;
    static const int nenu;
    static const int nenp;

    //node information
    ID connectedExternalNodes ;  //eight node numbers
    Node *nodePointers[20] ;      //pointers to eight nodes

    //material information
    NDMaterial **materialPointers; // pointer to the ND material objects

    //local nodal coordinates, three coordinates for each of twenty nodes
    //    static double xl[3][20] ;
    static double xl[3][20] ;
    double b[3];		// Body forces
    double appliedB[3]; // Body forces applied by load pattern, C.McGann, U.Washington
    int applyLoad;      // flag for body forces applied by load, C.McGann, U.Washington
    double rho;			// Fluid mass per unit volume
    double kc;   // combined bulk modulus
    double perm[3];  // permeability

    static double shgu[4][20][27];	// Stores shape functions and derivatives (overwritten)
    static double shgp[4][8][8];	// Stores shape functions and derivatives (overwritten)
    static double shgq[4][20][8];	// Stores shape functions and derivatives (overwritten)
    static double shlu[4][20][27];	// Stores shape functions and derivatives
    static double shlp[4][8][8];	// Stores shape functions and derivatives
    static double shlq[4][20][8];	// Stores shape functions and derivatives
    static double wu[27];		// Stores quadrature weights
    static double wp[8];		// Stores quadrature weights
    static double dvolu[27];  // Stores detJacobian (overwritten)
    static double dvolp[8];  // Stores detJacobian (overwritten)
    static double dvolq[8];  // Stores detJacobian (overwritten)
    //inertia terms
    void formInertiaTerms( int tangFlag ) ;
    //damping terms
    void formDampingTerms( int tangFlag ) ;
    // Mixture mass density at integration point i
    double mixtureRho(int ipt);

    //compute coordinate system
    void computeBasis( ) ;

    Vector *load;
    Matrix *Ki;

	// compute local shape functions
	void compuLocalShapeFunction();
	void Jacobian3d(int gaussPoint, double& xsj, int mode);
	const Matrix&  getStiff( int flag );
} ;

#endif
