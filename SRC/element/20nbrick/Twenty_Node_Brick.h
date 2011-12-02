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
// 20NodeBrick element
//
#ifndef TWENTY_NODE_BRICK_H
#define TWENTY_NODE_BRICK_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class Twenty_Node_Brick : public Element {

 public:
    //null constructor
    Twenty_Node_Brick( ) ;

    //full constructor
    Twenty_Node_Brick( int tag,
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
		NDMaterial &theMaterial,
		double b1 = 0.0, double b2 = 0.0, double b3 = 0.0) ;

    //destructor
    virtual ~Twenty_Node_Brick( ) ;

    const char *getClassType(void) const {return "Twenty_Node_Brick";};

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

    Response *setResponse(const char **argv, int argc, Information &eleInformation, OPS_Stream &s);
    int getResponse(int responseID, Information &eleInformation);

    //plotting
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

private :

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damp ;

    //quadrature data
    static const int nintu;
    static const int nenu;


    //node information
    ID connectedExternalNodes ;  //eight node numbers
    Node *nodePointers[20] ;      //pointers to eight nodes


    //material information
    NDMaterial **materialPointers; // pointer to the ND material objects

    //local nodal coordinates, three coordinates for each of twenty nodes
    //    static double xl[3][20] ;
    static double xl[3][20] ;
    double b[3];		// Body forces

    static double shgu[4][20][27];	// Stores shape functions and derivatives (overwritten)
    static double shlu[4][20][27];	// Stores shape functions and derivatives
    static double wu[27];		// Stores quadrature weights
    static double dvolu[27];  // Stores detJacobian (overwritten)

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
