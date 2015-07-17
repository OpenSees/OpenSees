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

///////////////////////////////////////////////////////////////////////////////
// Description: This file contains the class declaration for BBarBrickUP,    //
// an 8-node cubic element for solid-fluid fully coupled analysis.           //
// This implementation is a simplified u-p formulation of Biot theory        //
// (u - solid displacement, p - fluid pressure). Each element node has three //
// DOFs for u and 1 DOF for p.                                               //
// Constant volume/pressure integration (BBar method) is used for integration//
// of the volumetric component of solid phase and the fulid phase.           //
//                                                                           //
// Written by Zhaohui Yang	(September 2009)                                 //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// $Revision: 1.1 $
// $Date: 2009-10-07 20:02:23 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/UP-ucsd/BBarBrickUP.h,v $

// by Zhaohui Yang (Modified based on Ed "C++" Love's Brick element)
//
// Eight node BBarBrickUP element
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class BBarBrickUP : public Element {

  public :

    //null constructor
    BBarBrickUP( ) ;

    //full constructor
    BBarBrickUP( int tag,
			int node1,
			int node2,
		    int node3,
			int node4,
			int node5,
			int node6,
			int node7,
			int node8,
			NDMaterial &theMaterial,double bulk, double rhof,
			double perm1, double perm2, double perm3,
		   double b1 = 0.0, double b2 = 0.0, double b3 = 0.0) ;

    //destructor
    virtual ~BBarBrickUP( ) ;

    const char *getClassType(void) const {return "BBarBrickUP";};

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

    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);

    //plotting
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);


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
    static const double root3 ;
    static const double one_over_root3 ;
    static const double sg[2] ;
    static const double wg[8] ;


    //node information
    ID connectedExternalNodes ;  //eight node numbers
    Node *nodePointers[8] ;      //pointers to eight nodes


    //material information
    NDMaterial *materialPointers[8] ; //pointers to eight materials

    //local nodal coordinates, three coordinates for each of 8 nodes
    static double xl[4][8] ;
    double b[3];		// Body forces

	double appliedB[3]; // Body forces applied by load pattern, C.McGann, U.Washington
	int applyLoad;      // flag for body forces applied by load, C.McGann, U.Washington
	
    double rho;			// Fluid mass per unit volume
    double kc;   // combined bulk modulus
    double perm[3];  // permeability

    // [0,1,2=derivative wrt x,y,z;3=shape func][node][Gauss point]
    static double Shape[4][8][8]; // Stores shape functions and derivatives (overwritten)

    // [x,y,z][node]
    static double shpBar[3][8]; // Stores averaged shap functions (overwritten)

    // [row][col][node][Gauss point]
    static double BBar[6][3][8][8];  // Stores strain-displacement matrix (overwritten)

    // [col][node][Gauss point]  Note: there is only one row in Bp matrix
    static double BBarp[3][8][8]; // Stores strain-displacement matrix for fluid phase (overwritten)

    static double dvol[8];  // Stores detJacobian (overwritten)

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //damping terms
    void formDampingTerms( int tangFlag ) ;

    //form residual and tangent
    void formResidAndTangent( int tang_flag ) ;

	// Mixture mass density at integration point i
	double mixtureRho(int ipt);

    //compute coordinate system
    void computeBasis( ) ;

    //compute Bbar matrix
    void computeBBar() ;

    //compute B matrix
    const Matrix& computeB( int node, int Guass ) ;

    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    Vector *load;
    Matrix *Ki;
} ;
