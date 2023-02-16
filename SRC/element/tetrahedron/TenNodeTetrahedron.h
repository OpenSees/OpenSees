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

// ============================================================================
// 2022 By Jose Abell and Jose Larenas @ Universidad de los Andes, Chile
// www.joseabell.com | https://github.com/jaabell | jaabell@miuandes.cl
// ============================================================================
// Implements a standard 10-node tetrahedron element.
//
// This element has 4 Gauss points of integration.
//
//
// ============================================================================


#ifndef TenNodeTetrahedron_H
#define TenNodeTetrahedron_H


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

class TenNodeTetrahedron : public Element {

public :

    //null constructor
    TenNodeTetrahedron();

    //full constructor
    TenNodeTetrahedron(int tag,
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
                       NDMaterial &theMaterial,
                       double b1 = 0.0, double b2 = 0.0, double b3 = 0.0);

    //destructor
    virtual ~TenNodeTetrahedron( ) ;

    const char *getClassType(void) const {return "TenNodeTetrahedron";};

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

    // update
    int update(void);

    //print out element data
    void Print( OPS_Stream &s, int flag ) ;

    //return stiffness matrix
    const Matrix &getTangentStiff();
    const Matrix &getInitialStiff();
    const Matrix &getMass();

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
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);
    void onActivate();
    void onDeactivate();

private :

  //Number of Gauss-points
  enum {NumGaussPoints=4};
  enum {NumNodes=10};
  enum {NumDOFsPerNode=3};
  enum {NumStressComponents=6};
  enum {NumDOFsTotal=NumNodes*NumDOFsPerNode};
  
    // Routine to compute shape functions and their derivatives. These get stored as follows. 
    // for node n:
    //   shp[0][n] --> dN_n / d x, 
    //   shp[1][n] --> dN_n / d y
    //   shp[2][n] --> dN_n / d z
    //   shp[3][n] --> N_n, shape function n value at the z values 
    void shp3d( 
        const double zeta[4],  // Tetrahedral coordinates  (input)
        double &xsj,         // Jacobian determinant (output)
        double shp[4][NumNodes], // Shape function and derivatives values at the tetrahedral coordinates (output) 
        const double xl[3][NumNodes]   ); // Node coordinates (input)

    //
    // private attributes
    //

    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[NumNodes] ;      //pointers to eight nodes

    //material information
    NDMaterial *materialPointers[NumGaussPoints]; //pointers to eight materials

    double b[3];        // Body forces
    double appliedB[3];     // Body forces applied with load
    int applyLoad;

    Vector *load;
    Matrix *Ki;

    //
    // static attributes
    //

    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;
  static Matrix B;
  
    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;
    static const double alpha ;
    static const double beta ;
    static const double sg[4] ;
    static const double wg[1] ;

    //local nodal coordinates, three coordinates for each of four nodes
    static double xl[3][NumNodes] ;

    //
    // private methods
    //

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent
    void formResidAndTangent( int tang_flag ) ;

    //compute coordinate system
    void computeBasis( ) ;

    //compute B matrix
    const Matrix& computeB( int node, const double shp[4][NumNodes] ) ;

    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;
    Vector initDisp[NumNodes];

    int do_update;
} ;

#endif

