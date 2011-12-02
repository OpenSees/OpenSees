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
                                                                        
// $Revision: 1.9 $
// $Date: 2003-02-25 23:33:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/ShellMITC4.h,v $

// Ed "C++" Love
//
// B-bar four node shell element with membrane and drill
//

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <ID.h> 
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <R3vectors.h>

class ShellMITC4 : public Element {

 public:
  
  //null constructor
  ShellMITC4();
  
  //full constructor
  ShellMITC4(int tag, 
	     int node1,
	     int node2,
	     int node3,
	     int node4,
	     SectionForceDeformation &theMaterial ) ;
  
  //destructor 
  virtual ~ShellMITC4( ) ;

  //set domain because frank is a dumb ass 
  void setDomain( Domain *theDomain ) ;
  
  //get the number of external nodes
  int getNumExternalNodes( ) const ;
    
    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();

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
    const Matrix &getInitialStiff();
    const Matrix &getMass();

    // methods for applying loads
    void zeroLoad(void);	
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
      
    //plotting 
    int displaySelf(Renderer &theViewer, int displayMode, float fact);

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static double sg[4] ;
    static double tg[4] ;
    static double wg[4] ;

    //B-bar data
    static Matrix  **GammaB1pointer ;  
    static Matrix  **GammaD1pointer ;
    static Matrix  **GammaA2pointer ;
    static Matrix  **GammaC2pointer ;

    //Bhat data
    static Matrix **Bhat ;

    //node information
    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[4] ;      //pointers to four nodes

    //drilling stiffness
    double Ktt ;

    //material information
    SectionForceDeformation *materialPointers[4] ; //pointers to four materials
					  
    //local nodal coordinates, two coordinates for each of four nodes
    //static double xl[][4] ; 
    double xl[2][4] ; 

    //shell basis vectors
    double g1[3] ;
    double g2[3] ;
    double g3[3] ;

    //compute local coordinates and basis
    void computeBasis( ) ;
        
    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute Jacobian matrix and inverse at point {L1,L2}
    void  computeJacobian( double L1, double L2, 
			   const double x[2][4], 
                           Matrix &JJ, 
                           Matrix &JJinv ) ;

    //compute Bdrill matrix
    double* computeBdrill( int node, const double shp[3][4] ) ;

    //assemble a B matrix 
    const Matrix& assembleB( const Matrix &Bmembrane,
			     const Matrix &Bbend, 
			     const Matrix &Bshear ) ;
  
    //compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, const double shp[3][4] ) ;
  
    //compute Bbend matrix
    const Matrix& computeBbend( int node, const double shp[3][4] ) ;
  
    //compute standard Bshear matrix
    const Matrix&  computeBshear( int node, const double shp[3][4] ) ;

    //compute Bbar shear matrix
    const Matrix&  computeBbarShear( int node, double L1, double L2,
				     const Matrix& Jinv ) ;
  
			     
    //compute the gamma's
    void  computeGamma( const double xl[2][4], const Matrix &J ) ;

    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
		  const double x[2][4], 
		  double shp[3][4], 
		  double &xsj ) ;

    // vector for applying loads
    Vector *load;
    Matrix *Ki;
} ; 







