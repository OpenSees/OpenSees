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
// $Date: 2014/09/30 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/ShellNLDKGQ/ShellNLDKGQ.h,v $

//Written: Lisha Wang, Xinzheng Lu, Linlin Xie, Song Cen & Quan Gu
//
//  four node flat shell element with membrane and drill DOF
//  considering geometric nonlinear, form nonlinear shell element
//  using updated Lagrangian formula
// Ref: Plate Bending Part - DKQ, thin plate element
//      Membrane Part - GQ12, a membrane element with drilling DOF
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

class ShellNLDKGQ : public Element {

 public:
  
  //null constructor
  ShellNLDKGQ( );
  
  //full constructor
  ShellNLDKGQ( int tag, 
	          int node1,
	          int node2,
			  int node3,
			  int node4,
			  SectionForceDeformation &theMaterial ) ;
  
  //destructor 
  virtual ~ShellNLDKGQ( ) ;

  //set domain because frank is a dumb ass 
  void setDomain( Domain *theDomain ) ;
  
  //get the number of external nodes
  int getNumExternalNodes( ) const ;
    
    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs( );

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
    const Matrix &getInitialStiff( );
    const Matrix &getMass( );

    // methods for applying loads
    void zeroLoad( void );	
    int addLoad( ElementalLoad *theLoad, double loadFactor );
    int addInertiaLoadToUnbalance( const Vector &accel );

    //get residual
    const Vector &getResistingForce( ) ;
    
    //get residual with inertia terms
    const Vector &getResistingForceIncInertia( ) ;

    // public methods for element output
    int sendSelf ( int commitTag, Channel &theChannel );
    int recvSelf ( int commitTag, Channel &theChannel, FEM_ObjectBroker 
		           &theBroker );


    Response* setResponse( const char **argv, int argc, OPS_Stream &output );
    int getResponse( int responseID, Information &eleInfo );
      
    //plotting 
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

	//last resid
	// Vector CstrainGauss,TstrainGauss;
	Vector CstrainGauss,TstrainGauss;  //modify for geometric nonlinearity

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static double sg[4] ;
    static double tg[4] ;
    static double wg[4] ;

    //node information
    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[4] ;      //pointers to four nodes

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
//start Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
    void updateBasis( ) ;
//end Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
        
    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute Bdrill matrix
    //double* computeBdrill( int node, const double shp[3][4] ) ;

    //assemble a B matrix 
    const Matrix& assembleB( const Matrix &Bmembrane,
			                 const Matrix &Bbend, 
			                 const Matrix &Bshear ) ;
  
    //compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, const double shp[3][4],
		                         const double shpDrill[4][4]) ;
  
    //compute Bbend matrix
    const Matrix& computeBbend( int node, const double shpBend[6][12] ) ;

	//add for geometric nonlinearity
	const Matrix& computeBG(int node, const double shpBend[6][12]);
	const Vector& computeNLdstrain(const Matrix &BG,const Vector &dispIncLocalBend);
  
    //Matrix transpose
    //Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
		          const double x[2][4], 
		          double shp[3][4], 
		          double &xsj ,double sx[2][2]) ;

	//shape function routine for membrane
	void shapeDrill(double ss, double tt, 
		          const double x[2][4],
				  double sx[2][2], double shpDrill[4][4]);
	//shape function routine for bending
	void shapeBend(double ss, double tt, const double x[2][4],
		          double sx[2][2], double shpBend[6][12]);

    // vector for applying loads
    Vector *load;
    Matrix *Ki;
} ; 
