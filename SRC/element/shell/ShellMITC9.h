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
                                                                        
// Written: Leopoldo Tesser, Diego Talledo
// 9-node lagrandian shell element with membrane and drill
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

class ShellMITC9 : public Element
{
  public:

    //null constructor
    ShellMITC9();
    //full constructor
    ShellMITC9(int tag, 
	     int node1,
	     int node2,
	     int node3,
	     int node4,
	     int node5,
	     int node6,
	     int node7,
	     int node8,
	     int node9,
	     SectionForceDeformation &theMaterial ) ;
    //destructor 
    virtual ~ShellMITC9( ) ;

    //get the number of external nodes
    int getNumExternalNodes( ) const ;
    //return connected external nodes
    const ID &getExternalNodes( ) ;
    Node **getNodePtrs();
    //return number of dofs
    int getNumDOF( ) ;

    //set domain 
    void setDomain( Domain *theDomain ) ;
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


    Response* setResponse(const char **argv, int argc, OPS_Stream &output);
    int getResponse(int responseID, Information &eleInfo);
      
    //plotting 
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double root3_over_root5 ;
    static double sg[9] ;
    static double tg[9] ;
    static double wg[9] ;

    //node information
    ID connectedExternalNodes ;  //nine node numbers
    //pointers to nine nodes
    Node *nodePointers[9] ;

    //drilling stiffness
    double Ktt ;

    //material information: pointers to four materials
    SectionForceDeformation *materialPointers[9] ;
					  
    //local nodal coordinates, two coordinates for each of nine nodes
    //static double xl[][9] ; 
    double xl[2][9] ;

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
    //void  computeJacobian( double L1, double L2,const double x[2][9], 
    //                       Matrix &JJ,Matrix &JJinv ) ;

    //compute Bdrill matrix
    double* computeBdrill( int node, const double shp[3][9] ) ;

    //assemble a B matrix 
    const Matrix& assembleB( const Matrix &Bmembrane,
			     const Matrix &Bbend, 
			     const Matrix &Bshear ) ;
    
    //compute Bmembrane matrix
    const Matrix& computeBmembrane( int node, const double shp[3][9] ) ;
    
    //compute Bbend matrix
    const Matrix& computeBbend( int node, const double shp[3][9] ) ;
    
    //compute standard Bshear matrix
    const Matrix&  computeBshear( int node, const double shp[3][9] ) ;
    
    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    //shape function routine for four node quads
    void shape2d( double ss, double tt, 
		  const double x[2][9], 
		  double shp[3][9], 
		  double &xsj ) ;

    // vector for applying loads
    Vector *load;
    Matrix *Ki;
} ; 







