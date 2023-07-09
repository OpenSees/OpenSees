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
                                                                        
// $Revision: 1.7 $
// $Date: 2007-02-02 01:35:22 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/fourNodeQuad/NineNodeMixedQuad.h,v $

// Ed "C++" Love
//
// Constant Presssure/Volume Four Node Quadrilateral
// Plane Strain (NOT PLANE STRESS)
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

class NineNodeMixedQuad : public Element {

  public :
    
    //null constructor
    NineNodeMixedQuad( ) ;
  
    //full constructor
    NineNodeMixedQuad( int tag, 
		       int node1,
		       int node2,
		       int node3,
		       int node4,
		       int node5,
		       int node6,
		       int node7,
		       int node8,
		       int node9,
		       NDMaterial &theMaterial ) ;

    //destructor 
    ~NineNodeMixedQuad( ) ;

    const char *getClassType(void) const {return "NineNodeMixedQuad";};

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
    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);

    int getResponse(int responseID, Information &eleInformation);

    int sendSelf (int commitTag, Channel &theChannel);
    int recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker 
		  &theBroker);

    //plotting 
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
  
  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;
    
    
    //quadrature data
    static double root06;
    static double sg[3] ;
    static double wg[3] ;

    //node information
    ID connectedExternalNodes ;  //nine node numbers
    Node *nodePointers[9] ;      //pointers to nine nodes
					
    //material information
    NDMaterial *materialPointers[9] ; //pointers to nine materials
					  
    //nodal coordinates, two coordinates for each of nine nodes
    static double xl[][9] ; 
    
    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    const Matrix& computeBbar( int node, 
			       const double natCoor[2], 
			       const double shp[3][9], 
			       double shpBar[3][9][3] ) ;

    //shape function routine for four node quads
    void shape2dNine( double coor[2], 
		  const double x[2][9], 
		  double shp[3][9], 
		  double &xsj ) ;

    //nodal coordinates
    void computeBasis( ) ;

    //1d quadratic shape functions
    double shape1d( int code, int node, double xi ) ;

    Vector *load;
    Matrix *Ki;
} ; 
