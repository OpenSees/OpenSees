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
                                                                        
// $Revision: 1.12 $
// $Date: 2007-06-27 00:24:34 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/brick/BbarBrick.h,v $

// Ed "C++" Love
//
// Eight node BbarBrick element 
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

class BbarBrick : public Element {

  public :
    
    //null constructor
    BbarBrick( ) ;
  
    //full constructor
    BbarBrick( int tag, 
			int node1,
			int node2,
		        int node3,
			int node4,
			int node5,
			int node6,
			int node7,
			int node8,
			NDMaterial &theMaterial, 
			double b1 = 0.0, double b2 = 0.0, double b3 = 0.0 ) ;

    //destructor 
    virtual ~BbarBrick( ) ;

    const char *getClassType(void) const {return "BbarBrick";};

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

  private : 

    //static data
    static Matrix stiff ;
    static Vector resid ;
    static Matrix mass ;
    static Matrix damping ;

    //quadrature data
    static const double root3 ;
    static const double one_over_root3 ;    
    static const double sg[2] ;
    static const double wg[8] ;

  
    //node information
    ID connectedExternalNodes ;  //four node numbers
    Node *nodePointers[8] ;      //pointers to four nodes


    //material information
    NDMaterial *materialPointers[8] ; //pointers to eight materials
					  
    //local nodal coordinates, three coordinates for each of four nodes
    //    static double xl[3][8] ; 
    static double xl[][8] ; 

	double b[3];		// Body forces
	
	double appliedB[3]; // Body forces applied with load pattern, C.McGann, U.Washington
	int applyLoad;      // flag for body force in load, C.McGann, U.Washington

    //inertia terms
    void formInertiaTerms( int tangFlag ) ;

    //form residual and tangent					  
    void formResidAndTangent( int tang_flag ) ;

    //compute coordinate system
    void computeBasis( ) ;

    //compute Bbar matrix
    const Matrix& computeBbar( int node, 
			       const double shp[4][8], 
			       const double shpBar[4][8] ) ;
  
    //Matrix transpose
    Matrix transpose( int dim1, int dim2, const Matrix &M ) ;

    Vector *load;
    Matrix *Ki;
} ; 







