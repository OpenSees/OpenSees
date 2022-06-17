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

#ifndef IGAKLShell_h
#define IGAKLShell_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <NDMaterial.h>

// IGA stuff
#include <IGASurfacePatch.h>

// Van-Paepegem stuff
#include <Information.h>
#include <Parameter.h>

class IGAKLShell : public Element {

public:

    //null constructor
    IGAKLShell( );

    //full constructor
    IGAKLShell( int tag,
                IGASurfacePatch *myPatch_,
                const ID& nodes,
                int quadorder_,
                const Vector& xiE_,
                const Vector& etaE_,
                const ID& matTags) ;

    //destructor
    virtual ~IGAKLShell( ) ;

    void setDomain( Domain *theDomain ) ;

    const char *getClassType(void) const {return "IGAKLShell";};

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
    Response* emulateSectionSetResponse(const char **argv, int argc,
                                      OPS_Stream &output, int gaussPointNum, double xi, double eta);

    //plotting
    int displaySelf(Renderer &, int mode, float fact, const char **displayModes = 0, int numModes = 0);

    Matrix transpose( int dim1, int dim2, const Matrix &M );

    // Added by Felipe Elgueta
    void shellGeo(Matrix G, Matrix H, Vector& G3, double& dA, Vector& N, Matrix& Gab, Vector& Bv, Matrix& T_Gcon_E, Matrix& T_E_G, Matrix& T_G_E); // Get geometric quantities
    bool pointInElement(double xi, double eta) const;
private :
    //static data
    static Matrix* stiff ;
    static Vector* resid ;
    static Matrix* mass ;
    static Matrix* damping ;
    Vector* load; // For follower loads




    int quadorder, ngauss, nLayers;

    IGASurfacePatch *myPatch;

    Vector xiE;
    Vector etaE;

    int noFuncs;

    Matrix *quadPoint;
    Vector *quadWeight;

    ID connectedExternalNodes ;  // node numbers

    NDMaterial ***materialPointers ; //pointers to  materials
    Node **nodePointers ;      //pointers to  nodes


    int applyLoad;
    double appliedB[3];     // Body forces applied with load
    double pressure;        // Surface pressure.. positive is traction

    //form residual and tangent                   
    void formResidAndTangent( int tang_flag ) ;

    void formTangentNguyen();

    void formStrainsLinear();

    
    int setParameter(const char **argv, int argc, Parameter &param);
    int updateParameter(int parameterID, Information &info);


} ;



#endif
