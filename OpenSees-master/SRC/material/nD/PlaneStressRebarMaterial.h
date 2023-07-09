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


// written: fmk derived from code in PlateRebarMaterial from
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <UniaxialMaterial.h>
#include <NDMaterial.h>

class PlaneStressRebarMaterial: public NDMaterial{
  public : 
    PlaneStressRebarMaterial();
    PlaneStressRebarMaterial(int tag,
		       UniaxialMaterial &uniMat,
		       double ang );

    virtual ~PlaneStressRebarMaterial( );

    //make a clone of this material
    NDMaterial *getCopy( );
    NDMaterial *getCopy( const char *type );

    int getOrder( ) const ;
    const char *getType( ) const ;

    int commitState( ); 
    int revertToLastCommit( );
    int revertToStart( );

    int setTrialStrain( const Vector &strainFromElement );

    const Vector& getStrain( );
    const Vector& getStress( );
    const Matrix& getTangent( );
    const Matrix& getInitialTangent( );

    double getRho( );

    void Print( OPS_Stream &s, int flag );
    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private :
    UniaxialMaterial *theMat ;
    double angle, c, s;

    Vector strain ;
    static Vector stress ;
    static Matrix tangent ;

} ;





