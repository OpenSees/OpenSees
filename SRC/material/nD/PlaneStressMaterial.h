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
                                                                        
// $Revision: 1.4 $
// $Date: 2003-02-14 23:01:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressMaterial.h,v $

// Ed "C++" Love
//
// Generic Plane Stress Material
//


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>



class PlaneStressMaterial: public NDMaterial{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    PlaneStressMaterial( ) ;

    //full constructor
    PlaneStressMaterial(   int    tag, 
                           NDMaterial &the3DMaterial ) ;


    //destructor
    virtual ~PlaneStressMaterial( ) ;

    //make a clone of this material
    NDMaterial *getCopy( ) ;
    NDMaterial *getCopy( const char *type ) ;

    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const char *getType( ) const ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain 
    int setTrialStrain( const Vector &strainFromElement ) ;

    //send back the strain
    const Vector& getStrain( ) ;

    //send back the stress 
    const Vector& getStress( ) ;

    //send back the tangent 
    const Matrix& getTangent( ) ;
    const Matrix& getInitialTangent( ) ;

    //density
    double getRho( ) ;

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    int setParameter(const char **argv, int argc, Parameter &param);

    const Vector& getStressSensitivity(int gradIndex,
				       bool conditional);

  private :

    //out of plane strains .. trial and committed
    double Tstrain22 ;
    double Tgamma02 ;
    double Tgamma12 ; 
    double Cstrain22 ;
    double Cgamma02 ;
    double Cgamma12 ; 

    NDMaterial *theMaterial ;  //pointer to three dimensional material

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;
} ; //end of PlaneStressMaterial declarations





