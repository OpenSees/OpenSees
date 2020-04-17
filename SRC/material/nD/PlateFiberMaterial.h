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
                                                                        
// $Revision: 1.6 $
// $Date: 2006-08-04 18:18:38 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateFiberMaterial.h,v $

// Ed "C++" Love
//
// Generic Plate Fiber Material
//

#ifndef PlateFiberMaterial_h
#define PlateFiberMaterial_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>

class PlateFiberMaterial: public NDMaterial{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    PlateFiberMaterial( ) ;

    //full constructor
    PlateFiberMaterial(   int    tag, 
                           NDMaterial &the3DMaterial ) ;


    //destructor
    virtual ~PlateFiberMaterial( ) ;

    virtual const char *getClassType(void) const {return "PlateFiberMAterial";};

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
    const Matrix &getInitialTangent(void);

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

    //out of plane strain
    double Tstrain22 ;
    double Cstrain22 ;

    NDMaterial *theMaterial ;  //pointer to three dimensional material

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;
} ; //end of PlateFiberMaterial declarations





#endif
