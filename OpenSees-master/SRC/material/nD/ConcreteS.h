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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-05-25 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/ConcreteS.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Simple Concrete Model, plane stress
// 


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>



class ConcreteS: public NDMaterial{
  public : 
    ConcreteS( ) ;
    ConcreteS(int tag, double rE, double rnu, double rfc, double rft, double rEs) ;

    virtual ~ConcreteS( ) ;

    //make a clone of this material
    NDMaterial *getCopy( ) ;
    NDMaterial *getCopy( const char *type ) ;

    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const char *getType( ) const ;

    void setInitials() ;

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

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

private :
    Matrix tangent, eTangent ;
    Vector strain0, strain, stress0, stress, stressd;
    double E, nu, fc, ft, Es, EmEp1, Ep, beta, ftmin, cStrain0, cStrain;

} ;
