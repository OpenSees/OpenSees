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
// $Date: 2012-06-08 21:11:45 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlaneStressUserMaterial.h,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Plane Stress User Defined Material
// 


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h> 
#include <NDMaterial.h>



class PlaneStressUserMaterial: public NDMaterial{
  public : 
    PlaneStressUserMaterial( ) ;
    PlaneStressUserMaterial(int tag, int istatevs, int iprops, double *props) ;

    virtual ~PlaneStressUserMaterial( ) ;

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

    //cracking output - added by V.K. Papanikolaou [AUTh] - start
    const Vector& getCracking();
    //cracking output - added by V.K. Papanikolaou [AUTh] - end

    //set/getResponse - added by V.K. Papanikolaou [AUTh] - start
    Response* setResponse(const char** argv, int argc, OPS_Stream& s);
    int getResponse(int responseID, Information& matInformation);
    //set/getResponse - added by V.K. Papanikolaou [AUTh] - end

private :
    Vector strain0, strain, stress0, stress;
    Matrix tangent, eTangent;

    Vector *vprops, *statev0, *statev;
    double strain0data[3], straindata[3], dstraindata[3];
    double stress0data[3], stressdata[3], tangentdata[9];
    double *props, *statevdata;
    int nstatevs, nprops;

} ;
