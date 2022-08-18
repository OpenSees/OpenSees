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
// $Date: 2006-08-03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/MembranePlateFiberSection.h,v $

// Ed "C++" Love
//
// Generic Plate Section with membrane
//

#ifndef MembranePlateFiberSection_h
#define MembranePlateFiberSection_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>


class MembranePlateFiberSection : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    MembranePlateFiberSection( ) ;

    //full constructor
    MembranePlateFiberSection(   int    tag, 
                                 double thickness, 
                                 NDMaterial &Afiber,
				 int integrType = 0) ;


    const char *getClassType(void) const {return "MembranePlateFiberSection";};

    //destructor
    virtual ~MembranePlateFiberSection( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;

    //mass per unit area
    double getRho() ;

    //send back order of strain in vector form
    int getOrder( ) const ;

    //send back order of strain in vector form
    const ID& getType( ) ;

    //swap history variables
    int commitState( ) ; 

    //revert to last saved state
    int revertToLastCommit( ) ;

    //revert to start
    int revertToStart( ) ;

    //get the strain and integrate plasticity equations
    int setTrialSectionDeformation( const Vector &strain_from_element ) ;

    //send back the strain
    const Vector& getSectionDeformation( ) ;

    //send back the stress 
    const Vector& getStressResultant( ) ;

    //send back the tangent 
    const Matrix& getSectionTangent( ) ;

    //send back the initial tangent 
    const Matrix& getInitialTangent( ) {return this->getSectionTangent();}

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

    Response *setResponse(const char **argv, int argc, 
			  OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

  private :

    enum {numFibers = 5};

    //quadrature data
    static const double sgLobatto[numFibers] ;
    static const double wgLobatto[numFibers] ;
    static const double sgGauss[numFibers] ;
    static const double wgGauss[numFibers] ;  

    double h ; //plate thicknes

  int integrationType; // 0 = Lobatto, 1 = Gauss

    NDMaterial *theFibers[numFibers] ;  //pointers to five materials (fibers)

    static const double root56 ; // =sqrt(5/6) 

    Vector strainResultant ;

    static Vector stressResultant ;

    static Matrix tangent ;

    static ID array ;  

} ; //end of MembranePlateFiberSection declarations





#endif
