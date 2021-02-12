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
                                                                        
// $Revision: 1.9 $
// $Date: 2006-08-03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticMembranePlateSection.h,v $

// Ed "C++" Love
//
//  Elastic Plate Section with membrane
//

// Out-of-Plane stiffness modifier added by Pearl Ranchal
// Supported by Degenkolb Engineers

#ifndef ElasticMembranePlateSection_h
#define ElasticMembranePlateSection_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>


class ElasticMembranePlateSection : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    ElasticMembranePlateSection( ) ;

    //full constructor
    ElasticMembranePlateSection(int     tag, 
                                double  E,
                                double  nu,
                                double  h       = 1.0,
			                    double  rho     = 0.0,
                                double  Ep_mod  = 1.0);

    //destructor
    ~ElasticMembranePlateSection( ) ;

    //make a clone of this material
    SectionForceDeformation *getCopy( ) ;

    const char *getClassType(void) const {return "ElasticMembranePlate";};

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
    const Matrix& getInitialTangent( ) ;

    //print out data
    void Print( OPS_Stream &s, int flag ) ;

    //density per unit area
    double getRho() ;

    int sendSelf(int commitTag, Channel &theChannel);
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);


  private :

    double Em;  // elastic modulus for membrane (in-plane) action
    double Ep;  // elastic modulus for plate (out-of-plane) action
    double nu ; // poisson ratio
    double h  ; // MembranePlate thickness
    double rhoH ; //mass per unit 2D area

    static const double five6 ; // =5/6 = shear correction factor

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;

    static ID array ;  

} ; //end of ElasticMembranePlateSection declarations





#endif
