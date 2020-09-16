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
// $Date: 2020-09-16 16:45:00 $
// $Source: /OpenSees/SRC/material/section/ElasticMembranePlateSection2.h,v $

// Developed by Pearl Ranchal
// with support from Degenkolb Engineers
// Elastic shell section with decoupled in-plane and out-of-plane stiffness parameters

// Adapted from 'ElasticMembranePlateSection.h' by:
// Ed "C++" Love

#ifndef ElasticMembranePlateSection2_h
#define ElasticMembranePlateSection2_h

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>

#include <SectionForceDeformation.h>


class ElasticMembranePlateSection2 : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    ElasticMembranePlateSection2( ) ;

    //full constructor
    ElasticMembranePlateSection2(   int tag, 
                                    double Em,
                                    double Eb,
                                    double nu,
                                    double h = 1.0,
			                        double rho = 0.0);

    //destructor
    ~ElasticMembranePlateSection2( ) ;

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

    double Em ; // elastic modulus for membrane
    double Eb;  // elastic modulus for bending
    double nu ; // poisson ratio
    double h  ; // MembranePlate thickness
    double rhoH ; //mass per unit 2D area

    static const double five6 ; // =5/6 = shear correction factor

    Vector strain ;

    static Vector stress ;

    static Matrix tangent ;

    static ID array ;  

} ; //end of ElasticMembranePlateSection2 declarations





#endif
