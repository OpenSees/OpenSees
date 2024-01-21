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
// $Date: 2006/08/03 23:49:46 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/MembranePlateFiberSectionThermal.h,v $

// Ed "C++" Love
//
// Generic Plate Section with membrane
//


#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

#include <Vector.h>
#include <Matrix.h>
#include <ID.h>
#include <NDMaterial.h>

#include <SectionForceDeformation.h>


class MembranePlateFiberSectionThermal : public SectionForceDeformation{

//-------------------Declarations-------------------------------

  public : 

    //null constructor
    MembranePlateFiberSectionThermal( ) ;

    //full constructor
    MembranePlateFiberSectionThermal(   int    tag, 
                                 double thickness, 
                                 NDMaterial &Afiber ) ;


    const char *getClassType(void) const {return "MembranePlateFiberSectionThermal";};

    //destructor
    virtual ~MembranePlateFiberSectionThermal( ) ;

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

	const Vector &getTemperatureStress(const Vector&); //J.Jiang add to get Ft=EA*Elongation//

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

	  Response *setResponse(const char **argv, int argc, OPS_Stream &s);
    int getResponse(int responseID, Information &info);

    // parameters
    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int parameterID, Information& info);

  private :

    //quadrature data
    static const double sg[5] ;
    static const double wg[5] ;

    double h ; //plate thickness

    NDMaterial *theFibers[5] ;  //pointers to five materials (fibers)

    static const double root56 ; // =sqrt(5/6) 

    Vector strainResultant ;

    static Vector stressResultant ;

    static Matrix tangent ;

    static ID array ; 

	Vector  sT;  //  Pointer to sTData
	double  ThermalElongation[5]; // Temperature dependent elasticity modulus
	int countnGauss;
	double ThermalGradientShink;

} ; //end of MembranePlateFiberSectionThermal declarations





