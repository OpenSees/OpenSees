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
                                                                        
// $Revision: 1.8 $
// $Date: 2007-05-03 23:03:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/DoubleMembranePlateFiberSection.cpp,v $

// Ed "C++" Love
//
//  Generic Plate Section with membrane
//
// MHS - Modified for "back to back" sections separated by a specified distance
//


#include <DoubleMembranePlateFiberSection.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

void * OPS_ADD_RUNTIME_VPV(OPS_DoubleMembranePlateFiberSection)
{
    int numdata = OPS_GetNumRemainingInputArgs();
    if (numdata < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section DoublePlateFiber tag? matTag? h? d?" << endln;
	return 0;
    }
    
    int idata[2];
    numdata = 2;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
	opserr << "WARNING: invalid tags\n";
	return 0;
    }

    double h, d;
    numdata = 1;
    if (OPS_GetDoubleInput(&numdata, &h) < 0) {
	opserr << "WARNING: invalid h\n";
	return 0;
    }
    if (OPS_GetDoubleInput(&numdata, &d) < 0) {
	opserr << "WARNING: invalid d\n";
	return 0;
    }    

    NDMaterial *theMaterial = OPS_getNDMaterial(idata[1]);
    if (theMaterial == 0) {
	opserr << "WARNING nD material does not exist\n";
	opserr << "nD material: " << idata[1]; 
	opserr << "\nPlateFiber section: " << idata[0] << endln;
	return 0;
    }

    return new DoubleMembranePlateFiberSection(idata[0], h, d, *theMaterial);
}

//parameters
const double DoubleMembranePlateFiberSection::root56 = sqrt(5.0/6.0) ; //shear correction

//static vector and matrices
Vector  DoubleMembranePlateFiberSection::stressResultant(8) ;
Matrix  DoubleMembranePlateFiberSection::tangent(8,8) ;
ID      DoubleMembranePlateFiberSection::array(8) ;


const double  DoubleMembranePlateFiberSection::sg[] = { -1, 
						  -0.65465367, 
					           0, 
					           0.65465367, 
					           1 } ;
 
const double  DoubleMembranePlateFiberSection::wg[] = { 0.1, 
					          0.5444444444, 
						  0.7111111111, 
						  0.5444444444, 
						  0.1  };

/*      from Ham-O
        case 5:
         xi(0,0) = -1.;
         xi(1,0) = -0.65465367;
         xi(2,0) =  0.;
         xi(3,0) =  0.65465367;
         xi(4,0) =  1.;
      
         w(0) =  0.1;
         w(1) =  0.5444444444;
         w(2) =  0.7111111111;
         w(3) =  0.5444444444;
         w(4) =  0.1;
      break;
*/


//null constructor
DoubleMembranePlateFiberSection::DoubleMembranePlateFiberSection( ) : 
SectionForceDeformation( 0, SEC_TAG_DoubleMembranePlateFiberSection ), 
strainResultant(8) 
{ 
  for ( int i = 0; i < numFibers; i++ ) {
      theFibers[i] = 0 ;
      theFibers[i+numFibers] = 0 ;
  }
}



//full constructor
DoubleMembranePlateFiberSection::DoubleMembranePlateFiberSection(    
				   int tag, 
                                   double thickness,
				   double distance,
                                   NDMaterial &Afiber ) :
SectionForceDeformation( tag, SEC_TAG_DoubleMembranePlateFiberSection ),
strainResultant(8)
{
  this->h  = thickness ;
  this->d = distance;

  for (int i = 0; i < numFibers; i++ ) {
      theFibers[i] = Afiber.getCopy( "PlateFiber" ) ;
      theFibers[i+numFibers] = Afiber.getCopy( "PlateFiber" ) ;
  }

}



//destructor
DoubleMembranePlateFiberSection::~DoubleMembranePlateFiberSection( ) 
{ 
  for (int i = 0; i < numFibers; i++ ) {
     delete theFibers[i] ; 
     delete theFibers[i+numFibers] ; 
  }
} 



//make a clone of this material
SectionForceDeformation  *DoubleMembranePlateFiberSection::getCopy( ) 
{
  DoubleMembranePlateFiberSection *clone ;   //new instance of this class

  clone = new DoubleMembranePlateFiberSection( this->getTag(), 
					       this->h, this->d,
					       *theFibers[0] ) ; //make the copy
  return clone ;
}



//send back order of strainResultant in vector form
int DoubleMembranePlateFiberSection::getOrder( ) const
{
  return 8 ;
}


//send back order of strainResultant in vector form
const ID& DoubleMembranePlateFiberSection::getType( ) 
{
    static bool initialized = false;
    if (!initialized) {
        array(0) = SECTION_RESPONSE_FXX;
        array(1) = SECTION_RESPONSE_FYY;
        array(2) = SECTION_RESPONSE_FXY;
        array(3) = SECTION_RESPONSE_MXX;
        array(4) = SECTION_RESPONSE_MYY;
        array(5) = SECTION_RESPONSE_MXY;
        array(6) = SECTION_RESPONSE_VXZ;
        array(7) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return array;
}



//swap history variables
int DoubleMembranePlateFiberSection::commitState( ) 
{
  int success = 0 ;

  for (int i = 0; i < numFibers; i++ ) {
    success += theFibers[i]->commitState( ) ;
    success += theFibers[i+numFibers]->commitState( ) ;
  }

  return success ;
}



//revert to last saved state
int DoubleMembranePlateFiberSection::revertToLastCommit( )
{
  int success = 0 ;

  for (int i = 0; i < numFibers; i++ ) {
    success += theFibers[i]->revertToLastCommit( ) ;
    success += theFibers[i+numFibers]->revertToLastCommit( ) ;
  }

  return success ;
}

//revert to start
int DoubleMembranePlateFiberSection::revertToStart( )
{
  int success = 0 ;

  for (int i = 0; i < numFibers; i++ ) {
    success += theFibers[i]->revertToStart( ) ;
    success += theFibers[i+numFibers]->revertToStart( ) ;
  }

  return success ;
}


//mass per unit area
double
DoubleMembranePlateFiberSection::getRho( )
{

  double weight ;

  double rhoH = 0.0 ;

  for ( int i = 0; i < numFibers; i++ ) {
    
    weight = ( 0.5*h ) * wg[i] ;

    rhoH += ( theFibers[i]->getRho() ) * weight ;
    rhoH += ( theFibers[i+numFibers]->getRho() ) * weight ;

  }

  return rhoH ;

}


//receive the strainResultant 
int DoubleMembranePlateFiberSection ::
setTrialSectionDeformation( const Vector &strainResultant_from_element)
{
  this->strainResultant = strainResultant_from_element ;

  static Vector strain(numFibers) ;

  int success = 0 ;

  int i ;

  double z ;

  strain(3) =  root56*strainResultant(6) ;
  strain(4) =  root56*strainResultant(7) ;
  
  for ( i = 0; i < numFibers; i++ ) {

    z = 0.5*(d+h) + ( 0.5*h ) * sg[i] ;
  
    strain(0) =  strainResultant(0)  - z*strainResultant(3) ;
    
    strain(1) =  strainResultant(1)  - z*strainResultant(4) ;
    
    strain(2) =  strainResultant(2)  - z*strainResultant(5) ;
    
    success += theFibers[i]->setTrialStrain( strain ) ;


    z = -z;

    strain(0) =  strainResultant(0)  - z*strainResultant(3) ;
    
    strain(1) =  strainResultant(1)  - z*strainResultant(4) ;
    
    strain(2) =  strainResultant(2)  - z*strainResultant(5) ;
    
    success += theFibers[i+numFibers]->setTrialStrain( strain ) ;

  } //end for i

  return success ;
}


//send back the strainResultant
const Vector& DoubleMembranePlateFiberSection::getSectionDeformation( )
{
  return this->strainResultant ;
}


//send back the stressResultant 
const Vector&  DoubleMembranePlateFiberSection::getStressResultant( )
{

  static Vector stress(numFibers) ;

  int i ;

  double z, weight ;

  stressResultant.Zero( ) ;

  for ( i = 0; i < numFibers; i++ ) {

      weight = ( 0.5*h ) * wg[i] ;


      z = 0.5*(d+h) + ( 0.5*h ) * sg[i] ;

      stress = theFibers[i]->getStress( ) ;
  
      //membrane
      stressResultant(0)  +=  stress(0)*weight ;

      stressResultant(1)  +=  stress(1)*weight ;

      stressResultant(2)  +=  stress(2)*weight ;

      //bending moments
      stressResultant(3)  +=  ( z*stress(0) ) * weight ;

      stressResultant(4)  +=  ( z*stress(1) ) * weight ;

      stressResultant(5)  +=  ( z*stress(2) ) * weight ;

      //shear
      stressResultant(6)  += stress(3)*weight ;

      stressResultant(7)  += stress(4)*weight ;


      z = -z;

      stress = theFibers[i+numFibers]->getStress( ) ;
  
      //membrane
      stressResultant(0)  +=  stress(0)*weight ;

      stressResultant(1)  +=  stress(1)*weight ;

      stressResultant(2)  +=  stress(2)*weight ;

      //bending moments
      stressResultant(3)  +=  ( z*stress(0) ) * weight ;

      stressResultant(4)  +=  ( z*stress(1) ) * weight ;

      stressResultant(5)  +=  ( z*stress(2) ) * weight ;

      //shear
      stressResultant(6)  += stress(3)*weight ;

      stressResultant(7)  += stress(4)*weight ;
  
  } //end for i

   //modify shear 
   stressResultant(6) *= root56 ;  
   stressResultant(7) *= root56 ;

   return this->stressResultant ;
}


//send back the tangent 
const Matrix&  DoubleMembranePlateFiberSection::getSectionTangent( )
{
  static Matrix dd(5,5) ;

  static Matrix Aeps(5,8) ;

  static Matrix Asig(8,5) ;

  int i ;

  double z, weight ;

  tangent.Zero( ) ;

  for ( i = 0; i < numFibers; i++ ) {

      weight = (0.5*h) * wg[i] ;


      z = 0.5*(d+h) + ( 0.5*h ) * sg[i] ;

/*      //compute Aeps

      Aeps.Zero( ) ;

      Aeps(0,0) = 1.0 ;
      Aeps(0,3) = -z ;

      Aeps(1,1) = 1.0 ;
      Aeps(1,4) = -z ;

      Aeps(2,2) = 1.0 ;
      Aeps(2,5) = -z ;

      Aeps(3,6) = root56 ;
      Aeps(4,7) = root56 ;

      //compute Asig

      Asig.Zero( ) ;

      Asig(0,0) = 1.0 ;
      Asig(3,0) = z ;

      Asig(1,1) = 1.0 ;
      Asig(4,1) = z ;

      Asig(2,2) = 1.0 ;
      Asig(5,2) = z ;

      Asig(6,3) = root56 ;
      Asig(7,4) = root56 ;
*/

      //compute the tangent

      dd = theFibers[i]->getTangent( ) ;

      dd *= weight ;

      //tangent +=  ( Asig * dd * Aeps ) ;   

//from MATLAB : tangent = 
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
 
      //row 1
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
      tangent(0,0) +=  dd(0,0) ;
      tangent(0,1) +=  dd(0,1) ;
      tangent(0,2) +=  dd(0,2) ;      
      tangent(0,3) +=  -z*dd(0,0) ;      
      tangent(0,4) +=  -z*dd(0,1) ;
      tangent(0,5) +=  -z*dd(0,2) ;
      tangent(0,6) +=  root56*dd(0,3) ;
      tangent(0,7) +=  root56*dd(0,4) ;

      //row 2
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
      tangent(1,0) +=  dd(1,0) ;
      tangent(1,1) +=  dd(1,1) ;
      tangent(1,2) +=  dd(1,2) ;      
      tangent(1,3) +=  -z*dd(1,0) ;      
      tangent(1,4) +=  -z*dd(1,1) ;
      tangent(1,5) +=  -z*dd(1,2) ;
      tangent(1,6) +=  root56*dd(1,3) ;
      tangent(1,7) +=  root56*dd(1,4) ;

      //row 3
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
      tangent(2,0) +=  dd(2,0) ;
      tangent(2,1) +=  dd(2,1) ;
      tangent(2,2) +=  dd(2,2) ;      
      tangent(2,3) +=  -z*dd(2,0) ;      
      tangent(2,4) +=  -z*dd(2,1) ;
      tangent(2,5) +=  -z*dd(2,2) ;
      tangent(2,6) +=  root56*dd(2,3) ;
      tangent(2,7) +=  root56*dd(2,4) ;

      //row 4
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
      tangent(3,0) +=  z*dd(0,0) ;
      tangent(3,1) +=  z*dd(0,1) ;
      tangent(3,2) +=  z*dd(0,2) ;      
      tangent(3,3) +=  -z*z*dd(0,0) ;      
      tangent(3,4) +=  -z*z*dd(0,1) ;
      tangent(3,5) +=  -z*z*dd(0,2) ;
      tangent(3,6) +=  z*root56*dd(0,3) ;
      tangent(3,7) +=  z*root56*dd(0,4) ;

      //row 5
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
      tangent(4,0) +=  z*dd(1,0) ;
      tangent(4,1) +=  z*dd(1,1) ;
      tangent(4,2) +=  z*dd(1,2) ;      
      tangent(4,3) +=  -z*z*dd(1,0) ;      
      tangent(4,4) +=  -z*z*dd(1,1) ;
      tangent(4,5) +=  -z*z*dd(1,2) ;
      tangent(4,6) +=  z*root56*dd(1,3) ;
      tangent(4,7) +=  z*root56*dd(1,4) ;

      //row 6
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
      tangent(5,0) +=  z*dd(2,0) ;
      tangent(5,1) +=  z*dd(2,1) ;
      tangent(5,2) +=  z*dd(2,2) ;      
      tangent(5,3) +=  -z*z*dd(2,0) ;      
      tangent(5,4) +=  -z*z*dd(2,1) ;
      tangent(5,5) +=  -z*z*dd(2,2) ;
      tangent(5,6) +=  z*root56*dd(2,3) ;
      tangent(5,7) +=  z*root56*dd(2,4) ;

      //row 7
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
      tangent(6,0) +=  root56*dd(3,0) ;
      tangent(6,1) +=  root56*dd(3,1) ;
      tangent(6,2) +=  root56*dd(3,2) ;      
      tangent(6,3) +=  -root56*z*dd(3,0) ;      
      tangent(6,4) +=  -root56*z*dd(3,1) ;
      tangent(6,5) +=  -root56*z*dd(3,2) ;
      tangent(6,6) +=  root56*root56*dd(3,3) ;
      tangent(6,7) +=  root56*root56*dd(3,4) ;

      //row 8 
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
      tangent(7,0) +=  root56*dd(4,0) ;
      tangent(7,1) +=  root56*dd(4,1) ;
      tangent(7,2) +=  root56*dd(4,2) ;      
      tangent(7,3) +=  -root56*z*dd(4,0) ;      
      tangent(7,4) +=  -root56*z*dd(4,1) ;
      tangent(7,5) +=  -root56*z*dd(4,2) ;
      tangent(7,6) +=  root56*root56*dd(4,3) ;
      tangent(7,7) +=  root56*root56*dd(4,4) ;






      z = -z;

/*      //compute Aeps

      Aeps.Zero( ) ;

      Aeps(0,0) = 1.0 ;
      Aeps(0,3) = -z ;

      Aeps(1,1) = 1.0 ;
      Aeps(1,4) = -z ;

      Aeps(2,2) = 1.0 ;
      Aeps(2,5) = -z ;

      Aeps(3,6) = root56 ;
      Aeps(4,7) = root56 ;

      //compute Asig

      Asig.Zero( ) ;

      Asig(0,0) = 1.0 ;
      Asig(3,0) = z ;

      Asig(1,1) = 1.0 ;
      Asig(4,1) = z ;

      Asig(2,2) = 1.0 ;
      Asig(5,2) = z ;

      Asig(6,3) = root56 ;
      Asig(7,4) = root56 ;
*/

      //compute the tangent

      dd = theFibers[i+numFibers]->getTangent( ) ;

      dd *= weight ;

      //tangent +=  ( Asig * dd * Aeps ) ;   

//from MATLAB : tangent = 
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
 
      //row 1
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
      tangent(0,0) +=  dd(0,0) ;
      tangent(0,1) +=  dd(0,1) ;
      tangent(0,2) +=  dd(0,2) ;      
      tangent(0,3) +=  -z*dd(0,0) ;      
      tangent(0,4) +=  -z*dd(0,1) ;
      tangent(0,5) +=  -z*dd(0,2) ;
      tangent(0,6) +=  root56*dd(0,3) ;
      tangent(0,7) +=  root56*dd(0,4) ;

      //row 2
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
      tangent(1,0) +=  dd(1,0) ;
      tangent(1,1) +=  dd(1,1) ;
      tangent(1,2) +=  dd(1,2) ;      
      tangent(1,3) +=  -z*dd(1,0) ;      
      tangent(1,4) +=  -z*dd(1,1) ;
      tangent(1,5) +=  -z*dd(1,2) ;
      tangent(1,6) +=  root56*dd(1,3) ;
      tangent(1,7) +=  root56*dd(1,4) ;

      //row 3
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
      tangent(2,0) +=  dd(2,0) ;
      tangent(2,1) +=  dd(2,1) ;
      tangent(2,2) +=  dd(2,2) ;      
      tangent(2,3) +=  -z*dd(2,0) ;      
      tangent(2,4) +=  -z*dd(2,1) ;
      tangent(2,5) +=  -z*dd(2,2) ;
      tangent(2,6) +=  root56*dd(2,3) ;
      tangent(2,7) +=  root56*dd(2,4) ;

      //row 4
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
      tangent(3,0) +=  z*dd(0,0) ;
      tangent(3,1) +=  z*dd(0,1) ;
      tangent(3,2) +=  z*dd(0,2) ;      
      tangent(3,3) +=  -z*z*dd(0,0) ;      
      tangent(3,4) +=  -z*z*dd(0,1) ;
      tangent(3,5) +=  -z*z*dd(0,2) ;
      tangent(3,6) +=  z*root56*dd(0,3) ;
      tangent(3,7) +=  z*root56*dd(0,4) ;

      //row 5
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
      tangent(4,0) +=  z*dd(1,0) ;
      tangent(4,1) +=  z*dd(1,1) ;
      tangent(4,2) +=  z*dd(1,2) ;      
      tangent(4,3) +=  -z*z*dd(1,0) ;      
      tangent(4,4) +=  -z*z*dd(1,1) ;
      tangent(4,5) +=  -z*z*dd(1,2) ;
      tangent(4,6) +=  z*root56*dd(1,3) ;
      tangent(4,7) +=  z*root56*dd(1,4) ;

      //row 6
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
      tangent(5,0) +=  z*dd(2,0) ;
      tangent(5,1) +=  z*dd(2,1) ;
      tangent(5,2) +=  z*dd(2,2) ;      
      tangent(5,3) +=  -z*z*dd(2,0) ;      
      tangent(5,4) +=  -z*z*dd(2,1) ;
      tangent(5,5) +=  -z*z*dd(2,2) ;
      tangent(5,6) +=  z*root56*dd(2,3) ;
      tangent(5,7) +=  z*root56*dd(2,4) ;

      //row 7
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
      tangent(6,0) +=  root56*dd(3,0) ;
      tangent(6,1) +=  root56*dd(3,1) ;
      tangent(6,2) +=  root56*dd(3,2) ;      
      tangent(6,3) +=  -root56*z*dd(3,0) ;      
      tangent(6,4) +=  -root56*z*dd(3,1) ;
      tangent(6,5) +=  -root56*z*dd(3,2) ;
      tangent(6,6) +=  root56*root56*dd(3,3) ;
      tangent(6,7) +=  root56*root56*dd(3,4) ;

      //row 8 
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
      tangent(7,0) +=  root56*dd(4,0) ;
      tangent(7,1) +=  root56*dd(4,1) ;
      tangent(7,2) +=  root56*dd(4,2) ;      
      tangent(7,3) +=  -root56*z*dd(4,0) ;      
      tangent(7,4) +=  -root56*z*dd(4,1) ;
      tangent(7,5) +=  -root56*z*dd(4,2) ;
      tangent(7,6) +=  root56*root56*dd(4,3) ;
      tangent(7,7) +=  root56*root56*dd(4,4) ;

  } //end for i

  return this->tangent ;
}


//print out data
void  DoubleMembranePlateFiberSection::Print( OPS_Stream &s, int flag )
{
  s << "DoubleMembranePlateFiberSection: \n " ;
  s <<  "  Thickness, h = "        <<  h  <<  endln ;
  s <<  "  Distance, d = "        <<  d  <<  endln ;

  for (int i = 0; i < numFibers; i++) {
    theFibers[i]->Print( s, flag ) ;
	theFibers[i+numFibers]->Print( s, flag ) ;
  }

  return ;
}

int 
DoubleMembranePlateFiberSection::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(2*numFibers+1);
  
  int i;
  for (i = 0; i < numFibers; i++) {
    idData(i) = theFibers[i]->getClassTag();
    matDbTag = theFibers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theFibers[i]->setDbTag(matDbTag);
    }
    idData(i+numFibers) = matDbTag;
  }
  
  idData(2*numFibers) = this->getTag();

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING DoubleMembranePlateFiberSection::sendSelf() - " << this->getTag() << " failed to send ID\n";
			    
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < numFibers; i++) {
    res += theFibers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING DoubleMembranePlateFiberSection::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}


int 
DoubleMembranePlateFiberSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(2*numFibers+1);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING DoubleMembranePlateFiberSection::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(2*numFibers));

  int i;

  if (theFibers[0] == 0) {
    for (i = 0; i < numFibers; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+numFibers);
      // Allocate new material with the sent class tag
      theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theFibers[i] == 0) {
	opserr << "DoubleMembranePlateFiberSection::recvSelf() - " <<
	  "Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      theFibers[i]->setDbTag(matDbTag);
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "DoubleMembranePlateFiber::recvSelf() - material " << i << "failed to recv itself\n";
	  
	return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < numFibers; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+numFibers);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theFibers[i]->getClassTag() != matClassTag) {
	delete theFibers[i];
	theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (theFibers[i] == 0) {
	  opserr << "DoubleMembranePlateFiberSection::recvSelf() - " << 
	    "Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      theFibers[i]->setDbTag(matDbTag);
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "DoubleMembranePlateFiberSection::recvSelf() - material " << 
	  i << ", failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}
 


Response*
DoubleMembranePlateFiberSection::setResponse(const char **argv, int argc,
				       OPS_Stream &output)
{
  Response *theResponse =0;

  if (argc > 2 && (strcmp(argv[0], "fiber") == 0 || strcmp(argv[0], "Fiber") == 0)) {
    
    int passarg = 2;
    int key = atoi(argv[1]);    
    
    if (key > 0 && key <= 2*numFibers) {
      int quadrature_id = key - 1;
      if (key > numFibers)
          quadrature_id -= numFibers;
      double fiber_thickness = 0.5 * h * wg[quadrature_id];
      double fiber_location = 0.5 * (d + h) + (0.5 * h) * sg[quadrature_id];
      if (key > numFibers)
          fiber_location = -fiber_location;
      output.tag("FiberOutput");
      output.attr("number", key);
      output.attr("zLoc", fiber_location);
      output.attr("thickness", fiber_thickness);
      theResponse =  theFibers[key-1]->setResponse(&argv[passarg], argc-passarg, output);
      output.endTag();
    }
  }

  // If not a fiber response, call the base class method
  if (theResponse == 0)
      return SectionForceDeformation::setResponse(argv, argc, output);

  return theResponse;
}


int 
DoubleMembranePlateFiberSection::getResponse(int responseID, Information &sectInfo)
{
  // Just call the base class method ... don't need to define
  // this function, but keeping it here just for clarity
  return SectionForceDeformation::getResponse(responseID, sectInfo);
}
