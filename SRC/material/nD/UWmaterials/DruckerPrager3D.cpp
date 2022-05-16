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
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2010-02-04 20:50:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DruckerPrager3D.cpp,v $

// Written: K.Petek, U.Washington

#include <DruckerPrager3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>


//null constructor
DruckerPrager3D ::  DruckerPrager3D( ) : 
DruckerPrager( )
{  }


//full constructor
DruckerPrager3D::DruckerPrager3D(int tag, double bulk, double shear, double s_y,
							 double r, double r_bar, double Kinfinity, double Kinit, 
							 double d1, double d2, double H, double t, double mDen, double atm) : 
DruckerPrager(tag, ND_TAG_DruckerPrager3D, bulk, shear, s_y, r, r_bar, Kinfinity,
	                                       Kinit, d1, d2, H, t, mDen, atm)
{
}

   
//destructor
DruckerPrager3D :: ~DruckerPrager3D( ) 
{} 


//make a clone of this material
NDMaterial* DruckerPrager3D :: getCopy( ) 
{ 
  DruckerPrager3D  *clone;
  clone = new DruckerPrager3D( ) ;   //new instance of this class
  *clone = *this ;          //assignment to make copy
  return clone ;
}


//send back type of material
const char* DruckerPrager3D :: getType( ) const 
{
  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int DruckerPrager3D :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int DruckerPrager3D :: setTrialStrain( const Vector &strain_from_element) 
{
	mEpsilon = strain_from_element;
	this->plastic_integrator( ) ;
	return 0 ;
}


//unused trial strain functions
int DruckerPrager3D::setTrialStrain (const Vector &v, const Vector &r)
{
  opserr << "YOU SHOULD NOT SEE THIS: DruckerPrager::setTrialStrain (const Vector &v, const Vector &r)" << endln;
  return this->setTrialStrain (v);
}


//send back the strain
const Vector& DruckerPrager3D :: getStrain( ) 
{
  return mEpsilon ;
} 


//send back the stress 
const Vector& DruckerPrager3D :: getStress( ) 
{
  return mSigma ;
}

//send back the tangent 
const Matrix& DruckerPrager3D :: getTangent( ) 
{
  return mCep ;
} 

//send back the tangent 
const Matrix& DruckerPrager3D :: getInitialTangent( ) 
{
  return mCe ;
} 

