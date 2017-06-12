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
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/DruckerPrager3DThermal.cpp,v $

// Written: K.Petek, U.Washington
//Modified by Jian Jiang, Liming Jiang [http://openseesforfire.github.io]

#include <DruckerPrager3DThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>


//null constructor
DruckerPrager3DThermal ::  DruckerPrager3DThermal( ) : 
DruckerPragerThermal( )
{  }


//full constructor
DruckerPrager3DThermal::DruckerPrager3DThermal(int tag, double bulk, double shear, double s_y,
							 double r, double r_bar, double Kinfinity, double Kinit, 
							 double d1, double d2, double H, double t, double mDen, double atm) : 
DruckerPragerThermal(tag, ND_TAG_DruckerPrager3DThermal, bulk, shear, s_y, r, r_bar, Kinfinity,
	                                       Kinit, d1, d2, H, t, mDen, atm)
{
}

   
//destructor
DruckerPrager3DThermal :: ~DruckerPrager3DThermal( ) 
{ } 


//make a clone of this material
NDMaterial* DruckerPrager3DThermal :: getCopy( ) 
{ 
  DruckerPrager3DThermal  *clone;
  clone = new DruckerPrager3DThermal( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* DruckerPrager3DThermal :: getType( ) const 
{
  return "ThreeDimensionalThermal" ;
}


//send back order of strain in vector form
int DruckerPrager3DThermal :: getOrder( ) const 
{ 
  return 6 ; 
} 


//get the strain and integrate plasticity equations
int DruckerPrager3DThermal :: setTrialStrain( const Vector &strain_from_element) 
{
	//}
	mEpsilon = strain_from_element;
	this->plastic_integrator( ) ;
	return 0 ;
}


//unused trial strain functions
int DruckerPrager3DThermal::setTrialStrain (const Vector &v, const Vector &r)
{
  opserr << "YOU SHOULD NOT SEE THIS: DruckerPrager::setTrialStrain (const Vector &v, const Vector &r)" << endln;
  return this->setTrialStrain (v);
}


//send back the strain
const Vector& DruckerPrager3DThermal :: getStrain( ) 
{
  return mEpsilon ;
} 


//send back the stress 
const Vector& DruckerPrager3DThermal :: getStress( ) 
{
  return mSigma ;
}

//send back the tangent 
const Matrix& DruckerPrager3DThermal :: getTangent( ) 
{
  return mCep ;
} 

//send back the tangent 
const Matrix& DruckerPrager3DThermal :: getInitialTangent( ) 
{
  return mCe ;
} 

//send back TempAndElong(Liming,UoE)
const Vector& DruckerPrager3DThermal::getTempAndElong( ){

 return TempAndElong;
}