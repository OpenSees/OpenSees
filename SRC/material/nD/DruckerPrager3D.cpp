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

// Written: Ed "C++" Love

#include <DruckerPrager3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <Information.h>
#include <Parameter.h>

//#include <myDebug.h>

//static vectors and matrices
//Vector DruckerPrager3D :: strain_vec(6) ;
//Vector DruckerPrager3D :: stress_vec(6) ;
//Matrix DruckerPrager3D :: tangent_matrix(6,6) ;


//null constructor
DruckerPrager3D ::  DruckerPrager3D( ) : 
DruckerPrager( )
{  }


//full constructor
DruckerPrager3D::DruckerPrager3D(int tag, double bulk, double shear, double s_y,
							 double r, double r_bar, double Kinfinity, double Kinit, 
							 double d1, double d2, double H, double t, double atm) : 
DruckerPrager(tag, ND_TAG_DruckerPrager3D, bulk, shear, s_y, r, r_bar, Kinfinity,
	Kinit, d1, d2, H, t, atm)
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::DruckerPrager3D(...)" << endln;
	#endif
}

   
//destructor
DruckerPrager3D :: ~DruckerPrager3D( ) 
{ } 


//make a clone of this material
NDMaterial* DruckerPrager3D :: getCopy( ) 
{ 
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getCopy()" << endln;
	#endif

  DruckerPrager3D  *clone;
  clone = new DruckerPrager3D( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* DruckerPrager3D :: getType( ) const 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getType()" << endln;
	#endif

  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int DruckerPrager3D :: getOrder( ) const 
{ 
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getOrder()" << endln;
	#endif

  return 6 ; 
} 


//get the strain and integrate plasticity equations
int DruckerPrager3D :: setTrialStrain( const Vector &strain_from_element) 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::setTrialStrain(...)" << endln;
	#endif

	//}
	mEpsilon = strain_from_element;
	this->plastic_integrator( ) ;
	return 0 ;
}


//unused trial strain functions
int DruckerPrager3D::setTrialStrain (const Vector &v, const Vector &r)
{
#ifdef DEBUG
        opserr << "DruckerPrager::setTrialStrain(... ...)" << endln;
#endif

  opserr << "YOU SHOULD NOT SEE THIS: DruckerPrager::setTrialStrain (const Vector &v, const Vector &r)" << endln;
  return this->setTrialStrain (v);
}


//send back the strain
const Vector& DruckerPrager3D :: getStrain( ) 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getStrain()" << endln;
	#endif

  return mEpsilon ;
} 


//send back the stress 
const Vector& DruckerPrager3D :: getStress( ) 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getStress()" << endln;
	#endif

  return mSigma ;
}

//send back the tangent 
const Matrix& DruckerPrager3D :: getTangent( ) 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getTangent()" << endln;
	#endif

  return mCep ;
} 

//send back the tangent 
const Matrix& DruckerPrager3D :: getInitialTangent( ) 
{
	#ifdef DEBUG
        opserr << "DruckerPrager3D::getInitialTangent()" << endln;
	#endif
  return mCe ;
} 

//this is mike's problem
int DruckerPrager3D :: setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int DruckerPrager3D :: setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int DruckerPrager3D :: setTrialStrainIncr(const Tensor &v) 
{
  return -1 ;
}

int DruckerPrager3D :: setTrialStrainIncr(const Tensor &v, const Tensor &r) 
{
  return -1 ;
}
