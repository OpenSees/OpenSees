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
// $Date: 2010-03-04 19:10:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoundingCamClay3D.cpp,v $

// $Revision: 1.2 $
// $Date: 2010-03-04 19:10:52 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/BoundingCamClay3D.cpp,v $

// Written: kap	
// Created: 12/04

#include <BoundingCamClay3D.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//#include <myDebug.h>


//null constructor
BoundingCamClay3D ::  BoundingCamClay3D( ) : 
BoundingCamClay( )
{  }


//full constructor
BoundingCamClay3D::BoundingCamClay3D(int tag, double C, double r, double R, double p_o, double kappa,
							 double mu_o, double alpha, double lambda, double h, 
							 double m, double epsE_vo) : 
BoundingCamClay(tag, ND_TAG_BoundingCamClay3D, C, r, R, p_o, kappa, mu_o, alpha, lambda, 
				h, m, epsE_vo)
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::BoundingCamClay3D(...)" << endln;
	#endif
}
   


//destructor
BoundingCamClay3D :: ~BoundingCamClay3D( ) 
{ } 


//make a clone of this material
NDMaterial* BoundingCamClay3D :: getCopy( ) 
{ 
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getCopy()" << endln;
	#endif

  BoundingCamClay3D  *clone;
  clone = new BoundingCamClay3D( ) ;   //new instance of this class
  *clone = *this ;          //asignment to make copy
  return clone ;
}


//send back type of material
const char* BoundingCamClay3D :: getType( ) const 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getType()" << endln;
	#endif

  return "ThreeDimensional" ;
}


//send back order of strain in vector form
int BoundingCamClay3D :: getOrder( ) const 
{ 
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getOrder()" << endln;
	#endif

  return 6 ; 
} 


//get the strain and integrate plasticity equations
int BoundingCamClay3D :: setTrialStrain( const Vector &strain_from_element) 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::setTrialStrain(...)" << endln;
	#endif

	mEpsilon = strain_from_element;

    this->plastic_integrator( ) ;
	
	return 0 ;
}


//unused trial strain functions
int BoundingCamClay3D::setTrialStrain (const Vector &v, const Vector &r)
{
#ifdef DEBUG
        opserr << "BoundingCamClay::setTrialStrain(... ...)" << endln;
#endif

  opserr << "YOU SHOULD NOT SEE THIS: BoundingCamClay::setTrialStrain (const Vector &v, const Vector &r)" << endln;
  return this->setTrialStrain (v);
}




//send back the strain
const Vector& BoundingCamClay3D :: getStrain( ) 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getStrain()" << endln;
	#endif

  return mEpsilon ;
} 


//send back the stress 
const Vector& BoundingCamClay3D :: getStress( ) 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getStress()" << endln;
	#endif

  return mSigma ;
}

//send back the tangent 
const Matrix& BoundingCamClay3D :: getTangent( ) 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getTangent()" << endln;
	#endif

  return mCep;
} 

//send back the tangent 
const Matrix& BoundingCamClay3D :: getInitialTangent( ) 
{
	#ifdef DEBUG
        opserr << "BoundingCamClay3D::getInitialTangent()" << endln;
	#endif
  return mCep;
} 

//this is mike's problem
int BoundingCamClay3D :: setTrialStrain(const Tensor &v) 
{
  return -1 ;
}

int BoundingCamClay3D :: setTrialStrain(const Tensor &v, const Tensor &r)     
{
  return -1 ;
}

int BoundingCamClay3D :: setTrialStrainIncr(const Tensor &v) 
{
  return -1 ;
}

int BoundingCamClay3D :: setTrialStrainIncr(const Tensor &v, const Tensor &r) 
{
  return -1 ;
}







