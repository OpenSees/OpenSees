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
                                                                        
// Date: 2019-02-15
// /usr/local/cvs/OpenSees/SRC/material/section/ElasticMembranePlateSectionUncoupled.h,v $

// KK: Minor modification to ElasticMembranePlateSection to decouple membrane (in-plane) and bending (out-of-plane) action
//
//  Elastic Plate Section with membrane
//


#include <ElasticMembranePlateSectionUncoupled.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>

//parameters
const double ElasticMembranePlateSectionUncoupled::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  ElasticMembranePlateSectionUncoupled::stress(8) ;
Matrix  ElasticMembranePlateSectionUncoupled::tangent(8,8) ;
ID      ElasticMembranePlateSectionUncoupled::array(8) ;


//null constructor
ElasticMembranePlateSectionUncoupled::ElasticMembranePlateSectionUncoupled( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticMembranePlateSectionUncoupled),
strain(8) 
{ 

}



//full constructor
ElasticMembranePlateSectionUncoupled::ElasticMembranePlateSectionUncoupled(  
											int    tag, 
											double young,
											double poisson,
											double thickness,
											double r, 
											double tbend, 
											double tmembrane):
SectionForceDeformation( tag, SEC_TAG_ElasticMembranePlateSectionUncoupled),
strain(8)
{
  this->E  = young ;
  this->nu = poisson ;
  this->h  = thickness ;
  this->rhoH = r*thickness ;
  this->TBND = tbend ;
  this->TMEMB = tmembrane;
}



//destructor
ElasticMembranePlateSectionUncoupled::~ElasticMembranePlateSectionUncoupled( ) 
{ 

} 



//make a clone of this material
SectionForceDeformation*  ElasticMembranePlateSectionUncoupled::getCopy( ) 
{
  ElasticMembranePlateSectionUncoupled *clone ;   

  clone = new ElasticMembranePlateSectionUncoupled(this->getTag(), E, nu, h, rhoH, TBND, TMEMB) ; //new instance of this class

  //    *clone = *this ; //assignment to make copy
  clone->rhoH = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

//density per unit area
double
ElasticMembranePlateSectionUncoupled::getRho( )
{
  return rhoH ;
}


//send back order of strain in vector form
int ElasticMembranePlateSectionUncoupled::getOrder( ) const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticMembranePlateSectionUncoupled::getType( )
{
  return array ;
}



//swap history variables
int ElasticMembranePlateSectionUncoupled::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticMembranePlateSectionUncoupled::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticMembranePlateSectionUncoupled::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticMembranePlateSectionUncoupled ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticMembranePlateSectionUncoupled::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  ElasticMembranePlateSectionUncoupled::getStressResultant( )
{

  double Mm  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double Gm =  0.5 * E / ( 1.0 + nu ) ; //shear modulus for membrane
  double Gb = 0.5 * E / (1.0 + nu); //shear modulus for bending

  Mm *= h * TMEMB;

  Gm *= h * TMEMB;  //multiply by thickness
  Gb *= h * TBND;  //multiply by thickness

  //membrane resultants

  stress(0) =  Mm*strain(0) + (nu*Mm)*strain(1)  ;
 
  stress(1) =  (nu*Mm)*strain(0) +  Mm*strain(1)  ;

  stress(2) =  Gm*strain(2) ;
 
  //bending resultants
  Gb *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h*TBND*TBND*TBND) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  stress(3) = -( D*strain(3) + nu*D*strain(4) ) ;
 
  stress(4) = -( nu*D*strain(3) + D*strain(4) ) ;

  stress(5) = -0.5*D*( 1.0 - nu )*strain(5) ;

  stress(6) = Gb*strain(6) ;

  stress(7) = Gb*strain(7) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix&  ElasticMembranePlateSectionUncoupled::getSectionTangent( )
{
	
  double Mm = E / (1.0 - nu * nu); //membrane modulus

  double Gm = 0.5 * E / (1.0 + nu); //shear modulus for membrane
  double Gb = 0.5 * E / (1.0 + nu); //shear modulus for bending

  Mm *= h * TMEMB;

  Gm *= h * TMEMB;  //multiply by thickness
  Gb *= h * TBND;  //multiply by thickness

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = Mm ;
  tangent(1,1) = Mm ;

  tangent(0,1) = nu*Mm ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = Gm ;


  //bending tangent terms
  Gb *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h*TBND*TBND*TBND) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = Gb ;

  tangent(7,7) = Gb ;

  return this->tangent ;
}


//send back the initial tangent 
const Matrix&  ElasticMembranePlateSectionUncoupled::getInitialTangent( )
{

  double Mm = E / (1.0 - nu * nu); //membrane modulus

  double Gm = 0.5 * E / (1.0 + nu); //shear modulus for membrane
  double Gb = 0.5 * E / (1.0 + nu); //shear modulus for bending

  Mm *= h * TMEMB;

  Gm *= h * TMEMB;  //multiply by thickness
  Gb *= h * TBND;  //multiply by thickness

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = Mm ;
  tangent(1,1) = Mm ;

  tangent(0,1) = nu*Mm ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = Gm ;


  //bending tangent terms

  Gb *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h*TBND*TBND*TBND) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = Gb ;

  tangent(7,7) = Gb ;

  return this->tangent ;
}


//print out data
void  ElasticMembranePlateSectionUncoupled::Print( OPS_Stream &s, int flag )
{
  s << "ElasticMembranePlateSectionUncoupled: \n " ;
  s <<  "  Young's Modulus E = "  <<  E  <<  endln ;
  s <<  "  Poisson's Ratio nu = " <<  nu <<  endln ;
  s <<  "  Thickness h = "        <<  h  <<  endln ;
  s <<  "  Density rho = "        <<  (rhoH/h)  <<  endln ;
  s << "  Bending modifier = " << TBND << endln;
  s << "  Membrane modifier = " << TMEMB << endln;

  return ;
}

int 
ElasticMembranePlateSectionUncoupled::sendSelf(int cTag, Channel &theChannel) 
{
  int res = 0;
  static Vector data(7);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = h;
  data(4) = rhoH;
  data(5) = TBND;
  data(6) = TMEMB;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSectionUncoupled::sendSelf() - failed to send data\n";

  return res;
}


int 
ElasticMembranePlateSectionUncoupled::recvSelf(int cTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(7);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSectionUncoupled::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E    = data(1);
    nu   = data(2);
    h    = data(3);
    rhoH = data(4);
	TBND = data(5);
	TMEMB = data(6);
  }

  return res;
}
 
