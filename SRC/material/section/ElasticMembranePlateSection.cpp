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
// $Date: 2003-02-14 23:01:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticMembranePlateSection.cpp,v $

// Ed "C++" Love
//
//  Elastic Plate Section with membrane
//


#include <ElasticMembranePlateSection.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//parameters
const double ElasticMembranePlateSection::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  ElasticMembranePlateSection::stress(8) ;
Matrix  ElasticMembranePlateSection::tangent(8,8) ;
ID      ElasticMembranePlateSection::array(8) ;

void* OPS_ElasticMembranePlateSection()
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticMembranePlateSection tag? E? nu? h? <rho?>\n";
	return 0;
    }

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[4] = {0,0,0,0.0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 4) numdata = 4;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double values\n";
	return 0;
    }

    return new ElasticMembranePlateSection(tag,data[0],data[1],data[2],data[3]);
}

//null constructor
ElasticMembranePlateSection::ElasticMembranePlateSection( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticMembranePlateSection ), 
strain(8) 
{ 

}



//full constructor
ElasticMembranePlateSection::ElasticMembranePlateSection(  
                                           int    tag, 
                                           double young,
                                           double poisson,
                                           double thickness,
					   double r ) :
SectionForceDeformation( tag, SEC_TAG_ElasticMembranePlateSection ),
strain(8)
{
  this->E  = young ;
  this->nu = poisson ;
  this->h  = thickness ;
  this->rhoH = r*thickness ;
}



//destructor
ElasticMembranePlateSection::~ElasticMembranePlateSection( ) 
{ 

} 



//make a clone of this material
SectionForceDeformation*  ElasticMembranePlateSection::getCopy( ) 
{
  ElasticMembranePlateSection *clone ;   

  clone = new ElasticMembranePlateSection(this->getTag(), E, nu, h, rhoH) ; //new instance of this class

  //    *clone = *this ; //assignment to make copy
  clone->rhoH = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

//density per unit area
double
ElasticMembranePlateSection::getRho( )
{
  return rhoH ;
}


//send back order of strain in vector form
int ElasticMembranePlateSection::getOrder( ) const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticMembranePlateSection::getType( )
{
  return array ;
}



//swap history variables
int ElasticMembranePlateSection::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticMembranePlateSection::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticMembranePlateSection::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticMembranePlateSection ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticMembranePlateSection::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  ElasticMembranePlateSection::getStressResultant( )
{

  double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus
 
  G *= h ;  //multiply by thickness
  M *= h ;

  //membrane resultants

  stress(0) =  M*strain(0) + (nu*M)*strain(1)  ;
 
  stress(1) =  (nu*M)*strain(0) +  M*strain(1)  ;

  stress(2) =  G*strain(2) ;

 

  G *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending resultants

  stress(3) = -( D*strain(3) + nu*D*strain(4) ) ;
 
  stress(4) = -( nu*D*strain(3) + D*strain(4) ) ;

  stress(5) = -0.5*D*( 1.0 - nu )*strain(5) ;

  stress(6) = G*strain(6) ;

  stress(7) = G*strain(7) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix&  ElasticMembranePlateSection::getSectionTangent( )
{

  double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  G *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending tangent terms

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = G ;

  tangent(7,7) = G ;

  return this->tangent ;
}


//send back the initial tangent 
const Matrix&  ElasticMembranePlateSection::getInitialTangent( )
{

  double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  G *= five6 ;  //multiply by shear correction factor

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending tangent terms

  tangent(3,3) = -D ;
  tangent(4,4) = -D ;

  tangent(3,4) = -nu*D ;
  tangent(4,3) = tangent(3,4) ;

  tangent(5,5) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(6,6) = G ;

  tangent(7,7) = G ;

  return this->tangent ;
}


//print out data
void  ElasticMembranePlateSection::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticMembranePlateSection: \n ";
        s << "  Young's Modulus E = " << E << endln;
        s << "  Poisson's Ratio nu = " << nu << endln;
        s << "  Thickness h = " << h << endln;
        s << "  Density rho = " << (rhoH/h) << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticMembranePlateSection\", ";
        s << "\"E\": " << E << ", ";
        s << "\"nu\": " << nu << ", ";
        s << "\"thickness\": " << h << ", ";
        s << "\"masspervolume\": " << (rhoH/h) << "}";
    }
}


int ElasticMembranePlateSection::sendSelf(int cTag, Channel &theChannel) 
{
  int res = 0;
  static Vector data(5);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = h;
  data(4) = rhoH;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::sendSelf() - failed to send data\n";

  return res;
}


int ElasticMembranePlateSection::recvSelf(int cTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E    = data(1);
    nu   = data(2);
    h    = data(3);
    rhoH = data(4);
  }

  return res;
}
