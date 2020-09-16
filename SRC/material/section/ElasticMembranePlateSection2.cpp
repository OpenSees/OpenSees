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
// $Source: /OpenSees/SRC/material/section/ElasticMembranePlateSection2.cpp,v $

// Developed by Pearl Ranchal
// with support from Degenkolb Engineers
// Elastic shell section with decoupled in-plane and out-of-plane stiffness parameters

// Adapted from 'ElasticMembranePlateSection.cpp' by:
// Ed "C++" Love


#include <ElasticMembranePlateSection2.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

//parameters
const double ElasticMembranePlateSection2::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector  ElasticMembranePlateSection2::stress(8) ;
Matrix  ElasticMembranePlateSection2::tangent(8,8) ;
ID      ElasticMembranePlateSection2::array(8) ;

void* OPS_ElasticMembranePlateSection2()
{
    if (OPS_GetNumRemainingInputArgs() < 5) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticMembranePlateSection2 tag? Em? Eb? nu? h? <rho?>\n";
	return 0;
    }

    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid tag\n";
	return 0;
    }

    double data[5] = {0,0,0,0,0.0};
    numdata = OPS_GetNumRemainingInputArgs();
    if (numdata > 5) numdata = 5;
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid double values\n";
	return 0;
    }

    return new ElasticMembranePlateSection2(tag,data[0],data[1],data[2],data[3],data[4]);
}

//null constructor
ElasticMembranePlateSection2::ElasticMembranePlateSection2( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticMembranePlateSection2 ), 
strain(8) 
{ 

}



//full constructor
ElasticMembranePlateSection2::ElasticMembranePlateSection2( int    tag, 
                                                            double young_m,
                                                            double young_b,
                                                            double poisson,
                                                            double thickness,
					                                        double r ) :
SectionForceDeformation( tag, SEC_TAG_ElasticMembranePlateSection2 ),
strain(8)
{
  this->Em   = young_m ;
  this->Eb   = young_b;
  this->nu   = poisson ;
  this->h    = thickness ;
  this->rhoH = r*thickness ;
}



//destructor
ElasticMembranePlateSection2::~ElasticMembranePlateSection2( ) 
{ 

} 



//make a clone of this material
SectionForceDeformation*  ElasticMembranePlateSection2::getCopy( ) 
{
  ElasticMembranePlateSection2 *clone ;   

  clone = new ElasticMembranePlateSection2(this->getTag(), Em, Eb, nu, h, rhoH) ; //new instance of this class

  //    *clone = *this ; //assignment to make copy
  clone->rhoH   = this->rhoH ;
  clone->strain = this->strain;

  return clone ;
}

//density per unit area
double
ElasticMembranePlateSection2::getRho( )
{
  return rhoH ;
}


//send back order of strain in vector form
int ElasticMembranePlateSection2::getOrder( ) const
{
  return 8 ;
}


//send back order of strain in vector form
const ID& ElasticMembranePlateSection2::getType( )
{
  return array ;
}



//swap history variables
int ElasticMembranePlateSection2::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticMembranePlateSection2::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticMembranePlateSection2::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticMembranePlateSection2 ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticMembranePlateSection2::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector&  ElasticMembranePlateSection2::getStressResultant( )
{

  double M  = Em / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; //shear modulus
 
  G *= h ;  //multiply by thickness
  M *= h ;

  //membrane resultants

  stress(0) =  M*strain(0) + (nu*M)*strain(1)  ;
 
  stress(1) =  (nu*M)*strain(0) +  M*strain(1)  ;

  stress(2) =  G*strain(2) ;

 

  // G *= five6 ;  //multiply by shear correction factor
  G *= (Eb / Em);  //multiply by ratio of bending to membrane moduli

  double D  =  Eb * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

  //bending resultants

  stress(3) = -( D*strain(3) + nu*D*strain(4) ) ;
 
  stress(4) = -( nu*D*strain(3) + D*strain(4) ) ;

  stress(5) = -0.5*D*( 1.0 - nu )*strain(5) ;

  stress(6) = G*strain(6) ;

  stress(7) = G*strain(7) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix&  ElasticMembranePlateSection2::getSectionTangent( )
{

  double M  = Em / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



 // G *= five6 ;  //multiply by shear correction factor
  G *= (Eb/Em);  //multiply by ratio of bending to membrane moduli

  double D  =  Eb * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

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
const Matrix&  ElasticMembranePlateSection2::getInitialTangent( )
{

  double M  = Em / ( 1.0 - nu*nu ) ; //membrane modulus

  double G  =  0.5 * Em / ( 1.0 + nu ) ; //shear modulus

  G *= h ;  //multiply by thickness
  M *= h ;

  tangent.Zero() ;

  //membrane tangent terms

  tangent(0,0) = M ;
  tangent(1,1) = M ;

  tangent(0,1) = nu*M ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = G ;



  // G *= five6 ;  //multiply by shear correction factor
  G *= (Eb / Em);  //multiply by ratio of bending to membrane moduli

  double D  =  Eb * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;  //bending modulus

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
void  ElasticMembranePlateSection2::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticMembranePlateSection: \n ";
        s << "  Young's Modulus-Membrane Em = " << Em << endln;
        s << "  Young's Modulus-Bending Eb = " << Eb << endln;
        s << "  Poisson's Ratio nu = " << nu << endln;
        s << "  Thickness h = " << h << endln;
        s << "  Density rho = " << (rhoH/h) << endln;
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticMembranePlateSection2\", ";
        s << "\"Em\": " << Em << ", ";
        s << "\"Eb\": " << Eb << ", ";
        s << "\"nu\": " << nu << ", ";
        s << "\"thickness\": " << h << ", ";
        s << "\"masspervolume\": " << (rhoH/h) << "}";
    }
}


int ElasticMembranePlateSection2::sendSelf(int cTag, Channel &theChannel) 
{
  int res = 0;
  static Vector data(5);
  data(0) = this->getTag();
  data(1) = Em;
  data(2) = Eb;
  data(3) = nu;
  data(4) = h;
  data(5) = rhoH;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection2::sendSelf() - failed to send data\n";

  return res;
}


int ElasticMembranePlateSection2::recvSelf(int cTag, Channel &theChannel, 
				      FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(5);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticMembranePlateSection2::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    Em   = data(1);
    Eb   = data(1);
    nu   = data(2);
    h    = data(3);
    rhoH = data(4);
  }

  return res;
}
