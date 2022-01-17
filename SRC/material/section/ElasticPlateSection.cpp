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
// $Date: 2003-02-14 23:01:33 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/ElasticPlateSection.cpp,v $

// Ed "C++" Love
//
//  Elastic Plate Section
//


#include <ElasticPlateSection.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

void *
OPS_ADD_RUNTIME_VPV(OPS_ElasticPlateSection)
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr << "WARNING insufficient arguments\n";
	opserr << "Want: section ElasticPlateSection tag? E? nu? h? " << endln;
	return 0;
    }
	
    int tag;
    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid section ElasticPlateSection tag" << endln;
	return 0;
    }

    numdata = 3;
    double data[3];
    if (OPS_GetDoubleInput(&numdata, data) < 0) {
	opserr << "WARNING invalid section ElasticPlateSection double inputs" << endln;
	return 0;
    }
    double E = data[0];
    double nu = data[1];
    double h = data[2];

    return new ElasticPlateSection (tag, E, nu, h);
}

//parameters
const double ElasticPlateSection::five6 = 5.0/6.0 ; //shear correction

//static vector and matrices
Vector ElasticPlateSection::stress(5) ;
Matrix ElasticPlateSection::tangent(5,5) ;
ID     ElasticPlateSection::array(5) ;


//null constructor
ElasticPlateSection::ElasticPlateSection( ) : 
SectionForceDeformation( 0, SEC_TAG_ElasticPlateSection ), 
strain(5) 
{ }



//full constructor
ElasticPlateSection::ElasticPlateSection(  int    tag, 
                                           double young,
                                           double poisson,
                                           double thickness ) :
SectionForceDeformation( tag, SEC_TAG_ElasticPlateSection ),
strain(5)
{
  this->E  = young ;
  this->nu = poisson ;
  this->h  = thickness ;
}



//destructor
ElasticPlateSection::~ElasticPlateSection( ) 
{ } 



//make a clone of this material
SectionForceDeformation* ElasticPlateSection::getCopy( ) 
{
  ElasticPlateSection *clone ;   

  clone = new ElasticPlateSection( ) ; //new instance of this class

  *clone = *this ; //assignment to make copy

  return clone ;
}



//send back order of strain in vector form
int ElasticPlateSection::getOrder( ) const
{
  return 5 ;
}


//send back order of strain in vector form
const ID& ElasticPlateSection::getType( ) 
{
    static bool initialized = false;
    if (!initialized) {
        array(0) = SECTION_RESPONSE_MXX;
        array(1) = SECTION_RESPONSE_MYY;
        array(2) = SECTION_RESPONSE_MXY;
        array(3) = SECTION_RESPONSE_VXZ;
        array(4) = SECTION_RESPONSE_VYZ;
        initialized = true;
    }
    return array;
}



//swap history variables
int ElasticPlateSection::commitState( ) 
{
  return 0 ;
}



//revert to last saved state
int ElasticPlateSection::revertToLastCommit( )
{
  return 0 ;
}

//revert to start
int ElasticPlateSection::revertToStart( )
{
  return 0 ;
}


//get the strain 
int ElasticPlateSection ::
setTrialSectionDeformation( const Vector &strain_from_element)
{
  this->strain = strain_from_element ;

  return 0 ;
}


//send back the strain
const Vector& ElasticPlateSection::getSectionDeformation( )
{
  return this->strain ;
}


//send back the stress 
const Vector& ElasticPlateSection::getStressResultant( )
{
  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ; //bending modulus

  double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus
 
  G *= five6 ;
  G *= h ;


  stress(0) = -( D*strain(0) + nu*D*strain(1) ) ;
 
  stress(1) = -( nu*D*strain(0) + D*strain(1) ) ;

  stress(2) = -0.5*D*( 1.0 - nu )*strain(2) ;

  stress(3) = G*strain(3) ;

  stress(4) = G*strain(4) ;

 
  return this->stress ;
}


//send back the tangent 
const Matrix& ElasticPlateSection::getSectionTangent( )
{

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;

  double G  =  0.5 * E / ( 1.0 + nu ) ;


  tangent.Zero() ;

  tangent(0,0) = -D ;
  tangent(1,1) = -D ;

  tangent(0,1) = -nu*D ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(3,3) = five6*G*h ;

  tangent(4,4) = tangent(3,3) ;


  return this->tangent ;
}

//send back the initial tangent 
const Matrix& ElasticPlateSection::getInitialTangent( )
{

  double D  =  E * (h*h*h) / 12.0 / ( 1.0 - nu*nu ) ;

  double G  =  0.5 * E / ( 1.0 + nu ) ;


  tangent.Zero() ;

  tangent(0,0) = -D ;
  tangent(1,1) = -D ;

  tangent(0,1) = -nu*D ;
  tangent(1,0) = tangent(0,1) ;

  tangent(2,2) = -0.5 * D * ( 1.0 - nu ) ;

  tangent(3,3) = five6*G*h ;

  tangent(4,4) = tangent(3,3) ;


  return this->tangent ;
}

//print out data
void  ElasticPlateSection::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION) {
        s << "ElasticPlateSection: \n ";
        s << "  Young's Modulus E  = " << E << endln;
        s << "  Poisson's Ratio nu = " << nu << endln;
        s << "  Thickness h = " << h << endln;
    }

    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"ElasticPlateSection\", ";
        s << "\"E\": " << E << ", ";
        s << "\"nu\": " << nu << ", ";
        s << "\"thickness\": " << h << "}";
    }
}


int 
ElasticPlateSection::sendSelf(int cTag, Channel &theChannel) 
{
  int res = 0;
  static Vector data(4);
  data(0) = this->getTag();
  data(1) = E;
  data(2) = nu;
  data(3) = h;

  res = theChannel.sendVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPlateSection::sendSelf() - failed to send data\n";

  return res;
}


int 
ElasticPlateSection::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  static Vector data(4);
  res = theChannel.recvVector(this->getDbTag(), cTag, data);
  if (res < 0) 
    opserr << "ElasticPlateSection::recvSelf() - failed to recv data\n";
  else {
    this->setTag(data(0));
    E    = data(1);
    nu   = data(2);
    h    = data(3);
  }

  return res;
}
