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
// $Date: 2012-05-26 22:03:16 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/nD/PlateRebarMaterialThermal.cpp,v $

//
// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Generic Plate Rebar Material
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/

#include <PlateRebarMaterialThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>   //Antonios Vytiniotis used for the recorder
#include <math.h>

//static vector and matrices
Vector  PlateRebarMaterialThermal::stress(5) ;
Matrix  PlateRebarMaterialThermal::tangent(5,5) ;

//null constructor
PlateRebarMaterialThermal::PlateRebarMaterialThermal( ) : 
NDMaterial(0, ND_TAG_PlateRebarMaterialThermal ), 
strain(5) 
{ }


//full constructor
PlateRebarMaterialThermal::PlateRebarMaterialThermal(int tag,
                                       UniaxialMaterial &uniMat,
                                       double ang) :
NDMaterial( tag, ND_TAG_PlateRebarMaterialThermal ),
strain(5),angle(ang), temperature(0)
{ 
  theMat = uniMat.getCopy() ;
  double rang = ang * 0.0174532925;
  c = cos(rang);
  s = sin(rang);
}


//destructor
PlateRebarMaterialThermal::~PlateRebarMaterialThermal( ) 
{ 
  if (theMat != 0) delete theMat ;
} 


//make a clone of this material
NDMaterial*
PlateRebarMaterialThermal::getCopy( ) 
{
  PlateRebarMaterialThermal *clone ;   //new instance of this class

  clone = new PlateRebarMaterialThermal( this->getTag(), 
                                  *theMat,
                                  angle ) ; //make the copy

  return clone ;
}


//make a clone of this material
NDMaterial* 
PlateRebarMaterialThermal::getCopy( const char *type ) 
{
  if (strcmp(type,this->getType()) == 0)
    return this->getCopy( ) ;
  else
    return 0;
}


//send back order of strain in vector form
int 
PlateRebarMaterialThermal::getOrder( ) const
{
  return 5 ;
}


const char*
PlateRebarMaterialThermal::getType( ) const 
{
  return "PlateFiberThermal" ; 
}



//swap history variables
int 
PlateRebarMaterialThermal::commitState( ) 
{
  return theMat->commitState( ) ;
}



//revert to last saved state
int 
PlateRebarMaterialThermal::revertToLastCommit( )
{
  return theMat->revertToLastCommit( ) ;
}


//revert to start
int
PlateRebarMaterialThermal::revertToStart( )
{
  strain.Zero();
  return theMat->revertToStart( ) ;
}


//mass per unit volume
double
PlateRebarMaterialThermal::getRho( )
{
  return theMat->getRho( ) ;
}


//receive the strain
int 
PlateRebarMaterialThermal::setTrialStrain( const Vector &strainFromElement )
{
  strain(0) = strainFromElement(0) ;
  strain(1) = strainFromElement(1) ;
  strain(2) = strainFromElement(2) ;
  strain(3) = strainFromElement(3) ;
  strain(4) = strainFromElement(4) ;

#ifdef _DEBUG
 // opserr<<"PlateRebar"<< strain<<endln;
#endif

  return theMat->setTrialStrain(   strain(0) * c * c
                                 + strain(1) * s * s
                                 + strain(2) * c * s,temperature,
                                 0) ;
}

//Added temperature history
double 
PlateRebarMaterialThermal::getThermalTangentAndElongation(double &TempT, double&ET, double&Elong)
{
    temperature = TempT;
	double tangent =0.0;
	double ThermalElongation =0.0;
    static Vector tData(4);
    static Information iData(tData);
    tData(0) = temperature;
	tData(1) = tangent;
	tData(2) = ThermalElongation;
    tData(3) = temperature;
    iData.setVector(tData);
    theMat->getVariable("ElongTangent", iData);   //Actually here update initial tangent and thermalElongation corresponding  to current temperature
    tData = iData.getData();
    ET = tData(1);
    ThermalElongation = tData(2);
	Elong = ThermalElongation;
return 0;
}


//send back the strain
const Vector& 
PlateRebarMaterialThermal::getStrain( )
{
  return strain ;
}


//send back the stress 
const Vector&  
PlateRebarMaterialThermal::getStress( )
{
  double sig = theMat->getStress();
  stress(0) = sig * c * c;
  stress(1) = sig * s * s;
  stress(2) = sig * c * s;
  stress(3) = 0.0;
  stress(4) = 0.0;
  //opserr<<"PlaterRebarThermal : stress "<<stress<<endln;
  return stress ;
}


//send back the tangent 
const Matrix&  
PlateRebarMaterialThermal::getTangent( )
{
  double tan = theMat->getTangent( ) ;

  tangent(0,0) = tan * c * c * c * c ;
  tangent(0,1) = tan * c * c * c * s ;
  tangent(0,2) = tan * c * c * s * s ;
  tangent(1,0) = tangent(0,1) ;
  tangent(1,1) = tangent(0,2) ;
  tangent(1,2) = tan * c * s * s * s ;
  tangent(2,0) = tangent(0,2) ;
  tangent(2,1) = tangent(1,2) ;
  tangent(2,2) = tan * s * s * s * s ;

  return tangent ;
}

const Matrix&  
PlateRebarMaterialThermal::getInitialTangent
( )
{
  double tan = theMat->getInitialTangent( ) ;

  tangent(0,0) = tan * c * c * c * c ;
  tangent(0,1) = tan * c * c * c * s ;
  tangent(0,2) = tan * c * c * s * s ;
  tangent(1,0) = tangent(0,1) ;
  tangent(1,1) = tangent(0,2) ;
  tangent(1,2) = tan * c * s * s * s ;
  tangent(2,0) = tangent(0,2) ;
  tangent(2,1) = tangent(1,2) ;
  tangent(2,2) = tan * s * s * s * s ;

  return tangent ;
}


//print out data
void  
PlateRebarMaterialThermal::Print( OPS_Stream &s, int flag )
{
  s << "PlateRebar Material tag: " << this->getTag() << endln ; 
  s << "using uniaxialmaterials : " << endln ;

  theMat->Print( s, flag ) ;

  return ;
}


int 
PlateRebarMaterialThermal::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  int matDbTag;
  
  static ID idData(3);
  idData(0) = dataTag;
  idData(1) = theMat->getClassTag();
  matDbTag = theMat->getDbTag();
  if (matDbTag == 0) {
    matDbTag = theChannel.getDbTag();
    theMat->setDbTag(matDbTag);
  }
  idData(2) = matDbTag;

  res = theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to send data" << endln;
    return res;
  }

  static Vector vecData(1);
  vecData(0) = angle;

  res = theChannel.sendVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to send data" << endln;
    return res;
  }

  // now send the materials data
  res += theMat->sendSelf(commitTag, theChannel);
  if (res < 0) 
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to send material1" << endln;

  return res;
}

int 
PlateRebarMaterialThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;

  int dataTag = this->getDbTag();

  // recv an id containg the tag and associated materials class and db tags
  static ID idData(3);
  res = theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to receive id data" << endln;
    return res;
  }

  this->setTag(idData(0));
  int matClassTag = idData(1);
  if (theMat->getClassTag() != matClassTag) {
    if (theMat != 0) delete theMat;
    theMat = theBroker.getNewUniaxialMaterial(matClassTag);
    if (theMat == 0) {
      opserr << "PlateRebarMaterialThermal::recvSelf() - failed to get a material of type: " << matClassTag << endln;
      return -1;
    }
  }
  theMat->setDbTag(idData(2));

  static Vector vecData(1);
  res = theChannel.recvVector(dataTag, commitTag, vecData);
  if (res < 0) {
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to receive vector data" << endln;
    return res;
  }
  angle = vecData(0);
  double rang = angle * 0.0174532925;
  c = cos(rang);
  s = sin(rang);

  // now receive the materials data
  res = theMat->recvSelf(commitTag, theChannel, theBroker);
  if (res < 0) 
    opserr << "PlateRebarMaterialThermal::sendSelf() - failed to receive material1" << endln;
  
  return res;
}
 

Response*
PlateRebarMaterialThermal::setResponse (const char **argv, int argc, OPS_Stream &output)
{
	static Vector tempData(2);
	  static Information infoData(tempData);
	Response *theResponse =0;
	const char *matType = this->getType();

	output.tag("UniaxialMaterialOutput");
	output.attr("matType",this->getClassType());
	output.attr("matTag",this->getTag());

	if (strcmp(argv[0],"stress") == 0 ){
		output.tag("ResponseType", "sigma11");
		theResponse =  new MaterialResponse(theMat, 1, theMat->getStress());
	}
	else if (strcmp(argv[0],"strain") == 0 ){
		output.tag("ResponseType", "eps11");
    theResponse =  new MaterialResponse(theMat, 3, theMat->getStrain());
  }
	else if (strcmp(argv[0], "tangent") == 0){
    output.tag("ResponseType", "C11");
    theResponse =  new MaterialResponse(theMat, 2, theMat->getTangent());
  }
	else if (strcmp(argv[0], "TempAndElong") == 0){
	output.tag("ResponseType", "temp11");
	if((theMat->getVariable("TempAndElong", infoData))!=0){
			 opserr<<"Warning: invalid tag in uniaxialMaterial:getVariable"<<endln;
			 return 0;
	 }
	tempData = infoData.getData();
    theResponse =  new MaterialResponse(this, 4, tempData);
  }
	else
		return 0;

	return theResponse;
}

int PlateRebarMaterialThermal::getResponse (int responseID, Information &matInfo)
{
	  static Vector tempData(2);
	  static Information infoData(tempData);
	switch (responseID) {
		case -1:
			return -1;
		case 1:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = theMat->getStress();
			return 0;
		case 2:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = theMat->getStrain();
			return 0;
		case 3:
			if (matInfo.theVector != 0)
				*(matInfo.theVector) = theMat->getTangent();
			return 0;
		case 4:
			if((theMat->getVariable("TempAndElong", infoData))!=0){
			 opserr<<"Warning: invalid tag in uniaxialMaterial:getVariable"<<endln;
			 return -1;
			}
			tempData = infoData.getData();
			matInfo.setVector(tempData);
			return 0;
		default:
			return -1;
	}
}