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
// $Date: 2012-05-21 23:03:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/LayeredShellFiberSection.cpp,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Layered Shell Section
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/


#include <LayeredShellFiberSection.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>
#include <elementAPI.h>

void* OPS_LayeredShellFiberSection()
{
    if (OPS_GetNumRemainingInputArgs() < 4) {
	opserr << "WARNING insufficient arguments" << endln;
	opserr << "Want: section LayeredShell tag? nLayers? matTag1? h1? ... matTagn? hn? " << endln;
	return 0;
    }
      
    int tag, nLayers, matTag;
    double h, *thickness;
    NDMaterial **theMats;

    int numdata = 1;
    if (OPS_GetIntInput(&numdata, &tag) < 0) {
	opserr << "WARNING invalid section LayeredShell tag" << endln;
	return 0;
    }

    if (OPS_GetIntInput(&numdata, &nLayers) < 0) {
	opserr << "WARNING invalid nLayers" << endln;
	opserr << "LayeredShell section: " << tag << endln;	    	    
	return 0;
    }
	
    if (nLayers < 3) {
	opserr << "ERROR number of layers must be larger than 2" << endln;
	opserr << "LayeredShell section: " << tag << endln;	    	    
	return 0;
    }
      
    theMats   = new NDMaterial*[nLayers];
    thickness = new double[nLayers];
      
    for (int iLayer = 0; iLayer < nLayers; iLayer++) {
	if (OPS_GetNumRemainingInputArgs() < 2) {
	    opserr << "WARNING must provide "<<2*nLayers<<"inputs\n";
	    return 0;
	}
	if (OPS_GetIntInput(&numdata, &matTag) < 0) {
	    opserr << "WARNING invalid matTag" << endln;
	    opserr << "LayeredShell section: " << tag << endln;
	    return 0;
	}
	
	theMats[iLayer] = OPS_getNDMaterial(matTag);
	if (theMats[iLayer] == 0) {
	    opserr << "WARNING nD material does not exist" << endln;;
	    opserr << "nD material: " << matTag; 
	    opserr << "LayeredShell section: " << tag << endln;
	    return 0;
	}

	if (OPS_GetDoubleInput(&numdata, &h) < 0) {
	    opserr << "WARNING invalid h" << endln;
	    opserr << "LayeredShell section: " << tag << endln;	    	    
	    return 0;
	}
	
	if (h < 0) {
	    opserr << "WARNING invalid h" << endln;
	    opserr << "PlateFiber section: " << tag << endln;	    	    
	    return 0;
	}
	
	thickness[iLayer] = h;
    }
      
    SectionForceDeformation* theSection = new LayeredShellFiberSection(tag, nLayers, thickness, theMats);
    if (thickness != 0) delete [] thickness;
    if (theMats != 0) delete [] theMats;

    return theSection;
}

//static vector and matrices
Vector  LayeredShellFiberSection::stressResultant(8) ;
Matrix  LayeredShellFiberSection::tangent(8,8) ;
ID      LayeredShellFiberSection::array(8) ;

//null constructor
LayeredShellFiberSection::LayeredShellFiberSection( ) : 
SectionForceDeformation( 0, SEC_TAG_LayeredShellFiberSection ), 
strainResultant(8), nLayers(0)
{

}

//full constructor
LayeredShellFiberSection::LayeredShellFiberSection(    
                                   int tag, 
                                   int iLayers, 
                                   double *thickness, 
                                   NDMaterial **fibers ) :
SectionForceDeformation( tag, SEC_TAG_LayeredShellFiberSection ),
strainResultant(8)
{
  this->nLayers = iLayers;
  sg = new double[iLayers];
  wg = new double[iLayers];
  theFibers = new NDMaterial*[iLayers];

  h = 0.0;
  int i;
  for ( i = 0; i < iLayers; i++ )
  {
    h = h + thickness[i];
    theFibers[i] = fibers[i]->getCopy( "PlateFiber" ) ;
    if (theFibers[i]==0) {
      opserr << "LayeredShellFiberSection::ERROR: Could Not return a PlateFiber Material: ";
      opserr << fibers[i]->getTag() << endln;
      exit(-1);
    }
  }

  for ( i = 0; i < iLayers; i++ ) wg[i] = 2.0 * thickness[i] / h;
  double currLoc = 0.0;
  double h1 = 1.0 / h;
  for ( i = 0; i < iLayers; i++ )
  {
    currLoc = currLoc + thickness[i];
    sg[i] = currLoc * h1 - 1.0;
    currLoc = currLoc + thickness[i];
  }
}

//destructor
LayeredShellFiberSection::~LayeredShellFiberSection( ) 
{ 
  int i ;
  if (sg != 0) delete sg;
  if (wg != 0) delete wg;
  if (theFibers != 0)
  {
    for ( i = 0; i < nLayers; i++ )
    {
      if (theFibers[i] != 0) delete theFibers[i] ;
    }
    delete [] theFibers;
  }
} 

//make a clone of this material
SectionForceDeformation  *LayeredShellFiberSection::getCopy( ) 
{
  LayeredShellFiberSection *clone = 0;   //new instance of this class
  
  double *thickness = new double[nLayers];
  if (thickness != 0)
  {
    for (int i = 0; i < nLayers; i++ ) 
      thickness[i] = 0.5 * wg[i] * h;

    
    clone = new LayeredShellFiberSection(this->getTag(),
					 nLayers,
					 thickness,
					 theFibers ) ; //make the copy
    delete [] thickness;
  }
  
  return clone ;
}



//send back order of strainResultant in vector form
int LayeredShellFiberSection::getOrder( ) const
{
  return 8 ;
}


//send back order of strainResultant in vector form
const ID& LayeredShellFiberSection::getType( ) 
{
  return array ;
}



//swap history variables
int LayeredShellFiberSection::commitState( ) 
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->commitState( ) ;

  return success ;
}

//revert to last saved state
int LayeredShellFiberSection::revertToLastCommit( )
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->revertToLastCommit( ) ;

  return success ;
}

//revert to start
int LayeredShellFiberSection::revertToStart( )
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->revertToStart( ) ;

  return success ;
}

//mass per unit area
double
LayeredShellFiberSection::getRho( )
{

  double weight ;

  double rhoH = 0.0 ;

  for ( int i = 0; i < nLayers; i++ ) {
    
    weight = ( 0.5*h ) * wg[i] ;

    rhoH += ( theFibers[i]->getRho() ) * weight ;

  }

  return rhoH ;

}

Response*
LayeredShellFiberSection::setResponse(const char **argv, int argc,
                                      OPS_Stream &output)
{
  const ID &type = this->getType();
  int typeSize = this->getOrder();
  
  Response *theResponse =0;

  output.tag("SectionOutput");
  output.attr("secType", this->getClassType());
  output.attr("secTag", this->getTag());

  // deformations
  if (strcmp(argv[0],"deformations") == 0 || strcmp(argv[0],"deformation") == 0) {
    output.tag("ResponseType","eps11");
    output.tag("ResponseType","eps22");
    output.tag("ResponseType","gamma12");
    output.tag("ResponseType","theta11");
    output.tag("ResponseType","theta22");
    output.tag("ResponseType","theta33");
    output.tag("ResponseType","gamma13");
    output.tag("ResponseType","gamma23");
    theResponse =  new MaterialResponse(this, 1, this->getSectionDeformation());
  // forces
  } else if (strcmp(argv[0],"forces") == 0 || strcmp(argv[0],"force") == 0) {
    output.tag("ResponseType","p11");
    output.tag("ResponseType","p22");
    output.tag("ResponseType","p12");
    output.tag("ResponseType","m11");
    output.tag("ResponseType","m22");
    output.tag("ResponseType","m12");
    output.tag("ResponseType","q1");
    output.tag("ResponseType","q2");
    theResponse =  new MaterialResponse(this, 2, this->getStressResultant());
  
  // force and deformation
  } else if (strcmp(argv[0],"forceAndDeformation") == 0) { 
    output.tag("ResponseType","eps11");
    output.tag("ResponseType","eps22");
    output.tag("ResponseType","gamma12");
    output.tag("ResponseType","theta11");
    output.tag("ResponseType","theta22");
    output.tag("ResponseType","theta33");
    output.tag("ResponseType","gamma13");
    output.tag("ResponseType","gamma23");
    output.tag("ResponseType","p11");
    output.tag("ResponseType","p22");
    output.tag("ResponseType","p12");
    output.tag("ResponseType","m11");
    output.tag("ResponseType","m22");
    output.tag("ResponseType","m12");
    output.tag("ResponseType","q1");
    output.tag("ResponseType","q2");
    theResponse =  new MaterialResponse(this, 4, Vector(2*this->getOrder()));
  }  
  else if (strcmp(argv[0],"fiber") == 0 || strcmp(argv[0],"Fiber") == 0) {
    if (argc < 3) {
      opserr << "LayeredShellFiberSection::setResponse() - need to specify more data\n";
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= nLayers) {
      
      output.tag("FiberOutput");
      output.attr("number",pointNum);
      output.attr("zLoc",0.5*h*sg[pointNum-1]);
      output.attr("thickness",0.5*h*wg[pointNum-1]);
      
      theResponse =  theFibers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();
    }
  }
  output.endTag(); // SectionOutput
  return theResponse;
}

int 
LayeredShellFiberSection::getResponse(int responseID, Information &secInfo)
{
  switch (responseID) {
  case 1:
    return secInfo.setVector(this->getSectionDeformation());
    
  case 2:
    return secInfo.setVector(this->getStressResultant());
    
  case 4: {
    Vector &theVec = *(secInfo.theVector);
    const Vector &e = this->getSectionDeformation();
    const Vector &s = this->getStressResultant();
    for (int i = 0; i < 8; i++) {
      theVec(i) = e(i);
      theVec(i+8) = s(i);
    }
    
    return secInfo.setVector(theVec);
  }
  default:
    return -1;
  }
}


//receive the strainResultant 
int LayeredShellFiberSection ::
setTrialSectionDeformation( const Vector &strainResultant_from_element)
{
  this->strainResultant = strainResultant_from_element ;

  static Vector strain(5) ;

  int success = 0 ;

  int i ;

  double z ;

  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i] ;
  
      strain(0) =  strainResultant(0)  - z*strainResultant(3) ;

      strain(1) =  strainResultant(1)  - z*strainResultant(4) ;

      strain(2) =  strainResultant(2)  - z*strainResultant(5) ;

      strain(3) =  strainResultant(6) ;

      strain(4) =  strainResultant(7) ;
  
      success += theFibers[i]->setTrialStrain( strain ) ;

  } //end for i

  return success ;
}


//send back the strainResultant
const Vector& LayeredShellFiberSection::getSectionDeformation( )
{
  return this->strainResultant ;
}


//send back the stressResultant 
const Vector&  LayeredShellFiberSection::getStressResultant( )
{

  static Vector stress(5) ;

  int i ;

  double z, weight ;

  stressResultant.Zero( ) ;

  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i] ;

      weight = ( 0.5*h ) * wg[i] ;

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
  
  } //end for i

   return this->stressResultant ;
}


//send back the tangent 
const Matrix&  LayeredShellFiberSection::getSectionTangent( )
{
  static Matrix dd(5,5) ;

//  static Matrix Aeps(5,8) ;

//  static Matrix Asig(8,5) ;

  int i ;

  double z, weight ;

  tangent.Zero( ) ;

  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i] ;

      weight = (0.5*h) * wg[i] ;

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
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14,    d15]
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24,    d25]
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34,    d35]
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14,  z*d15]
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24,  z*d25]
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34,  z*d35]
//[       d41,           d42,           d43,        -d41*z,        -d42*z,        -d43*z,    d44,    d45]
//[       d51,           d52,           d53,        -d51*z,        -d52*z,        -d53*z,    d54,    d55]
 
      //row 1
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14,    d15]
      tangent(0,0) +=     dd(0,0) ;
      tangent(0,1) +=     dd(0,1) ;
      tangent(0,2) +=     dd(0,2) ;      
      tangent(0,3) +=  -z*dd(0,0) ;      
      tangent(0,4) +=  -z*dd(0,1) ;
      tangent(0,5) +=  -z*dd(0,2) ;
      tangent(0,6) +=     dd(0,3) ;
      tangent(0,7) +=     dd(0,4) ;

      //row 2
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24,    d25]
      tangent(1,0) +=     dd(1,0) ;
      tangent(1,1) +=     dd(1,1) ;
      tangent(1,2) +=     dd(1,2) ;      
      tangent(1,3) +=  -z*dd(1,0) ;      
      tangent(1,4) +=  -z*dd(1,1) ;
      tangent(1,5) +=  -z*dd(1,2) ;
      tangent(1,6) +=     dd(1,3) ;
      tangent(1,7) +=     dd(1,4) ;

      //row 3
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34,    d35]
      tangent(2,0) +=     dd(2,0) ;
      tangent(2,1) +=     dd(2,1) ;
      tangent(2,2) +=     dd(2,2) ;      
      tangent(2,3) +=  -z*dd(2,0) ;      
      tangent(2,4) +=  -z*dd(2,1) ;
      tangent(2,5) +=  -z*dd(2,2) ;
      tangent(2,6) +=     dd(2,3) ;
      tangent(2,7) +=     dd(2,4) ;

      //row 4
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14,  z*d15]
      tangent(3,0) +=     z*dd(0,0) ;
      tangent(3,1) +=     z*dd(0,1) ;
      tangent(3,2) +=     z*dd(0,2) ;      
      tangent(3,3) +=  -z*z*dd(0,0) ;      
      tangent(3,4) +=  -z*z*dd(0,1) ;
      tangent(3,5) +=  -z*z*dd(0,2) ;
      tangent(3,6) +=     z*dd(0,3) ;
      tangent(3,7) +=     z*dd(0,4) ;

      //row 5
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24,  z*d25]
      tangent(4,0) +=     z*dd(1,0) ;
      tangent(4,1) +=     z*dd(1,1) ;
      tangent(4,2) +=     z*dd(1,2) ;      
      tangent(4,3) +=  -z*z*dd(1,0) ;      
      tangent(4,4) +=  -z*z*dd(1,1) ;
      tangent(4,5) +=  -z*z*dd(1,2) ;
      tangent(4,6) +=     z*dd(1,3) ;
      tangent(4,7) +=     z*dd(1,4) ;

      //row 6
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34,  z*d35]
      tangent(5,0) +=     z*dd(2,0) ;
      tangent(5,1) +=     z*dd(2,1) ;
      tangent(5,2) +=     z*dd(2,2) ;      
      tangent(5,3) +=  -z*z*dd(2,0) ;      
      tangent(5,4) +=  -z*z*dd(2,1) ;
      tangent(5,5) +=  -z*z*dd(2,2) ;
      tangent(5,6) +=     z*dd(2,3) ;
      tangent(5,7) +=     z*dd(2,4) ;

      //row 7
//[  d41,    d42,    d43, -d41*z, -d42*z, -d43*z,  d44,  d45]
      tangent(6,0) +=     dd(3,0) ;
      tangent(6,1) +=     dd(3,1) ;
      tangent(6,2) +=     dd(3,2) ;      
      tangent(6,3) +=  -z*dd(3,0) ;      
      tangent(6,4) +=  -z*dd(3,1) ;
      tangent(6,5) +=  -z*dd(3,2) ;
      tangent(6,6) +=     dd(3,3) ;
      tangent(6,7) +=     dd(3,4) ;

      //row 8 
//[  d51,    d52,    d53, -d51*z, -d52*z, -d53*z,  d54,  d55]
      tangent(7,0) +=     dd(4,0) ;
      tangent(7,1) +=     dd(4,1) ;
      tangent(7,2) +=     dd(4,2) ;      
      tangent(7,3) +=  -z*dd(4,0) ;      
      tangent(7,4) +=  -z*dd(4,1) ;
      tangent(7,5) +=  -z*dd(4,2) ;
      tangent(7,6) +=     dd(4,3) ;
      tangent(7,7) +=     dd(4,4) ;

  } //end for i

  return this->tangent ;
}


//print out data
void  LayeredShellFiberSection::Print( OPS_Stream &s, int flag )
{
    if (flag == OPS_PRINT_PRINTMODEL_SECTION || flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
        s << "LayeredShellFiber Section tag: " << this->getTag() << endln;
        s << "Total thickness h = " << h << endln;
        for (int i = 0; i < nLayers; i++) {
            s << "Layer " << i + 1 << ", thickness h = " << 0.5 * wg[i] * h << endln;
            theFibers[i]->Print(s, flag);
            s << endln;
        }
    }
    
    if (flag == OPS_PRINT_PRINTMODEL_JSON) {
        s << "\t\t\t{";
        s << "\"name\": \"" << this->getTag() << "\", ";
        s << "\"type\": \"LayeredShellFiberSection\", ";
        s << "\"totalThickness\": " << h << ", ";
        s << "\"fibers\": [\n";
        for (int i = 0; i < nLayers; i++) {
            s << "\t\t\t\t{\"layer\": " << i+1 << ", ";
            s << "\"thickness\": " << 0.5*wg[i]*h << ", ";
            s << "\"material\": \"" << theFibers[i]->getTag() << "\"";
            if (i < nLayers - 1)
                s << "},\n";
            else
                s << "}\n";
        }
        s << "\t\t\t]}";
    }
}

int 
LayeredShellFiberSection::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID iData(3);
  iData(0) = this->getTag();
  iData(1) = nLayers;

  res += theChannel.sendID(dataTag, commitTag, iData);
  if (res < 0) {
    opserr << "WARNING LayeredShellFiberSection::sendSelf() - " << this->getTag() << " failed to send data" << endln;
    return res;
  }
  
  if (nLayers > 0)
  {
    Vector vecData(2*nLayers+1);
    int i;
    for (i = 0; i < nLayers; i++) {
      vecData(i)         = sg[i];
      vecData(i+nLayers) = wg[i];
    }
    vecData(2*nLayers) = h;
    res += theChannel.sendVector(dataTag, commitTag, vecData);
    if (res < 0) {
      opserr << "WARNING LayeredShellFiberSection::sendSelf() - " << this->getTag() << " failed to send data" << endln;
      return res;
    }
    
    // Send the ids of its materials
    
    int matDbTag;
    ID idData(nLayers*2);
    for (i = 0; i < nLayers; i++) {
      idData(i) = theFibers[i]->getClassTag();
      matDbTag = theFibers[i]->getDbTag();
      // ensure that the material has a database tag
      if (matDbTag == 0) {
        matDbTag = theChannel.getDbTag();
        if (matDbTag != 0) theFibers[i]->setDbTag(matDbTag);
      }
      idData(i+nLayers) = matDbTag;
    }
    
    res += theChannel.sendID(dataTag, commitTag, idData);
    if (res < 0) {
      opserr << "WARNING LayeredShellFiberSection::sendSelf() - " << this->getTag() << " failed to send ID" << endln;
      return res;
    }
    
    // Finally, quad asks its material objects to send themselves
    for (i = 0; i < nLayers; i++) {
      res += theFibers[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "WARNING LayeredShellFiberSection::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
        return res;
      }
    }
  }

  return res;
}


int 
LayeredShellFiberSection::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID iData(3);
  res += theChannel.recvID(dataTag, commitTag, iData);

  if (res < 0) {
    opserr << "WARNING LayeredShellFiberSection::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
   return res;
  } 

  this->setTag(iData(0));
  
  int i;
  if (nLayers != iData(1))
  {
    nLayers = iData(1);
    if (sg != 0) delete sg;
    sg = new double[nLayers];
    if (wg != 0) delete sg;
    wg = new double[nLayers];
    if (theFibers !=0)
    {
      for ( i = 0; i < nLayers; i++ )
      {
        if (theFibers[i] != 0) delete theFibers[i] ;
      }
      delete [] theFibers;
    }
    theFibers = new NDMaterial*[nLayers];
  }

  if (nLayers > 0)
  {
    Vector vecData(2*nLayers+1);
    res += theChannel.recvVector(dataTag, commitTag, vecData);
    if (res < 0) {
    opserr << "WARNING LayeredShellFiberSection::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
    return res;
    }  
    for (i = 0; i < nLayers; i++) {
      sg[i] = vecData[i];
      wg[i] = vecData[i+nLayers];
    }
    h = vecData[2*nLayers];
    ID idData(nLayers*2);
    res += theChannel.recvID(dataTag, commitTag, idData);
    if (res < 0) {
      opserr << "WARNING LayeredShellFiberSection::recvSelf() - " << this->getTag() << " failed to receive ID" << endln;
      return res;
    }

    for (i = 0; i < nLayers; i++) {
      int matClassTag = idData(i);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theFibers[i]->getClassTag() != matClassTag) {
        if (theFibers[i] != 0) delete theFibers[i];
        theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
        if (theFibers[i] == 0) {
          opserr << "LayeredShellFiberSection::recvSelf() - " << 
            "Broker could not create NDMaterial of class type" << matClassTag << endln;
         return -1;
        }
      }
      theFibers[i]->setDbTag(idData(i+nLayers));
      // Receive the material
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LayeredShellFiberSection::recvSelf() - material " << 
          i << ", failed to recv itself" << endln;
        return res;
      }
    }
  }
    
  return res;
}
 


