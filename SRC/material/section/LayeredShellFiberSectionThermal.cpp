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
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/LayeredShellFiberSectionThermal.cpp,v $

// Yuli Huang (yulihuang@gmail.com) & Xinzheng Lu (luxz@tsinghua.edu.cn)
//
// Layered Shell Section
//
/* Ref: Lu X, Lu XZ, Guan H, Ye LP, Collapse simulation of reinforced 
concrete high-rise building induced by extreme earthquakes, 
Earthquake Engineering & Structural Dynamics, 2013, 42(5): 705-723*/
// Modified for SIF modelling by Liming Jiang [http://openseesforfire.github.io] 


#include <LayeredShellFiberSectionThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>

//static vector and matrices
Vector  LayeredShellFiberSectionThermal::stressResultant(8) ;
Matrix  LayeredShellFiberSectionThermal::tangent(8,8) ;
ID      LayeredShellFiberSectionThermal::array(8) ;
//const double  LayeredShellFiberSectionThermal::root56 = sqrt(5.0/6.0) ; //shear correction
const double  LayeredShellFiberSectionThermal::root56 = 1 ; //shear correction
//null constructor
LayeredShellFiberSectionThermal::LayeredShellFiberSectionThermal( ) : 
SectionForceDeformation( 0, SEC_TAG_LayeredShellFiberSectionThermal ), 
strainResultant(8), nLayers(0),countnGauss(0),AverageThermalMomentP(0), AverageThermalForceP(0),sT(0), ThermalElongation(0)
{

}

//full constructor
LayeredShellFiberSectionThermal::LayeredShellFiberSectionThermal(    
                                   int tag, 
                                   int iLayers, 
                                   double *thickness, 
                                   NDMaterial **fibers ) :
SectionForceDeformation( tag, SEC_TAG_LayeredShellFiberSectionThermal ),
strainResultant(8), countnGauss(0), AverageThermalMomentP(0), AverageThermalForceP(0),sT(0), ThermalElongation(0)
{
  this->nLayers = iLayers;
  sg = new double[iLayers];
  wg = new double[iLayers];
  theFibers = new NDMaterial*[iLayers];
   ThermalElongation = new double[iLayers];

  h = 0.0;
  int i;
  for ( i = 0; i < iLayers; i++ )
  {
    h = h + thickness[i];
    theFibers[i] = fibers[i]->getCopy( "PlateFiberThermal" ) ;
  }
  for ( i = 0; i < iLayers; i++ ) 
	  wg[i] = 2.0 * thickness[i] / h;

  double currLoc = 0.0;
  double h1 = 1.0 / h;
  for ( i = 0; i < iLayers; i++ )
  {
    currLoc = currLoc + thickness[i];
    sg[i] = currLoc * h1 - 1.0;
    currLoc = currLoc + thickness[i];
	ThermalElongation[i] =0.0;  //Added  by LMJ
  }

  sT =new Vector(2);
  sT->Zero();

}

//destructor
LayeredShellFiberSectionThermal::~LayeredShellFiberSectionThermal( ) 
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
SectionForceDeformation  *LayeredShellFiberSectionThermal::getCopy( ) 
{
  LayeredShellFiberSectionThermal *clone = 0;   //new instance of this class
  
  double *thickness = new double[nLayers];
  if (thickness != 0)
  {
    for (int i = 0; i < nLayers; i++ ) 
      thickness[i] = 0.5 * wg[i] * h;

    clone = new LayeredShellFiberSectionThermal( this->getTag(),
					  nLayers,
					  thickness,
					  theFibers ) ; //make the copy
    delete [] thickness;
  }
  return clone ;
}



//send back order of strainResultant in vector form
int LayeredShellFiberSectionThermal::getOrder( ) const
{
  return 8 ;
}


//send back order of strainResultant in vector form
const ID& LayeredShellFiberSectionThermal::getType( ) 
{
  return array ;
}



//swap history variables
int LayeredShellFiberSectionThermal::commitState( ) 
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->commitState( ) ;

  return success ;
}

//revert to last saved state
int LayeredShellFiberSectionThermal::revertToLastCommit( )
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->revertToLastCommit( ) ;

  return success ;
}

//revert to start
int LayeredShellFiberSectionThermal::revertToStart( )
{
  int success = 0 ;

  for (int i = 0; i < nLayers; i++ )
    success += theFibers[i]->revertToStart( ) ;

  return success ;
}

//mass per unit area
double
LayeredShellFiberSectionThermal::getRho( )
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
LayeredShellFiberSectionThermal::setResponse(const char **argv, int argc,
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
      opserr << "LayeredShellFiberSectionThermal::setResponse() - need to specify more data\n";
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
LayeredShellFiberSectionThermal::getResponse(int responseID, Information &secInfo)
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
int LayeredShellFiberSectionThermal ::
setTrialSectionDeformation( const Vector &strainResultant_from_element)
{
 
  this->strainResultant = strainResultant_from_element ;

  static Vector strain(6) ;

  int success = 0 ;

  int i ;

  double z ;

 #ifdef _DEBUG
	  if (strainResultant(8) ==111)
	  opserr<< "Sec strain  "<<strainResultant<<endln;
	#endif

  for ( i = 0; i < nLayers; i++ ) {

      z = ( 0.5*h ) * sg[i] ;
  
      strain(0) =  strainResultant(0)  - z*strainResultant(3)-ThermalElongation[i] ;

      strain(1) =  strainResultant(1)  - z*strainResultant(4)-ThermalElongation[i] ;

      strain(2) =  strainResultant(2)  - z*strainResultant(5) ;

      strain(3) =  root56*strainResultant(6) ;

      strain(4) =  root56*strainResultant(7) ;
    
	#ifdef _SDEBUG
	  if (strainResultant(8) ==111&&i==0){
	  opserr<<"LayeredShellSection z : "<<z<< "strain  "<<strain<<endln;
	  strain(5)=1110;
	  }
	  else 
		strain(5)=1115;
	#endif

//#ifdef _DEBUG
//opserr<<" Layer"<< i <<endln;
//#endif
      success += theFibers[i]->setTrialStrain( strain ) ;
	  #ifdef _SDEBUG
	  if (strainResultant(8) ==111&&i==0){
	  opserr<<"LayeredShellSection z : "<<z<<"stress  "<<theFibers[i]->getStress()<<endln;
	  opserr<<" tangent "<<theFibers[i]->getTangent()<<endln;
		  }
	#endif
     
  } //end for i

  return success ;
}


//send back the strainResultant
const Vector& LayeredShellFiberSectionThermal::getSectionDeformation( )
{
  return this->strainResultant ;
}


const Vector&
LayeredShellFiberSectionThermal::getTemperatureStress(const Vector& dataMixed)
{
   this->countnGauss = 0;             
             
  double* ThermalTangent = new double[nLayers];
  //double ThermalElongation[5];
  for (int i = 0; i < nLayers; i++) {
          ThermalTangent[i]=0;
          ThermalElongation[i]=0;
  }
  double FiberTemperature = 0 ;
  double tangent, elongation;
  double averageThermalForce =0.0;
  double averageThermalMoment =0.0;
  //double *thickness = new double[nLayers];

  for (int i = 0; i < nLayers; i++) {
	  
	double thickness = 0.5*h*wg[i];

    double yi = ( 0.5*h ) * sg[i];

	double tangent, elongation;

	FiberTemperature = this->determineFiberTemperature( dataMixed, yi);

	theFibers[i]->getThermalTangentAndElongation(FiberTemperature, tangent, elongation);

	ThermalTangent[i]=tangent;
	ThermalElongation[i]=elongation;
	averageThermalForce += elongation*thickness*tangent;
	averageThermalMoment+= yi*thickness*elongation*tangent;
	}

      (*sT)(0) = averageThermalForce - AverageThermalForceP;

      (*sT)(1) = averageThermalMoment - AverageThermalMomentP;
	  AverageThermalForceP = averageThermalForce;
	  AverageThermalMomentP = averageThermalMoment;
      return *sT;

}


//send back the stressResultant 
const Vector&  LayeredShellFiberSectionThermal::getStressResultant( )
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

     //modify shear 
   stressResultant(6) *= root56 ;  
   stressResultant(7) *= root56 ;


   return this->stressResultant ;
}


//send back the tangent 
const Matrix&  LayeredShellFiberSectionThermal::getSectionTangent( )
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
void  LayeredShellFiberSectionThermal::Print( OPS_Stream &s, int flag )
{
  s << "LayeredShellFiber Section tag: " << this->getTag() << endln ; 
  s << "Total thickness h = " << h << endln ;

  for (int i = 0; i < nLayers; i++) {
  s << "Layer " << i+1 << ", thickness h = " << 0.5 * wg[i] * h << endln;
    theFibers[i]->Print( s, flag ) ;
  s << endln;
  }

}

int 
LayeredShellFiberSectionThermal::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  int dataTag = this->getDbTag();

  static ID iData(3);
  iData(0) = this->getTag();
  iData(1) = nLayers;

  res += theChannel.sendID(dataTag, commitTag, iData);
  if (res < 0) {
    opserr << "WARNING LayeredShellFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send data" << endln;
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
      opserr << "WARNING LayeredShellFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send data" << endln;
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
      opserr << "WARNING LayeredShellFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send ID" << endln;
      return res;
    }
    
    // Finally, quad asks its material objects to send themselves
    for (i = 0; i < nLayers; i++) {
      res += theFibers[i]->sendSelf(commitTag, theChannel);
      if (res < 0) {
        opserr << "WARNING LayeredShellFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send its Material" << endln;
        return res;
      }
    }
  }

  return res;
}


int 
LayeredShellFiberSectionThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID iData(3);
  res += theChannel.recvID(dataTag, commitTag, iData);

  if (res < 0) {
    opserr << "WARNING LayeredShellFiberSectionThermal::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
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
    opserr << "WARNING LayeredShellFiberSectionThermal::recvSelf() - " << this->getTag() << " failed to receive data" << endln;
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
      opserr << "WARNING LayeredShellFiberSectionThermal::recvSelf() - " << this->getTag() << " failed to receive ID" << endln;
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
          opserr << "LayeredShellFiberSectionThermal::recvSelf() - " << 
            "Broker could not create NDMaterial of class type" << matClassTag << endln;
         return -1;
        }
      }
      theFibers[i]->setDbTag(idData(i+nLayers));
      // Receive the material
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
        opserr << "LayeredShellFiberSectionThermal::recvSelf() - material " << 
          i << ", failed to recv itself" << endln;
        return res;
      }
    }
  }
    
  return res;
}
 



double
LayeredShellFiberSectionThermal::determineFiberTemperature(const Vector& DataMixed, double fiberLoc) 
{
		double FiberTemperature = 0;

		double dataTempe[18]; 
		for (int i = 0; i < 18; i++) { 
			dataTempe[i] = DataMixed(i);
		}

		if (  fiberLoc <= dataTempe[1])
		{
			opserr <<"FiberSection2dThermal::setTrialSectionDeformationTemperature -- fiber loc is out of the section";
		}
		else if (fiberLoc <= dataTempe[3])
		{
			FiberTemperature = dataTempe[0] - (dataTempe[1] - fiberLoc) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
		}
		else if (   fiberLoc <= dataTempe[5] )
		{
			FiberTemperature = dataTempe[2] - (dataTempe[3] - fiberLoc) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
			
		}
		else if ( fiberLoc <= dataTempe[7] )
		{
			FiberTemperature = dataTempe[4] - (dataTempe[5] - fiberLoc) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
		}
		else if ( fiberLoc <= dataTempe[9] )
		{
			FiberTemperature = dataTempe[6] - (dataTempe[7] - fiberLoc) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
		
		}
		else if (fiberLoc <= dataTempe[11] )
		{
			FiberTemperature = dataTempe[8] - (dataTempe[9] - fiberLoc) * (dataTempe[8] - dataTempe[10])/(dataTempe[9] - dataTempe[11]);

		}
		else if (fiberLoc <= dataTempe[13] )
		{
			FiberTemperature = dataTempe[10] - (dataTempe[11] - fiberLoc) * (dataTempe[10] - dataTempe[12])/(dataTempe[11] - dataTempe[13]);

		}
		else if (fiberLoc <= dataTempe[15] )
		{
			FiberTemperature = dataTempe[12] - (dataTempe[13] - fiberLoc) * (dataTempe[12] - dataTempe[14])/(dataTempe[13] - dataTempe[15]);

		}
		else if ( fiberLoc <= dataTempe[17] )
		{
			FiberTemperature = dataTempe[14] - (dataTempe[15] - fiberLoc) * (dataTempe[14] - dataTempe[16])/(dataTempe[15] - dataTempe[17]);

		}
		else 
		{
			opserr <<"FiberSection2dThermal::setTrialSectionDeformation -- fiber loc is out of the section";
		}

		return FiberTemperature;
}


