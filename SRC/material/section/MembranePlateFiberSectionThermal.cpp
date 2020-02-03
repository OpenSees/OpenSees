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
// $Date: 2007/05/03 23:03:01 $
// $Source: /usr/local/cvs/OpenSees/SRC/material/section/MembranePlateFiberSectionThermal.cpp,v $

// Ed "C++" Love
//
//  Generic Plate Section with membrane
// Modified for SIF modelling by Jian Jiang,Liming Jiang [http://openseesforfire.github.io] 

//


#include <MembranePlateFiberSectionThermal.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MaterialResponse.h>
#include <Information.h>


//parameters
//const double MembranePlateFiberSectionThermal::root56 = 1 ; //shear correction
const double MembranePlateFiberSectionThermal::root56 = sqrt(5.0/6.0) ; //shear correction

//static vector and matrices
Vector  MembranePlateFiberSectionThermal::stressResultant(8) ;
Matrix  MembranePlateFiberSectionThermal::tangent(8,8) ;
ID      MembranePlateFiberSectionThermal::array(8) ;


const double  MembranePlateFiberSectionThermal::sg[] = { -1, 
						  -0.65465367, 
					           0, 
					           0.65465367, 
					           1 } ;
 
const double  MembranePlateFiberSectionThermal::wg[] = { 0.1, 
					          0.5444444444, 
						  0.7111111111, 
						  0.5444444444, 
						  0.1  };

/*      from Ham-O
        case 5:
         xi(0,0) = -1.;
         xi(1,0) = -0.65465367;
         xi(2,0) =  0.;
         xi(3,0) =  0.65465367;
         xi(4,0) =  1.;
      
         w(0) =  0.1;
         w(1) =  0.5444444444;
         w(2) =  0.7111111111;
         w(3) =  0.5444444444;
         w(4) =  0.1;
      break;
*/


//null constructor
MembranePlateFiberSectionThermal::MembranePlateFiberSectionThermal( ) : 
SectionForceDeformation( 0, SEC_TAG_MembranePlateFiberSectionThermal ), 
strainResultant(8) 
{ 
  for ( int i = 0; i < 5; i++ )
      theFibers[i] = 0 ;
   this->countnGauss = 0;
}



//full constructor
MembranePlateFiberSectionThermal::MembranePlateFiberSectionThermal(    
				   int tag, 
                                   double thickness, 
                                   NDMaterial &Afiber ) :
SectionForceDeformation( tag, SEC_TAG_MembranePlateFiberSectionThermal ),
strainResultant(8)
{
  this->h  = thickness ;

  int i ;
  for ( i = 0; i < 5; i++ )
      theFibers[i] = Afiber.getCopy( "PlateFiberThermal" ) ;
 //J.Jiang add 
   sT = new Vector(sTData,2);   
   sTData[0] = 0.0;             
   sTData[1] = 0.0;  
   ThermalElongation[0] = 0.0;
   ThermalElongation[1] = 0.0;
   ThermalElongation[2] = 0.0;
   ThermalElongation[3] = 0.0;
   ThermalElongation[4] = 0.0;
   this->countnGauss = 0;
   ThermalGradientShink = 0.0;
}



//destructor
MembranePlateFiberSectionThermal::~MembranePlateFiberSectionThermal( ) 
{ 
  int i ;
  for ( i = 0; i < 5; i++ )
     delete theFibers[i] ;
} 



//make a clone of this material
SectionForceDeformation  *MembranePlateFiberSectionThermal::getCopy( ) 
{
  MembranePlateFiberSectionThermal *clone ;   //new instance of this class

  clone = new MembranePlateFiberSectionThermal( this->getTag(), 
                                         this->h,
                                         *theFibers[0] ) ; //make the copy
  return clone ;
}



//send back order of strainResultant in vector form
int MembranePlateFiberSectionThermal::getOrder( ) const
{
  return 8 ;
}


//send back order of strainResultant in vector form
const ID& MembranePlateFiberSectionThermal::getType( ) 
{
  return array ;
}



//swap history variables
int MembranePlateFiberSectionThermal::commitState( ) 
{
  int success = 0 ;

  int i ;
  for ( i = 0; i < 5; i++ )
    success += theFibers[i]->commitState( ) ;

  return success ;
}



//revert to last saved state
int MembranePlateFiberSectionThermal::revertToLastCommit( )
{
  int success = 0 ;

  int i ;
  for ( i = 0; i < 5; i++ )
    success += theFibers[i]->revertToLastCommit( ) ;

  return success ;
}

//revert to start
int MembranePlateFiberSectionThermal::revertToStart( )
{
  int success = 0 ;

  int i ;
  for ( i = 0; i < 5; i++ )
    success += theFibers[i]->revertToStart( ) ;

  return success ;
}


//mass per unit area
double
MembranePlateFiberSectionThermal::getRho( )
{

  double weight ;

  double rhoH = 0.0 ;

  for ( int i = 0; i < 5; i++ ) {
    
    weight = ( 0.5*h ) * wg[i] ;

    rhoH += ( theFibers[i]->getRho() ) * weight ;

  }

  return rhoH ;

}


//receive the strainResultant 
int MembranePlateFiberSectionThermal ::
setTrialSectionDeformation( const Vector &strainResultant_from_element)
{
  this->strainResultant = strainResultant_from_element ;
  //J.Jiang to see strainresultant
  double seestrain[8];
  seestrain[0] = strainResultant(0);
  seestrain[1] = strainResultant(1);
  seestrain[2] = strainResultant(2);
  seestrain[3] = strainResultant(3);
  seestrain[4] = strainResultant(4);
  seestrain[5] = strainResultant(5);
  seestrain[6] = strainResultant(6);
  seestrain[7] = strainResultant(7);

  static Vector strain(5) ;

  int success = 0 ;

  int i ;

  double z ;

  for ( i = 0; i < 5; i++ ) {

      z = ( 0.5*h ) * sg[i] ;
  
      strain(0) =  strainResultant(0)  - z*strainResultant(3) - ThermalElongation[i];

      strain(1) =  strainResultant(1)  - z*strainResultant(4) - ThermalElongation[i];

	  //strain(0) =  strainResultant(0)  - z*strainResultant(3) - ThermalElongation[i] + ThermalGradientShink;

      //strain(1) =  strainResultant(1)  - z*strainResultant(4) - ThermalElongation[i] + ThermalGradientShink;

      strain(2) =  strainResultant(2)  - z*strainResultant(5) ;

      strain(3) =  root56*strainResultant(6) ;

      strain(4) =  root56*strainResultant(7) ;
	  	  //J.Jiang add to see strain
	  double strain0, strain1,strain2, strain3, strain4;
	  strain0=strain(0);
	  strain1=strain(1);
	  strain2=strain(2);
	  strain3=strain(3);
	  strain4=strain(4);
	//opserr<<"Section "<<i<<" Membrane STRAIN "<< strain;
      success += theFibers[i]->setTrialStrain( strain ) ;

  } //end for i
  

 

  this->countnGauss++;
  return success ;
}


//send back the strainResultant
const Vector& MembranePlateFiberSectionThermal::getSectionDeformation( )
{
  return this->strainResultant ;
}


//send back the stressResultant 
const Vector&  MembranePlateFiberSectionThermal::getStressResultant( )
{

  static Vector stress(5) ;

  int i ;

  double z, weight ;

  stressResultant.Zero( ) ;

  for ( i = 0; i < 5; i++ ) {

      z = ( 0.5*h ) * sg[i] ;

      weight = ( 0.5*h ) * wg[i] ;

      stress = theFibers[i]->getStress( ) ;
#ifdef _SDEBUG
	  opserr<<"MembranePlateFiberSection:resultantstress "<<endln<<stress;
#endif
  
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
const Vector&
MembranePlateFiberSectionThermal::getTemperatureStress(const Vector& dataMixed)
{
   this->countnGauss = 0;             
  double dataTempe[18];
  for (int i = 0; i < 18; i++) {
	  dataTempe[i] = dataMixed(i);
  }
             
  double ThermalTangent[5];
  //double ThermalElongation[5];
  for (int i = 0; i < 5; i++) {
          ThermalTangent[i]=0;
          ThermalElongation[i]=0;
  }
  double FiberTemperature = 0 ;
  double tangent, elongation;
  for (int i = 0; i < 5; i++) {

    double yi = ( 0.5*h ) * sg[i] ;;
    
	 
	if ( fabs(dataTempe[1]) <= 1e-10 && fabs(dataTempe[17]) <= 1e-10 ) 
	{
		FiberTemperature = 0;
	}
	else
	{
	//calculate the fiber tempe, T=T1-(Y-Y1)*(T1-T2)/(Y1-Y2)
	if (  yi < dataTempe[1]) 
	{
		opserr <<"MembranePlateFiberSectionThermal::setTrialSectionDeformationup -- fiber loc is out of the section";
	}
	else if (yi <= dataTempe[3])
	{
		FiberTemperature = dataTempe[0] - (dataTempe[1] - yi) * (dataTempe[0] - dataTempe[2])/(dataTempe[1] - dataTempe[3]);
	}
	else if (   yi <= dataTempe[5] )
	{
		FiberTemperature = dataTempe[2] - (dataTempe[3] - yi) * (dataTempe[2] - dataTempe[4])/(dataTempe[3] - dataTempe[5]);
	}
	else if ( yi <= dataTempe[7] )
	{
		FiberTemperature = dataTempe[4] - (dataTempe[5] - yi) * (dataTempe[4] - dataTempe[6])/(dataTempe[5] - dataTempe[7]);
	}
	else if ( yi <= dataTempe[9] )
	{
		FiberTemperature = dataTempe[6] - (dataTempe[7] - yi) * (dataTempe[6] - dataTempe[8])/(dataTempe[7] - dataTempe[9]);
	}
	else if (yi <= dataTempe[11] )
	{
		FiberTemperature = dataTempe[8] - (dataTempe[9] - yi) * (dataTempe[8] - dataTempe[10])/(dataTempe[9] - dataTempe[11]);
	}
	else if (yi <= dataTempe[13] )
	{
		FiberTemperature = dataTempe[10] - (dataTempe[11] - yi) * (dataTempe[10] - dataTempe[12])/(dataTempe[11] - dataTempe[13]);
	}
	else if (yi <= dataTempe[15] )
	{
		FiberTemperature = dataTempe[12] - (dataTempe[13] - yi) * (dataTempe[12] - dataTempe[14])/(dataTempe[13] - dataTempe[15]);
	}
	else if ( yi <= dataTempe[17] )
	{
		FiberTemperature = dataTempe[14] - (dataTempe[15] - yi) * (dataTempe[14] - dataTempe[16])/(dataTempe[15] - dataTempe[17]);
	}

	else 
	{
		//opserr <<"MembranePlateFiberSectionThermal::setTrialSectionDeformationdown -- fiber loc is out of the section";
	}
	}

    // determine material strain and set it
    double tangent, elongation;

   theFibers[i]->getThermalTangentAndElongation(FiberTemperature, tangent, elongation);

  ThermalTangent[i]=tangent;
  ThermalElongation[i]=elongation;
  //ThermalGradientShink = 1-sin(2*12e-6*(dataTempe[0]-dataTempe[16])/0.1/2)/(2*12e-6*(dataTempe[0]-dataTempe[16])/0.1/2);
  }
  //update yBar = Ai*Ei*yi/(Ai*E*) 
  //calculate centroid of section yBar for composite section,i.e. yBar is related to tangent E
  
      sTData[0] = fabs(ThermalTangent[0])*h*12e-6*fabs(dataTempe[0]+dataTempe[16])/2;
      //sTData[0] = fabs(ThermalTangent[0])*h*12e-6*fabs(dataTempe[0]);
      //sTData[1] = fabs(dataTempe[0]-dataTempe[16])/h*ThermalTangent[0]*h*h*h/12*12e-6;
	  //J.Jiang second try to get Mt
      sTData[1] = fabs(dataTempe[0]-dataTempe[16])*12e-6*fabs(ThermalTangent[0])*h*h/12/0.7;
      return *sT;
}

//send back the tangent 
const Matrix&  MembranePlateFiberSectionThermal::getSectionTangent( )
{
  static Matrix dd(5,5) ;

  static Matrix Aeps(5,8) ;

  static Matrix Asig(8,5) ;

  int i ;

  double z, weight ;

  tangent.Zero( ) ;

  for ( i = 0; i < 5; i++ ) {

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
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
 
      //row 1
//[      d11,           d12,           d13,        -z*d11,        -z*d12,        -z*d13,    d14*root56,    d15*root56]
      tangent(0,0) +=  dd(0,0) ;
      tangent(0,1) +=  dd(0,1) ;
      tangent(0,2) +=  dd(0,2) ;      
      tangent(0,3) +=  -z*dd(0,0) ;      
      tangent(0,4) +=  -z*dd(0,1) ;
      tangent(0,5) +=  -z*dd(0,2) ;
      tangent(0,6) +=  root56*dd(0,3) ;
      tangent(0,7) +=  root56*dd(0,4) ;

      //row 2
//[      d21,           d22,           d23,        -z*d21,        -z*d22,        -z*d23,    d24*root56,    d25*root56]
      tangent(1,0) +=  dd(1,0) ;
      tangent(1,1) +=  dd(1,1) ;
      tangent(1,2) +=  dd(1,2) ;      
      tangent(1,3) +=  -z*dd(1,0) ;      
      tangent(1,4) +=  -z*dd(1,1) ;
      tangent(1,5) +=  -z*dd(1,2) ;
      tangent(1,6) +=  root56*dd(1,3) ;
      tangent(1,7) +=  root56*dd(1,4) ;

      //row 3
//[      d31,           d32,           d33,        -z*d31,        -z*d32,        -z*d33,    d34*root56,    d35*root56]
      tangent(2,0) +=  dd(2,0) ;
      tangent(2,1) +=  dd(2,1) ;
      tangent(2,2) +=  dd(2,2) ;      
      tangent(2,3) +=  -z*dd(2,0) ;      
      tangent(2,4) +=  -z*dd(2,1) ;
      tangent(2,5) +=  -z*dd(2,2) ;
      tangent(2,6) +=  root56*dd(2,3) ;
      tangent(2,7) +=  root56*dd(2,4) ;

      //row 4
//[     z*d11,         z*d12,         z*d13,      -z^2*d11,      -z^2*d12,      -z^2*d13,  z*d14*root56,  z*d15*root56]
      tangent(3,0) +=  z*dd(0,0) ;
      tangent(3,1) +=  z*dd(0,1) ;
      tangent(3,2) +=  z*dd(0,2) ;      
      tangent(3,3) +=  -z*z*dd(0,0) ;      
      tangent(3,4) +=  -z*z*dd(0,1) ;
      tangent(3,5) +=  -z*z*dd(0,2) ;
      tangent(3,6) +=  z*root56*dd(0,3) ;
      tangent(3,7) +=  z*root56*dd(0,4) ;

      //row 5
//[     z*d21,         z*d22,         z*d23,      -z^2*d21,      -z^2*d22,      -z^2*d23,  z*d24*root56,  z*d25*root56]
      tangent(4,0) +=  z*dd(1,0) ;
      tangent(4,1) +=  z*dd(1,1) ;
      tangent(4,2) +=  z*dd(1,2) ;      
      tangent(4,3) +=  -z*z*dd(1,0) ;      
      tangent(4,4) +=  -z*z*dd(1,1) ;
      tangent(4,5) +=  -z*z*dd(1,2) ;
      tangent(4,6) +=  z*root56*dd(1,3) ;
      tangent(4,7) +=  z*root56*dd(1,4) ;

      //row 6
//[     z*d31,         z*d32,         z*d33,      -z^2*d31,      -z^2*d32,      -z^2*d33,  z*d34*root56,  z*d35*root56]
      tangent(5,0) +=  z*dd(2,0) ;
      tangent(5,1) +=  z*dd(2,1) ;
      tangent(5,2) +=  z*dd(2,2) ;      
      tangent(5,3) +=  -z*z*dd(2,0) ;      
      tangent(5,4) +=  -z*z*dd(2,1) ;
      tangent(5,5) +=  -z*z*dd(2,2) ;
      tangent(5,6) +=  z*root56*dd(2,3) ;
      tangent(5,7) +=  z*root56*dd(2,4) ;

      //row 7
//[  root56*d41,    root56*d42,    root56*d43, -root56*d41*z, -root56*d42*z, -root56*d43*z,  root56^2*d44,  root56^2*d45]
      tangent(6,0) +=  root56*dd(3,0) ;
      tangent(6,1) +=  root56*dd(3,1) ;
      tangent(6,2) +=  root56*dd(3,2) ;      
      tangent(6,3) +=  -root56*z*dd(3,0) ;      
      tangent(6,4) +=  -root56*z*dd(3,1) ;
      tangent(6,5) +=  -root56*z*dd(3,2) ;
      tangent(6,6) +=  root56*root56*dd(3,3) ;
      tangent(6,7) +=  root56*root56*dd(3,4) ;

      //row 8 
//[  root56*d51,    root56*d52,    root56*d53, -root56*d51*z, -root56*d52*z, -root56*d53*z,  root56^2*d54,  root56^2*d55]
      tangent(7,0) +=  root56*dd(4,0) ;
      tangent(7,1) +=  root56*dd(4,1) ;
      tangent(7,2) +=  root56*dd(4,2) ;      
      tangent(7,3) +=  -root56*z*dd(4,0) ;      
      tangent(7,4) +=  -root56*z*dd(4,1) ;
      tangent(7,5) +=  -root56*z*dd(4,2) ;
      tangent(7,6) +=  root56*root56*dd(4,3) ;
      tangent(7,7) +=  root56*root56*dd(4,4) ;

  } //end for i


  //opserr<<"membrane tangent"<<tangent;
  return this->tangent ;
}


//print out data
void  MembranePlateFiberSectionThermal::Print( OPS_Stream &s, int flag )
{
  s << "MembranePlateFiberSectionThermal: \n " ;
  s <<  "  Thickness h = "        <<  h  <<  endln ;

  for (int i = 0; i < 5; i++) {
    theFibers[i]->Print( s, flag ) ;
  }

  return ;
}

int 
MembranePlateFiberSectionThermal::sendSelf(int commitTag, Channel &theChannel) 
{
  int res = 0;

  // note: we don't check for dataTag == 0 for Element
  // objects as that is taken care of in a commit by the Domain
  // object - don't want to have to do the check if sending data
  int dataTag = this->getDbTag();
  

  // Now quad sends the ids of its materials
  int matDbTag;
  
  static ID idData(11);
  
  int i;
  for (i = 0; i < 5; i++) {
    idData(i) = theFibers[i]->getClassTag();
    matDbTag = theFibers[i]->getDbTag();
    // NOTE: we do have to ensure that the material has a database
    // tag if we are sending to a database channel.
    if (matDbTag == 0) {
      matDbTag = theChannel.getDbTag();
			if (matDbTag != 0)
			  theFibers[i]->setDbTag(matDbTag);
    }
    idData(i+5) = matDbTag;
  }
  
  idData(10) = this->getTag();

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING MembranePlateFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send ID\n";
			    
    return res;
  }

  // Finally, quad asks its material objects to send themselves
  for (i = 0; i < 5; i++) {
    res += theFibers[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING MembranePlateFiberSectionThermal::sendSelf() - " << this->getTag() << " failed to send its Material\n";
      return res;
    }
  }

  return res;
}


int 
MembranePlateFiberSectionThermal::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  
  int dataTag = this->getDbTag();

  static ID idData(11);
  // Quad now receives the tags of its four external nodes
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING MembranePlateFiberSectionThermal::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(10));

  int i;

  if (theFibers[0] == 0) {
    for (i = 0; i < 5; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+5);
      // Allocate new material with the sent class tag
      theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
      if (theFibers[i] == 0) {
	opserr << "MembranePlateFiberSectionThermal::recvSelf() - " <<
	  "Broker could not create NDMaterial of class type " << matClassTag << endln;
	return -1;
      }
      // Now receive materials into the newly allocated space
      theFibers[i]->setDbTag(matDbTag);
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "MembranePlateFiber::recvSelf() - material " << i << "failed to recv itself\n";
	  
	return res;
      }
    }
  }
  // Number of materials is the same, receive materials into current space
  else {
    for (i = 0; i < 5; i++) {
      int matClassTag = idData(i);
      int matDbTag = idData(i+5);
      // Check that material is of the right type; if not,
      // delete it and create a new one of the right type
      if (theFibers[i]->getClassTag() != matClassTag) {
	delete theFibers[i];
	theFibers[i] = theBroker.getNewNDMaterial(matClassTag);
	if (theFibers[i] == 0) {
	  opserr << "MembranePlateFiberSectionThermal::recvSelf() - " << 
	    "Broker could not create NDMaterial of class type" << matClassTag << endln;
	  exit(-1);
	}
      }
      // Receive the material
      theFibers[i]->setDbTag(matDbTag);
      res += theFibers[i]->recvSelf(commitTag, theChannel, theBroker);
      if (res < 0) {
	opserr << "MembranePlateFiberSectionThermal::recvSelf() - material " << 
	  i << ", failed to recv itself\n";
	return res;
      }
    }
  }

  return res;
}
 


Response*
MembranePlateFiberSectionThermal::setResponse(const char **argv, int argc,
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
      opserr << "MembranePlateFiberSectionThermal::setResponse() - need to specify more data\n";
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= 5) {
      
      output.tag("FiberOutput");
      output.attr("number",pointNum);
      
      theResponse =  theFibers[pointNum-1]->setResponse(&argv[2], argc-2, output);
      
      output.endTag();
    }
  }
  output.endTag(); // SectionOutput
  return theResponse;
}

int 
MembranePlateFiberSectionThermal::getResponse(int responseID, Information &secInfo)
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

