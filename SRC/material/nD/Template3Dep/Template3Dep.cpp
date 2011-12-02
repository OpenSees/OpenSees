/*
################################################################################
# COPYRIGHT (C):     :-))                                                      #
# PROJECT:           Object Oriented Finite Element Program                    #
# PURPOSE:           General platform for elaso-plastic constitutive model     #
#                    implementation                                            #
# CLASS:             Template3Dep (the base class for all material point)     #
#                                                                              #
# VERSION:                                                                     #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
# TARGET OS:         DOS || UNIX || . . .                                      #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                               #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                               #
#                                                                              #
#                                                                              #
# DATE:              08-03-2000                                                #
# UPDATE HISTORY:    09-12-2000                                                #
#                                                                              #
#                                                                              #
#                                                                              #
#                                                                              #
# SHORT EXPLANATION: This file contains the class implementation for           #
#                    Template3Dep.                                             #
#                                                                              #
################################################################################
*/

#ifndef Template3Dep_CPP
#define Template3Dep_CPP

#define ITMAX 30
#define MAX_STEP_COUNT 30
#define	NUM_OF_SUB_INCR 30
#define KK 1000.0  //conversion between Pa and kPa, or N and kN 1 - kPa, kN; 1000 - Pa, N
//#define po 100.0 //Reference pressure Pa
//#include <string.h>

#include "Template3Dep.h"


//================================================================================
// Constructor
//================================================================================

Template3Dep::Template3Dep( int tag                       ,
                            YieldSurface     *YS_   ,         
                            PotentialSurface *PS_   , 
		  	    EPState          *EPS_  ,
	       	     	    EvolutionLaw_S   *ELS1_ , 
	       	     	    EvolutionLaw_S   *ELS2_ , 
	       	     	    EvolutionLaw_S   *ELS3_ , 
	       	     	    EvolutionLaw_S   *ELS4_ , 
	       	     	    EvolutionLaw_T   *ELT1_ ,
	       	     	    EvolutionLaw_T   *ELT2_ ,
	       	     	    EvolutionLaw_T   *ELT3_ ,
	       	     	    EvolutionLaw_T   *ELT4_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{      	      
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ ) 
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ ) 
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;

    if ( ELS3_ ) 
       ELS3 = ELS3_->newObj();
    else
       ELS3 = 0;

    if ( ELS4_ ) 
       ELS4 = ELS4_->newObj();
    else
       ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    if ( ELT2_ ) 
       ELT2 = ELT2_->newObj();
    else
       ELT2 = 0;

    if ( ELT3_ ) 
       ELT3 = ELT3_->newObj();
    else
       ELT3 = 0;

    if ( ELT4_ ) 
       ELT4 = ELT4_->newObj();
    else
       ELT4 = 0;

    //Initialze Eep using E-elastic
    tensor E  = ElasticStiffnessTensor();
    EPS->setEep( E );

}

//================================================================================
// Constructor 0
//================================================================================
Template3Dep::Template3Dep( int tag                     ,
                            YieldSurface     *YS_ ,        
                            PotentialSurface *PS_ ,
              	            EPState          *EPS_)
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    // Evolution laws
    ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}


//================================================================================
// Constructor 1
//================================================================================

Template3Dep::Template3Dep(   int tag                     ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	              EPState          *EPS_,
	       	              EvolutionLaw_S   *ELS1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{

    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}


//================================================================================
// Constructor 2
//================================================================================

Template3Dep::Template3Dep(   int tag                     ,
                              YieldSurface     *YS_ ,
                              PotentialSurface *PS_ ,
              	              EPState          *EPS_,
	       	              EvolutionLaw_T   *ELT1_ ) 
:NDMaterial(tag, ND_TAG_Template3Dep)
{

    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    // Evolution laws
    ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;
    
    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 3
//================================================================================

Template3Dep::Template3Dep(   int tag                     ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
	       	     	      EvolutionLaw_T   *ELT1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    // Evolution laws
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;
    ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 4
// Two scalar evolution laws and one tensorial evolution law are provided!
//================================================================================

Template3Dep::Template3Dep(   int tag                     ,
                              YieldSurface     *YS_ ,
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
     	       	              EvolutionLaw_S   *ELS2_, 
	       	     	      EvolutionLaw_T   *ELT1_ )
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ )
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;

    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;
    
    ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 5
// Two scalar evolution laws and two tensorial evolution law are provided!
//================================================================================

Template3Dep::Template3Dep(   int tag                     ,
                              YieldSurface     *YS_ ,        
                              PotentialSurface *PS_ ,
              	    	      EPState          *EPS_,
	       	     	      EvolutionLaw_S   *ELS1_, 
     	       	              EvolutionLaw_S   *ELS2_, 
	       	     	      EvolutionLaw_T   *ELT1_,
			      EvolutionLaw_T   *ELT2_)
:NDMaterial(tag, ND_TAG_Template3Dep)
{
    
    if ( YS_ ) 
       YS = YS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
       
    if ( PS_ ) 
       PS = PS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }

    if ( EPS_ ) 
       EPS = EPS_->newObj();
    else
    {
      opserr << "Template3Dep:: Template3Dep failed to construct the template3Dep\n";
       exit(-1);
    }
    
    if ( ELS1_ )
       ELS1 = ELS1_->newObj();
    else
       ELS1 = 0;

    if ( ELS2_ )
       ELS2 = ELS2_->newObj();
    else
       ELS2 = 0;
    ELS3 = 0;
    ELS4 = 0;

    if ( ELT1_ ) 
       ELT1 = ELT1_->newObj();
    else
       ELT1 = 0;

    if ( ELT2_ ) 
       ELT2 = ELT2_->newObj();
    else
       ELT2 = 0;
    ELT3 = 0;
    ELT4 = 0;
}

//================================================================================
// Constructor 6
// NO parameter is provided
//================================================================================

Template3Dep::Template3Dep()
:NDMaterial(0, ND_TAG_Template3Dep),ELS1(0), ELS2(0),ELS3(0), ELS4(0),
 ELT1(0), ELT2(0), ELT3(0), ELT4(0) 
{
    //YS = new DPYieldSurface();
    //PS = new DPPotentialSurface();
    //EPS = new EPState();
    YS  = 0;
    PS  = 0;
    EPS = 0;
}


//================================================================================
// Destructor 
//================================================================================

Template3Dep::~Template3Dep()
{   

    if (YS) 
       delete YS;

    if (PS) 
       delete PS;

     if (EPS) 
       delete EPS;

     if (ELS1) 
       delete ELS1;

     if (ELS2) 
       delete ELS2;
     
     if (ELS3) 
       delete ELS3;
     
     if (ELS4) 
       delete ELS4;

     if (ELT1) 
       delete ELT1;

     if (ELT2) 
       delete ELT2;
     
     if (ELT3) 
       delete ELT3;
     
     if (ELT4) 
       delete ELT4;
     

}

////================================================================================
////copy constructor
////================================================================================
//Template3Dep::Template3Dep(const  Template3Dep & rval) {   
//
//    YS = rval.YS->newObj();
//    PS = rval.PS->newObj();
//    EPS = rval.EPS->newObj();
//    
//    // Scalar Evolution Laws
//    if ( rval.getELS1() ) 
//       ELS1  = rval.getELS1()->newObj();
//    else
//       ELS1 = 0;
//
//    if ( !rval.getELS2() ) 
//       ELS2 = 0;
//    else
//       ELS2  = rval.getELS2()->newObj();
//    
//    if ( !rval.getELS3() ) 
//       ELS3 = 0;
//    else
//       ELS3  = rval.getELS3()->newObj();
//    
//    if ( !rval.getELS4() ) 
//       ELS4 = 0;
//    else
//       ELS4  = rval.getELS4()->newObj();
//    
//    // Tensorial Evolution Laws
//    if ( rval.getELT1() ) 
//       ELT1  = rval.getELT1()->newObj();
//    else
//       ELT1 = 0;
//
//    if ( !rval.getELT2() ) 
//       ELT2 = 0;
//    else
//       ELT2  = rval.getELT2()->newObj();
//
//    if ( !rval.getELT3() ) 
//       ELT3 = 0;
//    else
//       ELT3  = rval.getELT3()->newObj();
//
//    if ( !rval.getELT4() ) 
//       ELT4 = 0;
//    else
//       ELT4  = rval.getELT4()->newObj();
//
//}    	   
//    	   



//================================================================================
// Routine used to generate elastic compliance tensor D for this material point
//================================================================================
tensor Template3Dep::ElasticComplianceTensor(void) const
{
    tensor ret( 4, def_dim_4, 0.0 );

    int eflag = (this->EPS)->getElasticflag();    
    double Ey = this->EPS->getEo();
    double nuP =this->EPS->getnu();

    // Pressure dependent Isotropic
    if ( eflag == 1) 
      {
      //opserr << " Eo= " << Ey;
        if (Ey == 0)  
          {
            opserr << endln << "Ey = 0! Can't give you D!!" << endln;
            exit(1);
          }
      
      //Update E 
      stresstensor stc = (this->EPS)->getStress();
      double p = stc.p_hydrostatic();
      double po = EPS->getpo();
      
      //opserr << " p = " <<  p;
      
      //double po = 100.0; //kPa
      if (p <= 0.5000*KK) p = 0.500*KK;
      Ey = Ey * pow(p/po/KK, 0.5); //0.5
      //cerr << " Ec = " << Ey << endlnn;
      
      // Kronecker delta tensor
      tensor I2("I", 2, def_dim_2);
      
      tensor I_ijkl = I2("ij")*I2("kl");
      //I_ijkl.null_indices();
      tensor I_ikjl = I_ijkl.transpose0110();
      tensor I_iljk = I_ijkl.transpose0111();
      tensor I4s = (I_ikjl+I_iljk)*0.5;
      
      // Building compliance tensor
      ret = (I_ijkl * (-nuP/Ey)) + (I4s * ((1.0+nuP)/Ey));
      }

    // Pressure dependent Anisotropic
    else if (eflag == 3) 
      {
      // Form compliance matrix D
      double E = Ey;
      double nu = nuP;
      double Ev = this->EPS->getEv();
      double Ghv = this->EPS->getGhv();
      double nuhv = this->EPS->getnuhv();

      double A = 1.0/E;
      double B = 1.0/Ev;
      //opserr << "A " << A << " B " << B;

      Matrix D(6,6);
      D(0,0) = D(1,1) = A;
      D(2,2) = B;
      D(0,1) = D(1,0) = -nu*A;
      D(0,2) = D(2,0) = D(1,2) = D(2,1) = -nuhv*B;
      D(3,3) = (1.0+nu)*A;
      D(4,4) = D(5,5) = 0.5/Ghv;
      //opserr << " C " << D;

      //Convert Matric D to 4th order Elastic constants tensor ret;
      ret.val(1,1,1,1) = D(0,0); ret.val(1,1,2,2) = D(0,1); ret.val(1,1,3,3) = D(0,2); // --> Sigma_xx
      ret.val(1,2,1,2) = D(3,3); ret.val(1,2,2,1) = D(3,3); // --> Sigma_xy
      ret.val(1,3,1,3) = D(4,4); ret.val(1,3,3,1) = D(4,4); // --> Sigma_xz
      
      ret.val(2,1,1,2) = D(3,3); ret.val(2,1,2,1) = D(3,3); // --> Sigma_yx
      ret.val(2,2,1,1) = D(1,0); ret.val(2,2,2,2) = D(1,1); ret.val(2,2,3,3) = D(1,2); // --> Sigma_yy
      ret.val(2,3,2,3) = D(5,5); ret.val(2,3,3,2) = D(5,5); // --> Sigma_yz
      
      ret.val(3,1,1,3) = D(4,4); ret.val(3,1,3,1) = D(4,4); // --> Sigma_zx
      ret.val(3,2,2,3) = D(5,5); ret.val(3,2,3,2) = D(5,5); // --> Sigma_zy
      ret.val(3,3,1,1) = D(2,0); ret.val(3,3,2,2) = D(2,1); ret.val(3,3,3,3) = D(2,2); // --> Sigma_zz
      
      }
    else 
      {
      opserr << "Template3Dep::ElasticComplianceTensor failed to create elastic compliance tensor. Eflag must be 1 or 3: " <<  eflag << endln;
      exit(-1);
    
      }
            
    return ret;
}
  

//================================================================================
// Routine used to generate elastic stiffness tensor D for this material point
//================================================================================
tensor Template3Dep::ElasticStiffnessTensor(void) const
  {
    tensor ret( 4, def_dim_4, 0.0 );
	         	           
    int eflag = (this->EPS)->getElasticflag();
    
    double Ey = this->EPS->getEo();
    double nu =this->EPS->getnu();

    // Pressure dependent Isotropic
    if ( eflag == 1) {
       
       //Update E 
       stresstensor stc = (this->EPS)->getStress();
       double p = stc.p_hydrostatic();
       double po = EPS->getpo();
       //cerr << " p = " <<  p;
       
       //double po = 100.0; //kPa
       if (p <= 0.5000*KK) p = 0.5000*KK;
       double E = Ey * pow(p/po/KK, 0.5);//0.5
       //cerr << " Eo = " << Ey ;
       //cerr << " Ec = " << E << endlnn;
       
       				       
       // Kronecker delta tensor
       tensor I2("I", 2, def_dim_2);
       
       tensor I_ijkl = I2("ij")*I2("kl");
       
       
       //I_ijkl.null_indices();
       tensor I_ikjl = I_ijkl.transpose0110();
       tensor I_iljk = I_ijkl.transpose0111();
       tensor I4s = (I_ikjl+I_iljk)*0.5;
       
       //double x = I4s.trace();
       //opserr << "xxxxxx " << x << endlnn;
       
       //I4s.null_indices();
       
       // Building elasticity tensor
       ret = I_ijkl*( E*nu / ( (1.0+nu)*(1.0 - 2.0*nu) ) ) + I4s*( E / (1.0 + nu) );
       //ret.null_indices();
    }
    else if ( eflag == 3) { 
      // Form compliance matrix D
      double E = Ey;
      double Ev = this->EPS->getEv();
      double Ghv = this->EPS->getGhv();
      double nuhv = this->EPS->getnuhv();

      double A = 1.0/E;
      double B = 1.0/Ev;
      //opserr << "A " << A << " B " << B;

      Matrix D(6,6);
      D(0,0) = D(1,1) = A;
      D(2,2) = B;
      D(0,1) = D(1,0) = -nu*A;
      D(0,2) = D(2,0) = D(1,2) = D(2,1) = -nuhv*B;
      D(3,3) = (1.0+nu)*A;
      D(4,4) = D(5,5) = 0.5/Ghv;
      //opserr << " C " << D;

      // Do invertion once to get Elastic matrix D
      D.Invert( D );

      //Convert Matric D to 4th order Elastic constants tensor ret;
      ret.val(1,1,1,1) = D(0,0); ret.val(1,1,2,2) = D(0,1); ret.val(1,1,3,3) = D(0,2); // --> Sigma_xx
      ret.val(1,2,1,2) = D(3,3); ret.val(1,2,2,1) = D(3,3); // --> Sigma_xy
      ret.val(1,3,1,3) = D(4,4); ret.val(1,3,3,1) = D(4,4); // --> Sigma_xz
      
      ret.val(2,1,1,2) = D(3,3); ret.val(2,1,2,1) = D(3,3); // --> Sigma_yx
      ret.val(2,2,1,1) = D(1,0); ret.val(2,2,2,2) = D(1,1); ret.val(2,2,3,3) = D(1,2); // --> Sigma_yy
      ret.val(2,3,2,3) = D(5,5); ret.val(2,3,3,2) = D(5,5); // --> Sigma_yz
      
      ret.val(3,1,1,3) = D(4,4); ret.val(3,1,3,1) = D(4,4); // --> Sigma_zx
      ret.val(3,2,2,3) = D(5,5); ret.val(3,2,3,2) = D(5,5); // --> Sigma_zy
      ret.val(3,3,1,1) = D(2,0); ret.val(3,3,2,2) = D(2,1); ret.val(3,3,3,3) = D(2,2); // --> Sigma_zz
    
    }

    return ret;


}

//================================================================================
int Template3Dep::setTrialStrain(const Vector &v)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrain(const Vector &v, const Vector &r)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrainIncr(const Vector &v)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrainIncr(const Vector &v, const Vector &r)
{
//================================================================================
    // Not yet implemented
    return 0;
}

//================================================================================
const Matrix& Template3Dep::getTangent(void)
{
    // Not yet implemented
    Matrix *M = new Matrix();
    return *M;
}

const Matrix& Template3Dep::getInitialTangent(void)
{
  return this->getInitialTangent();
}

//================================================================================
const Vector& Template3Dep::getStress(void)
{
    // Not yet implemented
    Vector *V = new Vector();
    return *V;
}

//================================================================================
const Vector& Template3Dep::getStrain(void)
{
    // Not yet implemented
    Vector *V = new Vector();
    return *V;
}

//================================================================================
// what is the trial strain? Initial strain?
int Template3Dep::setTrialStrain(const Tensor &v)
{
    EPS->setStrain(v);
    return 0;
}


//================================================================================
int Template3Dep::setTrialStrain(const Tensor &v, const Tensor &r)
{
    EPS->setStrain(v);
    return 0;
}

//================================================================================

int Template3Dep::setTrialStrainIncr(const Tensor &v)
{
    
    //opserr << "\nBE: " << endlnn;
    EPState StartEPS( *(this->getEPS()) );
    stresstensor start_stress = StartEPS.getStress();
    //opserr << "start_stress 0 " << start_stress;
    
    //EPState tmp_EPS = BackwardEulerEPState(v);
    //if ( tmp_EPS.getConverged() ) {
    //     //setTrialEPS( tmp_EPS );
    //     setEPS( tmp_EPS );
    //	 //double p = (tmp_EPS.getStress()).p_hydrostatic();
    //	 //double ec = (tmp_EPS.getec()) - (tmp_EPS.getLam()) * log( p / (tmp_EPS.getpo()) );
    //	 //double	st = (tmp_EPS.getStrain()).Iinvariant1();
    //	 //double pl_s = (tmp_EPS.getPlasticStrain()).Iinvariant1();
    //	 //double dpl_s = (tmp_EPS.getdPlasticStrain()).Iinvariant1();
    //	 //cerr << "ec " << ec << " e " << tmp_EPS.gete() << " psi " << (tmp_EPS.gete() - ec) << " strain " << st << " t_pl " << pl_s << " d_pl " << dpl_s << "\n";
    //     return 0;
    //}
    
    //opserr << endlnn;    
    //setEPS( StartEPS );
    //int number_of_subincrements = 5;
    //Cascading subdividing in case that some incr_step is too big
    
    //int loop = 0;
    //while ( !tmp_EPS.getConverged()  && (loop < 1) ) {
    //
    //   setEPS( StartEPS );
    //   EPState startEPS( *(this->getEPS()) );
    //   stresstensor start_stress = startEPS.getStress();
    //   opserr << " Step Start Stress:" << start_stress << endln;
    //
    //   loop += 1;
    //   opserr << "\n "<< loop << "th Sub-BE, total subdivision: " << 10*loop*NUM_OF_SUB_INCR << endln;
    //   tmp_EPS = BESubIncrementation(v, 10*loop*NUM_OF_SUB_INCR);
    //   if ( tmp_EPS.getConverged() ) {
    //       //setTrialEPS( tmp_EPS );
    //       setEPS( tmp_EPS );
    //       return 0;
    //   }
    //   //else {
    //   //    opserr << "\n2nd Sub BE: " << 3*NUM_OF_SUB_INCR << endln;
    //   //	   tmp_EPS = BESubIncrementation(v, 3*NUM_OF_SUB_INCR);
    //   //	   
    //   //    if ( tmp_EPS.getConverged() ) {
    //   //       setEPS( tmp_EPS );
    //   //	      return 0;
    //   //	   }
    //   //	   else
    //   //	      return 1;
    //   //}
    //}
       
    //// for testing MD model only for no BE
    //EPState tmp_EPS = FESubIncrementation(v, NUM_OF_SUB_INCR);
    EPState tmp_EPS = ForwardEulerEPState(v);
    setEPS( tmp_EPS );
    //setEPS( StartEPS );
    return 0;
}

//================================================================================
int Template3Dep::setTrialStrainIncr(const Tensor &v, const Tensor &r)
{
    EPS->setStrain(v + EPS->getStrain() );
    return 0;
}

//================================================================================
const tensor& Template3Dep::getTangentTensor(void)
{
    //Tensor Eep = EPS->getEep();
    return EPS->Eep;
}

//================================================================================
const stresstensor  Template3Dep::getStressTensor(void)
{
    //opserr << *EPS;
    //stresstensor tmp;
    //tmp =  EPS->getStress();
    //opserr << "EPS->getStress() " << EPS->getStress() << endlnn;
    
    //Something funny!!! happened here when returning EPS->getStress()
    // This function will return wrong Stress.
    stresstensor ret = EPS->getStress();
    return ret;
}


//================================================================================
const straintensor Template3Dep::getStrainTensor(void)
{
    return EPS->getStrain();
}

//================================================================================
const straintensor Template3Dep::getPlasticStrainTensor(void)
{
    return EPS->getPlasticStrain();
}

//================================================================================
double Template3Dep::getpsi(void)
{
    return EPS->getpsi();
}

//================================================================================
int Template3Dep::commitState(void)
{
    int err;
    err = getEPS()->commitState();
    return err;
}

//================================================================================
int Template3Dep::revertToLastCommit(void)
{
	int err;
	err = EPS->revertToLastCommit();
	return err;
}
//================================================================================
int Template3Dep::revertToStart(void)
{
	int err;
	err = EPS->revertToStart();
	return err;
}

////================================================================================
//// Get one copy of the NDMaterial model
////NDMaterial *Template3Dep::getCopy(void)
//{     
//      Template3Dep *t3d;
//      return t3d;
//}

//================================================================================
// Get one copy of the NDMaterial model
NDMaterial * Template3Dep::getCopy(void)
{
    NDMaterial * tmp = 
            new Template3Dep( this->getTag()  ,
			      this->getYS()   ,
			      this->getPS()   ,
			      this->getEPS()  ,
			      this->getELS1() ,
			      this->getELS2() ,
			      this->getELS3() ,
			      this->getELS4() ,
			      this->getELT1() ,
			      this->getELT2() ,
			      this->getELT3() ,
			      this->getELT4() );
			      
    return tmp;

}


//================================================================================
NDMaterial * Template3Dep::getCopy(const char *code)
{
    if (strcmp(code,"Template3Dep") == 0)
    {
       Template3Dep * tmp = 
            new Template3Dep( this->getTag()  ,
			      this->getYS()   ,
			      this->getPS()   ,
			      this->getEPS()  ,
			      this->getELS1() ,
			      this->getELS2() ,
			      this->getELS3() ,
			      this->getELS4() ,
			      this->getELT1() ,
			      this->getELT2() ,
			      this->getELT3() ,
			      this->getELT4() );

       tmp->getEPS()->setInit();			      
       return tmp;
    }
    else
    {
      opserr << "Template3Dep::getCopy failed to get model: " <<  code << endln;
      exit(-1);
    }

}

//================================================================================
const char *Template3Dep::getType(void) const
{
    return "Template3Dep";
}

//================================================================================
//??What is the Order????????? might be the

int Template3Dep::getOrder(void) const
{
    return 6;
}

//================================================================================
int Template3Dep::sendSelf(int commitTag, Channel &theChannel)
{
    // Not yet implemented
    return 0;
}

//================================================================================
int Template3Dep::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    // Not yet implemented
    return 0;
}

//================================================================================
void
Template3Dep::Print(OPS_Stream &s, int flag)
{
     s << (*this);
}



//private utilities

//================================================================================
// Set new EPState
//================================================================================

void Template3Dep::setEPS( EPState & rval) 
{   
    //EPState eps = rval; older buggy one
    //EPS = rval.newObj();
/*
//EPS->setEo(rval.getEo());	                
EPS->setE(rval.getE());                 	 
//EPS->setnu(getnu());                 	 
//EPS->setrho(getrho());                	 
EPS->setStress(rval.getStress());             	 
EPS->setStrain(rval.getStrain());             	 
EPS->setElasticStrain(rval.getElasticStrain());      	 
EPS->setPlasticStrain(rval.getPlasticStrain());      	 
EPS->setdElasticStrain(rval.getdElasticStrain());     	 
EPS->setdPlasticStrain(rval.getdPlasticStrain());     	 
//EPS->setNScalarVar(rval.getNScalarVar());         	 
for (int i = 0; i <rval.getNScalarVar(); i++)
  EPS->setScalarVar(i, rval.getScalarVar(i));


//EPS->setNTensorVar(rval.getNTensorVar());         	 
EPS->setTensorVar(rval.getTensorVar());          	 
EPS->setEep(rval.getEep());           	 
EPS->Stress_commit=(rval.getStress_commit());      	 
EPS->Strain_commit=(rval.getStrain_commit());      	 

for (int i = 0; i <rval.getNScalarVar(); i++)
  EPS->ScalarVar_commit[i] = rval.getScalarVar_commit()[i];   	 

for (int i = 0; i <rval.getNTensorVar(); i++)
   EPS->TensorVar_commit[i] = (rval.getTensorVar_commit()[i]);   	 

EPS->Eep_commit = (rval.getEep_commit());         	 
//EPS->setStress_init(getStress_init());        	 
//EPS->setStrain_init(getStrain_init());        	 
//EPS->setScalarVar_init(getScalarVar_init());     	 
//EPS->setTensorVar_init(getTensorVar_init());     	 
//EPS->setEep_init(getEep_init());           	 
EPS->setConverged(rval.getConverged());           	 

*/   
    
    EPS->setElasticflag(rval.getElasticflag());

    EPS->setEo(rval.getEo());
    EPS->setE(rval.getE());
    EPS->setEv(rval.getEv());
    EPS->setGhv(rval.getGhv());
    EPS->setnu(rval.getnu());                 	 
    EPS->setnuhv(rval.getnuhv());

    EPS->setStress(rval.getStress());
    EPS->setStrain(rval.getStrain());
    EPS->setElasticStrain(rval.getElasticStrain());
    EPS->setPlasticStrain(rval.getPlasticStrain());
    EPS->setdElasticStrain(rval.getdElasticStrain());
    EPS->setdPlasticStrain(rval.getdPlasticStrain());
    EPS->setNScalarVar( rval.getNScalarVar() );

    int i;
    for (i = 0; i <rval.getNScalarVar(); i++)
      EPS->setScalarVar(i+1, rval.getScalarVar(i+1));
    
    
    EPS->setNTensorVar( rval.getNTensorVar() );
    EPS->setTensorVar(rval.getTensorVar());          	 
    EPS->setEep(rval.getEep());           	 
    EPS->setStress_commit(rval.getStress_commit());      	 
    EPS->setStrain_commit(rval.getStrain_commit());      	 
    
    for (i = 0; i <rval.getNScalarVar(); i++)
      EPS->setScalarVar_commit(i+1, rval.getScalarVar_commit(i+1));   	 
    
    for (i = 0; i <rval.getNTensorVar(); i++)
       EPS->setTensorVar_commit(i+1, rval.getTensorVar_commit(i+1));   	 
    
    EPS->Eep_commit = (rval.getEep_commit());
    EPS->Stress_init = rval.getStress_init();
    EPS->Strain_init = rval.getStrain_init();

    for (i = 0; i <rval.getNScalarVar(); i++)
       EPS->setScalarVar_init(i+1, rval.getScalarVar_init(i+1));

    for (i = 0; i <rval.getNTensorVar(); i++)
       EPS->setTensorVar_init(i+1, rval.getTensorVar_init(i+1));

    EPS->Eep_init = rval.getEep_init();
    EPS->setConverged(rval.getConverged());
    
    //Added Joey 02-13-03
    EPS->seteo(rval.geteo());
    EPS->setec(rval.getec());
    EPS->setLam(rval.getLam());
    EPS->setpo(rval.getpo());
    EPS->sete(rval.gete());
    EPS->setpsi(rval.getpsi());
    EPS->seta(rval.geta());
}    	   



//================================================================================
// Get the Yield Surface
//================================================================================
YieldSurface * Template3Dep::getYS() const 
{
    return YS;
}


//================================================================================
// Get the Potential Surface
//================================================================================
PotentialSurface * Template3Dep::getPS() const 
{
    return PS; 
}

//================================================================================
// Get the EPState
//================================================================================
EPState * Template3Dep::getEPS() const 
{ 
    return EPS; 
}


//================================================================================
// Get the 1st Salar evolution law
//================================================================================

EvolutionLaw_S * Template3Dep::getELS1() const 
{ 
    return ELS1; 
}

//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS2() const 
{ 
    return ELS2; 
}

//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS3() const 
{ 
    return ELS3; 
}
//================================================================================
// Get the 2ndst Salar evolution law
//================================================================================
EvolutionLaw_S * Template3Dep::getELS4() const 
{ 
    return ELS4; 
}


//================================================================================
// Get the 1st tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT1() const 
{ 
    return ELT1; 
}

//================================================================================
// Get the 2nd tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT2() const 
{ 
    return ELT2; 
}
//================================================================================
// Get the 3rd tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT3() const 
{ 
    return ELT3; 
}
//================================================================================
// Get the 4th tensorial evolution law
//================================================================================

EvolutionLaw_T * Template3Dep::getELT4() const 
{ 
    return ELT4; 
}


//================================================================================
// Predictor EPState by Forward, Backward, MidPoint Methods...
//================================================================================

EPState Template3Dep::PredictorEPState(straintensor & strain_increment)
{
      		    
    EPState Predictor = ForwardEulerEPState( strain_increment);
    //EPState Predictor = SemiBackwardEulerEPState( strain_increment, material_point);
    return Predictor;

}

//================================================================================
// New EPState using Forward Euler Algorithm
//================================================================================
EPState Template3Dep::ForwardEulerEPState( const straintensor &strain_increment)
{
    // Volumetric strain
    //double st_vol = strain_increment.p_hydrostatic();
    double st_vol = strain_increment.Iinvariant1(); //- = compressive
    //opserr << st_vol << " st_vol1 " << st_vol1 << "\n";

    //EPState forwardEPS( *(material_point->getEPS()) ); 
    EPState forwardEPS( *(this->getEPS()) ); 
    //opserr <<"start eps: " <<   forwardEPS;
    //opserr << "\nForwardEulerEPState  strain_increment " << strain_increment << endlnn;
    
    // Building elasticity tensor
    tensor E    = ElasticStiffnessTensor();
    //tensor Eep  = ElasticStiffnessTensor();
    tensor Eep  = E;
    tensor D    = ElasticComplianceTensor();
    E.null_indices();
    D.null_indices();
    
    //Checking E and D
    //tensor ed = E("ijpq") * D("pqkl");
    //double edd = ed.trace(); // = 3.
    
    straintensor strain_incr = strain_increment;
    strain_incr.null_indices();
    stresstensor stress_increment = E("ijpq") * strain_incr("pq");
    stress_increment.null_indices();
    //opserr << " stress_increment: " << stress_increment << endlnn;
    	 
    EPState startEPS( *(getEPS()) );
    stresstensor start_stress = startEPS.getStress();
    start_stress.null_indices();
    //opserr << "===== start_EPS =====: " << startEPS;
    
    double f_start = 0.0;
    double f_pred  = 0.0;
    
    EPState IntersectionEPS( startEPS );

    EPState ElasticPredictorEPS( startEPS );
    stresstensor elastic_predictor_stress = start_stress + stress_increment;
    ElasticPredictorEPS.setStress( elastic_predictor_stress );
    //opserr << " Elastic_predictor_stress: " << elastic_predictor_stress << endlnn;
    
    f_start = this->getYS()->f( &startEPS );  
    //::printf("\n##############  f_start = %.10e  ",f_start);
    //opserr << "\nFE  f_start = " << f_start;
    
    f_pred =  this->getYS()->f( &ElasticPredictorEPS );
    //::printf("##############  f_pred = %.10e\n\n",f_pred);
    //opserr << "  FE  f_pred = " << f_pred << "\n";
    
    stresstensor intersection_stress = start_stress; // added 20april2000 for forward euler
    stresstensor elpl_start_stress = start_stress;
    stresstensor true_stress_increment = stress_increment;
    straintensor El_strain_increment;
    
    if ( f_start <= 0 && f_pred <= 0 || f_start > f_pred )
      {
        //Updating elastic strain increment
        straintensor estrain = ElasticPredictorEPS.getElasticStrain();
        straintensor tstrain = ElasticPredictorEPS.getStrain();
        estrain = estrain + strain_incr;
        tstrain = tstrain + strain_incr;
        ElasticPredictorEPS.setElasticStrain( estrain );
        ElasticPredictorEPS.setStrain( tstrain );
        ElasticPredictorEPS.setdElasticStrain( strain_incr );
        
        //Evolve parameters like void ratio (e) according to elastic strain and elastic stress--for MD model specifically
        //double Delta_lambda = 0.0;
        //material_point.EL->UpdateVar( &ElasticPredictorEPS, 1);
        // Update E_Young and e according to current stress state before evaluate ElasticStiffnessTensor
        if ( getELT1() ) {
	    //getELT1()->updateEeDm(&ElasticPredictorEPS, st_vol, 0.0);
	    getELT1()->updateEeDm(&ElasticPredictorEPS, -st_vol, 0.0);
	    }
        
        //opserr <<" strain_increment.Iinvariant1() " << strain_increment.Iinvariant1() << endlnn;

        ElasticPredictorEPS.setEep(E);
        return ElasticPredictorEPS;
      }
    
    if ( f_start <= 0 && f_pred > 0 )
      {
        intersection_stress =
           yield_surface_cross( start_stress, elastic_predictor_stress);
        //opserr  << "    start_stress: " <<  start_stress << endlnn;
        //opserr  << "    Intersection_stress: " <<  intersection_stress << endlnn;
    
        IntersectionEPS.setStress( intersection_stress );
        //intersection_stress.reportshort("Intersection stress \n");
      
        elpl_start_stress = intersection_stress;
        //elpl_start_stress.reportshortpqtheta("elpl start stress AFTER \n");
      
        true_stress_increment = elastic_predictor_stress - elpl_start_stress;
        //true_stress_increment.null_indices();
 
	stresstensor EstressIncr = intersection_stress- start_stress;

        //forwardEPS.setStress( elpl_start_stress );
	//Should only count on that elastic portion, not st_vol...
        if ( getELT1() ) {
            El_strain_increment = D("ijpq") * EstressIncr("pq");
 	    double st_vol_El_incr = El_strain_increment.Iinvariant1();

	    //opserr << " FE crossing update... ";
	    getELT1()->updateEeDm(&IntersectionEPS, -st_vol_El_incr, 0.0);
	    //opserr << " FE crossing update... ";
	    getELT1()->updateEeDm(&forwardEPS, -st_vol_El_incr, 0.0);
	}
    
      }
    
    
    //forwardEPS.setStress( elpl_start_stress );
    
    //opserr <<"elpl start eps: " <<   forwardEPS;
    //double f_cross =  this->getYS()->f( &forwardEPS );
    //opserr << " #######  f_cross = " << f_cross << "\n";
    
    //set the initial value of D once the current stress hits the y.s. for Manzari-Dafalias Model
    //if ( f_start <= 0 && f_pred > 0 )
    //    material_point.EL->setInitD(&forwardEPS);
    //opserr << " inside ConstitutiveDriver after setInitD " << forwardEPS;
    
    
    //  pulling out some tensor and double definitions
    //tensor dFods( 2, def_dim_2, 0.0);
    //tensor dQods( 2, def_dim_2, 0.0);
    stresstensor dFods;
    stresstensor dQods;
    //  stresstensor s;  // deviator
    tensor H( 2, def_dim_2, 0.0);
    tensor temp1( 2, def_dim_2, 0.0);
    tensor temp2( 2, def_dim_2, 0.0);
    double lower = 0.0;
    tensor temp3( 2, def_dim_2, 0.0);
    
    double Delta_lambda = 0.0;
    double h_s[4]       = {0.0, 0.0, 0.0, 0.0};
    double xi_s[4]      = {0.0, 0.0, 0.0, 0.0};
    stresstensor h_t[4];
    stresstensor xi_t[4];
    double hardMod_     = 0.0;
    
    //double Dq_ast   = 0.0;
    //double q_ast_entry = 0.0;
    //double q_ast = 0.0;
    
    stresstensor plastic_stress;
    straintensor plastic_strain;
    stresstensor elastic_plastic_stress;
    // ::printf("\n北北北北北北...... felpred = %lf\n",felpred);
    
    if ( f_pred >= 0 ) {
        

        //dFods = getYS()->dFods( &forwardEPS );
        //dQods = getPS()->dQods( &forwardEPS );
        dFods = getYS()->dFods( &IntersectionEPS );
        dQods = getPS()->dQods( &IntersectionEPS );

        //opserr << "dF/ds" << dFods << endlnn;
        //opserr << "dQ/ds" << dQods << endlnn;
    
        // Tensor H_kl  ( eq. 5.209 ) W.F. Chen
        H = E("ijkl")*dQods("kl");       //E_ijkl * R_kl
        H.null_indices();
        temp1 = dFods("ij") * E("ijkl"); // L_ij * E_ijkl
        temp1.null_indices();
        temp2 = temp1("ij")*dQods("ij"); // L_ij * E_ijkl * R_kl
        temp2.null_indices();
        lower = temp2.trace();
        
        // Evaluating the hardening modulus: sum of  (df/dq*) * qbar
	
	hardMod_ = 0.0;
	//Of 1st scalar internal vars
	if ( getELS1() ) {
	   h_s[0]  = getELS1()->h_s(&IntersectionEPS, getPS());
           xi_s[0] = getYS()->xi_s1( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[0] * xi_s[0];
	}

	//Of 2nd scalar internal vars
	if ( getELS2() ) {
	   h_s[1]  = getELS2()->h_s( &IntersectionEPS, getPS());
           xi_s[1] = getYS()->xi_s2( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[1] * xi_s[1];
	}

	//Of 3rd scalar internal vars
	if ( getELS3() ) {
	   h_s[2]  = getELS3()->h_s( &IntersectionEPS, getPS());
           xi_s[2] = getYS()->xi_s3( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[2] * xi_s[2];
	}

	//Of 4th scalar internal vars
	if ( getELS4() ) {
	   h_s[3]  = getELS4()->h_s( &IntersectionEPS, getPS());
           xi_s[3] = getYS()->xi_s4( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[3] * xi_s[3];
	}
	      
	//Of tensorial internal var
	// 1st tensorial var
	if ( getELT1() ) {
	   h_t[0]  = getELT1()->h_t(&IntersectionEPS, getPS());
	   xi_t[0] = getYS()->xi_t1( &IntersectionEPS );
           tensor hm = (h_t[0])("ij") * (xi_t[0])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 2nd tensorial var
	if ( getELT2() ) {
	   h_t[1]  = getELT2()->h_t( &IntersectionEPS, getPS());
  	   xi_t[1] = getYS()->xi_t2( &IntersectionEPS );
           tensor hm = (h_t[1])("ij") * (xi_t[1])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 3rd tensorial var
	if ( getELT3() ) {
	   h_t[2]  = getELT3()->h_t( &IntersectionEPS, getPS());
	   xi_t[2] = getYS()->xi_t3( &IntersectionEPS );
           tensor hm = (h_t[2])("ij") * (xi_t[2])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 4th tensorial var
	if ( getELT4() ) {
	   h_t[3]  = getELT4()->h_t(&IntersectionEPS, getPS());
	   xi_t[3] = getYS()->xi_t4( &IntersectionEPS );
           tensor hm = (h_t[3])("ij") * (xi_t[3])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// Subtract accumulated hardMod_ from lower
        lower = lower - hardMod_;
        
        //Calculating Kp according to Kp = - (df/dq*) * qbar
        //double Kp = material_point.EL->getKp(&forwardEPS, norm_dQods);
        //Kp = 0.0;
        //opserr << endlnn << ">>>>>>>>>   Lower = " << lower << endlnn;           
        //lower = lower + Kp;
        //opserr << endlnn << ">>>>>>>>>    Kp = " << Kp << endlnn;           
        
        //opserr << " stress_increment "<< stress_increment << endlnn;
        //opserr << " true_stress_increment "<< true_stress_increment << endlnn;
    
        //temp3 = dFods("ij") * true_stress_increment("ij"); // L_ij * E_ijkl * d e_kl (true ep strain increment)
        //temp3.null_indices();
        //opserr << " temp3.trace() -- true_stress_incr " << temp3.trace() << endln;
	temp3 = temp1("ij")*strain_incr("ij");
        temp3.null_indices();
        //opserr << " temp3.trace() " << temp3.trace() << endlnn;
        Delta_lambda = (temp3.trace())/lower;
        //opserr << "Delta_lambda " <<  Delta_lambda << endln; 
        if (Delta_lambda<0.0) Delta_lambda=0.0;

        plastic_stress = H("kl") * Delta_lambda;
        plastic_strain = dQods("kl") * Delta_lambda; // plastic strain increment
        plastic_stress.null_indices();
        plastic_strain.null_indices(); 
        //opserr << " Delta_lambda " << Delta_lambda << "plastic_stress =   " << plastic_stress << endln;
        //opserr << "plastic_stress =   " << plastic_stress << endlnn;
        //opserr << "plastic_strain =   " << plastic_strain << endlnn;
        //opserr << "plastic_strain I1= " << plastic_strain.Iinvariant1() << endlnn;
        //opserr << "plastic_strain vol " << Delta_lambda * ( forwardEPS.getScalarVar( 2 ) )<< endlnn ; 
        //opserr << "  q=" << Delta_lambda * dQods.q_deviatoric()<< endlnn;
        //plastic_stress.reportshort("plastic stress (with delta_lambda)\n");
        
        elastic_plastic_stress = elastic_predictor_stress - plastic_stress;
        //elastic_plastic_stress.reportshortpqtheta("FE elastic plastic stress \n");
        	     
        //calculating elatic strain increment
        //stresstensor dstress_el = elastic_plastic_stress - start_stress;
        //straintensor elastic_strain = D("ijpq") * dstress_el("pq");
        straintensor elastic_strain = strain_incr - plastic_strain;  // elastic strain increment
        //opserr << "elastic_strain I1=" << elastic_strain.Iinvariant1() << endlnn;
        //opserr << "elastic_strain " << elastic_strain << endlnn;
        //opserr << "strain increment I1=" << strain_increment.Iinvariant1() << endlnn;
        //opserr << "strain increment    " << strain_increment << endlnn;
    
        straintensor estrain = forwardEPS.getElasticStrain(); //get old elastic strain
        straintensor pstrain = forwardEPS.getPlasticStrain(); //get old plastic strain 
    
        straintensor tstrain = forwardEPS.getStrain();        //get old total strain
        pstrain = pstrain + plastic_strain;
        estrain = estrain + elastic_strain;
        tstrain = tstrain + elastic_strain + plastic_strain;
        
        //Setting de_p, de_e, total plastic, elastic strain, and  total strain
        forwardEPS.setdPlasticStrain( plastic_strain );
        forwardEPS.setdElasticStrain( elastic_strain );
        forwardEPS.setPlasticStrain( pstrain );
        forwardEPS.setElasticStrain( estrain );
        forwardEPS.setStrain( tstrain );
        
        //================================================================
     	//Generating Eep using  dQods at the intersection point
        dFods = getYS()->dFods( &IntersectionEPS );
        dQods = getPS()->dQods( &IntersectionEPS );

	tensor upperE1 = E("pqkl")*dQods("kl");
        upperE1.null_indices();
	tensor upperE2 = dFods("ij")*E("ijmn");
        upperE2.null_indices();
	
	tensor upperE = upperE1("pq") * upperE1("mn");
        upperE.null_indices();

        /*//temp2 = upperE2("ij")*dQods("ij"); // L_ij * E_ijkl * R_kl
        temp2.null_indices();
        lower = temp2.trace();
        
        
	// Evaluating the hardening modulus: sum of  (df/dq*) * qbar
	
	hardMod_ = 0.0;
	//Of 1st scalar internal vars
	if ( getELS1() ) {
	   h_s[0]  = getELS1()->h_s( &IntersectionEPS, getPS());
           xi_s[0] = getYS()->xi_s1( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[0] * xi_s[0];
	}

	//Of 2nd scalar internal vars
	if ( getELS2() ) {
	   h_s[1]  = getELS2()->h_s( &IntersectionEPS, getPS());
           xi_s[1] = getYS()->xi_s2( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[1] * xi_s[1];
	}

	//Of 3rd scalar internal vars
	if ( getELS3() ) {
	   h_s[2]  = getELS3()->h_s( &IntersectionEPS, getPS());
           xi_s[2] = getYS()->xi_s3( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[2] * xi_s[2];
	}

	//Of 4th scalar internal vars
	if ( getELS4() ) {
	   h_s[3]  = getELS4()->h_s( &IntersectionEPS, getPS());
           xi_s[3] = getYS()->xi_s4( &IntersectionEPS );	   
   	   hardMod_ = hardMod_ + h_s[3] * xi_s[3];
	}
	      
	//Of tensorial internal var
	// 1st tensorial var
	if ( getELT1() ) {
	   h_t[0]  = getELT1()->h_t(&IntersectionEPS, getPS());
	   xi_t[0] = getYS()->xi_t1( &IntersectionEPS );
           tensor hm = (h_t[0])("ij") * (xi_t[0])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 2nd tensorial var
	if ( getELT2() ) {
	   h_t[1]  = getELT2()->h_t( &IntersectionEPS, getPS());
	   xi_t[1] = getYS()->xi_t2( &IntersectionEPS );
           tensor hm = (h_t[1])("ij") * (xi_t[1])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 3rd tensorial var
	if ( getELT3() ) {
	   h_t[2]  = getELT3()->h_t(&IntersectionEPS, getPS());
	   xi_t[2] = getYS()->xi_t3( &IntersectionEPS );
           tensor hm = (h_t[2])("ij") * (xi_t[2])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 4th tensorial var
	if ( getELT4() ) {
	   h_t[3]  = getELT4()->h_t( &IntersectionEPS, getPS());
	   xi_t[3] = getYS()->xi_t4( &IntersectionEPS );
           tensor hm = (h_t[3])("ij") * (xi_t[3])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// Subtract accumulated hardMod_ from lower
        
	lower = lower - hardMod_;
	*/

        tensor Ep = upperE*(1./lower);

	// elastoplastic constitutive tensor
	double h_L = 0.0; // Bug fixed Joey 07-21-02 added h(L) function
	if ( Delta_lambda > 0 ) h_L = 1.0;
	//opserr << " h_L = " << h_L << "\n";
        Eep =  Eep - Ep*h_L; 

     	//opserr <<" after calculation---Eep.rank()= " << Eep.rank() <<endlnn;
	//Eep.printshort(" IN template ");
         
        //--// before the surface is been updated !
        //--//        f_Final = Criterion.f(elastic_plastic_stress);
        //--
        //--        q_ast_entry = Criterion.kappa_get(elastic_plastic_stress);
        //--
        //--//        h_  = h(elastic_plastic_stress);
        //--        Dq_ast = Delta_lambda * h_ * just_this_PP;
        //--//        if (Dq_ast < 0.0 ) Dq_ast = 0.0;
        //--//        Dq_ast = fabs(Delta_lambda * h_ * just_this_PP); // because of softening response...
        //--//::printf(" h_=%.6e  q_ast_entry=%.6e  Dq_ast=%.6e   Delta_lambda=%.6e\n", h_, q_ast_entry, Dq_ast, Delta_lambda);
        //--
        //--        current_lambda_set(Delta_lambda);
        //--
        //--        q_ast = q_ast_entry + Dq_ast;
        //--        Criterion.kappa_set( elastic_plastic_stress, q_ast);
        //--//::fprintf(stdout," Criterion.kappa_get(elastic_plastic_stress)=%.8e\n",Criterion.kappa_get(elastic_plastic_stress));
        //--//::fprintf(stderr," Criterion.kappa_get(elastic_plastic_stress)=%.8e\n",Criterion.kappa_get(elastic_plastic_stress));
        //--
        //--
        //--//::fprintf(stdout," ######## predictor --> q_ast_entry=%.8e Dq_ast=%.8e q_ast=%.8e\n",q_ast_entry, Dq_ast, q_ast);
        //--//::fprintf(stderr," ######## predictor --> q_ast_entry=%.8e Dq_ast=%.8e q_ast=%.8e\n",q_ast_entry, Dq_ast, q_ast);
        //--
        //--//::fprintf(stdout,"ForwardEulerStress IN Criterion.kappa_get(start_stress)=%.8e\n",Criterion.kappa_get(start_stress));
        //--//::fprintf(stderr,"ForwardEulerStress IN Criterion.kappa_get(start_stress)=%.8e\n",Criterion.kappa_get(start_stress));
        //--
        
        //new EPState in terms of stress
        forwardEPS.setStress(elastic_plastic_stress); 
        //opserr <<"inside constitutive driver: forwardEPS "<< forwardEPS;

        forwardEPS.setEep(Eep); 
        //forwardEPS.getEep().printshort(" after set"); 
        
        //Before update all the internal vars
        double f_forward =  getYS()->f( &forwardEPS );
        //::printf("\n************  Before f_forward = %.10e\n",f_forward);
    
        //Evolve the surfaces and hardening vars
	int NS = forwardEPS.getNScalarVar();
	int NT = forwardEPS.getNTensorVar();

	double dS= 0.0;
	double S = 0.0;

	int ii;
	for (ii = 1; ii <= NS; ii++) {
              dS = Delta_lambda * h_s[ii-1] ;       // Increment to the scalar internal var
              S  = forwardEPS.getScalarVar(ii);     // Get the old value of the scalar internal var
              forwardEPS.setScalarVar(ii, S + dS ); // Update internal scalar var
	}
	
	stresstensor dT;
	stresstensor T;
	stresstensor new_T;

	for (ii = 1; ii <= NT; ii++) {
	      dT = h_t[ii-1]*Delta_lambda  ;       // Increment to the tensor internal var
              T  = forwardEPS.getTensorVar(ii);     // Get the old value of the tensor internal var
              new_T = T + dT;
              forwardEPS.setTensorVar(ii, new_T );
        }
       
        // Update E_Young and e according to current stress state
	//straintensor pl_strain_increment;
        int err = 0;
	if ( getELT1() ) {
	    //pl_strain_increment = strain_increment - El_strain_increment;
	    double st_vol_pl = plastic_strain.Iinvariant1();	    
	    //err = getELT1()->updateEeDm(&forwardEPS, st_vol_pl, Delta_lambda);
	    //cerr << "pl_vol= " << st_vol_pl << "|update before FE \n";
	    //opserr << " FE Pl update...";
	    // D > 0 compressive -> Iinv > 0  -> de < 0 correct!
	    err = getELT1()->updateEeDm(&forwardEPS, st_vol_pl, Delta_lambda);
        }
	   
        //tensor tempx  = plastic_strain("ij") * plastic_strain("ij");
        //double tempxd = tempx.trace();
        //double e_eq  = pow( 2.0 * tempxd / 3.0, 0.5 );
        ////opserr << "e_eq = " << e_eq << endlnn;
        //
        //double dalfa1 =  e_eq * 10;
        //double alfa1  = forwardEPS.getScalarVar(1);


        
    	//opserr << "UpdateAllVars " << forwardEPS<< endlnn;

        //After update all the internal vars
        f_forward =  getYS()->f( &forwardEPS );
        //opserr << "\n************  After f_forward = " <<  f_forward << "\n";
    
    
    }
    
    //::fprintf(stderr,"ForwardEulerStress EXIT Criterion.kappa_get(start_stress)=%.8e\n",Criterion.kappa_get(start_stress));
    //forwardEPS.setConverged(TRUE);

    //double p = (forwardEPS.getStress()).p_hydrostatic();
    //double ec = (forwardEPS.getec()) - (forwardEPS.getLam()) * log( p / (forwardEPS.getpo()) );
    //double st = (forwardEPS.getStrain()).Iinvariant1();
    //double pl_s = (forwardEPS.getPlasticStrain()).Iinvariant1();
    //double dpl_s = (forwardEPS.getdPlasticStrain()).Iinvariant1();
    //cerr << "P FE p=" << p << " ec " << ec << " e " << forwardEPS.gete() << " psi " << (forwardEPS.gete() - ec) << " strain " << st << " t_pl " << pl_s << " d_pl " << dpl_s << "\n";
 
    return forwardEPS;
}

//================================================================================
// Starting EPState using Semi Backward Euler Starting Point
//================================================================================
EPState Template3Dep::SemiBackwardEulerEPState( const straintensor &strain_increment)
{
    stresstensor start_stress;
    EPState SemibackwardEPS( *(this->getEPS()) ); 
    start_stress = SemibackwardEPS.getStress();
   
    // building elasticity tensor
    //tensor E = Criterion.StiffnessTensorE(Ey,nu);
    tensor E  = ElasticStiffnessTensor();
    // building compliance tensor
    //  tensor D = Criterion.ComplianceTensorD(Ey,nu);
    
    //  pulling out some tensor and double definitions
    tensor dFods(2, def_dim_2, 0.0);
    tensor dQods(2, def_dim_2, 0.0);
    tensor temp2(2, def_dim_2, 0.0);
    double lower = 0.0;
    double Delta_lambda = 0.0;
  
    EPState E_Pred_EPS( *(this->getEPS()) );

    straintensor strain_incr = strain_increment;
    stresstensor stress_increment = E("ijkl")*strain_incr("kl");
    stress_increment.null_indices();
    // stress_increment.reportshort("stress Increment\n");

  
    stresstensor plastic_stress;
    stresstensor elastic_predictor_stress;
    stresstensor elastic_plastic_stress;
    //..  double Felplpredictor = 0.0;
  
    double h_s[4]       = {0.0, 0.0, 0.0, 0.0};
    double xi_s[4]      = {0.0, 0.0, 0.0, 0.0};
    stresstensor h_t[4];
    stresstensor xi_t[4];
    double hardMod_ = 0.0;

    double S        = 0.0;
    double dS       = 0.0;
    stresstensor T;
    stresstensor dT;
    //double Dq_ast   = 0.0;
    //double q_ast_entry = 0.0;
    //double q_ast = 0.0;
  
    elastic_predictor_stress = start_stress + stress_increment;
    //..  elastic_predictor_stress.reportshort("ELASTIC PREDICTOR stress\n");
    E_Pred_EPS.setStress( elastic_predictor_stress );
  
    //  double f_start = Criterion.f(start_stress);
    //  ::printf("SEMI BE##############  f_start = %.10e\n",f_start);
    double f_pred =  getYS()->f( &E_Pred_EPS );
    //::printf("SEMI BE##############  f_pred = %.10e\n",f_pred);
  
    // second of alternative predictors as seen in MAC page 170-171
    if ( f_pred >= 0.0 )
    {
        //el_or_pl_range(1); // set to 1 ( plastic range )
        // PREDICTOR FASE
        //..     ::printf("\n\npredictor part  step_counter = %d\n\n", step_counter);
  
        dFods = getYS()->dFods( &E_Pred_EPS );
        dQods = getPS()->dQods( &E_Pred_EPS );
        //.... dFods.print("a","dF/ds  ");
        //.... dQods.print("a","dQ/ds  ");
  
        temp2 = (dFods("ij")*E("ijkl"))*dQods("kl");
        temp2.null_indices();
        lower = temp2.trace();
  
        // Evaluating the hardening modulus: sum of  (df/dq*) * qbar

 	//Of scalar internal var
	hardMod_ = 0.0;

	//Of 1st scalar internal vars
	if ( getELS1() ) {
	   h_s[0]  = getELS1()->h_s( &E_Pred_EPS, getPS());
           xi_s[0] = getYS()->xi_s1( &E_Pred_EPS );	   
   	   hardMod_ = hardMod_ + h_s[0] * xi_s[0];
	}

	//Of 2nd scalar internal vars
	if ( getELS2() ) {
	   h_s[1]  = getELS2()->h_s( &E_Pred_EPS, getPS());
           xi_s[1] = getYS()->xi_s2( &E_Pred_EPS );	   
   	   hardMod_ = hardMod_ + h_s[1] * xi_s[1];
	}

	//Of 3rd scalar internal vars
	if ( getELS3() ) {
	   h_s[2]  = getELS3()->h_s(&E_Pred_EPS, getPS());
           xi_s[2] = getYS()->xi_s3( &E_Pred_EPS );	   
   	   hardMod_ = hardMod_ + h_s[2] * xi_s[2];
	}

	//Of 4th scalar internal vars
	if ( getELS4() ) {
	   h_s[3]  = getELS4()->h_s( &E_Pred_EPS, getPS());
           xi_s[3] = getYS()->xi_s4( &E_Pred_EPS );	   
   	   hardMod_ = hardMod_ + h_s[3] * xi_s[3];
	}

	//Of tensorial internal var
	// 1st tensorial var
	if ( getELT1() ) {
	   h_t[0]  = getELT1()->h_t(&E_Pred_EPS, getPS());
	   xi_t[0] = getYS()->xi_t1( &E_Pred_EPS );
           tensor hm = (h_t[0])("ij") * (xi_t[0])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 2nd tensorial var
	if ( getELT2() ) {
	   h_t[1]  = getELT2()->h_t(&E_Pred_EPS, getPS());
	   xi_t[1] = getYS()->xi_t2( &E_Pred_EPS );
           tensor hm = (h_t[1])("ij") * (xi_t[1])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 3rd tensorial var
	if ( getELT3() ) {
	   h_t[2]  = getELT3()->h_t(&E_Pred_EPS, getPS());
	   xi_t[2] = getYS()->xi_t3( &E_Pred_EPS );
           tensor hm = (h_t[2])("ij") * (xi_t[2])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// 4th tensorial var
	if ( getELT4() ) {
	   h_t[3]  = getELT4()->h_t(&E_Pred_EPS, getPS());
	   xi_t[3] = getYS()->xi_t4( &E_Pred_EPS );
           tensor hm = (h_t[3])("ij") * (xi_t[3])("ij");
  	   hardMod_ = hardMod_ + hm.trace();
	}

	// Subtract accumulated hardMod_ from lower
        lower = lower - hardMod_;

        //h_s  = material_point.ELS1->h_s( &E_Pred_EPS, material_point.PS );
        //xi_s = material_point.YS->xi_s1( &E_Pred_EPS );
        //hardMod_ = h_s * xi_s;
        //lower = lower - hardMod_;

 	////Of tensorial internal var
	//h_t  = material_point.ELT1->h_t(&E_Pred_EPS, material_point.PS);
        //xi_t = material_point.YS->xi_t1( &E_Pred_EPS );
        //tensor hm = h_t("ij") * xi_t("ij");
	//hardMod_ = hm.trace();
        //lower = lower - hardMod_;
 

        Delta_lambda = f_pred/lower;
        if ( Delta_lambda < 0.0 )
          {
            //::fprintf(stderr,"\nP!\n");
          }
        plastic_stress = (E("ijkl")*dQods("kl"))*(-Delta_lambda);
        plastic_stress.null_indices();
        //.. plastic_stress.reportshort("plastic stress predictor II\n");
        //.. elastic_predictor_stress.reportshort("elastic predictor stress \n");
        elastic_plastic_stress = elastic_predictor_stress + plastic_stress;
        elastic_plastic_stress.null_indices();
  	
	SemibackwardEPS.setStress( elastic_plastic_stress );

        ////q_ast_entry = Criterion.kappa_get(elastic_plastic_stress);
        //S  = SemibackwardEPS.getScalarVar(1); // Get the old value of the internal var
        //h_s  = material_point.ELS1->h_s( &SemibackwardEPS, material_point.PS );
	//dS = Delta_lambda * h_s ;   // Increment to the internal scalar var
        //h_t  = material_point.ELT1->h_t( &SemibackwardEPS, material_point.PS );
	//dT = Delta_lambda * h_t ;   // Increment to the internal tensorial var

        //Evolve the surfaces and hardening vars
	int NS = SemibackwardEPS.getNScalarVar();
	int NT = SemibackwardEPS.getNTensorVar();

	int ii;
	for (ii = 1; ii <= NS; ii++) {
              dS = Delta_lambda * h_s[ii-1] ;       // Increment to the scalar internal var
              S  = SemibackwardEPS.getScalarVar(ii);     // Get the old value of the scalar internal var
              SemibackwardEPS.setScalarVar(ii, S + dS ); // Update internal scalar var
	}


        ////current_lambda_set(Delta_lambda);
        ////q_ast = q_ast_entry + Dq_ast;
        ////Criterion.kappa_set( elastic_plastic_stress, q_ast);
        //SemibackwardEPS.setScalarVar(1, S + dS );

	//stresstensor new_T = T + dT;
        //SemibackwardEPS.setTensorVar(1, new_T );

	stresstensor new_T;

	for (ii = 1; ii <= NT; ii++) {
	      dT = h_t[ii-1] * Delta_lambda;            // Increment to the tensor internal var
              T  = SemibackwardEPS.getTensorVar(ii);     // Get the old value of the tensor internal var
              new_T = T + dT;
              SemibackwardEPS.setTensorVar(ii, new_T );
        }

  
        //return elastic_plastic_stress;
        return SemibackwardEPS;
    }
    return E_Pred_EPS;
}
  
  
  

//================================================================================
// New EPState using Backward Euler Algorithm
//================================================================================
EPState Template3Dep::BackwardEulerEPState( const straintensor &strain_increment)

{
  // Temp matertial point
  //NDMaterial MP( material_point );
  //NDMaterial *MP = this->getCopy();

  // Volumetric strain
  //double st_vol = strain_increment.p_hydrostatic();
  double st_vol = strain_increment.Iinvariant1(); //- = compressive
  
  //EPState to be returned, it can be elastic or elastic-plastic EPState
  EPState backwardEPS( * (this->getEPS()) ); 
  
  EPState startEPS( *(this->getEPS()) );
  stresstensor start_stress = startEPS.getStress();

  //Output for plotting
  opserr.precision(5); 
  opserr.width(10);
  //opserr << " strain_increment " << strain_increment << "\n";
  
  opserr.precision(5); 
  opserr.width(10);
  //opserr << "start_stress " <<  start_stress;
      
  // Pulling out some tensor and double definitions
  tensor I2("I", 2, def_dim_2);
  tensor I_ijkl("I", 4, def_dim_4);
  I_ijkl = I2("ij")*I2("kl");
  I_ijkl.null_indices();
  tensor I_ikjl("I", 4, def_dim_4);
  I_ikjl = I_ijkl.transpose0110();
  
  //double Ey = Criterion.E();
  //double nu = Criterion.nu();
  //tensor E = StiffnessTensorE(Ey,nu);
  tensor E  = ElasticStiffnessTensor();
  // tensor D = ComplianceTensorD(Ey,nu);
  // stresstensor REAL_start_stress = start_stress;
  
  tensor dFods(2, def_dim_2, 0.0);
  tensor dQods(2, def_dim_2, 0.0);
  //  tensor dQodsextended(2, def_dim_2, 0.0);
  //  tensor d2Qodqast(2, def_dim_2, 0.0);
  tensor temp2(2, def_dim_2, 0.0);
  double lower = 0.0;  
  double Delta_lambda = 0.0; // Lambda
  double delta_lambda = 0.0; // dLambda

  double Felplpredictor    = 0.0;
  //Kai  double absFelplpredictor = 0.0;
  //  double Ftolerance = pow(d_macheps(),(1.0/2.0))*1000000.00; //FORWARD no iterations
  //double Ftolerance = pow( d_macheps(), 0.5)*1.00;
  
  double Ftolerance = pow( d_macheps(), 0.5)*1000*KK;  //Zhaohui UCD 10e6 for Pa, kg and m 1000 for kPa, ton and m
  //opserr << Ftolerance << endlnn;
  //  double Ftolerance = pow(d_macheps(),(1.0/2.0))*1.0;
  //  double entry_kappa_cone = Criterion.kappa_cone_get();
  //  double entry_kappa_cap  = Criterion.kappa_cap_get();

  tensor aC(2, def_dim_2, 0.0);
  stresstensor BEstress;
  stresstensor residual;
  tensor d2Qoverds2( 4, def_dim_4, 0.0);
  tensor T( 4, def_dim_4, 0.0);
  tensor Tinv( 4, def_dim_4, 0.0);

  double Fold = 0.0;
  tensor temp3lower;
  tensor temp5;
  double temp6 = 0.0;
  double upper = 0.0;

  stresstensor dsigma;
  //stresstensor Dsigma;
  stresstensor sigmaBack;
  straintensor dPlasticStrain; // delta plastic strain
  straintensor PlasticStrain;  // Total plastic strain
  
  //double dq_ast = 0.0;       // iterative change in internal variable (kappa in this case)
  //double Dq_ast = 0.0;       // incremental change in internal variable (kappa in this case)
  //double q_ast  = 0.0;       // internal variable (kappa in this case)
  //double q_ast_entry  = 0.0; //internal variable from previous increment (kappa in this case)
  
  int step_counter = 0;
  //int MAX_STEP_COUNT = ;
  //  int MAX_STEP_COUNT = 0;
  int flag = 0;

  straintensor strain_incr = strain_increment;
  strain_incr.null_indices();
  stresstensor stress_increment = E("ijkl")*strain_incr("kl");
  stress_increment.null_indices();

  stresstensor Return_stress; //  stress to be returned can be either predictor
                              // or elastic plastic stress.

  EPState ElasticPredictorEPS( startEPS );
  stresstensor elastic_predictor_stress = start_stress + stress_increment;

  ElasticPredictorEPS.setStress( elastic_predictor_stress );
  //  elastic_predictor_stress.reportshortpqtheta("\n . . . .  ELASTIC PREDICTOR stress");

  opserr.precision(5); 
  opserr.width(10);
  //opserr << "elastic predictor " <<  elastic_predictor_stress << endlnn;

  stresstensor elastic_plastic_predictor_stress;
  EPState EP_PredictorEPS( startEPS );

  //double f_start = material_point.YS->f( &startEPS );
  //opserr << " ************** f_start = " << f_start;
  //::fprintf(stdout,"tst##############  f_start = %.10e\n",f_start);
  // f_pred = Criterion.f(elastic_predictor_stress);
  //::fprintf(stdout,"tst##############  f_pred = %.10e\n",f_pred);
  //double f_start =  getYS()->f( &startEPS );
  //opserr << "\n  ************  Enter Backward   f_star **** " << f_start;

  double f_pred =  getYS()->f( &ElasticPredictorEPS );
  //opserr << "  **BE f_pred **** " << f_pred << endln;
  //int region = 5;

  //double h_s      = 0.0;
  //double xi_s     = 0.0;
  double hardMod_ = 0.0;
  
  double h_s[4]       = {0.0, 0.0, 0.0, 0.0};
  double xi_s[4]      = {0.0, 0.0, 0.0, 0.0};
  stresstensor h_t[4];
  stresstensor xi_t[4];

  //region
  //check for the region than proceede
  //region = Criterion.influence_region(elastic_predictor_stress);
  //if ( region == 1 )  // apex gray region
  //  {
  //    double pc_ = pc();
  //    elastic_plastic_predictor_stress =
  //      elastic_plastic_predictor_stress.pqtheta2stress(pc_, 0.0, 0.0);
  //    return  elastic_plastic_predictor_stress;
  //  }
  
  if ( f_pred <= Ftolerance  )
  {

      //Updating elastic strain increment
      straintensor estrain = ElasticPredictorEPS.getElasticStrain();
      straintensor tstrain = ElasticPredictorEPS.getStrain();
      estrain = estrain + strain_increment;
      tstrain = tstrain + strain_increment;
      ElasticPredictorEPS.setElasticStrain( estrain );
      ElasticPredictorEPS.setStrain( tstrain );
      ElasticPredictorEPS.setdElasticStrain( strain_increment );
      //opserr<< "Elastic:  Total strain" << tstrain << endln;
      
      if ( getELT1() ) {
         getELT1()->updateEeDm(&ElasticPredictorEPS, -st_vol, 0.0);
      }

      //Set Elasto-Plastic stiffness tensor
      ElasticPredictorEPS.setEep(E);
      ElasticPredictorEPS.setConverged(TRUE);
      //E.printshort(" BE -- Eep ");
       
      backwardEPS = ElasticPredictorEPS;
      //opserr <<  "\n backwardEPS e " <<  backwardEPS.gete();

      //double p = (backwardEPS.getStress()).p_hydrostatic();
      //double ec = (backwardEPS.getec()) - (backwardEPS.getLam()) * log( p / (backwardEPS.getpo()) );
      //double st = (backwardEPS.getStrain()).Iinvariant1();
      //double pl_s = (backwardEPS.getPlasticStrain()).Iinvariant1();
      //double dpl_s = (backwardEPS.getdPlasticStrain()).Iinvariant1();
      //cerr << "E ec " << ec << " e " << backwardEPS.gete() << " psi " << (backwardEPS.gete() - ec) << " strain " << st << " t_pl " << pl_s << " d_pl " << dpl_s << "\n";      
      return backwardEPS;   

      //opserr <<  "\nbackwardEPS" <<  backwardEPS;
      //opserr <<  "\nElasticPredictorEPS " <<  ElasticPredictorEPS;

  }
  if ( f_pred > 0.0 )
  {

      //Starting point by applying one Forward Euler step
      EP_PredictorEPS = PredictorEPState( strain_incr);
      //EP_PredictorEPS = ElasticPredictorEPS;
            
      //opserr << " ----------Predictor Stress" << EP_PredictorEPS.getStress();
      //Setting the starting EPState with the starting internal vars in EPState
      
      //MP->setEPS( EP_PredictorEPS );
      
      Felplpredictor =  getYS()->f(&EP_PredictorEPS);
      //opserr <<  " F_elplpredictor " << Felplpredictor << endlnn;


      //Kai     absFelplpredictor = fabs(Felplpredictor);
      if ( fabs(Felplpredictor) <= Ftolerance )
      {
	 //Forward Euler will do.
         backwardEPS = EP_PredictorEPS;

         //Return_stress = elastic_plastic_predictor_stress;
         flag = 1;
      }
      else {
        aC    = getPS()->dQods( &EP_PredictorEPS );	 
        dFods = getYS()->dFods( &EP_PredictorEPS );
        dQods = getPS()->dQods( &EP_PredictorEPS );
        
        //   aC    = Criterion.dQods(elastic_predictor_stress);
        //   dFods = Criterion.dFods(elastic_predictor_stress);
        //   dQods = Criterion.dQods(elastic_predictor_stress);
        

        temp2 = (dFods("ij")*E("ijkl"))*dQods("kl");
        temp2.null_indices();
        lower = temp2.trace();
        
        //     Delta_lambda = f_pred/lower; //?????????
        //::printf("  Delta_lambda = f_pred/lower = %.8e\n", Delta_lambda);
        ////     Delta_lambda = Felplpredictor/lower;
        ////::printf("  Delta_lambda = Felplpredictor/lower =%.8e \n", Delta_lambda);

	// Original segment
        //elastic_plastic_predictor_stress = elastic_predictor_stress - E("ijkl")*aC("kl")*Delta_lambda;
        //EP_PredictorEPS.setStress( elastic_plastic_predictor_stress );

	//Zhaohui modified, sometimes give much better convergence rate
        elastic_plastic_predictor_stress = EP_PredictorEPS.getStress();
	//opserr << "elastic_plastic_predictor_stress" << elastic_plastic_predictor_stress;
	
        //opserr.precision(5); 
        //opserr.width(10);
        //opserr << " " << EP_PredictorEPS.getStress().p_hydrostatic() << " ";
      
        //opserr.precision(5); 
        //opserr.width(10);
        //opserr << EP_PredictorEPS.getStress().q_deviatoric()<< " ";

        //opserr.precision(5); 
        //opserr.width(10);
        //opserr << Delta_lambda << endlnn;

	//elastic_plastic_predictor_stress.reportshort("......elastic_plastic_predictor_stress");
        //::printf("  F(elastic_plastic_predictor_stress) = %.8e\n", Criterion.f(elastic_plastic_predictor_stress));
        
        //h_s  = MP->ELS1->h_s( &EP_PredictorEPS, MP->PS );
        ////h_  = h(elastic_plastic_predictor_stress);
        ////  Dq_ast = Criterion.kappa_get(elastic_plastic_predictor_stress);
        
	//q_ast_entry = Criterion.kappa_get(elastic_plastic_predictor_stress);
        //Dq_ast = Delta_lambda * h_ * just_this_PP;
        //q_ast = q_ast_entry + Dq_ast;
	
	//Zhaohui comments: internal vars are alreary evolued in ForwardEulerEPS(...), not necessary here!!!
	//..dS = Delta_lambda * h_s ;   // Increment to the internal var
        //..S  = EP_PredictorEPS.getScalarVar(1); // Get the old value of the internal var
	//..new_S = S + dS;
	//..opserr << "Internal Var : " << new_S << endlnn;
	//..EP_PredictorEPS.setScalarVar(1, new_S); // Get the old value of the internal var
      	
	//::fprintf(stdout," ######## predictor --> Dq_ast=%.8e q_ast=%.8e\n", Dq_ast,        q_ast);
        //::fprintf(stderr," ######## predictor --> Dq_ast=%.8e q_ast=%.8e\n", Dq_ast,        q_ast);

        //Criterion.kappa_set( sigmaBack, q_ast);  //????
        //current_lambda_set(Delta_lambda);	   //????

        //::printf("  Delta_lambda = %.8e\n", Delta_lambda);
        //::printf("step = pre iteracija  #############################--   q_ast = %.10e \n", q_ast);
        //::printf("step = pre iteracija  posle predictora  ###########--   Dq_ast = %.10e \n",Dq_ast);
        //**********
        //**********
        //::printf("\nDelta_lambda  before BE = %.10e \n", Delta_lambda );
        }
        
        //========================== main part of iteration =======================
        //      while ( absFelplpredictor > Ftolerance &&

        while ( (fabs(Felplpredictor) > Ftolerance) && (step_counter < MAX_STEP_COUNT) ) // if more iterations than prescribed
        //out07may97      do
        {
          //opserr << "Iteration " << step_counter << " F " << Felplpredictor;
	  BEstress = elastic_predictor_stress - E("ijkl")*aC("kl")*Delta_lambda;
          //BEstress.reportshort("......BEstress ");
          /////          BEstress = elastic_plastic_predictor_stress - E("ijkl")*aC("kl")*Delta_lambda;
          BEstress.null_indices();
          //          Felplpredictor = Criterion.f(BEstress);
          //          ::printf("\nF_backward_Euler BE = %.10e \n", Felplpredictor);
          residual = elastic_plastic_predictor_stress - BEstress;
          //residual.reportshortpqtheta("\n......residual ");
          //          double ComplementaryEnergy = (residual("ij")*D("ijkl")*residual("ij")).trace();
          //::printf("\n Residual ComplementaryEnergy = %.16e\n", ComplementaryEnergy);

          /////          residual = elastic_predictor_stress - BEstress;
          
	  //d2Qoverds2 = Criterion.d2Qods2(elastic_plastic_predictor_stress);
          d2Qoverds2 = getPS()->d2Qods2( &EP_PredictorEPS );
          //d2Qoverds2.print();

	  T = I_ikjl + E("ijkl")*d2Qoverds2("klmn")*Delta_lambda;
          T.null_indices();

          Tinv = T.inverse();
          //dFods = Criterion.dFods(elastic_plastic_predictor_stress);
          //dQods = Criterion.dQods(elastic_plastic_predictor_stress);
          dFods = getYS()->dFods( &EP_PredictorEPS );
          dQods = getPS()->dQods( &EP_PredictorEPS );

          //Fold = Criterion.f(elastic_plastic_predictor_stress);
          Fold = getYS()->f( &EP_PredictorEPS );
          
	  lower = 0.0; // this is old temp variable used here again :-)
          //h_  = h(elastic_plastic_predictor_stress);
          //xi_ = xi(elastic_plastic_predictor_stress);
         
	  //h_s  = MP->ELS1->h_s( &EP_PredictorEPS, MP->PS );
          //xi_s = MP->YS->xi_s1( &EP_PredictorEPS );
          //hardMod_ = h_s * xi_s;
          
          // Evaluating the hardening modulus: sum of  (df/dq*) * qbar
	  
	  hardMod_ = 0.0;
	  //Of 1st scalar internal vars
	  if ( getELS1() ) {
	     h_s[0]  = getELS1()->h_s( &EP_PredictorEPS, getPS());
             xi_s[0] = getYS()->xi_s1( &EP_PredictorEPS );	   
   	     hardMod_ = hardMod_ + h_s[0] * xi_s[0];
	  }
	  
	  //Of 2nd scalar internal vars
	  if ( getELS2() ) {
	     h_s[1]  = getELS2()->h_s( &EP_PredictorEPS, getPS());
             xi_s[1] = getYS()->xi_s2( &EP_PredictorEPS );	   
   	     hardMod_ = hardMod_ + h_s[1] * xi_s[1];
	  }
	  
	  //Of 3rd scalar internal vars
	  if ( getELS3() ) {
	     h_s[2]  = getELS3()->h_s( &EP_PredictorEPS, getPS());
             xi_s[2] = getYS()->xi_s3( &EP_PredictorEPS );	   
   	     hardMod_ = hardMod_ + h_s[2] * xi_s[2];
	  }
	  
	  //Of 4th scalar internal vars
	  if ( getELS4() ) {
	     h_s[3]  = getELS4()->h_s( &EP_PredictorEPS, getPS());
             xi_s[3] = getYS()->xi_s4( &EP_PredictorEPS );	   
   	     hardMod_ = hardMod_ + h_s[3] * xi_s[3];
	  }
	        
	  //Of tensorial internal var
	  // 1st tensorial var
	  if ( getELT1() ) {
	     h_t[0]  = getELT1()->h_t( &EP_PredictorEPS, getPS());
	     xi_t[0] = getYS()->xi_t1( &EP_PredictorEPS );
             tensor hm = (h_t[0])("ij") * (xi_t[0])("ij");
  	     hardMod_ = hardMod_ + hm.trace();
	  }
	  
	  // 2nd tensorial var
	  if ( getELT2() ) {
	     h_t[1]  = getELT2()->h_t( &EP_PredictorEPS, getPS());
	     xi_t[1] = getYS()->xi_t2( &EP_PredictorEPS );
             tensor hm = (h_t[1])("ij") * (xi_t[1])("ij");
  	     hardMod_ = hardMod_ + hm.trace();
	  }
	  
	  // 3rd tensorial var
	  if ( getELT3() ) {
	     h_t[2]  = getELT3()->h_t( &EP_PredictorEPS, getPS());
	     xi_t[2] = getYS()->xi_t3( &EP_PredictorEPS );
             tensor hm = (h_t[2])("ij") * (xi_t[2])("ij");
  	     hardMod_ = hardMod_ + hm.trace();
	  }
	  
	  // 4th tensorial var
	  if ( getELT4() ) {
	     h_t[3]  = getELT4()->h_t( &EP_PredictorEPS, getPS());
	     xi_t[3] = getYS()->xi_t4( &EP_PredictorEPS );
             tensor hm = (h_t[3])("ij") * (xi_t[3])("ij");
  	     hardMod_ = hardMod_ + hm.trace();
	  }
	  
	  // Subtract accumulated hardMod_ from lower
          //lower = lower - hardMod_;
	  
	  //hardMod_ = hardMod_ * just_this_PP;
          //::printf("\n BackwardEulerStress ..  hardMod_ = %.10e \n", hardMod_ );
          //outfornow          d2Qodqast = d2Qoverdqast(elastic_plastic_predictor_stress);
          //outfornow          dQodsextended = dQods + d2Qodqast * Delta_lambda * h_;
          //outfornow          temp3lower = dFods("mn")*Tinv("ijmn")*E("ijkl")*dQodsextended("kl");
          temp3lower = dFods("mn")*Tinv("ijmn")*E("ijkl")*dQods("kl");
          temp3lower.null_indices();
          lower = temp3lower.trace();	  
          lower = lower - hardMod_;
	  //opserr << " 1st hardMod_ " <<  hardMod_ << "\n";

          temp5 = (residual("ij") * Tinv("ijmn")) * dFods("mn");
          temp6 = temp5.trace();
	  //The same as the above but more computation
          //temp5 = dFods("mn") * residual("ij") * Tinv("ijmn");
          //temp6 = temp5.trace();

          upper = Fold - temp6;

	  //================================================================================
	  //dlambda
          delta_lambda = upper / lower;
          Delta_lambda = Delta_lambda + delta_lambda;

	  // Zhaohui_____10-01-2000 not sure xxxxxxxxxxxxxxxxxxx
	  dPlasticStrain = dQods("kl") * delta_lambda;
	  PlasticStrain = PlasticStrain + dPlasticStrain;

          //::printf(" >> %d  Delta_lambda = %.8e", step_counter, Delta_lambda);
          // stari umesto dQodsextended za stari = dQods
          //outfornow          dsigma =
          //outfornow            ((residual("ij")*Tinv("ijmn"))+
          //outfornow            ((E("ijkl")*dQodsextended("kl"))*Tinv("ijmn")*Delta_lambda) )*-1.0;
          //::printf("    Delta_lambda = %.8e\n", Delta_lambda);
          dsigma = ( (residual("ij")*Tinv("ijmn") )+
                   ( (E("ijkl")*dQods("kl"))*Tinv("ijmn")*delta_lambda) )*(-1.0); //*-1.0???

          dsigma.null_indices();
	  //dsigma.reportshortpqtheta("\n......dsigma ");
          //::printf("  .........   in NR loop   delta_lambda = %.16e\n", delta_lambda);
          //::printf("  .........   in NR loop   Delta_lambda = %.16e\n", Delta_lambda);
          
	  //dq_ast = delta_lambda * h_ * just_this_PP;
          //Dq_ast += dq_ast;
          
	  //dS = delta_lambda * h_s ;   // Increment to the internal var
          //S  = EP_PredictorEPS.getScalarVar(1); // Get the old value of the internal var
	  //new_S = S + dS;
          //EP_PredictorEPS.setScalarVar(1, new_S); 
          
          //Evolve the surfaces and hardening vars
	  int NS = EP_PredictorEPS.getNScalarVar();
	  int NT = EP_PredictorEPS.getNTensorVar();
  	  
	  double dS = 0;
  	  double S  = 0;
  	  //double new_S = 0; 
  	  
  	  stresstensor dT;
  	  stresstensor T;
  	  stresstensor new_T;
  	  
	  //if ( delta_lambda < 0) delta_lambda = 0;	      
	  int ii;
	  for (ii = 1; ii <= NS; ii++) {
             dS = delta_lambda * h_s[ii-1] ;             // Increment to the scalar internal var
             S  = EP_PredictorEPS.getScalarVar(ii);      // Get the old value of the scalar internal var
             EP_PredictorEPS.setScalarVar(ii, S + dS );  // Update internal scalar var
	  }

	  for (ii = 1; ii <= NT; ii++) {
	     dT = h_t[ii-1] * delta_lambda;            // Increment to the tensor internal var
             T  = EP_PredictorEPS.getTensorVar(ii);     // Get the old value of the tensor internal var
             new_T = T + dT;
             EP_PredictorEPS.setTensorVar(ii, new_T );	// Update tensorial scalar var
          }	      
 	  
	  //=======          Dq_ast = Delta_lambda * h_ * just_this_PP;
          //q_ast = q_ast_entry + Dq_ast;
          
	  //::fprintf(stdout," ######## step = %3d --> Dq_ast=%.8e q_ast=%.8e\n",
          //             step_counter,         Dq_ast,        q_ast);
          //::fprintf(stderr," ######## step = %3d --> Dq_ast=%.8e q_ast=%.8e\n",
          //             step_counter,         Dq_ast,        q_ast);
          
	  //current_lambda_set(Delta_lambda);
          //....          elastic_plastic_predictor_stress.reportshort("elplpredstress");
          //....          dsigma.reportshort("dsigma");
          
          //sigmaBack.reportshortpqtheta("\n before======== SigmaBack");
	  sigmaBack = elastic_plastic_predictor_stress + dsigma;
          //sigmaBack.deviator().reportshort("\n after ======== SigmaBack");
          //sigmaBack.reportshortpqtheta("\n after ======== SigmaBack");
	  
	  //temp trick
       if  (sigmaBack.p_hydrostatic() > 0)
       {
          //======          sigmaBack = elastic_predictor_stress + Dsigma;
          //sigmaBack.reportshortpqtheta("BE................  NR sigmaBack   ");
          //sigmaBack.reportAnim();
          //::fprintf(stdout,"Anim BEpoint0%d   = {Sin[theta]*q, p, Cos[theta]*q} \n",step_counter+1);
          ////::fprintf(stdout,"Anim BEpoint0%dP = Point[BEpoint0%d] \n",step_counter+1,step_counter+1);
          //::fprintf(stdout,"Anim   \n");
          
	  //Criterion.kappa_set( sigmaBack, q_ast) ;
          EP_PredictorEPS.setStress( sigmaBack );
          
	  //Felplpredictor = Criterion.f(sigmaBack);
          //Kai          absFelplpredictor = fabs(Felplpredictor);
          //::printf("  F_bE=%.10e (%.10e)\n", Felplpredictor,Ftolerance);
	  Felplpredictor = getYS()->f( &EP_PredictorEPS );
          //::printf(" F_BE: step=%5d  F= %.10e (%.10e)\n", step_counter, Felplpredictor, Ftolerance);
       }
       else
       {   
          sigmaBack= sigmaBack.pqtheta2stress(0.1, 0.0, 0.0);
          Felplpredictor = 0;
          backwardEPS.setStress(sigmaBack);
          backwardEPS.setConverged(TRUE);
          return backwardEPS;
       }

          //double tempkappa1 = kappa_cone_get();
          //double tempdFodeta = dFoverdeta(sigmaBack);
          //double tempdetaodkappa = detaoverdkappa(tempkappa1);
          //::printf("    h_=%.6e  xi_=%.6e, dFodeta=%.6e, detaodkappa=%.6e, hardMod_=%.6e\n", 
          //     h_, xi_,tempdFodeta, tempdetaodkappa,  hardMod_);
          //::printf("   upper = %.6e    lower = %.6e\n", upper, lower);
          //::printf(" q_ast_entry=%.6e  Dq_ast=%.6e   Delta_lambda=%.6e\n", 
          //      q_ast_entry, Dq_ast, Delta_lambda);

          // now prepare new step
          elastic_plastic_predictor_stress = sigmaBack;
          
	  //Output for plotting
	  //opserr.precision(5); 
          //opserr.width(10);
          //opserr << " " << sigmaBack.p_hydrostatic() << " ";
          //opserr << " " << sigmaBack.p_hydrostatic() << " ";
	  
	  //opserr.precision(5); 
          //opserr.width(10);
          //opserr << sigmaBack.q_deviatoric() << " ";
          //cerr << sigmaBack.q_deviatoric() << " ";
	    
	  //opserr.precision(5); 
          //opserr.width(10);
          //opserr << " Felpl " << Felplpredictor;
          //cerr << " Felpl " << Felplpredictor;

	  //opserr.precision(5); 
          //opserr.width(10);
          //opserr << " tol " << Ftolerance << endlnn;
          //cerr << " tol " << Ftolerance << " " << step_counter << endlnn;

	  //::printf("         ...........................  end of step %d\n", step_counter);// getchar();
          step_counter++;

          //int err = 0;
          //if ( getELT1() ) {
	  //   double pl_st_vol = dPlasticStrain.Iinvariant1(); //Joey 02-17-03           Iinvariant1 > 0 for compression 
          //   opserr << " Pl update " << pl_st_vol << "\n";
	  //   err = getELT1()->updateEeDm(&EP_PredictorEPS, - pl_st_vol, delta_lambda); // D > 0 --> contraction
	  //}

        }

        //opserr << " " << sigmaBack.p_hydrostatic() << " ";
        //opserr << sigmaBack.q_deviatoric() << " ";
        //opserr << " Felpl " << Felplpredictor;
        //opserr << " tol " << Ftolerance << " " << step_counter << endln;
    
        //// Update E_Young and e according to current stress state before evaluate ElasticStiffnessTensor
        int err = 0;
        if ( getELT1() ) {
	   double pl_st_vol = PlasticStrain.Iinvariant1(); //Joey 02-17-03
	   // D > 0 compressive -> Iinv > 0  -> de < 0 correct!
           err = getELT1()->updateEeDm(&EP_PredictorEPS, pl_st_vol, Delta_lambda);
	}
	
	//out07may97      while ( absFelplpredictor > Ftolerance &&
        //out07may97              step_counter <= MAX_STEP_COUNT  ); // if more iterations than prescribed
     
        //**********
        //**********
        //**********
        //**********
        if ( step_counter >= MAX_STEP_COUNT  )
        {
           //g3ErrorHandler->warning("Template3Dep::BackwardEulerEPState   Step_counter > MAX_STEP_COUNT %d iterations", MAX_STEP_COUNT );
       	   EP_PredictorEPS.setConverged( false );
	   //::exit(1);
        }

        // already set everything
	// Need to genarate Eep and set strains and stresses
        //if ( ( flag !=1) && (step_counter < MAX_STEP_COUNT) ) 	 
        if ( ( flag !=1) ) {	 

           //Return_stress = elastic_plastic_predictor_stress;
           //Criterion.kappa_set( Return_stress, q_ast) ;
           
           // Generating Consistent Stiffness Tensor Eep
           tensor I2("I", 2, def_dim_2);
           tensor I_ijkl = I2("ij")*I2("kl");
           I_ijkl.null_indices();
           tensor I_ikjl = I_ijkl.transpose0110();
	   	 

           dQods = getPS()->dQods( &EP_PredictorEPS ); // this is m_ij
           tensor temp2 = E("ijkl")*dQods("kl");
           temp2.null_indices();
           dFods = getYS()->dFods( &EP_PredictorEPS ); // this is n_ij
           d2Qoverds2 = getPS()->d2Qods2( &EP_PredictorEPS );
           	   
           tensor T = I_ikjl + E("ijkl")*d2Qoverds2("klmn")*Delta_lambda;
  	   //tensor tt = E("ijkl")*d2Qoverds2("klmn")*Delta_lambda;
	   //tt.printshort("temp tt");
	   //T = I_ikjl + tt;	    
           T.null_indices();
           tensor Tinv = T.inverse();
           
           tensor R = Tinv("ijmn")*E("ijkl");
           R.null_indices();
	   
	   /*
           // Evaluating the hardening modulus: sum of  (df/dq*) * qbar at the final stress
	   hardMod_ = 0.0;
	   //Of 1st scalar internal vars
	   if ( getELS1() ) {
	      h_s[0]  = getELS1()->h_s( &EP_PredictorEPS, getPS());
              xi_s[0] = getYS()->xi_s1( &EP_PredictorEPS );	   
   	      hardMod_ = hardMod_ + h_s[0] * xi_s[0];
	   }
	   
	   //Of 2nd scalar internal vars
	   if ( getELS2() ) {
	      h_s[1]  = getELS2()->h_s( &EP_PredictorEPS, getPS());
              xi_s[1] = getYS()->xi_s2( &EP_PredictorEPS );	   
   	      hardMod_ = hardMod_ + h_s[1] * xi_s[1];
	   }
	   
	   //Of 3rd scalar internal vars
	   if ( getELS3() ) {
	      h_s[2]  = getELS3()->h_s( &EP_PredictorEPS, getPS());
              xi_s[2] = getYS()->xi_s3( &EP_PredictorEPS );	   
   	      hardMod_ = hardMod_ + h_s[2] * xi_s[2];
	   }
	   
	   //Of 4th scalar internal vars
	   if ( getELS4() ) {
	      h_s[3]  = getELS4()->h_s( &EP_PredictorEPS, getPS());
              xi_s[3] = getYS()->xi_s4( &EP_PredictorEPS );	   
   	      hardMod_ = hardMod_ + h_s[3] * xi_s[3];
	   }
	         
	   //Of tensorial internal var
	   // 1st tensorial var
	   if ( getELT1() ) {
	     h_t[0]  = getELT1()->h_t( &EP_PredictorEPS, getPS());
	     xi_t[0] = getYS()->xi_t1( &EP_PredictorEPS );
             tensor hm = (h_t[0])("ij") * (xi_t[0])("ij");
  	     hardMod_ = hardMod_ + hm.trace();
	   }
	   
	   // 2nd tensorial var
	   if ( getELT2() ) {
	      h_t[1]  = getELT2()->h_t( &EP_PredictorEPS, getPS());
	      xi_t[1] = getYS()->xi_t2( &EP_PredictorEPS );
              tensor hm = (h_t[1])("ij") * (xi_t[1])("ij");
  	      hardMod_ = hardMod_ + hm.trace();
	   }
	   
	   // 3rd tensorial var
	   if ( getELT3() ) {
	      h_t[2]  = getELT3()->h_t( &EP_PredictorEPS, getPS());
	      xi_t[2] = getYS()->xi_t3( &EP_PredictorEPS );
              tensor hm = (h_t[2])("ij") * (xi_t[2])("ij");
  	      hardMod_ = hardMod_ + hm.trace();
	   }
	   
	   // 4th tensorial var
	   if ( getELT4() ) {
	      h_t[3]  = getELT4()->h_t( &EP_PredictorEPS, getPS());
	      xi_t[3] = getYS()->xi_t4( &EP_PredictorEPS );
              tensor hm = (h_t[3])("ij") * (xi_t[3])("ij");
  	      hardMod_ = hardMod_ + hm.trace();
	   }
	   */ 	   
    	   tensor temp3lower = dFods("ot")*R("otpq")*dQods("pq");
    	   temp3lower.null_indices();
    	   
    	   double lower = temp3lower.trace();
    	   lower = lower - hardMod_;  // h
	   //opserr << " 2nd hardMod_ " <<  hardMod_ << "\n";
    	      
    	   //tensor upper = R("pqkl")*dQods("kl")*dFods("ij")*R("ijmn");
    	   tensor upper11 = R("pqkl")*dQods("kl");
    	   tensor upper22 = dFods("ij")*R("ijmn");
	   upper11.null_indices();
	   upper22.null_indices();
    	   tensor upper = upper11("pq")*upper22("mn");    	   
	   
	   upper.null_indices();
    	   tensor Ep = upper*(1./lower);
    	   tensor Eep =  R - Ep; // elastoplastic constitutive tensor

           //Set Elasto-Plastic stiffness tensor
           EP_PredictorEPS.setEep(Eep);
	   EP_PredictorEPS.setConverged(TRUE);

	   //set plastic strain and total strain
           straintensor elastic_strain = strain_increment - PlasticStrain;  // elastic strain increment
           straintensor estrain = EP_PredictorEPS.getElasticStrain(); //get old elastic strain
           straintensor pstrain = EP_PredictorEPS.getPlasticStrain(); //get old plastic strain 
    	   
           straintensor tstrain = EP_PredictorEPS.getStrain();        //get old total strain
           pstrain = pstrain + PlasticStrain;
           estrain = estrain + elastic_strain;
           tstrain = tstrain + strain_increment;
           //opserr<< "Plastic:  Total strain" << tstrain <<endln;
           
           //Setting de_p, de_e, total plastic, elastic strain, and  total strain
           EP_PredictorEPS.setdPlasticStrain( PlasticStrain );
           EP_PredictorEPS.setdElasticStrain( elastic_strain );
           EP_PredictorEPS.setPlasticStrain( pstrain );
           EP_PredictorEPS.setElasticStrain( estrain );
           EP_PredictorEPS.setStrain( tstrain );
		
	   backwardEPS = EP_PredictorEPS;
           
	   //double f_backward =  getYS()->f( &backwardEPS );
           //opserr << "\n************  Exit Backward = " <<  f_backward << "\n";

        }
  
  }
  //return Return_stress;
  backwardEPS = EP_PredictorEPS;

  //double p = (backwardEPS.getStress()).p_hydrostatic();
  //double ec = (backwardEPS.getec()) - (backwardEPS.getLam()) * log( p / (backwardEPS.getpo()) );
  //double st = (backwardEPS.getStrain()).Iinvariant1();
  //double pl_s = (backwardEPS.getPlasticStrain()).Iinvariant1();
  //double dpl_s = (backwardEPS.getdPlasticStrain()).Iinvariant1();
  //cerr << "P BE p=" << p << " ec " << ec << " e " << backwardEPS.gete() << " psi " << (backwardEPS.gete() - ec) << " strain " << st << " t_pl " << pl_s << " d_pl " << dpl_s << "\n";

  return backwardEPS;
}




//================================================================================
// New EPState using Forward Euler Subincrement Euler Algorithm
//================================================================================
EPState Template3Dep::FESubIncrementation( const straintensor & strain_increment,
                                           int number_of_subincrements)                                                 
{
    
    EPState old_EPS( *(this->getEPS()) ); 
    EPState FESI_EPS( *(this->getEPS()) ); 
    //NDMaterial MP( material_point );
    //NDMaterial *MP = this->getCopy();
    //opserr << "in FESubIncrementation MP " << MP;

    stresstensor back_stress;	       
    stresstensor begin_stress = this->getEPS()->getStress();
    //stresstensor begin_stress = start_stress;
    //::fprintf(stderr,"{");
  
    double sub = 1./( (double) number_of_subincrements );
    //stresstensor elastic_subincremental_stress = stress_increment * sub;

    straintensor tempp = strain_increment;
    straintensor elastic_subincremental_strain = tempp * sub;
    straintensor total_strain = elastic_subincremental_strain;
    //elastic_subincremental_stress.reportshort("SUB INCREMENT in stresses\n");
    //opserr << "INCREMENT strain " << strain_increment << endlnn ;
    //opserr << "SUB INCREMENT strain " << elastic_subincremental_strain << endlnn ;
   
    for( int steps=0 ; steps < number_of_subincrements ; steps++ ){

        //start_stress.reportshort("START stress\n");
        FESI_EPS = ForwardEulerEPState( elastic_subincremental_strain);
        
	// Update the EPState in Template3Dep
	this->setEPS( FESI_EPS );
	                      
        back_stress = FESI_EPS.getStress();
	//opserr.unsetf(ios::showpos);
	//opserr << setw(4);
        //opserr << "Step No. " << steps << "  ";

	//opserr.setPrecision(SCIENTIFIC);
	//opserr.precision(3);
	
	// opserr
	// opserr.setf(ios::showpos);
	// opserr.precision(3);

	//opserr << setw(7);
	//opserr << "p " << back_stress.p_hydrostatic() << "  "; 
	//opserr << setw(7);
	//opserr << "q " << back_stress.q_deviatoric() << "  "; 
	//opserr << setw(7);
	//opserr << " theta " << back_stress.theta() << "  "; 
	//opserr << setw(7);
	//opserr << "alfa1 " << FESI_EPS.getScalarVar(1) << "  "; 
	//opserr << setw(7);
	//opserr << "f = " << getYS()->f( &FESI_EPS ) << endlnn;
  
        //begin_stress = back_stress;  
        //total_strain = total_strain + elastic_subincremental_strain;

     }

     //    ::fprintf(stderr,"}");
     this->setEPS( old_EPS );
     return FESI_EPS;

}
  


//================================================================================
// New EPState using Backward Euler Subincrement Euler Algorithm
//================================================================================
EPState Template3Dep::BESubIncrementation( const straintensor & strain_increment,
                                           int number_of_subincrements)                                                 
{
    EPState old_EPS( *(this->getEPS()) ); 
    EPState BESI_EPS( *(this->getEPS()) ); 
    straintensor strain_incr = strain_increment;
    //NDMaterial MP( material_point );
    //NDMaterial *MP = getCopy();
    //opserr << "in FESubIncrementation MP " << MP;

    stresstensor back_stress;
    stresstensor begin_stress = old_EPS.getStress();
    //stresstensor begin_stress = getEPS()->getStress();
    //stresstensor begin_stress = start_stress;
    //::fprintf(stderr,"{");
  
    double sub = 1./( (double) number_of_subincrements );
    //stresstensor elastic_subincremental_stress = stress_increment * sub;

    straintensor elastic_subincremental_strain = strain_incr * sub;
    straintensor total_strain = elastic_subincremental_strain;
    //elastic_subincremental_stress.reportshort("SUB INCREMENT in stresses\n");
   
    for( int steps=0 ; steps < number_of_subincrements ; steps++ ){

        //start_stress.reportshort("START stress\n");
        BESI_EPS = BackwardEulerEPState( elastic_subincremental_strain);
        if ( BESI_EPS.getConverged() )
      	  this->setEPS( BESI_EPS );
	else {
      	  //this->setEPS( BESI_EPS );
          opserr << "Template3Dep::BESubIncrementation  failed to converge at " << steps << "th(of " 
		 << number_of_subincrements << "step sub-BackwardEuler Algor.\n";
	  //exit(1);
          //g3ErrorHandler->fatal("Template3Dep::BESubIncrementation  failed to converge using %d step sub-BackwardEuler Algor.", number_of_subincrements );
	  //exit(1);
	                      
        //back_stress = BESI_EPS.getStress();
	//opserr.unsetf(ios::showpos);
	//opserr << setw(4);
        //opserr << "Step No. " << steps << "  ";
	//
	//opserr.setf(ios::scientific);
	//opserr.setf(ios::showpos);
	//opserr.precision(3);
	//opserr << setw(7);
	//opserr << " back-stress  p " << back_stress.p_hydrostatic() << "  "; 
	//opserr << setw(7);
	//opserr << "q " << back_stress.q_deviatoric() << "  ";
	//opserr << setw(7);
	//opserr << "alfa1 " << BESI_EPS.getScalarVar(1) << "  ";
	//opserr << setw(7);
	//opserr << "f = " << MP->YS->f( &BESI_EPS ) << "  "<< endlnn;
  

	  //	opserr.setf(ios::scientific);
	  // opserr.setf(ios::showpos);
	  // opserr.precision(3);
	opserr.width(7);
	opserr << "\n intraStep: begin-stress  p " << begin_stress.p_hydrostatic() << "  "; 
	opserr.width(7);
	opserr << "q " << begin_stress.q_deviatoric() << "  ";
	opserr.width(7);
	opserr << "alfa1 " << old_EPS.getScalarVar(1) << "  ";
	//opserr.width(7);
	//opserr << "f = " << MP->YS->f( &BESI_EPS ) << "  "<< endlnn;

	opserr << "strain increment " << strain_incr << endln; 
	//fflush(opserr);
	break;
	//exit(1);
        //begin_stress = back_stress;  
        //total_strain = total_strain + elastic_subincremental_strain;
	}

     }

     //    ::fprintf(stderr,"}");
     //this->setEPS( old_EPS );
     return BESI_EPS;

}
  


////================================================================================
//// Routine used to generate elastic stiffness tensor E
////================================================================================
//tensor Template3Dep::ElasticStiffnessTensor( double E, double nu) const
//  {
//    tensor ret( 4, def_dim_4, 0.0 );
//
//    tensor I2("I", 2, def_dim_2);
//    tensor I_ijkl = I2("ij")*I2("kl");
//    I_ijkl.null_indices();
//    tensor I_ikjl = I_ijkl.transpose0110();
//    tensor I_iljk = I_ijkl.transpose0111();
//    tensor I4s = (I_ikjl+I_iljk)*0.5;
//
//    // Building elasticity tensor
//    ret = (I_ijkl*((E*nu*2.0)/(2.0*(1.0+nu)*(1-2.0*nu)))) + (I4s*(E/((1.0+nu))));
//
//    return ret;
//  }
//
//
////================================================================================
//// Routine used to generate elastic compliance tensor D
////================================================================================
//
//tensor Template3Dep::ElasticComplianceTensor( double E, double nu) const
//  {
//    if (E == 0) {
//      opserr << endln << "Ey = 0! Can't give you D!!" << endln;
//      exit(1);
//    }
//
//    tensor ret( 4, def_dim_4, 0.0 );
//    //tensor ret;
//    
//    tensor I2("I", 2, def_dim_2);
//    tensor I_ijkl = I2("ij")*I2("kl");
//    I_ijkl.null_indices();
//    tensor I_ikjl = I_ijkl.transpose0110();
//    tensor I_iljk = I_ijkl.transpose0111();
//    tensor I4s = (I_ikjl+I_iljk)*0.5;
//
//    // Building compliance tensor
//    ret = (I_ijkl * (-nu/E)) + (I4s * ((1.0+nu)/E));
//
//    return ret;
//  }
//

//================================================================================
// trying to find intersection point				  
// according to M. Crisfield's book				  
// "Non-linear Finite Element Analysis of Solids and Structures " 
// chapter 6.6.1 page 168.                                        
//================================================================================
stresstensor Template3Dep::yield_surface_cross(const stresstensor & start_stress,
                                               const stresstensor & end_stress)
{
  // Bounds
  double x1 = 0.0;
  double x2 = 1.0;
  
  // accuracy
  double const TOL = 1.0e-9;
  //opserr << "start_stress "<< start_stress;
  //opserr << "end_stress " << end_stress;
  //end_stress.reportshortpqtheta("end stress");
  
  double a = zbrentstress( start_stress, end_stress, x1, x2, TOL ); // Defined below 
  // ::printf("\n****\n a = %lf \n****\n",a);
  
  stresstensor delta_stress = end_stress - start_stress;
  stresstensor intersection_stress = start_stress + delta_stress * a;
  //***  intersection_stress.reportshort("FINAL Intersection stress\n");

  return intersection_stress;

}


//================================================================================
// Routine used by yield_surface_cross to 
// find the stresstensor at cross point   
//================================================================================

double Template3Dep::zbrentstress(const stresstensor & start_stress,
                                  const stresstensor & end_stress,
                                  double x1, double x2, double tol)
{
  double EPS = d_macheps();

  int iter;
  double a=x1;
  double b=x2;
  double c=0.0;
  double d=0.0;
  double e=0.0;
  double min1=0.0;
  double min2=0.0;
  double fa=func(start_stress, end_stress, a);
  double fb=func(start_stress, end_stress, b);
  //opserr << "fb = " << fb;

  double fc=0.0;
  double p=0.0;
  double q=0.0;
  double r=0.0;
  double s=0.0;
  double tol1=0.0;
  double xm=0.0;
  //::printf("\n############# zbrentstress iterations --\n");
  if (fb*fa > 0.0)
    {
      ::printf("\a\nRoot must be bracketed in ZBRENTstress\n");
      ::exit(1);
    }
  fc=fb;
  for ( iter=1 ; iter<=ITMAX ; iter++ )
  {
      //::printf("iter No. = %d  ;  b = %16.10lf\n", iter, b);
      if (fb*fc > 0.0)
        {
          c=a;
          fc=fa;
          e=d=b-a;
        }
      if (fabs(fc) < fabs(fb))
        {
          a=b;
          b=c;
          c=a;
          fa=fb;
          fb=fc;
          fc=fa;
        }
      tol1=2.0*EPS*fabs(b)+0.5*tol;
      xm=0.5*(c-b);
      if (fabs(xm) <= tol1 || fb == 0.0) return b;
      if (fabs(e) >= tol1 && fabs(fa) > fabs(fb))
        {
          s=fb/fa;
          if (a == c)
            {
              p=2.0*xm*s;
              q=1.0-s;
            }
          else
            {
              q=fa/fc;
              r=fb/fc;
              p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
              q=(q-1.0)*(r-1.0)*(s-1.0);
            }
          if (p > 0.0)  q = -q;
          p=fabs(p);
          min1=3.0*xm*q-fabs(tol1*q);
          min2=fabs(e*q);
          if (2.0*p < (min1 < min2 ? min1 : min2))
            {
              e=d;
              d=p/q;
            }
          else
            {
              d=xm;
              e=d;
            }
        }
      else
        {
          d=xm;
          e=d;
        }
      a=b;
      fa=fb;
      if (fabs(d) > tol1)
        b += d;
      else		    
        b += (xm > 0.0 ? fabs(tol1) : -fabs(tol1));
      fb=func(start_stress, end_stress, b);
  }
  //::printf("\a\nMaximum number of iterations exceeded in zbrentstress\n");
  return 0.0; // this just to full the compiler because of the warnings
}


//================================================================================
// routine used by zbrentstress, takes an alfa and returns the
// yield function value for that alfa
//================================================================================
double Template3Dep::func(const stresstensor & start_stress,
                          const stresstensor & end_stress,
                          double alfa )
{
   
   //EPState *tempEPS = getEPS()->newObj();
   EPState tempEPS( (*this->getEPS()) );
   stresstensor delta_stress = end_stress - start_stress;
   stresstensor intersection_stress = start_stress + delta_stress*alfa;

   tempEPS.setStress(intersection_stress); 
   
   //opserr << "*tempEPS" << *tempEPS;
   
   double f = getYS()->f( &tempEPS );
   return f; 
}

//================================================================================
// Overloaded Insertion Operator
// prints an EPState's contents 
//================================================================================
OPS_Stream& operator<< (OPS_Stream& os, const Template3Dep & MP)
{
    os << endln << "Template3Dep: " << endln;
    os << "\ttag: " << MP.getTag() << endln;
    os << "=================================================================" << endln;
    MP.getYS()->print();
    MP.getPS()->print();
    MP.getEPS()->print();

    os << endln << "Scalar Evolution Laws: " << endln; 
    if ( MP.ELS1 ){
       os << "\nFor 1st scalar var:\n";
       MP.ELS1->print();
    }
    
    if ( MP.ELS2 ){
       os << "\nFor 2nd scalar var:\n";
       MP.ELS2->print();
    }
    
    if ( MP.ELS3 ){
       os << "\nFor 3rd scalar var:\n";
       MP.ELS3->print();
    }
    
    if ( MP.ELS4 ){
       os << "\nFor 4th scalar var:\n";
       MP.ELS4->print();
    }
    

    os << endln << "Tensorial Evolution Laws: " << endln; 
    if ( MP.ELT1 ){
       os << "\nFor 1st tensorial var:\n";
       MP.ELT1->print();
    }
    if ( MP.ELT2 ){
       os << "\nFor 2nd tensorial var:\n";
       MP.ELT2->print();
    }
    if ( MP.ELT3 ){
       os << "\nFor 3rd tensorial var:\n";
       MP.ELT3->print();
    }
    if ( MP.ELT4 ){
       os << "\nFor 4th tensorial var:\n";
       MP.ELT4->print();
    }

    os << endln;           
    return os;
}  


#endif

