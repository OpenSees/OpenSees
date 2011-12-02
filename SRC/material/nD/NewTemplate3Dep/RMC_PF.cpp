///////////////////////////////////////////////////////////////////////////////
//   COPYLEFT (C): Woody's viral GPL-like license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              
// CLASS:             
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         
// DESIGNER:          Zhao Cheng, Boris Jeremic
// PROGRAMMER:        Zhao Cheng, 
// DATE:              Fall 2005
// UPDATE HISTORY:    
//
///////////////////////////////////////////////////////////////////////////////
//

#ifndef RMC_PF_CPP
#define RMC_PF_CPP

#include "RMC_PF.h"

straintensor RMC_PF::RMCm;

//================================================================================
RMC_PF::RMC_PF(int dilatant_which_in, int index_dilatant_in, int r_which_in, int index_r_in)
: dilatant_which(dilatant_which_in), index_dilatant(index_dilatant_in), 
  r_which(r_which_in), index_r(index_r_in)
{

}

//================================================================================
RMC_PF::~RMC_PF() 
{  

}

//================================================================================
PlasticFlow* RMC_PF::newObj() 
{  
     PlasticFlow  *new_PF = new RMC_PF(dilatant_which, index_dilatant, 
                                      r_which, index_r);
     
     return new_PF;
}

//================================================================================
const straintensor& RMC_PF::PlasticFlowTensor(const stresstensor &Stre, 
                                              const straintensor &Stra, 
                                              const MaterialParameter &MaterialParameter_in) const
{
	double d = getdilatant(MaterialParameter_in);
	double r = getr(MaterialParameter_in);
	
	double q = Stre.q_deviatoric();
	double theta = Stre.theta();
	
	// g = f1/f2, 1/g = f2/f1;
	// d(1/g) = d(f2/f1) = (df2*f1 - df1*f2)/(f1)^2;
	double f1 = RoundedFunctionf1(theta, r);
	double f2 = RoundedFunctionf2(theta, r);
	double df1 = RoundedFunctiondf1(theta, r);
	double df2 = RoundedFunctiondf2(theta, r);
	double dginv = (df2*f1 - df1*f2) /(f1*f1);
	
	double dfodp = -3.0*d;
	double dfodq = f1/(f2*sqrt(3.0));
	double dfodtheta = dginv*q/sqrt(3.0);
	
	RMCm = Stre.dpoverds() *dfodp + Stre.dqoverds() *dfodq + Stre.dthetaoverds() *dfodtheta;
		
    return RMCm;
}

//================================================================================
double RMC_PF::getdilatant(const MaterialParameter &MaterialParameter_in) const
{
	// to dilatant
	if ( dilatant_which == 0) {
		if ( index_dilatant <= MaterialParameter_in.getNum_Material_Parameter() && index_dilatant > 0)
			return MaterialParameter_in.getMaterial_Parameter(index_dilatant-1); 
		else {
			opserr << "RMC_PF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( dilatant_which == 1) {
		if ( index_dilatant <= MaterialParameter_in.getNum_Internal_Scalar() && index_dilatant > 0)
			return MaterialParameter_in.getInternal_Scalar(index_dilatant-1); 
		else {
			opserr << "RMC_PF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "RMC_PF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================   
double RMC_PF::getr(const MaterialParameter &MaterialParameter_in) const
{
	// to get r
	if ( r_which == 0) {
		if ( index_r <= MaterialParameter_in.getNum_Material_Parameter() && index_r > 0)
			return MaterialParameter_in.getMaterial_Parameter(index_r-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
	}
	else if ( r_which == 1) {
		if ( index_r <= MaterialParameter_in.getNum_Internal_Scalar() && index_r > 0)
			return MaterialParameter_in.getInternal_Scalar(index_r-1); 
		else {
			opserr << "RMC_YF: Invalid Input. " << endln;
			exit (1);
		}
    }
	else {
		opserr << "RMC_YF: Invalid Input. " << endln;
		exit(1);
	}
}

//================================================================================   
double RMC_PF::RoundedFunctionf1(double s, double r) const
{
	double t1 = sqrt(-4.0*r + 5.0*r*r + 4.0*(1.0 - r*r)*pow(cos(s),2));
	return 2.0*(1.0 - r*r)*cos(s) + (-1.0 + 2.0*r) * t1;
}

//================================================================================   
double RMC_PF::RoundedFunctionf2(double s, double r) const
{
	return pow(-1.0 + 2.0*r,2) + 4.0*(1.0 - r*r)*pow(cos(s),2);
}

//================================================================================   
double RMC_PF::RoundedFunctiondf1(double s, double r) const
{
	double t1 = sqrt(-4.0*r + 5.0*r*r + 4.0*(1 - r*r)*pow(cos(s),2));	
	return -2.0*(1.0 - r*r)*sin(s) - (4.0*(-1.0 + 2.0*r)*(1 - r*r)*cos(s)*sin(s))/t1;
}

//================================================================================   
double RMC_PF::RoundedFunctiondf2(double s, double r) const
{
	return -8.0*(1 - r*r)*cos(s)*sin(s);
}


#endif
