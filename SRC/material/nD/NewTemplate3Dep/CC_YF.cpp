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

#ifndef CC_YF_CPP
#define CC_YF_CPP

#include "CC_YF.h"

stresstensor CC_YF::CCst;

//================================================================================
CC_YF::CC_YF(int M_which_in, int index_M_in, 
             int p0_which_in, int index_p0_in)
: M_which(M_which_in), index_M(index_M_in), 
  p0_which(p0_which_in), index_p0(index_p0_in)
{

}

//================================================================================
CC_YF::~CC_YF() 
{

}

//================================================================================
YieldFunction* CC_YF::newObj() 
{
	YieldFunction  *new_YF = new CC_YF(M_which, index_M, p0_which, index_p0);

	return new_YF;
}

//================================================================================
double CC_YF::YieldFunctionValue( const stresstensor& Stre, 
                                  const MaterialParameter &MaterialParameter_in ) const
{
	// f = q*q - M*M*p*(po - p) = 0
	
	double M = getM(MaterialParameter_in);
	double p0 = getP0(MaterialParameter_in);
	double p = Stre.p_hydrostatic();
	double q = Stre.q_deviatoric();
	
	return  q*q - M*M*p*(p0 - p);
}

//================================================================================
const stresstensor& CC_YF::StressDerivative(const stresstensor& Stre, 
                                            const MaterialParameter &MaterialParameter_in) const
{
	double M = getM(MaterialParameter_in);
	double p0 = getP0(MaterialParameter_in);
	double p = Stre.p_hydrostatic();
	double q = Stre.q_deviatoric();
	double dFoverdp = -M*M*( p0 - 2.0*p );
	double dFoverdq = 2.0*q;
	BJtensor DpoDs = Stre.dpoverds();
	if (q != 0.0) {
		BJtensor DqoDs = Stre.dqoverds();
		CCst = DpoDs  *dFoverdp + DqoDs  *dFoverdq;
	}
	else
		CCst = DpoDs  *dFoverdp;
	
	return CCst;
}

//================================================================================
double CC_YF::InScalarDerivative(const stresstensor& Stre, 
                                 const MaterialParameter &MaterialParameter_in, 
                                 int which) const
{
	if (p0_which == 1 && which == index_p0) {
		double M = getM(MaterialParameter_in);
		double p = Stre.p_hydrostatic();
		return  (-1.0)*M*M*p;
	}
	else {
		cout << "Warning!! CC_YF: Invalid Input Parameter. " << endl;
		exit (1);
	}
}

//================================================================================
int CC_YF::getNumInternalScalar() const
{
	return 1;
}

//================================================================================
int CC_YF::getNumInternalTensor() const
{
	return 0;
}

//================================================================================   
int CC_YF::getYieldFunctionRank() const
{
	return 2;
}

//================================================================================   
double CC_YF::getM(const MaterialParameter &MaterialParameter_in) const
{
	// to get M
	if ( M_which == 0 && index_M <= MaterialParameter_in.getNum_Material_Parameter() && index_M > 0 )
		return MaterialParameter_in.getMaterial_Parameter(index_M-1);
	else {
		cout << "Warning!! CC_YF: Invalid Input (M). " << endl;
		exit (1);
	}
}

//================================================================================ 
double CC_YF::getP0(const MaterialParameter &MaterialParameter_in) const
{
	//to get P0
	if ( p0_which == 1 && index_p0 <= MaterialParameter_in.getNum_Internal_Scalar() & index_p0 > 0 )
		return MaterialParameter_in.getInternal_Scalar(index_p0-1);
	else {
		cout << "Warning!! CC_YF: Invalid Input (po). " << endl;
		exit (1);
	}
}


#endif

