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

// Ref: Dafalias and Manzari 2004: J. Eng. Mech. 130(6), pp 622-634
// Parameters:
//  1-e0:         initial void ratio at zero strain;
//  2- e_r:       reference void for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  3- lambda_c:  parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  4- xi:        parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  5- Pat:       atmospheris pressure for critical state line, ec = e0 - lambda_c*(pc/Pat)^xi;
//  6- m:         parameter in the yield function;
//  7- M:         critical state stress ration;
//  8- cc:        tension-compression strength ratio;
//  9- A0:        parameter;
// 10- nd         parameter;
// 11- alpha:     "back-stress" tensor in yield function; (the 1st tensorial internal variable);
// 12- z:         fabric dilatancy internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_PF_CPP
#define DM04_PF_CPP

#include "DM04_PF.h"

straintensor DM04_PF::DM04m;
stresstensor DM04_PF::DM04temp;

//================================================================================
DM04_PF::DM04_PF(int e0_which_in, int index_e0_in,
                 int e_r_which_in, int index_e_r_in, 
                 int lambda_c_which_in, int index_lambda_c_in,
                 int xi_which_in, int index_xi_in,
                 int Pat_which_in, int index_Pat_in,
                 int m_which_in, int index_m_in,
                 int M_cal_which_in, int index_M_cal_in,
                 int cc_which_in, int index_cc_in,
                 int A0_which_in, int index_A0_in,
                 int nd_which_in, int index_nd_in,
                 int alpha_which_in, int index_alpha_in,
                 int z_which_in, int index_z_in)
: e0_which(e0_which_in), index_e0(index_e0_in), 
  e_r_which(e_r_which_in), index_e_r(index_e_r_in), 
  lambda_c_which(lambda_c_which_in), index_lambda_c(index_lambda_c_in),
  xi_which(xi_which_in), index_xi(index_xi_in),
  Pat_which(Pat_which_in), index_Pat(index_Pat_in),
  m_which(m_which_in), index_m(index_m_in),
  M_cal_which(M_cal_which_in), index_M_cal(index_M_cal_in),
  cc_which(cc_which_in), index_cc(index_cc_in),
  A0_which(A0_which_in), index_A0(index_A0_in),
  nd_which(nd_which_in), index_nd(index_nd_in),
  alpha_which(alpha_which_in), index_alpha(index_alpha_in),
  z_which(z_which_in), index_z(index_z_in)
{

}

//================================================================================
DM04_PF::~DM04_PF() 
{  

}

//================================================================================
PlasticFlow* DM04_PF::newObj() 
{  
     PlasticFlow  *new_PF = new DM04_PF(e0_which, index_e0,
                                        e_r_which, index_e_r,
                                        lambda_c_which, index_lambda_c,
                                        xi_which, index_xi,
                                        Pat_which, index_Pat,
                                        m_which, index_m,
                                        M_cal_which, index_M_cal,
                                        cc_which, index_cc,
                                        A0_which, index_A0,
                                        nd_which, index_nd,
                                        alpha_which, index_alpha,
                                        z_which, index_z);
     
     return new_PF;
}

//================================================================================
const straintensor& DM04_PF::PlasticFlowTensor(const stresstensor& Stre, 
                                               const straintensor& Stra, 
                                               const MaterialParameter &MaterialParameter_in) const
{
	BJtensor KroneckerI("I", 2, def_dim_2);

    double e0 = gete0(MaterialParameter_in);
    double e_r = gete_r(MaterialParameter_in);
    double lambda_c = getlambda_c(MaterialParameter_in);
    double xi = getxi(MaterialParameter_in);
    double Pat = getPat(MaterialParameter_in);
    double m = getm(MaterialParameter_in);
    double M_cal = getM_cal(MaterialParameter_in);
    double cc = getcc(MaterialParameter_in);
    double A0 = getA0(MaterialParameter_in);
    double nd = getnd(MaterialParameter_in);    

    stresstensor alpha = getalpha(MaterialParameter_in);
    //alpha.print("a", " ");
    stresstensor z = getz(MaterialParameter_in);
    //z.print("z", " ");

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    stresstensor alpha_d;
    stresstensor alpha_d_alpha;
    double b0 = 0.0;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnd = 1.0;
    double ad = 0.0;
    double z_n = 0.0;
    double A_d = 0.0;
    double B_cal = 1.0;
    double C_cal = 0.0;
    double D_cal = 0.0;
    stresstensor s_bar;
    double _s_bar_ = 0.0;
    double epsilon_v = 0.0;
    double e = e0;
    double J3D;
    double cos3theta = 0.0;
        
    //if (p != 0.0 && m != 0.0)
    //    n = (s * (1.0/p) - alpha) *(1.0/(sqrt(2.0/3.0)*m));
    
    s_bar = Stre.deviator() - (alpha *p);
    _s_bar_ = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    if (_s_bar_ != 0.0)
        n = s_bar * (1.0/_s_bar_);
   
    //double J2D = n.Jinvariant2();
    //J3D = n.Jinvariant3();
    //if ( J2D != 0.0)
      //cos3theta = -1.5*sqrt(3.0) *J3D / pow(J2D, 1.5);  
    
    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;
    
    if (cos3theta > 1.0) cos3theta = 1.0;
    if (cos3theta < -1.0) cos3theta = -1.0;
    
    g = getg(cc, cos3theta);

    if ( (p/Pat) >= 0.0 )
      ec = getec(e_r, lambda_c, xi, Pat, p);
    
    epsilon_v = Stra.Iinvariant1();
    e = e0 + (1 + e0) *epsilon_v;
    
    stateParameter = e - ec;
    
    expnd = exp( nd *stateParameter );
    ad = sqrt(2.0/3.0) * ( g *M_cal *expnd - m);
    alpha_d  = n *ad;
    alpha_d_alpha = alpha_d - alpha;
    
    z_n = (z("ij")*n("ij")).trace();
    if (z_n < 0.0) 
      z_n = 0.0;
    A_d = A0 * (1.0 + z_n);

    D_cal = (alpha_d_alpha("ij")*n("ij")).trace() *A_d;

    B_cal = 1.0 + 1.5 *((1.0-cc)/cc) *g *cos3theta;
    C_cal = 3.0 *sqrt(1.5) *((1.0-cc)/cc) *g;
    
    stresstensor n_n = n("ik")*n("kj");
        n_n.null_indices();

    // note diffrent 'positive-negtive' since we assume extension (dilatant) positive 
    // which is diffrenet from the Ref. 
    DM04m = (n *B_cal) + (n_n *C_cal) + (KroneckerI *(( - C_cal - D_cal)/3.0));
    
                       
    return DM04m;
}

// to get e0
//================================================================================
double DM04_PF::gete0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e0_which, index_e0);
}

// to get e_r
//================================================================================
double DM04_PF::gete_r(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, e_r_which, index_e_r);
}

// to get lambda_c
//================================================================================
double DM04_PF::getlambda_c(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, lambda_c_which, index_lambda_c);
}

// to get xi
//================================================================================
double DM04_PF::getxi(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, xi_which, index_xi);

}

// to get Pat
//================================================================================
double DM04_PF::getPat(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, Pat_which, index_Pat);
}

// to get m
//================================================================================
double DM04_PF::getm(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, m_which, index_m);
}

// to get M
//================================================================================
double DM04_PF::getM_cal(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, M_cal_which, index_M_cal);
}

// to get c
//================================================================================
double DM04_PF::getcc(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, cc_which, index_cc);
}

// to get A0
//================================================================================
double DM04_PF::getA0(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, A0_which, index_A0);

}

// to get n_d
//================================================================================
double DM04_PF::getnd(const MaterialParameter &MaterialParameter_in) const
{
    return getParameters(MaterialParameter_in, nd_which, index_nd);
}

// to get alpha
//================================================================================
const stresstensor& DM04_PF::getalpha(const MaterialParameter &MaterialParameter_in) const
{
	if ( alpha_which == 2 && index_alpha <= MaterialParameter_in.getNum_Internal_Tensor() && index_alpha > 0) {
		DM04_PF::DM04temp = MaterialParameter_in.getInternal_Tensor(index_alpha-1);
		return DM04_PF::DM04temp;
	}
	else {
		cout << "DM04_PF: Invalid Input. " << endl;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& DM04_PF::getz(const MaterialParameter &MaterialParameter_in) const
{
	if ( z_which == 2 && index_z <= MaterialParameter_in.getNum_Internal_Tensor() && index_z > 0) {
		DM04_PF::DM04temp = MaterialParameter_in.getInternal_Tensor(index_z-1);
		return DM04_PF::DM04temp;
	}
	else {
		cout << "DM04_PF: Invalid Input. " << endl;
		exit (1);
	}
}


//================================================================================
double DM04_PF::getParameters(const MaterialParameter &MaterialParameter_in, int parIndex, int which) const
{
	if ( parIndex == 0 && which <= MaterialParameter_in.getNum_Material_Parameter() && which > 0)
		return MaterialParameter_in.getMaterial_Parameter(which-1);
	else if ( parIndex == 1 && which <= MaterialParameter_in.getNum_Internal_Scalar() && which > 0)
		return MaterialParameter_in.getInternal_Scalar(which-1);
	else {
		cout << "DM04_PF: Invalid Input. " << parIndex << " and " << which << endl;
		exit (1);
	}
} 


//================================================================================
double DM04_PF::getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda_c * pow(p_c/Pat, xi);
    else 
      cout << "Warning: DM04_PF - 'p_c/Pat' less than zero! " << endl;

    return ee;
}

//================================================================================
double DM04_PF::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}


#endif




