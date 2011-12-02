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
//  1- e0:        initial void ratio;
//  2- e_r:       reference void for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  3- lambda_c:  parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  4- xi:        parameter for critical state line, ec = e_r - lambda_c*(pc/Pat)^xi;
//  5- Pat:       atmospheris pressure for critical state line, ec = e0 - lambda_c*(pc/Pat)^xi;
//  6- m:         parameter in the yield function;
//  7- M:         critical state stress ration;
//  8- cc:        tension-compression strength ratio;
//  9- nb:        parameter;
// 10- h0:        parameter;
// 11- ch:        parameter;
// 12- G0:        parameter in the elastic part
// 13- alpha:     "back-stress ratio" tensor in yield function; (the 1st tensorial internal variable);
// 14- z:         fabric dilatancy internal tensor (the 2nd tensorial internal variable); 

#ifndef DM04_alpha_Eij_CPP
#define DM04_alpha_Eij_CPP

#include "DM04_alpha_Eij.h"

stresstensor DM04_alpha_Eij::DM04_alpha_t;

DM04_alpha_Eij::DM04_alpha_Eij(int e0_index_in,
                               int e_r_index_in,
                               int lambda_c_index_in,
                               int xi_index_in,
                               int Pat_index_in,
                               int m_index_in,
                               int M_cal_index_in,
                               int cc_index_in,
                               int nb_index_in,
                               int h0_index_in,
                               int ch_index_in,
                               int G0_index_in,
                               int alpha_index_in,
                               int z_index_in)
: e0_index(e0_index_in),
  e_r_index(e_r_index_in), 
  lambda_c_index(lambda_c_index_in),
  xi_index(xi_index_in),
  Pat_index(Pat_index_in),
  m_index(m_index_in),  
  M_cal_index(M_cal_index_in),
  cc_index(cc_index_in),
  nb_index(nb_index_in),
  h0_index(h0_index_in),   
  ch_index(ch_index_in),
  G0_index(G0_index_in),
  alpha_index(alpha_index_in),
  z_index(z_index_in) 
{
   a_index = 0;
   stresstensor zT;
   alpha_in = zT;
}

TensorEvolution* DM04_alpha_Eij::newObj()
{
    TensorEvolution* nObj = new DM04_alpha_Eij(this->e0_index,
                                               this->e_r_index,
                                               this->lambda_c_index,
                                               this->xi_index,
                                               this->Pat_index,
                                               this->m_index,
                                               this->M_cal_index,
                                               this->cc_index,
                                               this->nb_index,
                                               this->h0_index,
                                               this->ch_index,
                                               this->G0_index,
                                               this->alpha_index,
                                               this->z_index);
    return nObj;
}

const straintensor& DM04_alpha_Eij::Hij(const straintensor& plastic_flow, const stresstensor& Stre, 
                                        const straintensor& Stra, const MaterialParameter& material_parameter)
{
    stresstensor alpha_alpha_in;
    double a_in = 0.0;
    double h = 0.0;
    
    double e0 = gete0(material_parameter);
    double e_r = gete_r(material_parameter);
    double lambda_c = getlambda_c(material_parameter);
    double xi = getxi(material_parameter);
    double Pat = getPat(material_parameter);
    double m = getm(material_parameter);
    double M_cal = getM_cal(material_parameter);
    double cc = getcc(material_parameter);
    double nb = getnb(material_parameter);
    double h0 = geth0(material_parameter);
    double ch = getch(material_parameter);
    double G0 = getG0(material_parameter);        

    stresstensor alpha = getalpha(material_parameter);
    stresstensor z = getz(material_parameter);

    double p = Stre.p_hydrostatic();
    stresstensor s = Stre.deviator();

    stresstensor n;
    stresstensor alpha_b;
    stresstensor alpha_b_alpha;
    double b0 = 0.0;
    double g = 0.0;
    double ec = e_r;
    double stateParameter = 0.0;
    double expnb = 1.0;
    double ab = 0.0;
    stresstensor s_bar;
    double _s_bar_ = 0.0;

    double J3D;
    double cos3theta = 0.0;
    
    double e = e0 + (1.0 + e0) *Stra.Iinvariant1();

    //if (p != 0.0 && m != 0.0)
    //    n = (s * (1.0/p) - alpha) *(1.0/(sqrt(2.0/3.0)*m));
    
    s_bar = Stre.deviator() - (alpha *p);
    _s_bar_ = sqrt( (s_bar("ij")*s_bar("ij")).trace() );
    //cout << "_s_bar_ = " << _s_bar_ << endl;
    if (_s_bar_ > 0.0)
        n = s_bar * (1.0/_s_bar_);

    //double J2D = n.Jinvariant2();
    //J3D = n.Jinvariant3();
    //if ( J2D != 0.0)
      //cos3theta = -1.5*sqrt(3.0) *J3D / pow(J2D, 1.5);  
    
    J3D = n.Jinvariant3();
    cos3theta = -3.0*sqrt(6.0) *J3D;
    //cout << "cos3theta = " << cos3theta << endl;

    if (cos3theta > 1.0) cos3theta = 1.0;
    if (cos3theta < -1.0) cos3theta = -1.0;

    g = getg(cc, cos3theta);
    //cout << "g = " << g << endl;
    
    if ( (p/Pat) >= 0.0 )
      ec = getec(e_r, lambda_c, xi, Pat, p);

    stateParameter = e - ec;

    expnb = exp( -nb *stateParameter );
    ab = sqrt(2.0/3.0) * ( g *M_cal *expnb - m);
    alpha_b  = n *ab;
    alpha_b_alpha = alpha_b - alpha;
    
    if ( (p/Pat) >= 0.01 )
      b0 = G0 *h0 *(1.0-ch*e) *pow(p/Pat, -0.5);
    else
      b0 = G0 *h0 *(1.0-ch*e) *10.0;
    //cout << "b0 = " << b0 << endl;

    if (a_index == 0) {
        alpha_in = alpha;
        a_in = 0.0;
        a_index = 1;
        h = 1.0e10 * b0;
    }
    else {
        alpha_alpha_in = alpha - alpha_in;
        a_in = (alpha_alpha_in("ij")*n("ij")).trace();
        if ( a_in < 0.0 )
            a_index = 0;
        if ( a_in < 1.0e-10 )
            a_in = 1.0e-10;
        h = b0 /a_in;
    }
   
    TensorEvolution::TensorEvolutionHij = alpha_b_alpha *(h*2.0/3.0);
     
    return TensorEvolution::TensorEvolutionHij;
}

// to get e0
//================================================================================
double DM04_alpha_Eij::gete0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e0_index);
}

// to get e_r
//================================================================================
double DM04_alpha_Eij::gete_r(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, e_r_index);
}

// to get lambda_c
//================================================================================
double DM04_alpha_Eij::getlambda_c(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, lambda_c_index);
}

// to get xi
//================================================================================
double DM04_alpha_Eij::getxi(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, xi_index);

}

// to get Pat
//================================================================================
double DM04_alpha_Eij::getPat(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, Pat_index);
}

// to get m
//================================================================================
double DM04_alpha_Eij::getm(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, m_index);
}

// to get M
//================================================================================
double DM04_alpha_Eij::getM_cal(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, M_cal_index);
}

// to get c
//================================================================================
double DM04_alpha_Eij::getcc(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, cc_index);
}

// to get n_d
//================================================================================
double DM04_alpha_Eij::getnb(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, nb_index);
}


// to get h0
//================================================================================
double DM04_alpha_Eij::geth0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, h0_index);

}


// to get ch
//================================================================================
double DM04_alpha_Eij::getch(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, ch_index);

}

// to get G0
//================================================================================
double DM04_alpha_Eij::getG0(const MaterialParameter& material_parameter) const
{
    return getParameters(material_parameter, G0_index);

}

// to get alpha
//================================================================================
const stresstensor& DM04_alpha_Eij::getalpha(const MaterialParameter& material_parameter) const
{
	if ( alpha_index <= material_parameter.getNum_Internal_Tensor() && alpha_index > 0) {
		DM04_alpha_Eij::DM04_alpha_t = material_parameter.getInternal_Tensor(alpha_index-1);
		return DM04_alpha_Eij::DM04_alpha_t;
	}
	else {
		cout << "DM04_alpha: Invalid Input (alpha) " << endl;
		exit (1);
	}
}

// to get z
//================================================================================
const stresstensor& DM04_alpha_Eij::getz(const MaterialParameter& material_parameter) const
{
    if ( z_index <= material_parameter.getNum_Internal_Tensor() && z_index > 0) {
		DM04_alpha_Eij::DM04_alpha_t = material_parameter.getInternal_Tensor(z_index-1);
		return DM04_alpha_Eij::DM04_alpha_t;
	}
	else {
		cout << "DM04_alpha: Invalid Input (z) " << endl;
		exit (1);
	}
}

//================================================================================
double DM04_alpha_Eij::getParameters(const MaterialParameter& material_parameter, int which) const
{
	if ( which <= material_parameter.getNum_Material_Parameter() && which > 0)
		return material_parameter.getMaterial_Parameter(which-1);
	else {
		cout << "DM04_alpha: Invalid Input - #" << which << endl;
		exit (1);
	}
} 


//================================================================================
double DM04_alpha_Eij::getec(double e_r, double lambda_c, double xi, double Pat, double p_c) const
{
    double ee = e_r; 
    
    if ( p_c/Pat >= 0.0 )
      ee = e_r - lambda_c * pow(p_c/Pat, xi);
    else 
      cout << "Warning: DM04_alpha_Eij - 'p_c/Pat' less than zero! " << endl;

    return ee;
}

//================================================================================
double DM04_alpha_Eij::getg(double c, double cos3theta) const
{
    return 2.0 * c / ( (1.0+c) - (1.0-c)*cos3theta );
}


#endif

