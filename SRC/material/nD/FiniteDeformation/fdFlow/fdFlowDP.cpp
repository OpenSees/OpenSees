//===============================================================================
//# COPYRIGHT (C): Woody's license (by BJ):
//                 ``This    source  code is Copyrighted in
//                 U.S.,  for  an  indefinite  period,  and anybody
//                 caught  using it without our permission, will be
//                 mighty good friends of ourn, cause we don't give
//                 a  darn.  Hack it. Compile it. Debug it. Run it.
//                 Yodel  it.  Enjoy it. We wrote it, that's all we
//                 wanted to do.''
//
//# PROJECT:           Object Oriented Finite Element Program
//# PURPOSE:           Finite Deformation Hyper-Elastic classes
//# CLASS:
//#
//# VERSION:           0.6_(1803398874989) (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic
//#
//#
//# DATE:              July 2004
//# UPDATE HISTORY:
//#
//===============================================================================

#ifndef fdFlowDP_CPP
#define fdFlowDP_CPP

#include "fdFlowDP.h"

//--------------------------------------------------------------------
fdFlowDP::fdFlowDP(double DilatedAngle_in, double Cohesion_in, int ConeIndex_in)
:DilatedAngle(DilatedAngle_in), Cohesion(Cohesion_in), ConeIndex(ConeIndex_in)
{
   double pipi = 3.14159265358979323846;
   double Root3 = sqrt(3.0);
   double Root1o3 = 1.0/Root3;
   double archAngle = DilatedAngle*pipi/180.0; 

   switch (ConeIndex) {

     case 0: { // Compressive (Outer) Cone
       k1 = 2.0*Root1o3*sin(archAngle)/(3.0-sin(archAngle));
       k2 = 6.0*Root1o3*cos(archAngle)/(3.0-sin(archAngle));
       break;
     }

     case 1: { // Tensile (Inner) Cone
       k1 = 2.0*Root1o3*sin(archAngle)/(3.0+sin(archAngle));
       k2 = 6.0*Root1o3*cos(archAngle)/(3.0+sin(archAngle));
       break;
     }

     case 2: { // Mean  Cone
       k1 = Root3*sin(archAngle)/(9.0-sin(archAngle)*sin(archAngle));
       k2 = 2.0*Root3*cos(archAngle)/(9.0-sin(archAngle)*sin(archAngle));
       break;
     }

     case 3: { // Inner-tangent  Cone
       k1 = tan(archAngle)/sqrt(9.0+12.0*tan(archAngle)*tan(archAngle));
       k2 = 3.0/sqrt(9.0+12.0*tan(archAngle)*tan(archAngle));
       break;
     }   

     default: { // Compressive (Outer) Cone
       k1 = 2.0*Root1o3*sin(archAngle)/(3.0-sin(archAngle));
       k2 = 6.0*Root1o3*cos(archAngle)/(3.0-sin(archAngle));
     }
     
   }
}

//--------------------------------------------------------------------
fdFlow * fdFlowDP::newObj()
{
     fdFlow *newfdyd = new fdFlowDP(DilatedAngle, Cohesion, ConeIndex);
     return newfdyd;
}

//-------------------------------------------------------------------
// Q = k1*I1 + sqrt(0.5*Sij*Sji) - k2*(c+q) = 0, 
// Note here NumRank = 1, No Kinematic Hardening
// Yd = (3/2)(sij-p_bar*alpha_ij)(sij-p_bar*alpha_ij) - M^2*(p_bar)^2 = 0, p_bar = p + c0/tan() + c
// Note here NumRank = 2: with Kinematic hardening
//-------------------------------------------------------------------

//--------------------------------------------------------------------
stresstensor fdFlowDP::dFods(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    // NumRank=1, No Ki Hardeing
    //tensor tI2("I", 2 , def_dim_2);
    //stresstensor dev = sts.deviator();
    //tensor st = dev("ij")*dev("ij");  
    //  st.null_indices();
    //double x = st.trace();
    //return tI2*k1 + dev *(1.0/sqrt(2.0*x));

    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    tensor tI2("I", 2, def_dim_2);
    double root3 = sqrt(3.0);
    double M = 3.0*root3*k1;
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    tensor aijaij = a("ij")*a("ij");
      aijaij.null_indices();
    double aa = aijaij.trace();
    stresstensor w  = dev - a *pb;
    tensor wijaij = w("ij")*a("ij");
      wijaij.null_indices();
    double wa = wijaij.trace();
    return w*3.0 + tI2*(pb*aa +wa + (2.0/3.0)*M*M*pb);
}

//--------------------------------------------------------------------
double fdFlowDP::dFodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=1, No Ki Hardeing
    //return -k2;

    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    double root3 = sqrt(3.0);
    double M = 3.0*root3*k1;
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    stresstensor w  = dev - a *pb;
    tensor wijaij = w("ij")*a("ij");
      wijaij.null_indices();
    double wa = wijaij.trace();

    return wa*(-3.0) + pb*(-2.0*M*M);
}

//--------------------------------------------------------------------
stresstensor fdFlowDP::dFoda(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    stresstensor w  = dev - a*pb;

    return w*(-3.0*pb);
}

//--------------------------------------------------------------------
tensor fdFlowDP::d2Fodsds(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=1, No Ki Hardeing
    //tensor tI2("I", 2, def_dim_2);
    //stresstensor dev = sts.deviator();
    //tensor st = dev("ij")*dev("ij")*2.0;  
    //  st.null_indices();
    //double x = st.trace();
    //tensor I4 = tI2("ij")*tI2("kl"); 
    //  I4.null_indices();
    //I4 = I4 - I4*(1.0/3.0);  
    //tensor st4 = dev("ij")*dev("kl");
    //  st4.null_indices();
    //
    //return I4 *(1.0/sqrt(x)) - st4 *(2.0/x/sqrt(x));

    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    tensor tI2("I", 2, def_dim_2);
    double root3 = sqrt(3.0);
    double M = 3.0*root3*k1;
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    stresstensor w  = dev - a *pb;
    double akk = a.Iinvariant1();
    tensor aijaij = a("ij")*a("ij");
      aijaij.null_indices();
    double aa = aijaij.trace();
    tensor tI40 = tI2("ij")*tI2("kl");
      tI40.null_indices();	
    tensor tI41 = tI40.transpose0110(); // for non-symmetric
    tensor atI = a("ij")*tI2("kl");
      atI.null_indices();
    tensor tIa = tI2("ij")*a("kl");
      tIa.null_indices();

    return tI41*3.0 + tI40*((aa- akk)/3.0 - M*M*2.0/9.0 -1.0) + atI + tIa;    
}

//--------------------------------------------------------------------
stresstensor fdFlowDP::d2Fodsdq(const stresstensor &sts, const FDEPState &fdepstate ) const
{
    // NumRank=2, with Ki Hardeing
    tensor tI2("I", 2, def_dim_2);
    double root3 = sqrt(3.0);
    double M = 3.0*root3*k1;
    stresstensor a = fdepstate.getStressLikeKiVar();
    double akk = a.Iinvariant1();
    tensor aijaij = a("ij")*a("ij");
      aijaij.null_indices();
    double aa = aijaij.trace();

    return a*(-3.0) + tI2*(akk - aa + M*M*2.0/3.0); 
}

//--------------------------------------------------------------------
tensor fdFlowDP::d2Fodsda(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    tensor tI2("I", 2, def_dim_2);
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    stresstensor w  = dev - a *pb;
    tensor tI40 = tI2("ij")*tI2("kl");
      tI40.null_indices();	
    tensor tI41 = tI40.transpose0110(); // for non-symmetric
    tensor wtI = w("ij")*tI2("kl");
      wtI.null_indices();
    tensor atI = a("ij")*tI2("kl");
      atI.null_indices();

    return tI41*(-3.0*pb) + tI40*pb + atI*(-pb) + wtI;    
}
    
//--------------------------------------------------------------------
double fdFlowDP::d2Fodqdq(const stresstensor &sts, const FDEPState &fdepstate ) const
{
    // NumRank=2, with Ki Hardeing
    double root3 = sqrt(3.0);
    double M = 3.0*root3*k1;
    stresstensor a = fdepstate.getStressLikeKiVar();
    tensor aijaij = a("ij")*a("ij");
      aijaij.null_indices();
    double aa = aijaij.trace();

    return aa*3.0 - M*M*2.0; 
}

//--------------------------------------------------------------------
stresstensor fdFlowDP::d2Fodqda(const stresstensor &sts, const FDEPState &fdepstate ) const
{
    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    tensor tI2("I", 2, def_dim_2);
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    stresstensor w  = dev - a*pb;

    return a*(3.0*pb) - w*3.0; 
}

tensor fdFlowDP::d2Fodada(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=2, with Ki Hardeing
    double pipi = 3.14159265358979323846;
    double archAngle = DilatedAngle*pipi/180.0;
    tensor tI2("I", 2, def_dim_2);
    double cv = fdepstate.getStressLikeInVar();
    if ( tan(archAngle) < 1.0e-4 ) {
      opserr << "fdFlowDP-- (close) zero or invalid dilation angle, try VM flow rule\n";
      exit (-1);      
    }
    double pa = Cohesion/tan(archAngle) + cv;
    stresstensor a = fdepstate.getStressLikeKiVar();
    stresstensor dev = sts.deviator();
    double p = sts.p_hydrostatic();
    double pb = p + pa;
    tensor tI40 = tI2("ij")*tI2("kl");
      tI40.null_indices();	
    tensor tI41 = tI40.transpose0110(); // for non-symmetric

    return tI41*(3.0*pb*pb);
}

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdFlowDP &fdfdDP)
{
    os << "fdFlowDP Parameters: " << "\n";
    os << "DilatedAngle: " << fdfdDP.DilatedAngle << "\n";
    os << "Cohesion: " << fdfdDP.Cohesion << "\n";
    os << "ConeIndex: " << fdfdDP.ConeIndex << "\n";
    return os;
}


#endif

