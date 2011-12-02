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
fdFlowDP::fdFlowDP(double DilatedAngle_in, int ConeIndex_in)
:DilatedAngle(DilatedAngle_in), ConeIndex(ConeIndex_in)
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
     fdFlow *newfdyd = new fdFlowDP(DilatedAngle, ConeIndex);
     return newfdyd;
}

//-------------------------------------------------------------------
// Q = k1*I1 + sqrt(0.5*Sij*Sji) - k2*(c+q) = 0, 
// Note here NumRank = 1, No Kinematic Hardening
//-------------------------------------------------------------------

//--------------------------------------------------------------------
stresstensor fdFlowDP::dFods(const stresstensor &sts, const FDEPState &fdepstate) const
{    
    // NumRank=1, No Ki Hardeing
    tensor tI2("I", 2 , def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij");  
      st.null_indices();
    double x = st.trace();
    return tI2*k1 + dev *(1.0/sqrt(2.0*x));
}

//--------------------------------------------------------------------
double fdFlowDP::dFodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=1, No Ki Hardeing
    return -k2;
}

//--------------------------------------------------------------------
tensor fdFlowDP::d2Fodsds(const stresstensor &sts, const FDEPState &fdepstate ) const
{  
    // NumRank=1, No Ki Hardeing
    tensor tI2("I", 2, def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij")*2.0;  
      st.null_indices();
    double x = st.trace();
    tensor I4 = tI2("ij")*tI2("kl"); 
      I4.null_indices();
    I4 = I4 - I4*(1.0/3.0);  
    tensor st4 = dev("ij")*dev("kl");
      st4.null_indices();

    return I4 *(1.0/sqrt(x)) - st4 *(2.0/x/sqrt(x));
    
}

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdFlowDP &fdfdDP)
{
    os << "fdFlowDP Parameters: " << "\n";
    os << "DilatedAngle: " << fdfdDP.DilatedAngle << "\n";
    os << "ConeIndex: " << fdfdDP.ConeIndex << "\n";
    return os;
}


#endif

