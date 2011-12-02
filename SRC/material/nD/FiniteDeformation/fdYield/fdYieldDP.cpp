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

#ifndef fdYieldDP_CPP
#define fdYieldDP_CPP

#include "fdYieldDP.h"

//--------------------------------------------------------------------
fdYieldDP::fdYieldDP(double FricAngle_in, double Cohension_in, int ConeIndex_in)
 : FricAngle(FricAngle_in), Cohension(Cohension_in), ConeIndex(ConeIndex_in)
{   
   double pipi = 3.14159265358979323846;
   double Root3 = sqrt(3.0);
   double Root1o3 = 1.0/Root3;
   double archAngle = FricAngle*pipi/180.0;
   
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
fdYield * fdYieldDP::newObj()
{
    fdYield *newfdyd = new fdYieldDP(FricAngle, Cohension, ConeIndex);
    return newfdyd;
}

int fdYieldDP::getNumRank()
{
    return 1;
}

double fdYieldDP::getTolerance()
{
    return 1.0e-8*Cohension > 1.0e-8?  1.0e-8*Cohension : 1.0e-8;
}

//--------------------------------------------------------------------------------------------
// Yd =  k1*I1 + sqrt(0.5*Sij*Sji) - k2*(c+q) = 0, 
// Note here NumRank = 1: No Kinematic hardening
//--------------------------------------------------------------------------------------------

double fdYieldDP::Yd(const stresstensor &sts, const FDEPState &fdepstate ) const
{    
    // NumRank=1, No Ki Hardeing
    double q = fdepstate.getStressLikeInVar();
    stresstensor dev = sts.deviator();
    double I1 = sts.Iinvariant1();
    tensor st = dev("ij")*dev("ij");
      st.null_indices();
    double x = st.trace();
    return k1*I1 + sqrt(0.5*x) - k2*(Cohension+q);
}

//--------------------------------------------------------------------
stresstensor fdYieldDP::dYods(const stresstensor &sts, const FDEPState &fdepstate ) const
{   
    // NumRank=1, No Ki Hardeing
    tensor tI2("I", 2, def_dim_2);
    stresstensor dev = sts.deviator();
    tensor st = dev("ij")*dev("ij");  
      st.null_indices();
    double x = st.trace();
    return tI2*k1 + dev*(1.0/sqrt(2.0*x));
}					      

//--------------------------------------------------------------------
double fdYieldDP::dYodq(const stresstensor &sts, const FDEPState &fdepstate ) const
{      
    // NumRank=1, No Ki Hardeing
    return -k2;
}

//--------------------------------------------------------------------
OPS_Stream& operator<<(OPS_Stream& os, const fdYieldDP &fdydDP)
{
    os << "fdYieldDP Parameters: " << "\n";
    os << "FricAngle: " << fdydDP.FricAngle << "\n";
    os << "Cohesion: " << fdydDP.Cohension << "\n";
    os << "ConeIndex: " << fdydDP.ConeIndex << "\n";
    return os;
}


#endif

