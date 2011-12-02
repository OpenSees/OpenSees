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

#ifndef fdYield_CPP
#define fdYield_CPP

#include "fdYield.h"

fdYield::fdYield ()
{ 

}

int fdYield::getNumRank( )
{ 
     return 0;
}

double fdYield::getTolerance( )
{ 
     return 1.0e-8;
}

double fdYield::dYodq(const stresstensor &sts,  const FDEPState &fdepstate ) const 
{ 
     return 0.0;
}

stresstensor fdYield::dYoda(const stresstensor &sts,  const FDEPState &fdepstate ) const 
{ 
     tensor tZ2(2, def_dim_2, 0.0);
     return tZ2;
}

OPS_Stream& operator<<(OPS_Stream& os, const fdYield & fdyd)
{
     os << "fdYield Parameters: " << "\n";
     return os;
}


#endif

