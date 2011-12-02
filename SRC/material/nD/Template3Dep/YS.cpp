/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elaso-plastic constitutive model       #
#                    implementation                                              #
# CLASS:             YS (the base class for all Yield Surface)                   #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                                 #
#                                                                                #
#                                                                                #
# DATE:              08-03-2000                                                  #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is the base class of all yield surface. Here defined   #
#                    are initial values of the first derivatives of dF/ds_i 	 #
#                    and dF/ot_i                                                 # 
#                                                                                #
//================================================================================
*/

#ifndef YS_CPP
#define YS_CPP

#include "YS.h"


//================================================================================
// First derivative of F over scalar internal variables 
// (at most 4 scalar internal vars allowed currently)
//================================================================================

//================================================================================
//First derivative of F over the 1st scalar var
//================================================================================

double YieldSurface::xi_s1( const EPState *EPS ) const
{ 
     return 0.0;
}

//================================================================================
//First derivative of F over the 2nd scalar var
//================================================================================

double YieldSurface::xi_s2( const EPState *EPS ) const
{ 
     return 0.0;
}

//================================================================================
//First derivative of F over the 3rd scalar var
//================================================================================

double YieldSurface::xi_s3( const EPState *EPS ) const
{ 
     return 0.0;
}

//================================================================================
//First derivative of F over the 4th scalar var
//================================================================================

double YieldSurface::xi_s4( const EPState *EPS ) const
{ 
     return 0.0;
}

//================================================================================
// First derivative of F over scalar internal variables 
// (at most 4 tensorial internal vars allowed currently)
//================================================================================

//================================================================================
//First derivative of F over the 1st scalar var
//================================================================================
tensor YieldSurface::xi_t1( const EPState *EPS ) const
{ 
     stresstensor temp;
  //   cout << "inside YieldSurface::xi_t1( const EPState *EPS ) const  " << temp;
     return temp;
}
//================================================================================
//First derivative of F over the 2nd scalar var
//================================================================================
tensor YieldSurface::xi_t2( const EPState *EPS ) const
{ 
     stresstensor temp;
//     cout << " inside YieldSurface::xi_t2 "   <<temp;
     return temp;
}

//================================================================================
//First derivative of F over the 3rd scalar var
//================================================================================
tensor YieldSurface::xi_t3( const EPState *EPS ) const
{ 
     stresstensor temp;
//     cout << "inside  YieldSurface::xi_t3 " <<temp;
     return temp;
}

//================================================================================
//First derivative of F over the 4th scalar var
//================================================================================
tensor YieldSurface::xi_t4( const EPState *EPS ) const
{ 
     stresstensor temp;
    // cout << " inside YieldSurface::xi_t4 "  << temp;
     return temp;
}

//================================================================================
// friend OPS_Stream functions for output
//================================================================================
OPS_Stream& operator<<(OPS_Stream& os, const YieldSurface & YS)
{
       os << "Yield Surface Parameters: " << endln;
       return os;
}


#endif

