/*
//================================================================================
# COPYRIGHT (C):     :-))                                                        #
# PROJECT:           Object Oriented Finite Element Program                      #
# PURPOSE:           General platform for elasto-plastic constitutive model      #
#                    implementation                                              #
# CLASS:             PS (the base class for all Potential Surfaces)              #
#                                                                                #
# VERSION:                                                                       #
# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )    #
# TARGET OS:         DOS || UNIX || . . .                                        #
# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                                 #
# PROGRAMMER(S):     Boris Jeremic, Zhao Cheng                                   #
#                                                                                #
# Date:           Jan 2004                                                       #
# UPDATE HISTORY:                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
#                                                                                #
# SHORT EXPLANATION: This is the base class of all Potential surfaces.           #
#                                                                                #
//================================================================================
*/

#ifndef PS_CPP
#define PS_CPP

#include "PS.h"

// At least 4 scalar and/or tensor internal variables are allowed at current time

//================================================================================
// The d(dQ/dsigma_ij)/ds1
//================================================================================
tensor PotentialSurface::d2Qodsds1( const EPState *EPS ) const
{ 
     tensor temp(2,def_dim_2,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/ds2
//================================================================================
tensor PotentialSurface::d2Qodsds2( const EPState *EPS ) const
{ 
     tensor temp(2,def_dim_2,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/ds3
//================================================================================
tensor PotentialSurface::d2Qodsds3( const EPState *EPS ) const
{ 
     tensor temp(2,def_dim_2,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/ds4
//================================================================================
tensor PotentialSurface::d2Qodsds4( const EPState *EPS ) const
{ 
     tensor temp(2,def_dim_2,0.0);
     return temp;
}
          
//================================================================================
// The d(dQ/dsigma_ij)/dt1_mn
//================================================================================
tensor PotentialSurface::d2Qodsdt1( const EPState *EPS ) const
{ 
     tensor temp(4,def_dim_4,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/dt2_mn
//================================================================================
tensor PotentialSurface::d2Qodsdt2( const EPState *EPS ) const
{ 
     tensor temp(4,def_dim_4,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/dt3_mn
//================================================================================
tensor PotentialSurface::d2Qodsdt3( const EPState *EPS ) const
{ 
     tensor temp(4,def_dim_4,0.0);
     return temp;
}

//================================================================================
// The d(dQ/dsigma_ij)/dt4_mn
//================================================================================
tensor PotentialSurface::d2Qodsdt4( const EPState *EPS ) const
{ 
     tensor temp(4,def_dim_4,0.0);
     return temp;
}

#endif
