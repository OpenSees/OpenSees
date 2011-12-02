                                                                        
// $Revision: 1.2 $                                                              
// $Date: 2000-10-20 03:51:47 $                                                                  
// $Source: /usr/local/cvs/OpenSees/SRC/matrix/BJvector.h,v $                                                                
                                                                        
//#############################################################################
//                                                                            #
//                                                                            #
//             /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/~~\                #
//            |                                          |____|               #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |        B A S E   C L A S S E S           |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |          C + +     H E A D E R           |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//            |                                          |                    #
//         /~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~/   |                    #
//         \_________________________________________\__/                     #
//                                                                            #
//                                                                            #
//#############################################################################
//#############################################################################
///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:             BJvector                                                    #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              November '92                                              #
//# UPDATE HISTORY:    05 - __ avgust '93.  redefined as derived class from      #
//#                                 nDarray class                                #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/


#ifndef VECTOR_HH
#define VECTOR_HH

#include "BJmatrix.h"

// All of this inheritance idioms are after
// Jim Coplien : "Advanced C++ programing styles and idioms".
// I tried to understand idioms and I think I succeded.

class BJvector : virtual public BJmatrix
  {
    public:
      BJvector(int order_n = 1, double initvalue = 0.0);  // default constructor

      BJvector(int order_n, double *initval);

      BJvector(const nDarray & x); // copy-initializer

//....       ~BJvector( );

      BJvector& operator=( const BJvector & x ); // BJvector assignment
//..      BJvector& operator=( const BJmatrix & x ); // BJvector assignment
//..      BJvector& operator=( const nDarray & x ); // BJvector assignment

//#######        BJmatrix operator*( BJvector &); // BJvector multiplication

//....     BJvector operator*( double arg); // scalar multiplication 
// this  ellipsis at the end are just to prevent the compiler
// from issuing a warning on hiding function from base class nDarray . . . 
     double & val(int subscript, ... );
     double cval(int subscript, ... ) const; // const
// THE ROW COUNTER STARTS FROM 1 ( NOT FROM 0 )

};
#endif
