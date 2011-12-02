
// $Revision: 1.1 $
// $Date: 2001-08-23 16:45:49 $
// $Source: /usr/local/cvs/OpenSees/SRC/nDarray/basics.h,v $

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
//##
///*
//################################################################################
//# COPY-YES  (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:                                                                     #
//# CLASS:                                                                       #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.10, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #                                                                              #
//# DATE:              November '92                                              #
//# UPDATE HISTORY:    05 - __ avgust '93.  redefined as derived class from      #
//#                                 nDarray class                                #
//#                    january 06 '93  added matrix2BJtensor_1, matrix2BJtensor_2    #
//#                                   matrix2BJtensor_3                            #
//#                    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//################################################################################
//*/
#ifndef BASICS_HH
#define BASICS_HH

#include <math.h>
//#include <values.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>


// Define BJmatrix as matrix type
#ifndef matrix
#define matrix BJmatrix
#endif
// Define BJvector as vector type
#ifndef vector
#define vector BJvector
#endif
// Define BJtensor as tensor type
#ifndef tensor
#define tensor BJtensor
#endif

#ifndef tensor
#define Tensor BJtensor
#endif

// redefine Pi from math.h M_PI
#ifndef PI
#define PI 3.14159265358979323846
#endif
//##############################################################################
#ifndef TWOOVERTHREE
#define TWOOVERTHREE 0.6666666666667
#endif
//##############################################################################
#ifndef ONEOVERTHREE
#define ONEOVERTHREE 0.3333333333333
#endif
//##############################################################################
//##############################################################################
// usefull arrays for constructors . . .
//#ifdef SASA
 static const int def_dim_4_2[]={2,2,2,2}; //  Sasa jan - 99
//#endif

#ifndef DEF_DIM
#define DEF_DIM
  static const int def_dim_1[] = {3};
  static const int def_dim_2[] = {3, 3}; // static-> see ARM pp289-290
  static const int def_dim_3[] = {3, 3, 3}; // static-> see ARM pp289-290
  static const int def_dim_4[] = {3, 3, 3, 3}; // static-> see ARM pp289-290
#endif

#ifndef TEST
#define TEST
  const int tst = 3;
#endif


#ifndef DZERO
#define DZERO 0.0
//  double ZERO = 0.0;
#endif

float       f_macheps();
double      d_macheps();
long double ld_macheps();

//double min(double , double );
//double max(double , double );
//
//int min(int , int );
//int max(int , int );

#endif

