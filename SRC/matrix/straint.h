///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           strain tensor with all necessery functions                #
//# CLASS:             straintensor                                              #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              July 25 '93                                               #
//# UPDATE HISTORY:    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC,        #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/
//
#ifndef STRAINTENSOR_HH
#define STRAINTENSOR_HH

#include "BJtensor.h"

class straintensor : public tensor
{
  public: // just send appropriate arguments to the base constructor

//    straintensor (int rank_of_tensor=2, double initval=0.00000003141528);
    straintensor (int rank_of_tensor=2, double initval=0.0);
// default constructor           // this is just PI/10^8 to check default constructor

    straintensor ( double *values );

    straintensor ( double initvalue );

    straintensor(const straintensor & x );
    straintensor(const tensor & x); // copy-initializer
    straintensor(const nDarray & x); // copy-initializer

//....     ~straintensor( );

    straintensor operator=(const straintensor & rval); // straintensor assignment
    straintensor operator=(const tensor & rval);// tensor assignment to straintensor
    straintensor operator=(const nDarray & rval);// nDarray assignment to straintensor

    straintensor deep_copy(void);
//..    straintensor * p_deep_copy(void);

//ini  // use "from" and initialize already allocated strain tensor from "from" values
//ini      void Initialize( const straintensor & from );

//___// operator() overloading for 3D Gauss points!
//___    straintensor & operator()(short ir, short is, short it,
//___                              short tr, short ts, short tt  );
 
    double Iinvariant1( ) const;
    double Iinvariant2( ) const;
    double Iinvariant3( ) const;

    double Jinvariant1( ) const;
    double Jinvariant2( ) const;
    double Jinvariant3( ) const;

    straintensor deviator( ) const;
    straintensor principal( ) const;

    double sigma_octahedral( ) const;
    double tau_octahedral( ) const;

    double ksi( ) const;
    double ro( ) const;
    double theta( ) const;
    double thetaPI( ) const;

    double p_hydrostatic( ) const;
    double q_deviatoric( ) const;


    straintensor pqtheta2strain( double, double, double );
    straintensor evoleq2strain( double, double );

    void report(char *) const;
    void reportshort(char *) const;

//..// polinomial root solver friend functions definitions
//..public:
//..friend void laguer(complex *, int , complex *, double , int );
//..friend void zroots(complex *, int , complex *, int );
//..
};

#endif

