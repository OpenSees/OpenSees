///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           strain tensor with all necessery functions                #
//# CLASS:             Cosseratstraintensor                                      #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Alireza Tabarrei, Boris Jeremic                           #
//# PROGRAMMER(S):     Alireza Tabarrei,Boris Jeremic                            #
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
#ifndef COSSERATSTRAINTENSOR_HH
#define COSSERATSTRAINTENSOR_HH

#include "BJtensor.h"

class Cosseratstraintensor : public tensor
{
  public: // just send appropriate arguments to the base constructor

//    Cosseratstraintensor (int rank_of_tensor=2, double initval=0.00000003141528);
    Cosseratstraintensor (int rank_of_tensor=2, double initval=0.0);
// default constructor           // this is just PI/10^8 to check default constructor

    Cosseratstraintensor ( double *values );

    Cosseratstraintensor ( double initvalue );

    Cosseratstraintensor(const Cosseratstraintensor & x );
    Cosseratstraintensor(const tensor & x); // copy-initializer
    Cosseratstraintensor(const nDarray & x); // copy-initializer

    //~Cosseratstraintensor( );
    
    Cosseratstraintensor operator=(const Cosseratstraintensor & rval); // Cosseratstraintensor assignment
    Cosseratstraintensor operator=(const tensor & rval);// tensor assignment to Cosseratstraintensor
    Cosseratstraintensor operator=(const nDarray & rval);// nDarray assignment to Cosseratstraintensor

    Cosseratstraintensor deep_copy(void);
//..    Cosseratstraintensor * p_deep_copy(void);

//ini  // use "from" and initialize already allocated strain tensor from "from" values
//ini      void Initialize( const Cosseratstraintensor & from );

//___// operator() overloading for 3D Gauss points!
//___    Cosseratstraintensor & operator()(short ir, short is, short it,
//___                              short tr, short ts, short tt  );
 
    double Iinvariant1( ) const;
    double Iinvariant2( ) const;
    double Iinvariant3( ) const;

    double Jinvariant1( ) const;
    double Jinvariant2( ) const;
    double Jinvariant3( ) const;

    double equivalent( ) const;	  //Zhaohui added 09-02-2000

    Cosseratstraintensor deviator( ) const;
    Cosseratstraintensor principal( ) const;

    double sigma_octahedral( ) const;
    double tau_octahedral( ) const;

    double ksi( ) const;
    double ro( ) const;
    double theta( ) const;
    double thetaPI( ) const;

    double p_hydrostatic( ) const;
    double q_deviatoric( ) const;


    Cosseratstraintensor pqtheta2strain( double, double, double );
    Cosseratstraintensor evoleq2strain( double, double );

    void report(char *) const;
    void reportshort(char *) const;

//..// polinomial root solver friend functions definitions
//..public:
//..friend void laguer(complex *, int , complex *, double , int );
//..friend void zroots(complex *, int , complex *, int );
//..
};

#endif

