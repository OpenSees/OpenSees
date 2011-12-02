///*
//################################################################################
//# COPYRIGHT (C):     :-))                                                      #
//# PROJECT:           Object Oriented Finite Element Program                    #
//# PURPOSE:           stress tensor with all necessery functions                #
//# CLASS:             stresstensor                                              #
//#                                                                              #
//# VERSION:                                                                     #
//# LANGUAGE:          C++.ver >= 2.0 ( Borland C++ ver=3.00, SUN C++ ver=2.1 )  #
//# TARGET OS:         DOS || UNIX || . . .                                      #
//# DESIGNER(S):       Boris Jeremic                                             #
//# PROGRAMMER(S):     Boris Jeremic                                             #
//#                                                                              #
//#                                                                              #
//# DATE:              July 22 '93                                               #
//# UPDATE HISTORY:    August 22-29 '94 choped to separate files and worked on   #
//#                                   const and & issues                         #
//#                    August 30-31 '94 added use_def_dim to full the CC         #
//#                                   resolved problem with temoraries for       #
//#                                   operators + and - ( +=, -= )               #
//#                    13 septembar '96 added reportAnim  :)                     #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef STRESSTENSOR_H
#define STRESSTENSOR_H

#include <OPS_Globals.h>
#include "BJtensor.h"
class Material_Model;


class stresstensor : public BJtensor
{
  public:
    friend class Material_Model;

  public:
    // just send appropriate arguments to the base constructor
    stresstensor(int rank_of_tensor=2, double initval=0.0); // default constructor
    stresstensor( double *values );
    stresstensor( double initvalue );

    stresstensor(const stresstensor & x );
    stresstensor(const BJtensor & x); // copy-initializer
    stresstensor(const nDarray & x); // copy-initializer

    //~stresstensor( );
    

    stresstensor operator=(const stresstensor & rval);// stresstensor assignment
    stresstensor operator=(const BJtensor & rval);// tensor assignment to stresstensor
    stresstensor operator=(const nDarray & rval);// nDarray assignment to stresstensor

    stresstensor deep_copy(void);
    //..    stresstensor * p_deep_copy(void);

    //ini  // use "from" and initialize already allocated stress tensor from "from" values
    //ini      void Initialize( const stresstensor & from );

    //___// operator() overloading for 3D Gauss points!
    //___    stresstensor & operator()(short ir, short is, short it,
    //___                              short tr, short ts, short tt  );
    

    double Iinvariant1( ) const;
    double Iinvariant2( ) const;
    double Iinvariant3( ) const;

    double Jinvariant1( ) const;
    double Jinvariant2( ) const;
    double Jinvariant3( ) const;

    stresstensor deviator( ) const;
    stresstensor principal( ) const;

    double sigma_octahedral( ) const;
    double tau_octahedral( ) const;

    double ksi( )     const;
    double xi( )      const;
    double ro( )      const;
    double rho( )      const;
    double theta()   const;
    double thetaPI( ) const;

    double p_hydrostatic( ) const;
    double q_deviatoric( ) const;

    tensor dpoverds( void  ) const;
    tensor dqoverds( void ) const;
    tensor dthetaoverds( void ) const;
    tensor d2poverds2( void ) const;
    tensor d2qoverds2( void  ) const;
    tensor d2thetaoverds2( void ) const;
	     	          


    //--    stresstensor yield_surface_cross(stresstensor & end_stress,
    //--                                     Material_Model & YC);

    stresstensor pqtheta2stress( double, double, double );

    void report(char *) const;
    void reportshort(char *) const;
    void reportshortpqtheta(char *) const;
    void reportSHORTpqtheta(char *) const;
    void reportSHORTs1s2s3(char *) const;
    void reportKLOTpqtheta(char *) const;
    void reportshortI1J2J3(char *) const;
    void reportAnim(void) const;
    void reportTensor(char *) const;

    //================================================================================
    // Overloaded Insertion Operator	  ZHaohui Added Aug. 13, 2000
    // prints an stresstensor's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const stresstensor & rhs);

    //  // routine used by root finder, takes an alfa and returns the
    //  // yield function value for that alfa
    //    public:
    //      double func( stresstensor & start_stress,
    //                   stresstensor & end_stress,
    //                   Material_Model & YC,
    //                   double alfa );
    //  
    //  
    //  //..// polinomial root solver friend functions definitions
    //  //..public:
    //  //..friend void laguer(complex *, int , complex *, double , int );
    //  //..friend void zroots(complex *, int , complex *, int );
    //  //..
    //  
    // zero of function
    friend double zbrentstress(stresstensor   & start_stress,
                             stresstensor   & end_stress,
                             Material_Model & YC,
                             double x1, double x2, double tol);
  
    //  friend double zbrent(double x1, double x2, double tol);
    //  
    //  
};

#endif

