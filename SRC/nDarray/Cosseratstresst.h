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
//# DESIGNER(S):       Alireza Tabarrei                                          #
//# PROGRAMMER(S):     Alireza Tabarrei                                          #
//#                                                                              #
//#                                                                              #
//# DATE:              June 2004																																																 #
//# UPDATE HISTORY:    																																																									 #
//#                                                                              #
//#                                                                              #
//#                                                                              #
//################################################################################
//*/

#ifndef COSSERATSTRESSTENSOR_H
#define COSSERATSTRESSTENSOR_H

#include <OPS_Globals.h>
#include "BJtensor.h"
class Material_Model;

class Cosseratstresstensor : public BJtensor
{
  public:
    friend class Material_Model;

  public:
    // just send appropriate arguments to the base constructor
    Cosseratstresstensor(int rank_of_tensor=2, double initval=0.0); // default constructor
    Cosseratstresstensor( double *values );
    Cosseratstresstensor( double initvalue );

    Cosseratstresstensor(const Cosseratstresstensor & x );
    Cosseratstresstensor(const BJtensor & x); // copy-initializer
    Cosseratstresstensor(const nDarray & x); // copy-initializer

    //~Cosseratstresstensor( );
    

    Cosseratstresstensor operator=(const Cosseratstresstensor & rval);// Cosseratstresstensor assignment
    Cosseratstresstensor operator=(const BJtensor & rval);// tensor assignment to Cosseratstresstensor
    Cosseratstresstensor operator=(const nDarray & rval);// nDarray assignment to Cosseratstresstensor

    Cosseratstresstensor deep_copy(void);
    //..    Cosseratstresstensor * p_deep_copy(void);

    //ini  // use "from" and initialize already allocated stress tensor from "from" values
    //ini      void Initialize( const Cosseratstresstensor & from );

    //___// operator() overloading for 3D Gauss points!
    //___    Cosseratstresstensor & operator()(short ir, short is, short it,
    //___                              short tr, short ts, short tt  );
    

    double Iinvariant1( ) const;
    double Iinvariant2( ) const;
    double Iinvariant3( ) const;

    double Jinvariant1( ) const;
    double Jinvariant2( ) const;
    double Jinvariant3( ) const;

    Cosseratstresstensor deviator( ) const;
    Cosseratstresstensor principal( ) const;

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
	     	          


    //--    Cosseratstresstensor yield_surface_cross(Cosseratstresstensor & end_stress,
    //--                                     Material_Model & YC);

    Cosseratstresstensor pqtheta2stress( double, double, double );

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
    // prints an Cosseratstresstensor's contents 
    //================================================================================
    friend OPS_Stream& operator<< (OPS_Stream& os, const Cosseratstresstensor & rhs);

    //  // routine used by root finder, takes an alfa and returns the
    //  // yield function value for that alfa
    //    public:
    //      double func( Cosseratstresstensor & start_stress,
    //                   Cosseratstresstensor & end_stress,
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
    friend double zbrentstress(Cosseratstresstensor   & start_stress,
                             Cosseratstresstensor   & end_stress,
                             Material_Model & YC,
                             double x1, double x2, double tol);
  
    //  friend double zbrent(double x1, double x2, double tol);
    //  
    //  
};

#endif

